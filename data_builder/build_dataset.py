# -*- coding: utf-8 -*-
"""
Created on Tue Aug 08 11:32:50 2017

@author: csce
"""

import os
import argparse
import logging
import cPickle
from glob import glob

import numpy as np
import h5py

from collections import defaultdict
from dicom_parser import DicomParse
from annotation_parser import LIDC_XML_Parser

class DatasetBuilder():
    def __init__(self,dataset_name,dataset_dir,output_dir):
        self.dataset_name = dataset_name
        self.dataset_dir = dataset_dir
        self.output_dir = output_dir
        self.dicom_parser = DicomParse()
        if dataset_name =="LIDC":
            self.ann_parser = LIDC_XML_Parser()
        
        
        
    def build_spie_dataset(self):
        if not os.path.exists(self.output_dir):
            os.mkdir(self.output_dir)
        
        vol_dataset = h5py.File(os.path.join(self.output_dir,"vol.hdf5"),'x')
        spie2uid = dict()
        spie2img = dict()
        
        for patient_dir in os.listdir(self.dataset_dir):
            for study_dir in glob('%s/%s/*' % (self.dataset_dir, patient_dir)):
                for series_name in os.listdir(study_dir):
                    series_dir = os.path.join(study_dir, series_name)
                    
                    vol, sop, loc, spacing, slice_thickness, series_uid = self.dicom_parser.parse_dicom(series_dir)
                    spie2img[patient_dir] = sop
                    vol_data = vol_dataset.create_dataset(patient_dir, data=vol, dtype=np.int16,
                                                      chunks=(512, 512, 1), compression='gzip')
                    vol_data.attrs['spacing'] = spacing
                if slice_thickness is not None:
                    vol_data.attrs['slice_thickness'] = slice_thickness

                spie2uid[patient_dir] = series_uid
    
        with open(r"C:\Users\csce\Desktop\nova\lung_tumer\data\spie2uid.pkl",'a') as f:
            cPickle.dump(spie2uid,f)
        with open(r"C:\Users\csce\Desktop\nova\lung_tumer\data\spie2img.pkl",'a') as f:
            cPickle.dump(spie2img,f)
        vol_dataset.close()
        
    def build_lidc_dataset(self,lidc_dir, output_dir, version='v9'):
        """
        Build the dataset for lidc. Distinct version differs in the order of dicoms.
        :param lidc_dir: lidc original data, which could be downloaded free
        :param output_dir: target output directory to save the result
        :param version: 'v9' ordered by
        :return: None
        """
        if version == 'v9':
            order = 'Location'
        elif version == 'v15':
            order = 'InstanceNumber'
        else:
            raise ValueError('Unknown version')
    
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
    
        uid2lidc = np.load(r'C:\Users\csce\Desktop\nova\lung_tumer\data\LIDC-IDRI\uid2lidc.pkl')
        lidc2img = dict()
        label_dict = dict()
        vol_dataset = h5py.File(os.path.join(output_dir, 'vol.hdf5'), 'x')
    
        # main
        patient_dirs = os.listdir(lidc_dir)
        patient_dirs.sort()
        for patient_dir in patient_dirs:
            for study_dir in glob('%s/%s/*' % (lidc_dir, patient_dir)):
                for series_name in os.listdir(study_dir):
                    series_dir = os.path.join(study_dir, series_name)
                    if not self.dicom_parser.check_dcm_dir(series_dir):
                        continue
    
                    scan_name = uid2lidc.get(series_name) # scans of interest
                    if scan_name is None:
                        continue
                    else:
                        print(scan_name)
                    
                    vol,sop, loc, spacing, slice_thickness, series_uid = self.dicom_parser.parse_dicom(series_dir, order=order)
                    #parse xml file to get roi info
                    xml_files = glob("%s/*.xml" % series_dir)
                    #assert (len(xml_files)==1),"One dicom directory must only have one xml file %s" % series_dir
                    xml_file = xml_files[0]
                    nodule_slices = self.ann_parser.parse(xml_file)
                    vols = self.build_vol(nodule_slices,sop)
                    #assert(len(nodule_slices)==len(vols)), "num of nodules and num of vols dont match %d vs %d" %(len(nodule_slices),len(vols))
                    label_dict[scan_name] = vols
                                
                    lidc2img[scan_name] = sop
                    error_info = self.dicom_parser.check_location(loc, spacing[-1])
                    if error_info is not None:
                        logging.warning(error_info)
    
                    vol_data = vol_dataset.create_dataset(scan_name, data=vol, dtype=np.int16,
                                                          chunks=(512, 512, 1), compression='gzip')
                    vol_data.attrs['spacing'] = spacing
                    if slice_thickness is not None:
                        vol_data.attrs['slice_thickness'] = slice_thickness
    
        vol_dataset.close()
        with open(r'C:\Users\csce\Desktop\nova\lung_tumer\data\LIDC-IDRI\label_dict.pkl','a') as f:
              cPickle.dump(label_dict,f)
    def build_label(self,lidc_dir):
        label_dict = dict()
        uid2lidc = np.load(r'C:\Users\csce\Desktop\nova\lung_tumer\data\LIDC-IDRI\uid2lidc.pkl')
        uid2img = np.load(r'C:\Users\csce\Desktop\nova\lung_tumer\data\LIDC-IDRI\lidc2img.pkl')
        patient_dirs = os.listdir(lidc_dir)
        patient_dirs.sort()
        for patient_dir in patient_dirs:
            for study_dir in glob('%s/%s/*' % (lidc_dir, patient_dir)):
                for series_name in os.listdir(study_dir):
                    series_dir = os.path.join(study_dir, series_name)
                    if not self.dicom_parser.check_dcm_dir(series_dir):
                        continue
    
                    scan_name = uid2lidc.get(series_name) # scans of interest
                    if scan_name is None:
                        continue
                    else:
                        print(scan_name)
                    sop = uid2img.get(scan_name)
                    #parse xml file to get roi info
                    xml_files = glob("%s/*.xml" % series_dir)
                    #assert (len(xml_files)==1),"One dicom directory must only have one xml file %s" % series_dir
                    xml_file = xml_files[0]
                    nodule_slices = self.ann_parser.parse(xml_file)
                    vols = self.build_vol(nodule_slices,sop)
                    #assert(len(nodule_slices)==len(vols)), "num of nodules and num of vols dont match %d vs %d" %(len(nodule_slices),len(vols))
                    label_dict[scan_name] = vols
        with open(r'C:\Users\csce\Desktop\nova\lung_tumer\data\LIDC-IDRI\label_dict.pkl','a') as f:
              cPickle.dump(label_dict,f)
                    
        
    def build_vol(self,nodule_slices,sop):
        vols = list()
        for slices in nodule_slices:
            tmp = list()
            for idx,sop_uid in enumerate(sop):
                rect = slices.get(sop_uid)
                if (rect  is not None):
                    #print('slice match')
                
                    rect = np.append(rect,idx)
                    tmp.append(rect)
            tmp = np.array(tmp,dtype=np.int64)
            if tmp.shape[0] ==0:
                logging.warning("No nodule deteced")
                vol = np.array([0,0,0,0,0],dtype=np.int64)
            else:
                x_min = min(np.min(tmp[:,0]),np.min(tmp[:,2]))
                y_min = min(np.min(tmp[:,1]),np.min(tmp[:,3]))
                x_max = max(np.max(tmp[:,0]),np.max(tmp[:,2]))
                y_max = max(np.max(tmp[:,1]),np.max(tmp[:,3]))
                z_min = np.min(tmp[:,4])
                z_max = np.max(tmp[:,4])
                vol = np.array([x_min,y_min,z_min,x_max,y_max,z_max])
                vols.append(vol)
        return vols
                    
        
                
                
    def uid2dir(self,dataset_dir,dataset_name,output_dir):
          print("making user id to dataset map")
          uid2ds = dict()
          for patient_dir in os.listdir(dataset_dir):
              j = 0
              for study_dir in glob('%s/%s/*' % (dataset_dir, patient_dir)):
                  i = 0
                  j += 1
                  for series_name in os.listdir(study_dir):
                      uid2ds[series_name] = patient_dir + '_' +  str(j) + '_' +str(i)
                      i += 1
          output_name = "uid2{}.pkl".format(dataset_name)
          with open(r'C:\Users\csce\Desktop\nova\lung_tumer\data\LIDC-IDRI\uid2lidc.pkl','a') as f:
              cPickle.dump(uid2ds,f)
            
        

if __name__ =="__main__":
    dataset_name = 'LIDC'
    dataset_path = "/home/genonova/data/LIDC-IDRI/DOI"
    output_dir = "/home/genonova/dataser/LIDC-IDRI"
    Builder = DatasetBuilder(dataset_name,dataset_path,output_dir)
    #Builder.build_lidc_dataset(dataset_path,output_dir)
    #Builder.uid2dir(dataset_path,"lidc",output_dir)
    Builder.build_label(dataset_path)
    