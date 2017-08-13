# -*- coding: utf-8 -*-
"""
Created on Mon Aug 07 21:06:43 2017

@author: csce
"""
import logging
import dicom
import numpy as np
from glob import glob



class DicomParse():
        
    def parse_dicom(self,dicom_dir,order="InstanceNumber"):
        """parse info of single patient
        param:
            dicom_dir:path to sigle patient
        """
        dicom_names = glob("%s/*.dcm" % dicom_dir)
        slices_list = list()
        
        for dicom_name in dicom_names:
            slices_list.append(self.parse_single_dicom(dicom_name))
        if order=='InstanceNumber':
            slices_list.sort(key=(lambda x: int(x[1].InstanceNumber)))
        elif order=="Location":
            slices_list.sort(key=(lambda x: self.get_location(x[1])))
            mean_loaction = np.mean([self.get_location(x_) for _,x_ in slices_list])
            reverse_flag = bool(mean_loaction<=0)
            if reverse_flag:
                slices_list = slices_list[::-1]
                logging.warning("non positive mean location; then reverse")
        else:
            raise ValueError("unknown order")
        print("load dicom file done.")    
        vol = list()
        sop = list()
        loc = list()
        series_uid = None
        spacing = None
        slice_thickness = None
        
        for slice,dicom_ in slices_list:
            vol.append(slice)
            sop.append(str(dicom_.SOPInstanceUID))
            loc.append(self.get_location(dicom_))
            # check series_uid
        if series_uid == str(dicom_.SeriesInstanceUID):
            pass
        elif series_uid is None:
            series_uid = str(dicom_.SeriesInstanceUID)
        else:
            logging.error('series_uid inconsistent: %s and %s' %(series_uid, str(dicom_.SeriesInstanceUID)))
        # check spacing
        if spacing is None:
            spacing = np.asarray(dicom_.PixelSpacing)
        elif np.max(spacing - np.asarray(dicom.PixelSpacing)) > 0.1:
            logging.error('spacing inconsistent')

        # try to get slice_thickness
        if slice_thickness is None:
            slice_thickness = dicom_.get('SliceThickness')
        elif slice_thickness != dicom_.get('SliceThickness'):
            slice_thickness = -1

        vol = np.transpose(np.asarray(vol), (1, 2, 0)).astype(np.int16)
        emp_diff = np.fabs(loc[2] - loc[3])
        sum_diff = []
        for ind in range(len(loc) - 1):
            diff = np.fabs(loc[ind] - loc[ind + 1])
            sum_diff.append(diff)

        mean_diff = np.mean(sum_diff)
        if slice_thickness is not None and np.fabs(slice_thickness - mean_diff) < 0.1:
            spacing = np.append(spacing, slice_thickness)
        else:
            spacing = np.append(spacing, emp_diff)

        return vol, sop, loc, spacing, slice_thickness, series_uid
    def get_location(self,dicom_):
        location = None
        
        try:
            location = np.float(dicom_.SliceLocation)
        except AttributeError:
            logging.warning('dicom miss location info.')
        try:
            location = np.float(dicom_.ImagePositionPatient[2])
        except AttributeError:
            logging.warning('dicom miss image position.')

        return location        

    def parse_single_dicom(self,dicom_name):
        """parse dicom file
        params:
            dicom_name:data name of single dicom file
        return:
        """
        dicom_ = dicom.read_file(dicom_name)
        slice = dicom_.pixel_array
        try:
            intercept = dicom_.RescaleIntercept
            slope = dicom.RescaleSlope
        except AttributeError:
            intercept = 0
            slope = 1
            logging.warning("rescale info missing")
            
        finally:
            slice = slice * slope + intercept
        return slice,dicom_
    def check_dcm_dir(self,dicom_dir):
        dicom_names = glob('%s/*.dcm' % dicom_dir)

        if not dicom_names:
            logging.debug('%s contains no dicoms.' % dicom_dir)
            return False

        dicom_ = dicom.read_file(dicom_names[0], stop_before_pixels=True)
        
        if dicom_.Modality != 'CT':
            logging.debug('%s is not a CT directory.' % dicom_dir)
            return False
        
        return True
    def check_location(self,location, slice_thickness):
        max_diff = 0.5 * slice_thickness
        min_diff = 1e-2
    
        for ind in range(len(location)-1):
            diff = np.fabs(location[ind] - location[ind+1])
            if diff < min_diff:
                return 'too small diff: (%f, %f). slice_thickness = %f' \
                       % (location[ind], location[ind + 1], slice_thickness)
            elif np.fabs(diff-slice_thickness) >= max_diff:
                return 'too large diff: (%f, %f). slice_thickness = %f' \
                       % (location[ind], location[ind + 1], slice_thickness)
    
        return None

                    


if __name__=="__main__":
    data_path = "E:\NBIA\DOI"
    parser = DicomParse(data_path)
    _,sop,loc,spacing,slice_thickness,series_uid=parser.parse_dicom(r"E:\NBIA\DOI\CT-Training-BE001\1.2.840.113704.1.111.2112.1167842143.1\1.2.840.113704.1.111.2112.1167842347.17")