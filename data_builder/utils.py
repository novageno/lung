# -*- coding: utf-8 -*-
"""
Created on Wed Aug 09 10:35:49 2017
utils file for dataset creation
@author: csce
"""
import os
from glob import glob
import dicom
import matplotlib.pyplot as plt
def find_all_files(root, suffix=None):
    res = []
    for root, _, files in os.walk(root):
        for f in files:
            if suffix is not None and not f.endswith(suffix):
                continue
            res.append(os.path.join(root, f))
    return res

def vis_dcm_files(dicom_dir):
    dicom_names = glob("%s/*.dcm" % dicom_dir)
    for dicom_name in dicom_names:
        ax = plt.gca()
        plt.figure()
        dicom_ = dicom.read_file(dicom_name)
        img = dicom_.pixel_array
        ax.imshow(img)
        plt.show()
    
if __name__ =="__main__":
    dicom_dir = r"E:\LIDC-IDRI\DOI\LIDC-IDRI-0001\1.3.6.1.4.1.14519.5.2.1.6279.6001.298806137288633453246975630178\1.3.6.1.4.1.14519.5.2.1.6279.6001.179049373636438705059720603192"
    vis_dcm_files(dicom_dir)
        
        
        
                
    