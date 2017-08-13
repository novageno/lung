# -*- coding: utf-8 -*-
"""
Created on Sat Aug 12 08:21:59 2017

@author: csce
"""

import dicom
import h5py
import cPickle
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle
def check_img_vis(data,label):
    
    for key in data.keys():
        vol = data[key]
        labels = label[key]
        for lab in labels:
            begin = lab[2]
            end = lab[5]
            for i in range(begin,end):
                vis_single_img(vol[:,:,i],lab)

def vis_single_img(data,label):
    ax = plt.gca()
    plt.figure()
    ax.imshow(data)
    print(label[0],label[1],label[3] - label[0] + 1, label[4] - label[1]+1)
    gt_box_plot = Rectangle((label[0], label[1]), label[3]-label[0]+1, label[4]-label[1]+1, fill=False, edgecolor='red', linewidth=0.5)
    ax.add_patch(gt_box_plot)
    plt.show()
    
if __name__=="__main__":
    data = h5py.File("D:data/LIDC-IDRI/vol.hdf5")
    with open("C:/Users/csce/Desktop/nova/lung_tumer/data/LIDC-IDRI/label_dict.pkl") as f:
        label = cPickle.load(f)
    check_img_vis(data,label)