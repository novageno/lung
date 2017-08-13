# -*- coding: utf-8 -*-
"""
Created Faster RCNN data
"""

import numpy as np
import h5py,os
import array
import scipy.io as sio
import cPickle

#### create all nodule as aim




def cr_data(main_dir):

    if not os.path.exists(main_dir+'/data/faster_rcnn_data/v9/JPEGImages'):
        os.makedirs(main_dir+'/data/faster_rcnn_data/v9/JPEGImages')
    if not os.path.exists(main_dir+'/data/faster_rcnn_data/v9/Annotations'):
        os.makedirs(main_dir+'/data/faster_rcnn_data/v9/Annotations')
    if not os.path.exists(main_dir+'/data/faster_rcnn_data/v9/ImageSets'):
        os.makedirs(main_dir+'/data/faster_rcnn_data/v9/ImageSets')
    if not os.path.exists(main_dir+'/data/faster_rcnn_data/v9/ImageSets/fold1'):
        os.makedirs(main_dir+'/data/faster_rcnn_data/v9/ImageSets/fold1')
    if not os.path.exists(main_dir+'/data/faster_rcnn_data/v9/ImageSets/fold2'):
        os.makedirs(main_dir+'/data/faster_rcnn_data/v9/ImageSets/fold2')
    if not os.path.exists(main_dir+'/data/faster_rcnn_data/v9/ImageSets/fold3'):
        os.makedirs(main_dir+'/data/faster_rcnn_data/v9/ImageSets/fold3')        
    if not os.path.exists(main_dir+'/data/faster_rcnn_data/v9/ImageSets/fold4'):
        os.makedirs(main_dir+'/data/faster_rcnn_data/v9/ImageSets/fold4')
    if not os.path.exists(main_dir+'/data/faster_rcnn_data/v9/ImageSets/fold5'):
        os.makedirs(main_dir+'/data/faster_rcnn_data/v9/ImageSets/fold5')

    #f1=h5py.File(main_dir+'/data/lidc_v9/vol.hdf5','r')
    #f2=h5py.File(main_dir+'/data/kaggle/vol.hdf5','r')
    f1 = h5py.File('/home/genonova/dataset/LIDC-IDRI/vol.hdf5')
    #l=np.load(main_dir+'/data/LIDC-IDRI/label_dict.npy').item()
    with open('/home/genonova/nova/lung_tumer/data/LIDC-IDRI/label_dict.pkl','rb') as f:
        l = cPickle.load(f)
    #fold_list = np.load(main_dir+'/code/cr_frcnn_data/fold_list.npy')
    print(l)
    names = f1.keys()

    for name in names:
        lab=l[name]
        vol = f1[name][:]
            

        vol = vol + 1350.0
        vol [vol < 0] = 0
        vol [vol > 1500] = 1500
        vol = vol * 255.0 /1500
        vol = np.floor(vol)
        vol [vol < 0] = 0   

        labf = [[] for _ in xrange(len(lab))]
        labl = [[] for _ in xrange(len(lab))]

        for ilab in xrange (len(lab)):
            h=lab[ilab][5]-lab[ilab][2]
            labf[ilab]=np.round(lab[ilab][2]+0.2*h)
            labl[ilab]=np.round(lab[ilab][2]+0.8*h)

        for islice in xrange(1,vol.shape[2]-1):
            cr_img = 0
            cr_xml = 0
            for j in xrange (len(lab)):
                if (islice >=labf[j]) and (islice<=labl[j]):       
                    if cr_img == 0:
                        img=np.dstack((vol[:,:,islice-1],vol[:,:,islice],vol[:,:,islice+1]))
                        sio.savemat(main_dir+'/data/faster_rcnn_data/v9/JPEGImages/'+name+'_'+str(islice)+'.mat',{'img':img},do_compression=True)
                        # namelist.append(name+'_'+str(islice))
                        if cmp(name,'LIDC-IDRI-0671_01')==0 and islice==50:
                            sio.savemat(main_dir+'/data/faster_rcnn_data/v9/JPEGImages/'+name+'_'+str(islice)+'.mat',{'img':img})
                        
                        cr_img=1
                    if cr_xml==0:
                        fid=open(main_dir+'/data/faster_rcnn_data/v9/Annotations/'+name+'_'+str(islice)+'.xml','w')
                        print >> fid,'<annotation>\n'
                        print >> fid,'\t<size>\n'
                        print >>fid,'\t\t<width>512</width>\n'
                        print >> fid,'\t\t<height>512</height>\n'
                        print >> fid,'\t\t<depth>3</depth>\n'
                        print >>fid,'\t</size>\n'
                        print >>fid,'\t<object>\n'
                        print >>fid,'\t\t<name>nodule</name>\n'
                        print >>fid,'\t\t<difficult>0</difficult>\n'
                        print >>fid,'\t\t<bndbox>\n'
                        print >>fid,'\t\t\t<xmin>', int(np.round(lab[j][0])), '</xmin>\n'
                        print >>fid,'\t\t\t<ymin>', int(np.round(lab[j][1])), '</ymin>\n'
                        print >>fid,'\t\t\t<xmax>', int(np.round(lab[j][3])), '</xmax>\n'
                        print >>fid,'\t\t\t<ymax>', int(np.round(lab[j][4])), '</ymax>\n'
                        print>>fid,'\t\t</bndbox>\n'
                        print>>fid,'\t</object>\n'
                        cr_xml=1
                    else:
                        print >>fid,'\t<object>\n'
                        print >>fid,'\t\t<name>nodule</name>\n'
                        print >>fid,'\t\t<difficult>0</difficult>\n'
                        print >>fid,'\t\t<bndbox>\n'
                        print >>fid,'\t\t\t<xmin>', int(np.round(lab[j][0])), '</xmin>\n'
                        print >>fid,'\t\t\t<ymin>', int(np.round(lab[j][1])), '</ymin>\n'
                        print >>fid,'\t\t\t<xmax>', int(np.round(lab[j][3])), '</xmax>\n'
                        print >>fid,'\t\t\t<ymax>', int(np.round(lab[j][4])), '</ymax>\n'
                        print >>fid,'\t\t</bndbox>\n'
                        print >>fid,'\t</object>\n'
            if cr_xml:
                print >>fid,'</annotation>\n'
                fid.close()

def cr_data_fold(main_dir,k):
    """split dataset into k flod for k-fold cross evalution"""
    num_files = len(os.listdir(main_dir+'data/faster_rcn
    for i in xrange(k):
        dir_name = "fold{}".format(k)
        if not os.path.exists(main_dir+'/data/faster_rcnn_data/v9/ImageSets/'+dir_name):
            os.makedirs(main_dir+'/data/faster_rcnn_data/v9/ImageSets/' + dir_name)
    



if __name__ == '__main__':
    main_dir = "/home/genonova/nova/lung_tumer"
    if not os.path.exists(main_dir+'/data/faster_rcnn_data/v9/ImageSets'):
        os.makedirs(main_dir+'/data/faster_rcnn_data/v9/ImageSets')
    
    cr_data(main_dir)
