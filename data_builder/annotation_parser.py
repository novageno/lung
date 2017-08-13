# -*- coding: utf-8 -*-

import os
import xml.etree.ElementTree as etree
import numpy as np
import logging
import cPickle
import glob
from nodule import RadAnnotation, SmallNodule, NormalNodule, \
    NoduleRoi, NonNodule, AnnotationHeader
from utils import find_all_files

class LIDC_XML_Parser():
    def __init__(self):
        self.NS = {'nih': 'http://www.nih.gov'}

    def parse_dir(self,dirname,flatten=True,pickle=True):
    
        if not flatten:
            return self.parse_original_xmls(dirname, pickle)
        pickle_file = os.path.join(dirname,"annotaion_flatten.pkl")
            
        if os.path.isfile(pickle_file):
            logging.info("Loading annotations from file %s" % pickle_file)
            with open(pickle_file, 'r') as f:
                annotations = cPickle.load(f)
                logging.info("Load annotations complete")
                return annotations
        annotations = self.parse_original_xmls(dirname, pickle)
        annotations = self.flatten_annotation(annotations)
        if pickle:
            logging.info("Saving annotations to file %s" % pickle_file)
            with open(pickle_file, 'w') as f:
                cPickle.dump(annotations, f)
        return annotations

    def parse_original_xmls(self,dirname,pickle=True):
        pickle_file = pickle and os.path.join(dirname, 'annotation.pkl') or None
        if pickle and os.path.isfile(pickle_file):
            logging.info("Loading annotations from file %s" % pickle_file)
            with open(pickle_file, 'r') as f:
                annotations = cPickle.load(f)
            logging.info("Load annotations complete")
        else:
            logging.info("Reading annotations")
            annotations = []
            xml_files = find_all_files(dir_name,'.xml')
            print(len(xml_files))
            for f in xml_files:
                annotations.append(self.parse(f))
            if pickle and not os.path.isfile(pickle_file):
                logging.info("Saving annotations to file %s" % pickle_file)
                with open(pickle_file, 'w') as f:
                    cPickle.dump(annotations, f)
            return annotations
 

    def parse(self,xml_filename):
        logging.info("Parsing %s" % xml_filename)
        annotations = []
        # ET is the library we use to parse xml data
        tree = etree.parse(xml_filename)
        root = tree.getroot()
        # header = parse_header(root)
        # readingSession-> holds radiologist's annotation info
        for read_session in root.findall('nih:readingSession',self.NS):
            # to hold each radiologists annotation
            # i.e. readingSession in xml file
            rad_annotation = RadAnnotation()
            rad_annotation.version = \
                read_session.find('nih:annotationVersion',self.NS).text
            rad_annotation.id = \
                read_session.find('nih:servicingRadiologistID', self.NS).text

            # nodules
            nodule_nodes = read_session.findall('nih:unblindedReadNodule', self.NS)
            for node in nodule_nodes:
                nodule = self.parse_nodule(node)
                if nodule.is_small:
                    rad_annotation.small_nodules.append(nodule)
                else:
                    rad_annotation.nodules.append(nodule)
                    
                    # non-nodules
            non_nodule = read_session.findall('nih:nonNodule', self.NS)
            for node in non_nodule:
                nodule = self.parse_non_nodule(node)
                rad_annotation.non_nodules.append(nodule)
            annotations.append(rad_annotation)
            vols = self._filter(annotations)
        return vols
    def _filter(self,anns):
        """filter useless info"""
        ann = anns[-1]
        nodules = ann.nodules
        
        nodule_slices = list()
        for nodule in nodules:
            vol = dict()
            rois = nodule.rois
            for roi in rois:
               vol[roi.sop_uid] = np.array(roi.roi_rect,dtype=np.int64)
            nodule_slices.append(vol)
        assert(len(nodules)==len(nodule_slices)), "dim not match"
        return nodule_slices
        
    def parse_header(self,root):
        header = AnnotationHeader()
        print root.findall('nih:*', self.NS)
        resp_hdr = root.findall('nih:ResponseHeader', self.NS)[0]
        header.version = resp_hdr.find('nih:Version', self.NS).text
        header.message_id = resp_hdr.find('nih:MessageId', self.NS).text
        header.date_request = resp_hdr.find('nih:DateRequest', self.NS).text
        header.time_request = resp_hdr.find('nih:TimeRequest', self.NS).text
        header.task_desc = resp_hdr.find('nih:TaskDescription', self.NS).text
        header.series_instance_uid = resp_hdr.find('nih:SeriesInstanceUid', self.NS).text
        date_service = resp_hdr.find('nih:DateService', self.NS)
        if date_service is not None:
            header.date_service = date_service.text
            time_service = resp_hdr.find('nih:TimeService', self.NS)
        if time_service is not None:
            header.time_service = time_service.text
        header.study_instance_uid = resp_hdr.find('nih:StudyInstanceUID',self.NS).text
        return header

    def parse_nodule(self,xml_node):
        char_node = xml_node.find('nih:characteristics', self.NS)
        # if no characteristics, it is smallnodule  i.e. is_small=TRUE
        is_small = (char_node is None or len(char_node) == 0)
        nodule = is_small and SmallNodule() or NormalNodule()
        nodule.id = xml_node.find('nih:noduleID', self.NS).text
        if not is_small:
            subtlety = char_node.find('nih:subtlety', self.NS)
            nodule.characteristics.subtlety = int(subtlety.text)
            nodule.characteristics.internal_struct = \
                int(char_node.find('nih:internalStructure', self.NS).text)
            nodule.characteristics.calcification = \
                int(char_node.find('nih:calcification',self.NS).text)
            nodule.characteristics.sphericity = \
                int(char_node.find('nih:sphericity', self.NS).text)
            nodule.characteristics.margin = \
                int(char_node.find('nih:margin', self.NS).text)
            nodule.characteristics.lobulation = \
                int(char_node.find('nih:lobulation', self.NS).text)
            nodule.characteristics.spiculation = \
                int(char_node.find('nih:spiculation', self.NS).text)
            nodule.characteristics.texture = \
                int(char_node.find('nih:texture', self.NS).text)
            nodule.characteristics.malignancy = \
                int(char_node.find('nih:malignancy', self.NS).text)
        xml_rois = xml_node.findall('nih:roi', self.NS)
        for xml_roi in xml_rois:
            roi = NoduleRoi()
            roi.z = float(xml_roi.find('nih:imageZposition', self.NS).text)
            roi.sop_uid = xml_roi.find('nih:imageSOP_UID', self.NS).text
            # when inclusion = TRUE ->roi includes the whole nodule
            # when inclusion = FALSE ->roi is drown twice for one nodule
            # 1.ouside the nodule
            # 2.inside the nodule -> to indicate that the nodule has donut
            # hole(the inside hole is
            # not part of the nodule) but by forcing inclusion to be TRUE,
            # this situation is ignored
            roi.inclusion = (xml_roi.find('nih:inclusion', self.NS).text == "TRUE")
            edge_maps = xml_roi.findall('nih:edgeMap', self.NS)
            for edge_map in edge_maps:
                x = int(edge_map.find('nih:xCoord', self.NS).text)
                y = int(edge_map.find('nih:yCoord', self.NS).text)
                roi.roi_xy.append([x, y])
            xmax = np.array(roi.roi_xy)[:, 0].max()
            xmin = np.array(roi.roi_xy)[:, 0].min()
            ymax = np.array(roi.roi_xy)[:, 1].max()
            ymin = np.array(roi.roi_xy)[:, 1].min()
            if not is_small:  # only for normalNodules
                roi.roi_rect = (xmin, ymin, xmax, ymax)
                roi.roi_centroid = (
                    (xmax + xmin) / 2., (ymin + ymax) / 2.)  # center point
            nodule.rois.append(roi)
        return nodule  # is equivalent to unblindedRead
    def parse_non_nodule(self,xml_node):  # xml_node is one nonNodule
        nodule = NonNodule()
        nodule.id = xml_node.find('nih:nonNoduleID', self.NS).text
        roi = NoduleRoi()
        roi.z = float(xml_node.find('nih:imageZposition', self.NS).text)
        roi.sop_uid = xml_node.find('nih:imageSOP_UID', self.NS).text
        loci = xml_node.findall('nih:locus', self.NS)
        for locus in loci:
            x = int(locus.find('nih:xCoord', self.NS).text)
            y = int(locus.find('nih:yCoord', self.NS).text)
            roi.roi_xy.append((x, y))
        nodule.rois.append(roi)
        return nodule  # is equivalent to nonNodule(xml element)
    def flatten_annotation(self,annotation_dict):
        logging.info("Start flatten")
        res = {}
        for annotations in annotation_dict:
            # annotations in each file
            for anno in annotations:
                self.flatten_nodule(anno.nodules, 'nodules', res)
                self.flatten_nodule(anno.small_nodules, 'small_nodules', res)
                self.flatten_nodule(anno.non_nodules, 'non_nodules', res)
        logging.info("Flatten complete")
        return res
    def flatten_nodule(self,nodules, type, result):
        for nodule in nodules:
            for roi in nodule.rois:
                # logging.info(roi)
                sop_uid = roi.sop_uid
                # logging.info(sop_uid)
                # logging.info(result)
                if not result.has_key(sop_uid):
                    result[sop_uid] = {
                        'nodules': [], 'small_nodules': [], 'non_nodules': []
                    }
                centroid = type == 'nodules' and roi.roi_centroid or roi.roi_xy[0]
                point = {'centroid': centroid, 'pixels': roi.roi_xy, 'field': roi.roi_rect}
                result[sop_uid][type].append(point)            
if __name__ == "__main__":
    dir_name = r"E:\LIDC-IDRI\DOI\LIDC-IDRI-0016"
    LIDC_parser = LIDC_XML_Parser()
    anns = LIDC_parser.parse_dir(dir_name)