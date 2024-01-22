'''
Created on 20.08.2020

Author:
Michael Diedenhofen
Max Planck Institute for Metabolism Research, Cologne
'''

from __future__ import print_function

try:
    zrange = xrange
except NameError:
    zrange = range

import os
import sys

import numpy as np
import nibabel as nib

import proc_tools as pt

def create_rois_1(iatlas, labels, labels_hdr, labels_data, datatype=None, preserve=False):
    if datatype == 2:
        labels_dtype = np.uint8
    elif datatype == 4:
        labels_dtype = np.int16
    elif datatype == 8:
        labels_dtype = np.int32
    elif datatype == 16:
        labels_dtype = np.float32
    else:
        labels_dtype = labels_hdr[0].get_data_dtype()
    labels_shape = labels_hdr[0].get_data_shape()
    rois = np.zeros(labels_shape + (len(iatlas),), dtype=labels_dtype)
    if preserve:
        for k, index in enumerate(iatlas):
            data = labels_data[index-1]
            for label in labels[k]:
                rois[:, :, :, k][data==label] = label
    else:
        for k, index in enumerate(iatlas):
            data = labels_data[index-1]
            for label in labels[k]:
                rois[:, :, :, k][data==label] = 1

    return rois

def create_rois_2(iatlas, labels, labels_hdr, labels_data, datatype=None, preserve=False):
    if datatype == 2:
        labels_dtype = np.uint8
    elif datatype == 4:
        labels_dtype = np.int16
    elif datatype == 8:
        labels_dtype = np.int32
    elif datatype == 16:
        labels_dtype = np.float32
    else:
        labels_dtype = labels_hdr[0].get_data_dtype()
    labels_shape = labels_hdr[0].get_data_shape()
    rois = np.zeros(labels_shape + (len(iatlas),), dtype=labels_dtype)
    if preserve:
        for k, index in enumerate(iatlas):
            ires = []
            data = labels_data[index-1]
            for label in labels[k]:
                ires.append(np.where(data == label))
            indices = tuple(np.hstack(ires))
            rois[:, :, :, k][indices] = data[indices]
    else:
        for k, index in enumerate(iatlas):
            ires = []
            data = labels_data[index-1]
            for label in labels[k]:
                ires.append(np.where(data == label))
            indices = tuple(np.hstack(ires))
            rois[:, :, :, k][indices] = 1

    return rois

def create_rois_3(iatlas, labels, labels_hdr, labels_data, datatype=None, preserve=False):
    if datatype == 2:
        labels_dtype = np.uint8
    elif datatype == 4:
        labels_dtype = np.int16
    elif datatype == 8:
        labels_dtype = np.int32
    elif datatype == 16:
        labels_dtype = np.float32
    else:
        labels_dtype = labels_hdr[0].get_data_dtype()
    labels_shape = labels_hdr[0].get_data_shape()
    mask = np.zeros(labels_shape, dtype=np.bool)
    rois = np.zeros(labels_shape + (len(iatlas),), dtype=labels_dtype)
    if preserve:
        for k, index in enumerate(iatlas):
            data = labels_data[index-1]
            for label in labels[k]:
                mask = np.logical_or(mask, data == label)
            rois[:, :, :, k] = data * mask
            mask[:] = False
    else:
        for k, index in enumerate(iatlas):
            data = labels_data[index-1]
            for label in labels[k]:
                mask = np.logical_or(mask, data == label)
            rois[:, :, :, k] = mask
            mask[:] = False    

    return rois

def create_rois(path_labels, list_atlas, datatype=None, preserve=False):
    # read labels text file
    iatlas, labels = pt.read_labels(path_labels)

    # read 3D atlas labels files (NIfTI)
    labels_img = []
    labels_hdr = []
    labels_data = []
    labels_shape = []
    for k, path_atlas in enumerate(list_atlas):
        #print("Atlas%d:" % (k + 1,), path_atlas)
        labels_img.append(nib.load(path_atlas))
        labels_data.append(labels_img[k].get_data())
        #print("labels_data[%d].dtype:" % (k,), labels_data[k].dtype)
        #print("labels_data[%d].shape:" % (k,), labels_data[k].shape)
        labels_hdr.append(labels_img[k].get_header())
        labels_shape.append(labels_hdr[k].get_data_shape())
        #print("labels_shape[%d]:" % (k,), labels_shape[k])
        if len(labels_shape[k]) != 3:
            sys.exit("Error: Atlas%d labels %s don't have three dimensions." % (k, str(labels_shape[k])))

    for k in zrange(1, len(labels_shape)):
        if labels_shape[0] != labels_shape[k]:
            sys.exit("Error: Atlas1 labels %s and Atlas%d labels %s don't have the same shape." % (str(labels_shape[0]), k, str(labels_shape[k])))

    # create atlas labels hyperstack (4D)
    rois = create_rois_1(iatlas, labels, labels_hdr, labels_data, datatype=datatype, preserve=preserve)

    return (labels_hdr, rois)

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Create atlas seed ROIs.')
    parser.add_argument('in_labels', help='input labels text file name')
    parser.add_argument('in_atlas', nargs='+', help='input 3D atlas labels file names (NIfTI)')
    parser.add_argument('-o', '--out_rois', help='output 4D seed ROIs file name')
    parser.add_argument('-p', '--preserve', action='store_true', help='preserve label values')
    parser.add_argument('-t', '--datatype', type=int, choices=[2, 4, 8, 16], help='data type (2: char, 4: short, 8: int, 16: float)')
    args = parser.parse_args()

    ext_text = '.txt'
    ext_nifti = '.nii.gz'
    file_name = 'Seed_ROIs'

    # input labels text file with atlas index and seed regions (labels) in each line
    # Atlas (1 or 2), Label 1, Label 2, ...
    if len(args.in_labels) > 0: path_labels = args.in_labels
    if not os.path.isfile(path_labels):
        sys.exit("Error: '%s' is not a regular file." % (path_labels,))

    # input atlas labels files (NIfTI)
    if len(args.in_atlas) > 0: list_atlas = args.in_atlas
    for path_atlas in list_atlas:
        if not os.path.isfile(path_atlas):
            sys.exit("Error: '%s' is not a regular file." % (path_atlas,))

    # output seed ROIs file
    path_rois = os.path.join(os.getcwd(), args.out_rois) if args.out_rois != None else os.path.join(os.path.dirname(list_atlas[0]), file_name)

    # get date and time
    #print(pt.get_date())

    # create atlas labels hyperstack (4D)
    labels_hdr, rois = create_rois(path_labels, path_atlas, datatype=args.datatype, preserve=args.preserve)

    # save atlas labels file
    voxel_dims = labels_hdr[0].get_zooms()
    pt.save_data(rois, voxel_dims, path_rois, dtype=None)
