"""
Created on 10/08/2017

@author: Niklas Pallast
Neuroimaging & Neuroengineering
Department of Neurology
University Hospital Cologne

"""

from __future__ import print_function

import argparse
import os
import sys
import numpy as np
import nibabel as nib
from datetime import datetime

def startSeedPoint(in_labels,in_atlas):


    ext_text = '.txt'
    ext_nifti = '.nii.gz'
    sPathLabels = None
    PathAtlas = None
    preserve = 0
    datatype = None
    # input labels text file with atlas index and seed regions (labels) in each line
    # Atlas (1 or 2), Label 1, Label 2, ...
    sPathLabels = in_labels
    if not os.path.isfile(sPathLabels):
        sys.exit("Error: '%s' is not a regular file." % (sPathLabels,))

    # input atlas labels files (NIfTI)
    PathAtlas = list([in_atlas])
    #sPathAtlas = in_atlas
    for sPathAtlas in PathAtlas:
        if not os.path.isfile(sPathAtlas):
            sys.exit("Error: '%s' is not a regular file." % (sPathAtlas,))

    # output seed ROIs file
    sPathROIs = os.path.join(os.path.dirname(PathAtlas[0]), 'Seed_ROIs.nii.gz')

    # get date and time
    # print(get_date())

    # read labels text file
    iatlas = []
    labels = []

    for row in read_csv(sPathLabels):
        iatlas.append(int(row.split(',\t')[0]))
        labels.append(int(row.split(',\t')[1]))
    # print("iatlas:", iatlas)
    # print("labels:", labels)

    # read 3D atlas labels files (NIfTI)
    labels_img = []
    labels_hdr = []
    labels_data = []
    labels_shape = []
    # read 3D atlas labels files (NIfTI)
    for k, sPathAtlas in enumerate(PathAtlas):
        # print("Atlas%d:" % (k + 1,), sPathAtlas)
        labels_img.append(nib.load(sPathAtlas))
        labels_data.append(labels_img[k].get_data())
        # print("labels_data[%d].dtype:" % (k,), labels_data[k].dtype)
        # print("labels_data[%d].shape:" % (k,), labels_data[k].shape)
        labels_hdr.append(labels_img[k].get_header())
        labels_shape.append(labels_hdr[k].get_data_shape())
        # print("labels_shape[%d]:" % (k,), labels_shape[k])
        if len(labels_shape[k]) != 3:
            sys.exit("Error: Atlas%d labels %s don't have three dimensions." % (k, str(labels_shape[k])))

    for k in range(1, len(labels_shape)):
        if labels_shape[0] != labels_shape[k]:
            sys.exit("Error: Atlas1 labels %s and Atlas%d labels %s don't have the same shape." % (
            str(labels_shape[0]), k, str(labels_shape[k])))

    # create atlas labels hyperstack (4D)
    rois = create_rois_1(iatlas, labels, labels_hdr, labels_data, datatype=datatype, preserve=preserve)

    # save atlas labels file
    dataOrg = nib.load(sPathAtlas)

    niiData = nib.Nifti1Image(rois, dataOrg.affine)
    hdrIn = niiData.header
    hdrIn.set_xyzt_units('mm')
    scaledNiiData = nib.as_closest_canonical(niiData)
    nib.save(niiData, sPathROIs)

    print("output:", sPathROIs)
    return sPathROIs


def get_date():
    now = datetime.now()
    pvDate = now.strftime("%a %d %b %Y")
    pvTime = now.strftime("%H:%M:%S")
    return pvDate + ' ' + pvTime

def read_csv(sFilename):
    fid = open(sFilename)
    holeDataset = fid.read()
    rowsInDataset = holeDataset.split('\n')
    return rowsInDataset[1::]

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
                rois[:,:,:,k][data==label] = label
    else:
        for k, index in enumerate(iatlas):
            data = labels_data[index-1]
            rois[:,:,:,k][data==labels[k]] = 1
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
            rois[:,:,:,k][indices] = data[indices]
    else:
        for k, index in enumerate(iatlas):
            ires = []
            data = labels_data[index-1]
            for label in labels[k]:
                ires.append(np.where(data == label))
            indices = tuple(np.hstack(ires))
            rois[:,:,:,k][indices] = 1


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
            rois[:,:,:,k] = data * mask
            mask[:] = False
    else:
        for k, index in enumerate(iatlas):
            data = labels_data[index-1]
            for label in labels[k]:
                mask = np.logical_or(mask, data == label)
            rois[:,:,:,k] = mask
            mask[:] = False    
    return rois

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create atlas seed ROIs.')
    parser.add_argument('-i', '--in_atlas', nargs='+', help='input 3D atlas labels file names (NIfTI)')
    parser.add_argument('-l','--in_labels', help='input labels text file name',default='/Volumes/AG_Aswendt_Share/Scratch/Asw_fMRI2AllenBrain_Data/annotation_50CHANGEDanno_label_IDs.txt')
    parser.add_argument('-o', '--out_rois', help='output 4D seed ROIs file name')
    parser.add_argument('-p', '--preserve', action='store_true', help='preserve label values')
    parser.add_argument('-t', '--datatype', type=int, choices=[2, 4, 8, 16], help='data type (2: char, 4: short, 8: int, 16: float)')
    args = parser.parse_args()

    ext_text = '.txt'
    ext_nifti = '.nii.gz'

    # input labels text file with atlas index and seed regions (labels) in each line
    # Atlas (1 or 2), Label 1, Label 2, ...
    if len(args.in_labels) > 0: sPathLabels = args.in_labels
    if not os.path.isfile(sPathLabels):
        sys.exit("Error: '%s' is not a regular file." % (sPathLabels,))

    # input atlas labels files (NIfTI)
    if len(args.in_atlas) > 0: PathAtlas = args.in_atlas
    for sPathAtlas in PathAtlas:
        if not os.path.isfile(sPathAtlas):
            sys.exit("Error: '%s' is not a regular file." % (sPathAtlas,))

    # output seed ROIs file
    sPathROIs = os.path.join(os.getcwd(), args.out_rois) if args.out_rois != None else os.path.join(os.path.dirname(PathAtlas[0]), 'Seed_ROIs')

    # get date and time
    #print(get_date())

    # read labels text file
    iatlas = []
    labels = []

    for row in read_csv(sPathLabels):
        iatlas.append(int(row.split(',\t')[0]))
        labels.append(int(row.split(',\t')[1]))
    #print("iatlas:", iatlas)
    #print("labels:", labels)

    # read 3D atlas labels files (NIfTI)
    labels_img = []
    labels_hdr = []
    labels_data = []
    labels_shape = []
    # read 3D atlas labels files (NIfTI)
    for k, sPathAtlas in enumerate(PathAtlas):
        #print("Atlas%d:" % (k + 1,), sPathAtlas)
        labels_img.append(nib.load(sPathAtlas))
        labels_data.append(labels_img[k].get_data())
        #print("labels_data[%d].dtype:" % (k,), labels_data[k].dtype)
        #print("labels_data[%d].shape:" % (k,), labels_data[k].shape)
        labels_hdr.append(labels_img[k].get_header())
        labels_shape.append(labels_hdr[k].get_data_shape())
        #print("labels_shape[%d]:" % (k,), labels_shape[k])
        if len(labels_shape[k]) != 3:
            sys.exit("Error: Atlas%d labels %s don't have three dimensions." % (k, str(labels_shape[k])))

    for k in range(1, len(labels_shape)):
        if labels_shape[0] != labels_shape[k]:
            sys.exit("Error: Atlas1 labels %s and Atlas%d labels %s don't have the same shape." % (str(labels_shape[0]), k, str(labels_shape[k])))

    # create atlas labels hyperstack (4D)
    rois = create_rois_1(iatlas, labels, labels_hdr, labels_data, datatype=args.datatype, preserve=args.preserve)

    # save atlas labels file
    dataOrg = nib.load(sPathAtlas)


    niiData = nib.Nifti1Image(rois, dataOrg.affine)
    hdrIn = niiData.header
    hdrIn.set_xyzt_units('mm')
    scaledNiiData = nib.as_closest_canonical(niiData)
    nib.save(niiData, sPathROIs+ext_nifti)

    print("output:", sPathROIs+ext_nifti)
