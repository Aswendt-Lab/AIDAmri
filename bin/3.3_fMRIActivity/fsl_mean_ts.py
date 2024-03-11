"""
Created on 07.12.2015

@author: michaeld
"""

from __future__ import print_function

import argparse

import os
import sys

import numpy as np
import nibabel as nib
import scipy.io as io
import correlate_matrix

from datetime import datetime
def start_fsl_mean_ts(sPathData,sPathMask,labelNames,postTxt):
    # input data

    # Read 4D data file (NIfTI )
    data_img = nib.load(sPathData)
    data = data_img.get_data()
    data_hdr = data_img.header
    data_dtype = data_hdr.get_data_dtype()
    data_shape = data_hdr.get_data_shape()

    # output data
    sPathOut = os.path.abspath(os.path.join(sPathData, os.pardir, postTxt + os.path.basename(sPathData).split('_')[0]))
    sPathOut = sPathOut + '.txt'
    
    if os.path.basename(labelNames) == "annoVolume.nii.txt":
        PcorrR_matrix_path = os.path.abspath(os.path.join(sPathData, os.pardir,  'Matrix_PcorrR.' + os.path.basename(sPathData).split('_')[0])) + ".mat"
        PcorrP_matrix_path = os.path.abspath(os.path.join(sPathData, os.pardir,  'Matrix_PcorrP.' + os.path.basename(sPathData).split('_')[0])) + ".mat"
        PcorrZ_matirx_path = os.path.abspath(os.path.join(sPathData, os.pardir,  'Matrix_PcorrZ.' + os.path.basename(sPathData).split('_')[0])) + ".mat"
    elif os.path.basename(labelNames) == "annoVolume+2000_rsfMRI.nii.txt" :
        PcorrR_matrix_path = os.path.abspath(os.path.join(sPathData, os.pardir,  'Matrix_PcorrR_Split.' + os.path.basename(sPathData).split('_')[0])) + ".mat"
        PcorrP_matrix_path = os.path.abspath(os.path.join(sPathData, os.pardir,  'Matrix_PcorrP_Split.' + os.path.basename(sPathData).split('_')[0])) + ".mat"
        PcorrZ_matirx_path = os.path.abspath(os.path.join(sPathData, os.pardir,  'Matrix_PcorrZ_Split.' + os.path.basename(sPathData).split('_')[0])) + ".mat"
    elif os.path.basename(labelNames) == "SIGMA_InVivo_Anatomical_Brain_Atlas_Labels.txt":
        PcorrR_matrix_path = os.path.abspath(os.path.join(sPathData, os.pardir,  'Matrix_PcorrR.' + os.path.basename(sPathData).split('_')[0])) + ".mat"
        PcorrP_matrix_path = os.path.abspath(os.path.join(sPathData, os.pardir,  'Matrix_PcorrP.' + os.path.basename(sPathData).split('_')[0])) + ".mat"
        PcorrZ_matirx_path = os.path.abspath(os.path.join(sPathData, os.pardir,  'Matrix_PcorrZ.' + os.path.basename(sPathData).split('_')[0])) + ".mat"
    
    
    pcorr_paths = [PcorrR_matrix_path, PcorrP_matrix_path, PcorrZ_matirx_path]


    if len(data_shape) != 4:
        sys.exit("Error: data %s has no 4D shape." % (str(data_shape),))

    # Read 4D mask file (NIfTI)
    mask_img = nib.load(sPathMask)
    mask = mask_img.get_data()
    mask_hdr = mask_img.header
    #mask_dtype = mask_hdr.get_data_dtype()
    mask_shape = mask_hdr.get_data_shape()

    if len(mask_shape) != 4:
        sys.exit("Error: mask %s has no 4D shape." % (str(mask_shape),))

    if data_shape[:3] != mask_shape[:3]:
        sys.exit("Error: data %s and mask %s are not the same shape." % (str(data_shape[:3]), str(mask_shape[:3])))

    m = np.zeros((mask_shape[3], data_shape[3]), dtype=data_dtype)
    for k in range(mask_shape[3]):
        msk = np.array(mask[:, :, :, k]) > 0
        maskedData = data[msk, :]
        if maskedData.size > 0:
            m[k] = np.mean(data[msk, :], 0)

    fileNames = open(labelNames, 'r');
    lines = fileNames.readlines()
    mT = np.transpose(m)
    np.savetxt(sPathOut, mT, fmt='%.4f', delimiter=' ')
    matPathOut = os.path.join(os.path.dirname(sPathOut), os.path.basename(sPathOut) + '.mat')
    io.savemat(matPathOut, dict([('matrix', mT),('label',lines)]))


    correlate_matrix.calculate_p_corr_matrix(mT, lines, pcorr_paths)
        
    return sPathOut


def get_date():
    now = datetime.now()
    pvDate = now.strftime("%a %d %b %Y")
    pvTime = now.strftime("%H:%M:%S")
    return pvDate + ' ' + pvTime


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('in_data', nargs='?', default='', help='input 4D data (x, y, slc, rep)')
    parser.add_argument('in_mask', nargs='?', default='', help='input 4D mask (x, y, slc, msk)')
    #parser.add_argument('in_data', help='input 4D data (x, y, slc, rep)')
    #parser.add_argument('in_mask', help='input 4D mask (x, y, slc, msk)')
    parser.add_argument('-o', '--out_text', help='output text matrix')
    args = parser.parse_args()

    sPathData = None
    sPathMask = None

    # input data
    if len(args.in_data) > 0: sPathData = args.in_data
    if not os.path.isfile(sPathData):
        sys.exit("Error: '%s' is not a regular file." % (sPathData,))

    # input mask
    if len(args.in_mask) > 0: sPathMask = args.in_mask
    if not os.path.isfile(sPathMask):
        sys.exit("Error: '%s' is not a regular file." % (sPathMask,))

    # output data
    sPathOut = args.out_text if args.out_text is not None else os.path.abspath(os.path.join(sPathData, os.pardir, 'MasksTCs.' +  os.path.basename(sPathData).split('_')[0]))
    sPathOut = sPathOut + '.txt'

    #print(get_date())

    # Read 4D data file (NIfTI )
    #print(sPathData)
    data_img = nib.load(sPathData)
    data = data_img.get_data()
    #data = np.squeeze(data_img.get_data())
    #data = np.cast[np.float32](data_img.get_data())
    #print("data.dtype:", data.dtype)
    #print("data.shape:", data.shape)
    data_hdr = data_img.header
    data_dtype = data_hdr.get_data_dtype()
    data_shape = data_hdr.get_data_shape()
    #print("data_dtype:", data_dtype)
    #print("data_shape:", data_shape)
    if len(data_shape) != 4:
        sys.exit("Error: data %s has no 4D shape." % (str(data_shape),))

    # Read 4D mask file (NIfTI)
    #print(sPathMask)
    mask_img = nib.load(sPathMask)
    mask = mask_img.get_data()
    #mask = np.squeeze(mask_img.get_data())
    #mask = np.cast[np.float32](mask_img.get_data())
    #print("mask.dtype:", mask.dtype)
    #print("mask.shape:", mask.shape)
    mask_hdr = mask_img.header
    mask_dtype = mask_hdr.get_data_dtype()
    mask_shape = mask_hdr.get_data_shape()
    #print("mask_dtype:", mask_dtype)
    #print("mask_shape:", mask_shape)
    if len(mask_shape) != 4:
        sys.exit("Error: mask %s has no 4D shape." % (str(mask_shape),))

    if data_shape[:3] != mask_shape[:3]:
        sys.exit("Error: data %s and mask %s are not the same shape." % (str(data_shape[:3]), str(mask_shape[:3])))

    m = np.zeros((mask_shape[3], data_shape[3]), dtype=data_dtype)
    for k in range(mask_shape[3]):
        msk = np.array(mask[:,:,:,k]) > 0
        maskedData = data[msk,:]
        if maskedData.size > 0:
            m[k] = np.mean(data[msk,:], 0)

    #s = [['%.4f' % (x,) for x in line] for line in m.T.tolist()]
    #s = [map(lambda x: '%.4f' % (x,), line) for line in m.T.tolist()]



    mT = np.transpose(m)
    np.savetxt(sPathOut, mT , fmt='%.4f',delimiter=' ')
    matPathOut = os.path.join(os.path.dirname(sPathOut), os.path.basename(sPathOut) + '.mat')

    io.savemat(matPathOut, dict([('matrix', mT)]))
    #print(sPathOut)
    #save_csv(sPathOut, s)
    #save_data(sPathOut, s)
