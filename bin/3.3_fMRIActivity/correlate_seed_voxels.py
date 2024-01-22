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
import logging

def get_date():
    now = datetime.now()
    pvDate = now.strftime("%a %d %b %Y")
    pvTime = now.strftime("%H:%M:%S")
    return pvDate + ' ' + pvTime

def save_csv(sFilename, data):

    thefile = open(sFilename, 'w')
    for item in data:
        thefile.write("%s\n" % item)
    thefile.close()

def save_nifti(sFilename, data, index, ext_nii):

    # save matrix (NIfTI)
    image = nib.Nifti1Image(data, None)
    header = image.get_header()
    header.set_xyzt_units(xyz=None, t=None)
    image.to_filename(sFilename + '_%03d' % (index + 1,) + ext_nii)
    logging.info("Output:", image.get_filename())

def get_seed_stat(sPathMatrix, sPathTS, data, seed, ext_nii, r_to_z=False, save_mat=False, ignore_nan=False):
    seed_stat = np.zeros((5, seed.shape[3]), dtype=np.float64)
    for k in range(seed.shape[3]):
        msk = seed[:,:,:,k] > 0
        maskData = data[msk, :]
        if maskData.size == 0:
            matrix[0] = np.nan
        else:
            matrix = np.corrcoef(maskData, rowvar=True)
        if r_to_z:
            matrix = np.arctanh(matrix)
        if save_mat:
            #pos = np.where(seed[:,:,:,k] > 0)
            #labels = [', '.join(str(v) for v in t) for t in zip(pos[0], pos[1], pos[2])]
            save_nifti(sPathTS, data[msk,:].T, k, ext_nii)
            save_nifti(sPathMatrix, matrix, k, ext_nii)
        # get upper-triangle of matrix as list and set inf to nan
        triu_cc = matrix[np.triu_indices_from(matrix, k=1)]
        #triu_cc[np.where(np.isinf(triu_cc))] = np.nan
        triu_cc[np.isinf(triu_cc)] = np.nan
        seed_stat[0,k] = matrix.shape[0]
        if ignore_nan:
            seed_stat[1,k] = np.nanmin(triu_cc)
            seed_stat[2,k] = np.nanmax(triu_cc)
            seed_stat[3,k] = np.nanmean(triu_cc)
            seed_stat[4,k] = np.nanstd(triu_cc)
        else:
            seed_stat[1,k] = np.amin(triu_cc)
            seed_stat[2,k] = np.amax(triu_cc)
            seed_stat[3,k] = np.mean(triu_cc)
            seed_stat[4,k] = np.std(triu_cc)
    return seed_stat

def make_text_stat(sPathData, sPathSeed, seed_stat):
    text_mean = []
    line = ['Data', 'Seed_ROIs', 'ROI_Index', 'Voxels', 'Min', 'Max', 'Mean', 'StdDev']
    text_mean.append(line)
    for k in range(seed_stat.shape[1]):
        stat = seed_stat[:,k]
        line = [os.path.basename(sPathData), os.path.basename(sPathSeed), str(k + 1), '%d' % (stat[0],), '%.8f' % (stat[1],), '%.8f' % (stat[2],), '%.8f' % (stat[3],), '%.8f' % (stat[4],)]
        text_mean.append(line)
    return text_mean

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create correlation matrix of seed voxels.')
    parser.add_argument('in_data', help='input 4D EPI data file name (NIfTI)')
    parser.add_argument('in_seed', help='input 4D seed ROIs file name (NIfTI)')
    parser.add_argument('-m', '--out_matrix', help='output seed ROIs matrix file name (w/o ext.)')
    parser.add_argument('-o', '--out_stat', help='output seed ROIs statistic text file name')
    args = parser.parse_args()

    sPathData = None
    sPathSeed = None

    ext_nii  = '.nii'
    ext_gz   = '.gz' 
    ext_text = '.txt'

    # input 4D EPI data file (NIfTI)
    if args.in_data is not None: sPathData = args.in_data
    if not os.path.isfile(sPathData):
        sys.exit("Error: '%s' is not a regular file." % (sPathData,))

    # input 4D seed ROIs file (NIfTI)
    if args.in_seed is not None: sPathSeed = args.in_seed
    if not os.path.isfile(sPathSeed):
        sys.exit("Error: '%s' is not a regular file." % (sPathSeed,))

    sData = os.path.basename(sPathData)
    sData = sData[:-len(ext_nii)] if sData.endswith(ext_nii) else sData[:-len(ext_nii+ext_gz)]

    # output seed ROIs matrix file
    sPathMatrix = os.path.join(os.getcwd(), args.out_matrix) if args.out_matrix is not None else os.path.join(os.path.dirname(sPathData), 'Matrix.' + sData)

    # output seed ROIs statistic text file
    sPathStat = os.path.join(os.getcwd(), args.out_stat) if args.out_stat is not None else os.path.join(os.path.dirname(sPathData), 'Stat.' + sData + ext_text)

    # output seed ROIs time series file
    sPathTS = os.path.join(os.path.dirname(sPathData), 'TS.' + sData)

    # get date and time
    print(get_date())

    # read 3D data file (NIfTI)
    print("Data:", sPathData)
    data_img = nib.load(sPathData)
    data_data = data_img.get_data()
    #print("data_data.dtype:", data_data.dtype)
    #print("data_data.shape:", data_data.shape)
    data_hdr = data_img.get_header()
    data_shape = data_hdr.get_data_shape()
    #print("data_shape:", data_shape)
    if len(data_shape) != 4:
        sys.exit("Error: EPI data %s do not have four dimensions." % (str(data_shape),))

    # read 4D seed ROIs file (NIfTI)
    print("Seed:", sPathSeed)
    seed_img = nib.load(sPathSeed)
    seed_data = seed_img.get_data()
    #print("seed_data.dtype:", seed_data.dtype)
    #print("seed_data.shape:", seed_data.shape)
    seed_hdr = seed_img.get_header()
    seed_shape = seed_hdr.get_data_shape()
    #print("seed_shape:", seed_shape)
    if len(seed_shape) != 4:
        sys.exit("Error: Seed ROIs %s do not have four dimensions." % (str(seed_shape),))

    seed_stat = get_seed_stat(sPathMatrix, sPathTS, data_data, seed_data, ext_nii, r_to_z=True, save_mat=False if args.out_matrix is None else True, ignore_nan=True)
    text_stat = make_text_stat(sPathData, sPathSeed, seed_stat)
    save_csv(sPathStat, text_stat)
    print('Stat:', sPathStat)
