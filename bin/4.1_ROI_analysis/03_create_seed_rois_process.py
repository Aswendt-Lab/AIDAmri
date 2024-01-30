'''
Created on 19.10.2020

Author:
Michael Diedenhofen
Max Planck Institute for Metabolism Research, Cologne

Description:
Pre-requisite: 02_apply_xfm_process.py
Result:
rsfMRI - a Matlab file which contains two text files: 1) for each region one column with the averaged rsfMRI time series and 2) the atlas labels names.
DTI - atlas labels file modified to include individually shaped peri-infarct brain regions which replace the original regions

The text file annotation_50CHANGEDanno_label_IDs+2000.txt contains all atlas labels and another text file contains selected cortical peri-infarct atlas labels.
rsfMRI:
1. For all time points and for each subject of the two groups a hyperstack with modified selected cortical regions is created.
  - The peri-infarct mask <subject>_T2w_peri_mask_rsfMRI.nii.gz and the atlas labels <subject>_T2w_Anno_rsfMRI.nii.gz are located in the fMRI subfolder.
  - The output hyperstack Seed_ROIs_all_mod_peri.nii.gz is stored in the fMRI subfolder.
  - Create a first hyperstack with all regions (each region is one volume) from the atlas labels.
  - Create a second hyperstack with selected cortical regions from the atlas labels and apply the rsfMRI peri-infarct mask.
  - Create a third hyperstack with all regions but replaced selected cortical regions from the second hyperstack.
2. For each region of the modified hyperstack an averaged rsfMRI time series is computed and a text file with one column for each region is created.
  - The atlas labels names are listed in the annoVolume+2000_rsfMRI.nii.txt text file.
  - The input rsfMRI file <subject>_mcf_f_SFRGR.nii.gz is located in the fMRI/regr subfolder.
  - The resulting text file MasksTCsSplit_GV_all_mod_peri.txt and Matlab file MasksTCsSplit_GV_all_mod_peri.txt.mat are stored in the fMRI/regr subfolder.
  - The modified hyperstack is used with the rsfMRI data and the averaged time series of a region is computed from the voxels of the rsfMRI data which belong to this region.
  - The resulting text file contains for each hyperstack region one column with the averaged rsfMRI time series.
  - This text file is combined with the atlas labels names and stored as a Matlab file.
DTI:
1. For all time points and for each subject of the two groups a hyperstack and atlas labels with modified selected cortical regions are created.
  - The peri-infarct mask <subject>_T2w_peri_mask_DTI.nii.gz and the atlas labels <subject>_T2w_Anno_DTI.nii.gz are located in the DTI subfolder.
  - The output hyperstack Seed_ROIs_all_mod_peri.nii.gz is stored in the DTI subfolder.
  - The output atlas labels file <subject>_T2w_Anno_DTI_mod_peri_scaled.nii.gz is stored in the DTI/DSI_studio subfolder.
  - Create a first hyperstack with all regions (each region is one volume) from the atlas labels.
  - Create a second hyperstack with selected cortical regions from the atlas labels and apply the DTI peri-infarct mask.
  - Create a third hyperstack with all regions but replaced selected cortical regions from the second hyperstack.
  - A maximum intensity projection of the third hyperstack generates atlas labels with all regions but replaced selected cortical regions.
  - The voxel dimensions of the output atlas labels file are scaled by a factor of 10 which is required for DSI Studio.
'''

from __future__ import print_function

try:
    zrange = xrange
except NameError:
    zrange = range

import os
import sys

import shutil

import numpy as np

import create_seed_rois as csr
import fsl_mean_ts as mts
import proc_tools as pt

def read_mask(path_in_mask):
    if not os.path.isfile(path_in_mask):
        sys.exit("Error: '%s' is not a regular file." % (path_in_mask,))

    # read mask file (NIfTI)
    mask, voxel_dims = pt.read_data(path_in_mask)

    return (mask, voxel_dims)

def create_rois_1(path_labels, path_atlas, path_rois=None, mask=None):
    if not os.path.isfile(path_atlas):
        sys.exit("Error: '%s' is not a regular file." % (path_atlas,))

    # create atlas labels ROIs
    labels_hdr, rois = csr.create_rois(path_labels, [path_atlas], datatype=16, preserve=True)
    voxel_dims = labels_hdr[0].get_zooms()

    rois = np.squeeze(rois)

    # apply mask to ROIs
    if mask is not None:
        rois = np.multiply(rois, mask)

    # save ROIs file (NIfTI)
    if path_rois is not None:
        pt.save_data(rois, voxel_dims, path_rois, dtype=None)

    return (rois, voxel_dims)

def create_rois_2(path_labels, path_atlas, path_rois=None, mask=None, preserve=False):
    if not os.path.isfile(path_atlas):
        sys.exit("Error: '%s' is not a regular file." % (path_atlas,))

    # create atlas labels ROIs hyperstack (4D)
    labels_hdr, rois = csr.create_rois(path_labels, [path_atlas], datatype=16, preserve=preserve)
    voxel_dims = labels_hdr[0].get_zooms()

    # apply mask to each ROI
    if mask is not None:
        for k in zrange(rois.shape[3]):
            rois[:, :, :, k] = np.multiply(rois[:, :, :, k], mask)

    # save ROIs file (NIfTI)
    if path_rois is not None:
        pt.save_data(rois, voxel_dims, path_rois, dtype=None)

    return (rois, voxel_dims)

def replace_rois(labels_1, labels_2, rois_1, rois_2):
    for index_2, label in enumerate(labels_2):
        if label in labels_1:
            index_1 = labels_1.index(label)
            rois_1[:, :, :, index_1] = rois_2[:, :, :, index_2]

def create_rois_rsfMRI(timepoint, group, subject, labels, labels_2, label_names_2000, label_names_peri):
    in_dir = os.path.join(pt.proc_in_dir, timepoint, group, subject, 'fMRI')
    out_dir = os.path.join(pt.proc_out_dir, timepoint, group, subject, 'fMRI')
    regr_dir = os.path.join(pt.proc_out_dir, timepoint, group, subject, 'fMRI', 'regr')

    if not os.path.isdir(in_dir):
        sys.exit("Error: '%s' is not an existing directory." % (in_dir,))

    if not os.path.isdir(out_dir):
        sys.exit("Error: '%s' is not an existing directory." % (out_dir,))

    if not os.path.exists(regr_dir):
        os.makedirs(regr_dir)

    # input atlas labels file (NIfTI)
    #path_atlas = os.path.join(in_dir, subject + 'SmoothBet_AnnoSplit_rsfMRI.nii.gz')
    path_atlas = os.path.join(out_dir, subject + '_T2w_Anno_rsfMRI.nii.gz')

    # input peri-infarct mask file (NIfTI)
    path_mask = os.path.join(out_dir, subject + '_T2w_peri_mask_rsfMRI.nii.gz')
    peri_mask, _ = read_mask(path_mask)

    # input rsfMRI file (NIfTI)
    path_rsfMRI = os.path.join(in_dir, 'regr', subject + '_mcf_f_SFRGR.nii.gz')
    if not os.path.isfile(path_rsfMRI):
        sys.exit("Error: '%s' is not a regular file." % (path_rsfMRI,))

    # output ROIs files (NIfTI)
    #path_out_rois = os.path.join(out_dir, subject + '_seed_rois.nii.gz')
    #path_out_rois_1 = os.path.join(out_dir, subject + '_cortex_rois_1.nii.gz')
    #path_out_rois_2 = os.path.join(out_dir, subject + '_cortex_rois_2.nii.gz')
    #path_out_rois_x = os.path.join(out_dir, subject + '_seed_rois_mod.nii.gz')
    path_out_rois = os.path.join(out_dir, 'Seed_ROIs_all.nii.gz')
    path_out_rois_1 = os.path.join(out_dir, subject + '_rois_peri.nii.gz')
    path_out_rois_2 = os.path.join(out_dir, 'Seed_ROIs_peri.nii.gz')
    path_out_rois_x = os.path.join(out_dir, 'Seed_ROIs_all_mod_peri.nii.gz')

    # output time series text files
    #path_out_ts = os.path.join(out_dir, subject + '_ts.txt')
    #path_out_ts_2 = os.path.join(out_dir, subject + '_ts_cortex.txt')
    #path_out_ts_x = os.path.join(out_dir, subject + '_ts_mod.txt')
    path_out_ts = os.path.join(regr_dir, 'MasksTCsSplit_GV_all.txt')
    path_out_ts_2 = os.path.join(regr_dir, 'MasksTCsSplit_GV_peri.txt')
    path_out_ts_x = os.path.join(regr_dir, 'MasksTCsSplit_GV_all_mod_peri.txt')

    rois, rois_dims = create_rois_2(pt.path_labels, path_atlas, path_rois=path_out_rois)
    _, _ = create_rois_1(pt.path_labels_1, path_atlas, path_rois=path_out_rois_1, mask=peri_mask)
    rois_2, _ = create_rois_2(pt.path_labels_2, path_atlas, path_rois=path_out_rois_2, mask=peri_mask)

    rois_x = np.copy(rois)
    replace_rois(labels, labels_2, rois_x, rois_2)
    pt.save_data(rois_x, rois_dims, path_out_rois_x, dtype=None)

    mts.mean_ts(path_rsfMRI, path_out_rois, path_out_ts, label_names_2000)
    mts.mean_ts(path_rsfMRI, path_out_rois_2, path_out_ts_2, label_names_peri)
    mts.mean_ts(path_rsfMRI, path_out_rois_x, path_out_ts_x, label_names_2000)

def create_rois_DTI(timepoint, group, subject, labels, labels_2, scale=10.0):
    in_dir = os.path.join(pt.proc_in_dir, timepoint, group, subject, 'DTI')
    out_dir = os.path.join(pt.proc_out_dir, timepoint, group, subject, 'DTI')
    dti_dir = os.path.join(pt.proc_out_dir, timepoint, group, subject, 'DTI', 'DSI_studio')

    if not os.path.isdir(in_dir):
        #sys.exit("Error: '%s' is not an existing directory." % (in_dir,))
        return

    if not os.path.isdir(out_dir):
        sys.exit("Error: '%s' is not an existing directory." % (out_dir,))

    if not os.path.exists(dti_dir):
        os.makedirs(dti_dir)

    # input atlas labels file (NIfTI)
    #path_atlas = os.path.join(in_dir, subject + 'SmoothMicoBet_AnnoSplit_rsfMRI.nii.gz')
    path_atlas = os.path.join(out_dir, subject + '_T2w_Anno_DTI.nii.gz')
    if not os.path.isfile(path_atlas):
        #sys.exit("Error: '%s' is not a regular file." % (path_atlas,))
        return

    # input peri-infarct mask file (NIfTI)
    path_mask = os.path.join(out_dir, subject + '_T2w_peri_mask_DTI.nii.gz')
    peri_mask, _ = read_mask(path_mask)

    # output ROIs files (NIfTI)
    #path_out_rois = os.path.join(out_dir, subject + '_seed_rois.nii.gz')
    #path_out_rois_1 = os.path.join(out_dir, subject + '_cortex_rois_1.nii.gz')
    #path_out_rois_2 = os.path.join(out_dir, subject + '_cortex_rois_2.nii.gz')
    #path_out_rois_x = os.path.join(out_dir, subject + '_seed_rois_mod.nii.gz')
    path_out_rois = os.path.join(out_dir, 'Seed_ROIs_all.nii.gz')
    path_out_rois_1 = os.path.join(out_dir, subject + '_rois_peri.nii.gz')
    path_out_rois_2 = os.path.join(out_dir, 'Seed_ROIs_peri.nii.gz')
    path_out_rois_x = os.path.join(out_dir, 'Seed_ROIs_all_mod_peri.nii.gz')

    # output maximum intensity projection files (NIfTI)
    path_out_mip_rois = os.path.join(dti_dir, subject + '_T2w_Anno_DTI_scaled.nii.gz')
    path_out_mip_rois_x = os.path.join(dti_dir, subject + '_T2w_Anno_DTI_mod_peri_scaled.nii.gz')

    # output maximum intensity projection label names text files
    path_out_label_names = os.path.join(dti_dir, subject + '_T2w_Anno_DTI_scaled.nii.txt')
    path_out_label_names_x = os.path.join(dti_dir, subject + '_T2w_Anno_DTI_mod_peri_scaled.nii.txt')

    rois, rois_dims = create_rois_2(pt.path_labels, path_atlas, path_rois=path_out_rois, preserve=True)
    _, _ = create_rois_1(pt.path_labels_1, path_atlas, path_rois=path_out_rois_1, mask=peri_mask)
    rois_2, _ = create_rois_2(pt.path_labels_2, path_atlas, path_rois=path_out_rois_2, mask=peri_mask, preserve=True)

    rois_x = np.copy(rois)
    replace_rois(labels, labels_2, rois_x, rois_2)
    pt.save_data(rois_x, rois_dims, path_out_rois_x, dtype=None)

    mip_rois = np.max(rois, axis=3)
    mip_rois_x = np.max(rois_x, axis=3)
    mip_rois_dims = tuple((x * scale) for x in rois_dims)
    pt.save_data(mip_rois, mip_rois_dims, path_out_mip_rois, dtype=None)
    pt.save_data(mip_rois_x, mip_rois_dims, path_out_mip_rois_x, dtype=None)

    shutil.copyfile(pt.path_label_names_2000, path_out_label_names)
    shutil.copyfile(pt.path_label_names_2000, path_out_label_names_x)
    print(path_out_label_names)
    print(path_out_label_names_x)

def main():
    # read labels
    _, labels = pt.read_labels(pt.path_labels)
    _, labels_2 = pt.read_labels(pt.path_labels_2)
    labels = [x for list_x in labels for x in list_x]
    labels_2 = [x for list_x in labels_2 for x in list_x]

    label_names_2000 = pt.read_text(pt.path_label_names_2000)
    labels_2000 = [int(x.split('\t')[0]) for x in label_names_2000]
    label_names_peri = [label_names_2000[x] for x in list(np.where(np.in1d(labels_2000, labels_2))[0])]

    for index_t, timepoint in enumerate(pt.timepoints):
        for index_g, group in enumerate(pt.groups):
            for subject in pt.study[index_t][index_g]:
                if subject is not None:
                    create_rois_rsfMRI(timepoint, group, subject, labels, labels_2, label_names_2000, label_names_peri)
                    create_rois_DTI(timepoint, group, subject, labels, labels_2)

if __name__ == '__main__':
    main()
