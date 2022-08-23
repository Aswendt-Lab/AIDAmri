'''
Created on 19.10.2020

Author:
Michael Diedenhofen
Max Planck Institute for Metabolism Research, Cologne

Description:
Pre-requisits: stroke mask was defined at post stroke day 7 (P7)
Result: for all time points the peri-infarct mask is created aligned to the individual T2w MRI data

1. Time point P7: For each subject of the two groups a peri-infarct mask is generated from the stroke mask.
  - The input stroke mask <subject>Stroke_mask.nii.gz is located in the T2w subfolder.
  - The output peri-infarct mask <subject>_peri_mask_m3_n15.nii.gz is stored in the T2w subfolder.
  - The SciPy function binary_dilation() is called for each slice of the stroke mask with a parameter struct if at least one pixel value is greater than 0.
  - struct is a mask of size [2*R+1, 2*R+1] filled with a circular disk (radius R=15 pixels).
  - In order to obtain the peri-infarct mask the original stroke mask is subtracted from the dilated stroke mask.

2. Time point P7: For each subject of the two groups a non-rigid transformation (from template to T2w MRI) is inverted with NiftyReg.
  - The non-rigid transformation <subject>BiasBetMatrixBspline.nii (-invNrr filename1) is inverted with NiftyReg.
    NiftyReg command: reg_transform -ref <filename> -invNrr <filename1> <filename2> <filename3>

3. All time points except P7: For each subject of the two groups a peri-infarct mask from time point P7 is transformed with NiftyReg.
  - Compose the non-rigid transformation (template -> T2w) with the inverse non-rigid transformation from time point P7 (step 2).
    NiftyReg command: reg_transform -ref <filename> -ref2 <filename> -comp <filename1> <filename2> <filename3>
  - Apply the combined non-rigid transformation to the peri-infarct mask from time point P7.
    NiftyReg command: reg_resample -ref <filename> -flo <filename> -res <filename> -trans <filename> -inter 0
'''


from __future__ import print_function

try:
    zrange = xrange
except NameError:
    zrange = range

import os
import sys

import numpy as np

import create_seed_rois as csr
import dilate_mask as dm
import proc_tools as pt

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

def create_peri_mask(timepoint, group, subject, na=[15]):
    in_dir = os.path.join(pt.proc_in_dir, timepoint, group, subject, 'T2w')
    out_dir = os.path.join(pt.proc_out_dir, timepoint, group, subject, 'T2w')

    if not os.path.isdir(in_dir):
        sys.exit("Error: '%s' is not an existing directory." % (in_dir,))

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # input atlas labels file (NIfTI)
    #path_atlas = os.path.join(in_dir, subject + 'BiasBet_AnnorsfMRI.nii.gz')

    # input mask file (NIfTI)
    path_in_mask = os.path.join(in_dir, subject + 'Stroke_mask.nii.gz')
    if not os.path.isfile(path_in_mask):
        sys.exit("Error: '%s' is not a regular file." % (path_in_mask,))

    # output cortex ROIs file (NIfTI)
    #path_out_rois = os.path.join(out_dir, subject + '_cortex_rois_1.nii.gz')

    mask, mask_dims = pt.read_data(path_in_mask)

    #rois, rois_dims = create_rois_1(pt.path_labels_1, path_atlas, path_rois=path_out_rois)

    for model in range(3, 4):
        for n in na:
            # output mask file (NIfTI)
            path_out_mask = os.path.join(out_dir, subject + '_peri_mask_m%d_n%d.nii.gz' % (model, n))

            # output masked cortex ROIs file (NIfTI)
            #path_out_rois = os.path.join(out_dir, subject + '_cortex_rois_1_m%d_n%d.nii.gz' % (model, n))

            peri = np.copy(mask)

            if (model == 1) or (model == 2):
                for k in zrange(mask.shape[2]):
                    image = mask[:, :, k]
                    if np.any(image.astype(np.bool)):
                        peri[:, :, k] = dm.dilate_repeat(image, connectivity=model, n=n)
            else:
                struct = dm.circle_mask(n=n)
                for k in zrange(mask.shape[2]):
                    image = mask[:, :, k]
                    if np.any(image.astype(np.bool)):
                        peri[:, :, k] = dm.dilate_struct(image, struct)

            pt.save_data(peri.astype(np.float32), mask_dims, path_out_mask, dtype=None)

            #rois_masked = np.multiply(rois, peri.astype(np.float32))
            #pt.save_data(rois_masked, rois_dims, path_out_rois, dtype=None)

def xfm_inv(timepoint, group, subject):
    in_dir = os.path.join(pt.proc_in_dir, timepoint, group, subject, 'T2w')
    out_dir = os.path.join(pt.proc_out_dir, timepoint, group, subject, 'T2w')

    if not os.path.isdir(in_dir):
        sys.exit("Error: '%s' is not an existing directory." % (in_dir,))

    if not os.path.isdir(out_dir):
        sys.exit("Error: '%s' is not an existing directory." % (out_dir,))
    
    brain_template = os.path.join(pt.lib_in_dir, 'NP_template_sc0.nii.gz')
    input_volume = os.path.join(in_dir, subject + 'BiasBet.nii.gz')
    output_aff = os.path.join(in_dir, subject + 'BiasBetMatrixAff.txt')
    output_aff_inv = os.path.join(out_dir, subject + 'BiasBetMatrixAff_inv.txt')
    output_cpp = os.path.join(in_dir, subject + 'BiasBetMatrixBspline.nii')
    output_cpp_inv = os.path.join(out_dir, subject + 'BiasBetMatrixBspline_inv.nii.gz')

    if not os.path.isfile(brain_template):
        sys.exit("Error: '%s' is not a regular file." % (brain_template,))

    if not os.path.isfile(input_volume):
        sys.exit("Error: '%s' is not a regular file." % (input_volume,))

    if not os.path.isfile(output_cpp):
        sys.exit("Error: '%s' is not a regular file." % (output_cpp,))

    # inverse affine transformation
    command = 'reg_transform -invAff %s %s' % (output_aff, output_aff_inv)
    os.system(command)
    print(output_aff_inv)

    # inverse transformation
    command = 'reg_transform -ref %s -invNrr %s %s %s' % (input_volume, output_cpp, brain_template, output_cpp_inv)
    os.system(command)
    print(output_cpp_inv)

def xfm_peri_mask(timepoint_P7, timepoint, group, subject_P7, subject):
    in_dir_P7 = os.path.join(pt.proc_in_dir, timepoint_P7, group, subject_P7, 'T2w')
    in_dir = os.path.join(pt.proc_in_dir, timepoint, group, subject, 'T2w')
    out_dir_P7 = os.path.join(pt.proc_out_dir, timepoint_P7, group, subject_P7, 'T2w')
    out_dir = os.path.join(pt.proc_out_dir, timepoint, group, subject, 'T2w')

    if not os.path.isdir(in_dir_P7):
        sys.exit("Error: '%s' is not an existing directory." % (in_dir_P7,))

    if not os.path.isdir(in_dir):
        sys.exit("Error: '%s' is not an existing directory." % (in_dir,))

    if not os.path.isdir(out_dir_P7):
        sys.exit("Error: '%s' is not an existing directory." % (out_dir_P7,))

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    brain_template = os.path.join(pt.lib_in_dir, 'NP_template_sc0.nii.gz')
    input_volume = os.path.join(in_dir, subject + 'BiasBet.nii.gz')
    output_cpp = os.path.join(in_dir, subject + 'BiasBetMatrixBspline.nii')
    output_cpp_inv = os.path.join(out_dir_P7, subject_P7 + 'BiasBetMatrixBspline_inv.nii.gz')
    output_cpp_comp = os.path.join(out_dir, subject + 'BiasBetMatrixBspline_comp.nii.gz')
    mask_P7 = os.path.join(in_dir_P7, subject_P7 + 'Stroke_mask.nii.gz')
    mask = os.path.join(out_dir, subject + 'Stroke_mask.nii.gz')
    peri_P7 = os.path.join(out_dir_P7, subject_P7 + '_peri_mask_m3_n15.nii.gz')
    peri = os.path.join(out_dir, subject + '_peri_mask_m3_n15.nii.gz')

    if not os.path.isfile(brain_template):
        sys.exit("Error: '%s' is not a regular file." % (brain_template,))

    if not os.path.isfile(input_volume):
        sys.exit("Error: '%s' is not a regular file." % (input_volume,))

    if not os.path.isfile(output_cpp):
        sys.exit("Error: '%s' is not a regular file." % (output_cpp,))

    if not os.path.isfile(mask_P7):
        sys.exit("Error: '%s' is not a regular file." % (mask_P7,))

    if not os.path.isfile(peri_P7):
        sys.exit("Error: '%s' is not a regular file." % (peri_P7,))

    # compose transformations
    command = 'reg_transform -ref %s -ref2 %s -comp %s %s %s' % (input_volume, brain_template, output_cpp, output_cpp_inv, output_cpp_comp)
    os.system(command)

    # resample stroke mask
    command = 'reg_resample -ref %s -flo %s -res %s -trans %s -inter 0' % (input_volume, mask_P7, mask, output_cpp_comp)
    os.system(command)

    # resample peri-infarct mask
    command = 'reg_resample -ref %s -flo %s -res %s -trans %s -inter 0' % (input_volume, peri_P7, peri, output_cpp_comp)
    os.system(command)

def main():
    timepoint_P7 = pt.timepoints[1]

    # timepoint P7
    for index_g, group in enumerate(pt.groups):
        for subject in pt.study[1][index_g]:
            if subject is not None:
                create_peri_mask(timepoint_P7, group, subject)
                xfm_inv(timepoint_P7, group, subject)

    # all timepoints except P7
    for index_t, timepoint in enumerate(pt.timepoints):
        if index_t != 1:
            for index_g, group in enumerate(pt.groups):
                for index_s, subject in enumerate(pt.study[index_t][index_g]):
                    if subject is not None:
                        xfm_peri_mask(timepoint_P7, timepoint, group, pt.study[1][index_g][index_s], subject)

if __name__ == '__main__':
    main()
