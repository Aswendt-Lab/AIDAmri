'''
Created on 19.10.2020

Author:
Michael Diedenhofen
Max Planck Institute for Metabolism Research, Cologne

Description:
Pre-requisite: 01_dilate_mask_process.py
Result: for all time points the peri-infarct masks will be aligned in the rsfMRI and DTI space

Two scans of the same session can be aligned to each other without image registration if there was only very limited movement (otherwise image registration is necessary, e.g. between T2w and DTI).
The ParaVision visu_pars file provides a mapping (VisuCoreOrientation, VisuCorePosition) from subject (LPS) into the image coordinate system.
rsfMRI:
- For all time points and for each subject of the two groups a T2w peri-infarct mask is transformed to rsfMRI.
  - The input peri-infarct mask <subject>_peri_mask_m3_n15.nii.gz and atlas labels <subject>BiasBet_AnnorsfMRI.nii.gz are located in the T2w subfolder.
  - The transformed peri-infarct mask <subject>_T2w_peri_mask_rsfMRI.nii.gz and atlas labels <subject>_T2w_Anno_rsfMRI.nii.gz are stored in the fMRI subfolder.
  - Create the inverse of a rigid transformation matrix from the parameters of the T2w raw data visu_pars file.
  - Create a rigid transformation matrix from the parameters of the rsfMRI raw data visu_pars file.
  - Compose the inverse T2w matrix and the rsfMRI matrix to a rigid transformation matrix (rsfMRI -> T2w).
  - The function xfm_serial() (in apply_xfm.py) applies the inverse of the matrix (rsfMRI -> T2w) to the T2w peri-infarct mask.
  - The T2w peri-infarct mask is first flipped in x- and z-axis then transformed and then flipped back in x- and z-axis.
DTI:
- For all time points and for each subject of the two groups a T2w peri-infarct mask is transformed to DTI.
  - The transformed peri-infarct mask <subject>_T2w_peri_mask_DTI.nii.gz and atlas labels <subject>_T2w_Anno_DTI.nii.gz are stored in the DTI subfolder.
  - Create the inverse of a rigid transformation matrix from the parameters of the T2w raw data visu_pars file.
  - Create a rigid transformation matrix from the parameters of the DTI raw data visu_pars file.
  - Compose the inverse T2w matrix and the DTI matrix to a rigid transformation matrix (DTI -> T2w).
  - The function xfm_serial() applies the inverse of the matrix (DTI -> T2w) to the T2w peri-infarct mask.
  - The T2w peri-infarct mask is first flipped in x- and z-axis then transformed and then flipped back in x- and z-axis.
'''

from __future__ import print_function

import os
import sys

import numpy as np

import pv_reader as pvr
import proc_tools as pt
import apply_xfm as ax

def xfm_T2w_rsfMRI(raw_dir, timepoint_P7, timepoint, group, subject, expno_T2w, expno_rsfMRI, procno_T2w, procno_rsfMRI):
    if (expno_T2w is None) or (expno_rsfMRI is None) or (procno_T2w is None) or (procno_rsfMRI is None):
        return

    in_dir = os.path.join(pt.proc_in_dir, timepoint, group, subject, 'T2w')
    mask_dir = os.path.join(pt.proc_out_dir, timepoint, group, subject, 'T2w')
    out_dir = os.path.join(pt.proc_out_dir, timepoint, group, subject, 'fMRI')

    if not os.path.isdir(in_dir):
        sys.exit("Error: '%s' is not an existing directory." % (in_dir,))

    if not os.path.isdir(mask_dir):
        sys.exit("Error: '%s' is not an existing directory." % (mask_dir,))

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # input T2w atlas labels file
    path_in_anno = os.path.join(in_dir, subject + 'BiasBet_AnnorsfMRI.nii.gz')
    if not os.path.isfile(path_in_anno):
        sys.exit("Error: '%s' is not a regular file." % (path_in_anno,))

    # input T2w stroke mask file
    if timepoint == timepoint_P7:
        path_in_mask = os.path.join(in_dir, subject + 'Stroke_mask.nii.gz')
    else:
        path_in_mask = os.path.join(mask_dir, subject + 'Stroke_mask.nii.gz')

    # input T2w peri-infarct mask file
    path_in_peri = os.path.join(mask_dir, subject + '_peri_mask_m3_n15.nii.gz')
    if not os.path.isfile(path_in_peri):
        sys.exit("Error: '%s' is not a regular file." % (path_in_peri,))

    # output rsfMRI file
    #path_rsfMRI = os.path.join(out_dir, subject + '_rsfMRI.nii.gz')

    # output transformed T2w file
    path_T2w_rsfMRI = os.path.join(out_dir, subject + '_T2w_rsfMRI.nii.gz')

    # output transformed atlas labels file
    path_anno_rsfMRI = os.path.join(out_dir, subject + '_T2w_Anno_rsfMRI.nii.gz')

    # output transformed stroke mask file
    path_mask_rsfMRI = os.path.join(out_dir, subject + '_T2w_Stroke_mask_rsfMRI.nii.gz')

    # output transformed peri-infarct mask file
    path_peri_rsfMRI = os.path.join(out_dir, subject + '_T2w_peri_mask_rsfMRI.nii.gz')

    pvr.check_args(pt.proc_out_dir, raw_dir, subject, expno_T2w, procno_T2w)
    pvr.check_args(pt.proc_out_dir, raw_dir, subject, expno_rsfMRI, procno_rsfMRI)

    # T2w data
    pv = pvr.ParaVision(os.path.join(pt.proc_out_dir, timepoint, group), raw_dir, subject, expno_T2w, procno_T2w)
    pv.read_2dseq(map_raw=False, map_pv6=False, roll_fg=False, squeeze=False, compact=False, swap_vd=False, scale=1.0)
    #pv.save_nifti(ftype='NIFTI_GZ')
    matrix_T2w, matrix_T2w_inv = pv.get_matrix()
    data_T2w = pv.nifti_image.get_data()
    #data_dims_T2w = pv.data_dims[:3]
    #data_type_T2w = pv.data_type
    voxel_dims_T2w = pv.voxel_dims[:3]
    #voxel_unit_T2w = pv.voxel_unit

    # rsfMRI data
    pv = pvr.ParaVision(os.path.join(pt.proc_out_dir, timepoint, group), raw_dir, subject, expno_rsfMRI, procno_rsfMRI)
    pv.read_2dseq(map_raw=False, map_pv6=False, roll_fg=False, squeeze=False, compact=False, swap_vd=False, scale=1.0)
    #pv.save_nifti(ftype='NIFTI_GZ')
    matrix_rsfMRI, matrix_rsfMRI_inv = pv.get_matrix()
    #data_rsfMRI = np.mean(pv.nifti_image.get_data(), axis=3)
    data_dims_rsfMRI = pv.data_dims[:3]
    #data_type_rsfMRI = pv.data_type
    voxel_dims_rsfMRI = pv.voxel_dims[:3]
    #voxel_unit_rsfMRI = pv.voxel_unit

    # transformation matrix
    matrix_T2w_rsfMRI = np.dot(matrix_rsfMRI_inv, matrix_T2w)
    matrix_rsfMRI_T2w = np.dot(matrix_T2w_inv, matrix_rsfMRI)
    pt.save_matrix(os.path.join(out_dir, subject + '_T2w_rsfMRI.mat'), matrix_T2w_rsfMRI)
    pt.save_matrix(os.path.join(out_dir, subject + '_rsfMRI_T2w.mat'), matrix_rsfMRI_T2w)

    # save rsfMRI data as NIfTI file
    #pt.save_data(np.rot90(data_rsfMRI, k=2, axes=(0, 2)), voxel_dims_rsfMRI, path_rsfMRI, dtype=None)

    # save transformed T2w data as NIfTI file
    data_T2w_rsfMRI = ax.xfm_serial(data_T2w, matrix_rsfMRI_T2w, data_dims_rsfMRI, voxel_dims_rsfMRI, voxel_dims_T2w, interp=1, inverse=True)
    pt.save_data(np.rot90(data_T2w_rsfMRI, k=2, axes=(0, 2)), voxel_dims_rsfMRI, path_T2w_rsfMRI, dtype=None)

    # save transformed T2w atlas labels as NIfTI file
    data_anno, voxel_dims_anno = pt.read_data(path_in_anno)
    data_anno_rsfMRI = ax.xfm_serial(np.rot90(data_anno, k=2, axes=(0, 2)), matrix_rsfMRI_T2w, data_dims_rsfMRI, voxel_dims_rsfMRI, voxel_dims_anno, interp=0, inverse=True)
    pt.save_data(np.rot90(data_anno_rsfMRI, k=2, axes=(0, 2)), voxel_dims_rsfMRI, path_anno_rsfMRI, dtype=None)

    # save transformed T2w stroke mask as NIfTI file
    if os.path.isfile(path_in_mask):
        data_mask, voxel_dims_mask = pt.read_data(path_in_mask)
        data_mask_rsfMRI = ax.xfm_serial(np.rot90(data_mask, k=2, axes=(0, 2)), matrix_rsfMRI_T2w, data_dims_rsfMRI, voxel_dims_rsfMRI, voxel_dims_mask, interp=0, inverse=True)
        pt.save_data(np.rot90(data_mask_rsfMRI, k=2, axes=(0, 2)), voxel_dims_rsfMRI, path_mask_rsfMRI, dtype=None)

    # save transformed T2w peri-infarct mask as NIfTI file
    data_peri, voxel_dims_peri = pt.read_data(path_in_peri)
    data_peri_rsfMRI = ax.xfm_serial(np.rot90(data_peri, k=2, axes=(0, 2)), matrix_rsfMRI_T2w, data_dims_rsfMRI, voxel_dims_rsfMRI, voxel_dims_peri, interp=0, inverse=True)
    pt.save_data(np.rot90(data_peri_rsfMRI, k=2, axes=(0, 2)), voxel_dims_rsfMRI, path_peri_rsfMRI, dtype=None)

def xfm_T2w_DTI(raw_dir, timepoint_P7, timepoint, group, subject, expno_T2w, expno_DTI, procno_T2w, procno_DTI):
    if (expno_T2w is None) or (expno_DTI is None) or (procno_T2w is None) or (procno_DTI is None):
        return

    in_dir = os.path.join(pt.proc_in_dir, timepoint, group, subject, 'T2w')
    mask_dir = os.path.join(pt.proc_out_dir, timepoint, group, subject, 'T2w')
    out_dir = os.path.join(pt.proc_out_dir, timepoint, group, subject, 'DTI')

    if not os.path.isdir(in_dir):
        sys.exit("Error: '%s' is not an existing directory." % (in_dir,))

    if not os.path.isdir(mask_dir):
        sys.exit("Error: '%s' is not an existing directory." % (mask_dir,))

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # input T2w atlas labels file
    path_in_anno = os.path.join(in_dir, subject + 'BiasBet_AnnorsfMRI.nii.gz')
    if not os.path.isfile(path_in_anno):
        sys.exit("Error: '%s' is not a regular file." % (path_in_anno,))

    # input T2w stroke mask file
    if timepoint == timepoint_P7:
        path_in_mask = os.path.join(in_dir, subject + 'Stroke_mask.nii.gz')
    else:
        path_in_mask = os.path.join(mask_dir, subject + 'Stroke_mask.nii.gz')

    # input T2w peri-infarct mask file
    path_in_peri = os.path.join(mask_dir, subject + '_peri_mask_m3_n15.nii.gz')
    if not os.path.isfile(path_in_peri):
        sys.exit("Error: '%s' is not a regular file." % (path_in_peri,))

    # output DTI file
    #path_DTI = os.path.join(out_dir, subject + '_DTI.nii.gz')

    # output transformed T2w file
    path_T2w_DTI = os.path.join(out_dir, subject + '_T2w_DTI.nii.gz')

    # output transformed atlas labels file
    path_anno_DTI = os.path.join(out_dir, subject + '_T2w_Anno_DTI.nii.gz')

    # output transformed stroke mask file
    path_mask_DTI = os.path.join(out_dir, subject + '_T2w_Stroke_mask_DTI.nii.gz')

    # output transformed peri-infarct mask file
    path_peri_DTI = os.path.join(out_dir, subject + '_T2w_peri_mask_DTI.nii.gz')

    pvr.check_args(pt.proc_out_dir, raw_dir, subject, expno_T2w, procno_T2w)
    pvr.check_args(pt.proc_out_dir, raw_dir, subject, expno_DTI, procno_DTI)

    # T2w data
    pv = pvr.ParaVision(os.path.join(pt.proc_out_dir, timepoint, group), raw_dir, subject, expno_T2w, procno_T2w)
    pv.read_2dseq(map_raw=False, map_pv6=False, roll_fg=False, squeeze=False, compact=False, swap_vd=False, scale=1.0)
    #pv.save_nifti(ftype='NIFTI_GZ')
    matrix_T2w, matrix_T2w_inv = pv.get_matrix()
    data_T2w = pv.nifti_image.get_data()
    #data_dims_T2w = pv.data_dims[:3]
    #data_type_T2w = pv.data_type
    voxel_dims_T2w = pv.voxel_dims[:3]
    #voxel_unit_T2w = pv.voxel_unit

    # DTI data
    pv = pvr.ParaVision(os.path.join(pt.proc_out_dir, timepoint, group), raw_dir, subject, expno_DTI, procno_DTI)
    pv.read_2dseq(map_raw=False, map_pv6=False, roll_fg=False, squeeze=False, compact=False, swap_vd=False, scale=1.0)
    #pv.save_nifti(ftype='NIFTI_GZ')
    matrix_DTI, matrix_DTI_inv = pv.get_matrix()
    #data_DTI = np.mean(pv.nifti_image.get_data(), axis=3)
    data_dims_DTI = pv.data_dims[:3]
    #data_type_DTI = pv.data_type
    voxel_dims_DTI = pv.voxel_dims[:3]
    #voxel_unit_DTI = pv.voxel_unit

    # transformation matrix
    matrix_T2w_DTI = np.dot(matrix_DTI_inv, matrix_T2w)
    matrix_DTI_T2w = np.dot(matrix_T2w_inv, matrix_DTI)
    pt.save_matrix(os.path.join(out_dir, subject + '_T2w_DTI.mat'), matrix_T2w_DTI)
    pt.save_matrix(os.path.join(out_dir, subject + '_DTI_T2w.mat'), matrix_DTI_T2w)

    # save DTI data as NIfTI file
    #pt.save_data(np.rot90(data_DTI, k=2, axes=(0, 2)), voxel_dims_DTI, path_DTI, dtype=None)

    # save transformed T2w data as NIfTI file
    data_T2w_DTI = ax.xfm_serial(data_T2w, matrix_DTI_T2w, data_dims_DTI, voxel_dims_DTI, voxel_dims_T2w, interp=1, inverse=True)
    pt.save_data(np.rot90(data_T2w_DTI, k=2, axes=(0, 2)), voxel_dims_DTI, path_T2w_DTI, dtype=None)

    # save transformed T2w atlas labels as NIfTI file
    data_anno, voxel_dims_anno = pt.read_data(path_in_anno)
    data_anno_DTI = ax.xfm_serial(np.rot90(data_anno, k=2, axes=(0, 2)), matrix_DTI_T2w, data_dims_DTI, voxel_dims_DTI, voxel_dims_anno, interp=0, inverse=True)
    pt.save_data(np.rot90(data_anno_DTI, k=2, axes=(0, 2)), voxel_dims_DTI, path_anno_DTI, dtype=None)

    # save transformed T2w stroke mask as NIfTI file
    if os.path.isfile(path_in_mask):
        data_mask, voxel_dims_mask = pt.read_data(path_in_mask)
        data_mask_DTI = ax.xfm_serial(np.rot90(data_mask, k=2, axes=(0, 2)), matrix_DTI_T2w, data_dims_DTI, voxel_dims_DTI, voxel_dims_mask, interp=0, inverse=True)
        pt.save_data(np.rot90(data_mask_DTI, k=2, axes=(0, 2)), voxel_dims_DTI, path_mask_DTI, dtype=None)

    # save transformed T2w peri-infarct mask as NIfTI file
    data_peri, voxel_dims_peri = pt.read_data(path_in_peri)
    data_peri_DTI = ax.xfm_serial(np.rot90(data_peri, k=2, axes=(0, 2)), matrix_DTI_T2w, data_dims_DTI, voxel_dims_DTI, voxel_dims_peri, interp=0, inverse=True)
    pt.save_data(np.rot90(data_peri_DTI, k=2, axes=(0, 2)), voxel_dims_DTI, path_peri_DTI, dtype=None)

def xfm_T2w_DTI_reg(timepoint_P7, timepoint, group, subject, expno_T2w, expno_DTI, procno_T2w, procno_DTI):
    if (expno_T2w is None) or (expno_DTI is None) or (procno_T2w is None) or (procno_DTI is None):
        return

    in_dir = os.path.join(pt.proc_in_dir, timepoint, group, subject, 'T2w')
    dti_dir = os.path.join(pt.proc_in_dir, timepoint, group, subject, 'DTI')
    mask_dir = os.path.join(pt.proc_out_dir, timepoint, group, subject, 'T2w')
    out_dir = os.path.join(pt.proc_out_dir, timepoint, group, subject, 'DTI')

    if not os.path.isdir(in_dir):
        sys.exit("Error: '%s' is not an existing directory." % (in_dir,))

    if not os.path.isdir(dti_dir):
        sys.exit("Error: '%s' is not an existing directory." % (dti_dir,))

    if not os.path.isdir(mask_dir):
        sys.exit("Error: '%s' is not an existing directory." % (mask_dir,))

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # input T2w file
    path_in_T2w = os.path.join(in_dir, subject + 'BiasBet.nii.gz')
    if not os.path.isfile(path_in_T2w):
        sys.exit("Error: '%s' is not a regular file." % (path_in_T2w,))

    # input T2w atlas labels file
    path_in_anno = os.path.join(in_dir, subject + 'BiasBet_AnnorsfMRI.nii.gz')
    if not os.path.isfile(path_in_anno):
        sys.exit("Error: '%s' is not a regular file." % (path_in_anno,))

    # input T2w stroke mask file
    if timepoint == timepoint_P7:
        path_in_mask = os.path.join(in_dir, subject + 'Stroke_mask.nii.gz')
    else:
        path_in_mask = os.path.join(mask_dir, subject + 'Stroke_mask.nii.gz')

    # input T2w peri-infarct mask file
    path_in_peri = os.path.join(mask_dir, subject + '_peri_mask_m3_n15.nii.gz')
    if not os.path.isfile(path_in_peri):
        sys.exit("Error: '%s' is not a regular file." % (path_in_peri,))

    # DTI reference file
    path_ref = os.path.join(dti_dir, subject + 'SmoothMicoBet.nii.gz')
    if not os.path.isfile(path_ref):
        sys.exit("Error: '%s' is not a regular file." % (path_ref,))

    # transformation matrix
    path_xfm = os.path.join(dti_dir, subject + 'SmoothMicoBettransMatrixAff.txt')
    if not os.path.isfile(path_xfm):
        sys.exit("Error: '%s' is not a regular file." % (path_xfm,))

    # output transformed T2w file
    path_T2w_DTI = os.path.join(out_dir, subject + '_T2w_DTI.nii.gz')

    # output transformed atlas labels file
    path_anno_DTI = os.path.join(out_dir, subject + '_T2w_Anno_DTI.nii.gz')

    # output transformed stroke mask file
    path_mask_DTI = os.path.join(out_dir, subject + '_T2w_Stroke_mask_DTI.nii.gz')

    # output transformed peri-infarct mask file
    path_peri_DTI = os.path.join(out_dir, subject + '_T2w_peri_mask_DTI.nii.gz')

    # resample T2w
    command = 'reg_resample -ref %s -flo %s -res %s -trans %s -inter 1' % (path_ref, path_in_T2w, path_T2w_DTI, path_xfm)
    os.system(command)

    # resample atlas labels
    command = 'reg_resample -ref %s -flo %s -res %s -trans %s -inter 0' % (path_ref, path_in_anno, path_anno_DTI, path_xfm)
    os.system(command)

    # resample stroke mask
    command = 'reg_resample -ref %s -flo %s -res %s -trans %s -inter 0' % (path_ref, path_in_mask, path_mask_DTI, path_xfm)
    os.system(command)

    # resample peri-infarct mask
    command = 'reg_resample -ref %s -flo %s -res %s -trans %s -inter 0' % (path_ref, path_in_peri, path_peri_DTI, path_xfm)
    os.system(command)

def main():
    timepoint_P7 = pt.timepoints[1]
    procno = pt.procno
    for index_t, timepoint in enumerate(pt.timepoints):
        for index_g, group in enumerate(pt.groups):
            # raw data directory
            raw_dir = os.path.join(pt.raw_in_dir, group, timepoint)
            for index_s, subject in enumerate(pt.study[index_t][index_g]):
                if subject is not None:
                    expno_T2w = pt.expno_T2w[index_t][index_g][index_s]
                    expno_rsfMRI = pt.expno_rsfMRI[index_t][index_g][index_s]
                    expno_DTI = pt.expno_DTI[index_t][index_g][index_s]
                    xfm_T2w_rsfMRI(raw_dir, timepoint_P7, timepoint, group, subject, expno_T2w, expno_rsfMRI, procno, procno)
                    xfm_T2w_DTI(raw_dir, timepoint_P7, timepoint, group, subject, expno_T2w, expno_DTI, procno, procno)
                    #xfm_T2w_DTI_reg(timepoint_P7, timepoint, group, subject, expno_T2w, expno_DTI, procno, procno)

if __name__ == '__main__':
    main()
