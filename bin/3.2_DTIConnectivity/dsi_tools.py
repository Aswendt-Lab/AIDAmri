"""
Created on 10/08/2017

@author: Niklas Pallast
Neuroimaging & Neuroengineering
Department of Neurology
University Hospital Cologne

Documentation preface, added 23/05/09 by Victor Vera Frazao:
This document is currently in revision for improvement and fixing.
Specifically changes are made to allow compatibility of the pipeline with Ubuntu 18.04 systems 
and Ubuntu 18.04 Docker base images, respectively, as well as adapting to appearent changes of 
DSI-Studio that were applied since the AIDAmri v.1.1 release. As to date the DSI-Studio version 
used is the 2022/08/03 Ubuntu 18.04 release.
All changes and additional documentations within this script carry a signature with the writer's 
initials (e.g. VVF for Victor Vera Frazao) and the date at application, denoted after '//' at 
the end of the comment line. If code segments need clearance the comment line will be prefaced 
by '#?'. Changes are prefaced by '#>' and other comments are prefaced ordinalrily 
by '#'.

Updated (August 2025) by Paul B Camacho
Biomedical Imaging Center
Beckman Institute for Advanced Science & Technology
University of Illinois at Urbana Champaign
Changes:
Expanded options, CLI, N4BiasFieldCorrection support
Compatibility with Ubuntu 22.04, Python 3.10, DSI-Studio release 2025/04/16 (more efficient storage 
using '.fz' and '.sz' file formats)
"""


from __future__ import print_function

import os
import re
import sys
import time
import glob
import nibabel as nii
import numpy as np
import nipype.interfaces.fsl as fsl
import shutil
import subprocess
import pandas as pd

def scaleBy10(input_path, inv):
    data = nii.load(input_path)
    imgTemp = data.get_fdata()
    if inv is False:
        # scale = np.eye(4) * 10
        # scale[3][3] = 1
        # scaledNiiData = nii.Nifti1Image(imgTemp, data.affine * scale)
        scaledNiiData = nii.Nifti1Image(imgTemp, data.affine)
        # overwrite old nifti
        fslPath = os.path.join(os.path.dirname(input_path), 'fslScaleTemp.nii.gz')
        nii.save(scaledNiiData, fslPath)
        return fslPath
    elif inv is True:
        # scale = np.eye(4) / 10
        # scale[3][3] = 1
        # unscaledNiiData = nii.Nifti1Image(imgTemp, data.affine * scale)
        unscaledNiiData = nii.Nifti1Image(imgTemp, data.affine)
        hdrOut = unscaledNiiData.header
        hdrOut.set_xyzt_units('mm')

        nii.save(unscaledNiiData, input_path)
        return input_path
    else:
        sys.exit("Error: inv - parameter should be a boolean.")


def findSlicesData(path, pre):
    regMR_list = []
    fileALL = glob.iglob(path + '/' + pre + '*.nii.gz', recursive=True)
    for filename in fileALL:
        regMR_list.append(filename)
    regMR_list.sort()
    return regMR_list


def fsl_SeparateSliceMoCo(input_file, par_folder):
    # scale Nifti data by factor 10
    dataName = os.path.basename(input_file).split('.')[0]
    fslPath = scaleBy10(input_file, inv=False)

    aidamri_dir = os.getcwd()
    temp_dir = os.path.join(os.path.dirname(input_file), "temp")
    if not os.path.exists(temp_dir):
        os.mkdir(temp_dir)

    os.chdir(temp_dir)
    mySplit = fsl.Split(in_file=fslPath, dimension='z', out_base_name=dataName)
    mySplit.run()
    os.remove(fslPath)

    # separate ref and src volume in slices
    sliceFiles = findSlicesData(os.getcwd(), dataName)

    # start to correct motions slice by slice
    for i in range(len(sliceFiles)):
        slc = sliceFiles[i]
        output_file = os.path.join(par_folder, os.path.basename(slc))
        myMCFLIRT = fsl.preprocess.MCFLIRT(in_file=slc, out_file=output_file, save_plots=True, terminal_output='none')
        myMCFLIRT.run()
        os.remove(slc)

    # merge slices to a single volume
    mcf_sliceFiles = findSlicesData(par_folder, dataName)
    output_file = os.path.join(os.path.dirname(input_file),
                               os.path.basename(input_file).split('.')[0]) + '_mcf.nii.gz'
    myMerge = fsl.Merge(in_files=mcf_sliceFiles, dimension='z', merged_file=output_file)
    myMerge.run()

    for slc in mcf_sliceFiles: 
        os.remove(slc)

    # unscale result data by factor 10**(-1)
    output_file = scaleBy10(output_file, inv=True)
    
    os.chdir(aidamri_dir)

    return output_file


def fsl_eddy_correct(input_file, outputPath):
    """
    Correct eddy currents in DWI data using FSL's eddy command.
    """
    mask_file = os.path.join(os.path.dirname(input_file), os.path.basename(input_file).split('.')[0] + 'Bet_mask.nii.gz')
    index = os.path.join(os.path.dirname(input_file), os.path.basename(input_file).split('.')[0] + 'index.txt')
    acqp = os.path.join(os.path.dirname(input_file), os.path.basename(input_file).split('.')[0] + 'acqp.txt')
    bvec = os.path.join(os.path.dirname(input_file), os.path.basename(input_file).split('.')[0] + '.bvec')
    bval = os.path.join(os.path.dirname(input_file), os.path.basename(input_file).split('.')[0] + '.bval')
    myEddy = fsl.preprocess.Eddy(in_file=input_file,
                                  out_file=outputPath,
                                  in_mask=mask_file,
                                  in_index=index,
                                  in_acqp=acqp,
                                  in_bvec=bvec,
                                  in_bval=bval)
    myEddy.run()
    print("Eddy current correction completed")
    return outputPath


def fsl_topup(input_file, outputPath):
    mask_file = os.path.join(os.path.dirname(input_file), os.path.basename(input_file).split('.')[0] + 'Bet_mask.nii.gz')
    myTopup = fsl.preprocess.Topup(in_file=input_file, out_file=outputPath, in_mask=mask_file)
    myTopup.run()
    print("Topup completed")
    return outputPath


def make_dir(dir_out, dir_sub):
    """
    Creates new directory.
    """
    dir_out = os.path.normpath(os.path.join(dir_out, dir_sub))
    if not os.path.exists(dir_out):
        os.mkdir(dir_out)
        time.sleep(1.0)
        if not os.path.exists(dir_out):
            sys.exit("Could not create directory \"%s\"" % (dir_out,))
    return dir_out

def move_files(dir_in, dir_out, pattern):
    time.sleep(1.0)
    file_list = glob.glob(dir_in+pattern)
    file_list.sort()

    time.sleep(1.0)
    for file_mv in file_list: # move files from input to output directory
        file_in = os.path.join(dir_in, file_mv)
        shutil.copy(file_in, dir_out)

    for file_mv in file_list: # remove files in output directory
        file_out = os.path.join(dir_out, file_mv)
        if os.path.isfile(file_out):
            os.remove(file_out)


def erode_mask(input_file, outputPath, n_voxels=1):
    """
    Erodes the mask by n voxels (default = 1).
    """
    myErode = fsl.preprocess.ErodeImage(in_file=input_file, out_file=outputPath, kernel_shape='sphere', kernel_size=n_voxels)
    myErode.run()
    print("Mask erosion completed")
    return outputPath


def get_min_voxel_size_mm(nifti_path):
    nii_file = nii.load(nifti_path)
    header = nii_file.header
    mm_voxel_size_arr = header.get_zooms()
    spatial_dims = mm_voxel_size_arr[:3]
    min_vox_size_mm = str(min(spatial_dims))
    return min_vox_size_mm


def strip_nifti_suffix(path):
    """
    Return a filename stem without .nii or .nii.gz so DSI output patterns can
    be matched and renamed reliably.
    """
    name = os.path.basename(path)
    if name.endswith('.nii.gz'):
        return name[:-7]
    if name.endswith('.nii'):
        return name[:-4]
    return os.path.splitext(name)[0]


def find_matching_gradient_pair(folder_path, preferred_stem=None):
    """
    Find the .bval/.bvec pair that belongs to a DWI series.

    The exact DWI stem is preferred when it is known. Otherwise we fall back to
    the common *_dwi naming convention or to the only remaining unique pair.
    """
    bval_files = sorted(glob.glob(os.path.join(folder_path, '*.bval')))
    bvec_files = sorted(glob.glob(os.path.join(folder_path, '*.bvec')))

    bval_map = {
        os.path.basename(path).replace('.bval', ''): path
        for path in bval_files
    }
    bvec_map = {
        os.path.basename(path).replace('.bvec', ''): path
        for path in bvec_files
    }

    common_stems = sorted(set(bval_map) & set(bvec_map))
    if not common_stems:
        return None, "Both bval and bvec files must be present in the folder."

    if preferred_stem:
        preferred_stem = strip_nifti_suffix(os.path.basename(preferred_stem))
        if preferred_stem in common_stems:
            return {
                'stem': preferred_stem,
                'bval': bval_map[preferred_stem],
                'bvec': bvec_map[preferred_stem],
            }, None

    preferred_dwi_stems = [stem for stem in common_stems if stem.endswith('_dwi')]
    if len(preferred_dwi_stems) == 1:
        stem = preferred_dwi_stems[0]
        return {
            'stem': stem,
            'bval': bval_map[stem],
            'bvec': bvec_map[stem],
        }, None

    if len(common_stems) == 1:
        stem = common_stems[0]
        return {
            'stem': stem,
            'bval': bval_map[stem],
            'bvec': bvec_map[stem],
        }, None

    return None, (
        "Multiple matching bval/bvec pairs were found. "
        "Please keep only one pair in the folder."
    )


def connectivity(dsi_studio, dir_in, dir_seeds, dir_out, dir_con, make_isotropic=0, flip_image_y=False, legacy=False):
    """
    Calculates connectivity data (types: pass and end).
    """
    if not os.path.exists(dir_in):
        sys.exit("Input directory \"%s\" does not exist." % (dir_in,))

    dir_seeds = os.path.normpath(os.path.join(dir_in, dir_seeds))
    if not os.path.exists(dir_seeds):
        sys.exit("Seeds directory \"%s\" does not exist." % (dir_seeds,))

    if not os.path.exists(dir_out):
        sys.exit("Output directory \"%s\" does not exist." % (dir_out,))

    dir_con = make_dir(dir_out, dir_con)

    # change to input directory
    os.chdir(os.path.dirname(dir_in))
    cmd_ana = r'%s --action=%s --source=%s --tract=%s --connectivity=%s --connectivity_value=%s --connectivity_type=%s'

    if legacy:
        filename = glob.glob(dir_in+'/*.fib.gz')[0]
    else:
        filename = glob.glob(dir_in + '/*.fz')[0]
    file_trk = glob.glob(dir_in+'/*trk.gz')[0]
    file_seeds = dir_seeds

    # Dev note: if we resample the diffusion image, we need to resample the file_seeds here
    iso_value = None
    if make_isotropic == "auto":
        # Match srcgen(): when "auto" is requested, use the native voxel size
        # of the ROI/seed image for the optional connectivity resampling step.
        iso_value = float(get_min_voxel_size_mm(file_seeds))
    else:
        iso_value = float(make_isotropic)
    if iso_value is not None and iso_value > 0:
        # resample seeds image to isotropic voxel size using AFNI 3dresample
        resampled_seeds_path = os.path.join(dir_con, os.path.basename(file_seeds).replace('.nii', '_resampled.nii.gz'))

        cmd_resample = f"flirt -in {file_seeds} -ref {file_seeds} -applyisoxfm {iso_value} -nosearch -interp trilinear -out {resampled_seeds_path}"
        print(f'Resampling seeds image to {iso_value} mm isotropic voxel size')
        subprocess.run(cmd_resample, shell=True, check=True)
        # create a quality control image for the resampling, with the original image and the resampled image side by side
        # Generate PNG images for original and resampled seeds
        qc_orig_png = os.path.join(dir_con, 'qc_seeds_orig.png')
        qc_resampled_png = os.path.join(dir_con, 'qc_seeds_resampled.png')
        subprocess.run(f"slicer {file_seeds} -L -a {qc_orig_png}", shell=True, check=True)
        subprocess.run(f"slicer {resampled_seeds_path} -L -a {qc_resampled_png}", shell=True, check=True)
        # Combine the two PNGs using pngappend
        qc_combined_png = os.path.join(dir_con, 'qc_resampled_seeds_combined.png')
        cmd_qc = f"pngappend {qc_orig_png} - {qc_resampled_png} {qc_combined_png}"
        print(f'Creating quality control image for resampled seeds image')
        subprocess.run(cmd_qc, shell=True, check=True)
        # update file_seeds to the resampled path
        file_seeds = resampled_seeds_path
        # # command needs to use the --t1t2 t2_rare.nii.gz to align the ROI files
        # cmd_ana = r'%s --action=%s --source=%s --tract=%s --connectivity=%s --connectivity_value=%s --connectivity_type=%s --t1t2=%s'
        # # Inverse scale the file_seeds by 10
        # # file_seeds = scaleBy10(file_seeds, inv=True)
    # Performs analysis on every connectivity value within the list ('qa' may not be necessary; might be removed in the future.)
    connect_vals = ['qa', 'count']
    for i in connect_vals:
        parameters = (dsi_studio, 'ana', filename, file_trk, file_seeds, i, 'pass,end')
        # if make_isotropic != 0:
        #     parameters += (t2rare)
        os.system(cmd_ana % parameters)

        # DSI Studio reuses the same default connectivity filename regardless of
        # connectivity_value. Rename each result immediately so qa/count do not
        # overwrite each other before the files are moved to dir_con.
        tract_dir = os.path.dirname(file_trk)
        tract_base = os.path.basename(file_trk)
        roi_base = strip_nifti_suffix(file_seeds)
        connectivity_outputs = sorted(glob.glob(os.path.join(
            tract_dir, f'{tract_base}.{roi_base}.connectivity.*'
        )))
        for output_path in connectivity_outputs:
            output_dir = os.path.dirname(output_path)
            output_name = os.path.basename(output_path)
            if f'.{i}.connectivity.' in output_name:
                continue
            renamed_name = output_name.replace('.connectivity.', f'.{i}.connectivity.', 1)
            os.replace(output_path, os.path.join(output_dir, renamed_name))

    #move_files(dir_in, dir_con, re.escape(filename) + '\.' + re.escape(pre_seeds) + '.*(?:\.pass\.|\.end\.)')
    move_files(os.path.dirname(file_trk), dir_con, '/*.txt')
    move_files(os.path.dirname(file_trk), dir_con, '/*.mat')

def mapsgen(dsi_studio, dir_in, dir_msk, b_table, pattern_in, pattern_fib):
    """
    FUNCTION DEPRECATED. REMOVAL PENDING.
    """
    pre_msk = 'bet.bin.'

    ext_src = '.src.gz'
    ext_nii = '.nii.gz'

    if not os.path.exists(dir_in):
        sys.exit("Input directory \"%s\" does not exist." % (dir_in,))

    dir_msk = os.path.normpath(os.path.join(dir_in, dir_msk))
    if not os.path.exists(dir_msk):
        sys.exit("Masks directory \"%s\" does not exist." % (dir_msk,))

    b_table = os.path.join(dir_in, b_table)
    if not os.path.isfile(b_table):
        sys.exit("File \"%s\" does not exist." % (b_table,))

    # change to input directory
    os.chdir(dir_in)

    cmd_src = r'%s --action=%s --source=%s --output=%s --b_table=%s'
    # method: 0:DSI, 1:DTI, 4:GQI 7:QSDR, param0: 1.25 (in vivo) diffusion sampling lenth ratio for GQI and QSDR reconstruction, --thread_count: number of multi-threads used to conduct reconstruction 
    cmd_rec = r'%s --action=%s --source=%s --mask=%s --method=%d --param0=%s --thread_count=%d --check_btable=%d'

    file_list = [x for x in os.listdir(dir_in) if os.path.isfile(os.path.join(dir_in, x)) and re.match(pattern_in, x)]
    file_list.sort()

    for index, filename in enumerate(file_list):
        # create source files
        pos = filename.rfind('_')

        file_src = filename[:pos] + ext_src
        parameters = (dsi_studio, 'src', filename, file_src, b_table)
        subprocess.call(cmd_src % parameters)

        # create fib files
        file_msk = os.path.join(dir_msk, pre_msk + filename[:pos] + ext_nii)
        parameters = (dsi_studio, 'rec', file_src, file_msk, 3, '1.25', 2, 0)
        subprocess.call(cmd_rec % parameters)

    # extracts maps: 2 ways:
    cmd_exp = r'%s --action=%s --source=%s --export=%s'

    file_list = [x for x in os.listdir(dir_in) if os.path.isfile(os.path.join(dir_in, x)) and re.match(pattern_fib, x)]
    file_list.sort()

    for index, filename in enumerate(file_list):
        #file_fib = os.path.join(dir_in, filename)
        #parameters = (dsi_studio, 'exp', file_fib, 'fa')
        parameters = (dsi_studio, 'exp', filename, 'fa')
        print("%d of %d:" % (index + 1, len(file_list)), cmd_exp % parameters)
        subprocess.call(cmd_exp % parameters)

def srcgen(dsi_studio, dir_in, dir_msk, dir_out, b_table, recon_method='dti', vivo='in_vivo', make_isotropic=0, flip_image_y=False, template=6, legacy=False, gradient_pair=None):
    """
    Sources and creates fib files. Diffusivity and anisotropy metrics are exported from data.
    """
    dir_src = r'src'
    dir_fib = r'fib_map'
    dir_qa  = r'DSI_studio'
    dir_con = r'connectivity'
    # Support for backwards compatibility with pre-2024 DSI Studio (AIDAmri <= v2.0)
    if legacy:
        ext_src = '.src.gz'
        ext_fib = '.fib.gz'
    else:
        ext_src = '.sz'
        ext_fib = '.fz'

    if not os.path.exists(dir_in):
        sys.exit("Input directory \"%s\" does not exist." % (dir_in,))

    dir_msk = os.path.normpath(os.path.join(dir_in, dir_msk))
    if not os.path.exists(dir_msk):
        sys.exit("Masks directory \"%s\" does not exist." % (dir_msk,))

    if not os.path.exists(dir_out):
        sys.exit("Output directory \"%s\" does not exist." % (dir_out,))

    input_dir = os.path.dirname(os.path.abspath(dir_in))

    if gradient_pair is None:
        if not os.path.isabs(b_table):
            b_table = os.path.join(input_dir, b_table)
        if not os.path.isfile(b_table):
            sys.exit("File \"%s\" does not exist." % (b_table,))
    else:
        if not os.path.isfile(gradient_pair['bval']):
            sys.exit("File \"%s\" does not exist." % (gradient_pair['bval'],))
        if not os.path.isfile(gradient_pair['bvec']):
            sys.exit("File \"%s\" does not exist." % (gradient_pair['bvec'],))

    dir_src = make_dir(os.path.dirname(dir_out), dir_src)
    dir_fib = make_dir(os.path.dirname(dir_out), dir_fib)
    dir_qa  = make_dir(os.path.dirname(dir_out), dir_qa)

    # change to input directory
    os.chdir(input_dir)

    # Prefer explicit gradient files when they are available. This keeps the
    # source image and the original gradients coupled even after motion
    # correction changed the NIfTI filename to *_mcf.nii.gz.
    cmd_src_btable = r'%s --action=%s --source=%s --output=%s --b_table=%s'
    cmd_src_gradients = r'%s --action=%s --source=%s --output=%s --bval=%s --bvec=%s'
    # method: 0:DSI, 1:DTI, 4:GQI 7:QSDR, param0: 1.25 (in vivo) diffusion sampling lenth ratio for GQI and QSDR reconstruction, 
    # check_btable: Set –check_btable=1 to test b-table orientation and apply automatic flippin, thread_count: number of multi-threads used to conduct reconstruction
    cmd_rec = r'%s --action=%s --source=%s --mask=%s --method=%d --other_output=all --output=%s --check_btable=%d --cmd=%s' # Dev note: if not using slice-wise motion correction, --motion_correction 1

    # create source files
    filename = os.path.basename(dir_in)
    # pos = filename.rfind('.')
    # file_src = os.path.join(dir_src, filename[:pos] + ext_src)
    filename_base = filename.split('.')[0]
    file_src = os.path.join(dir_src, filename_base + ext_src)
    if gradient_pair is not None:
        parameters = (
            dsi_studio,
            'src',
            filename,
            file_src,
            gradient_pair['bval'],
            gradient_pair['bvec'],
        )
        os.system(cmd_src_gradients % parameters)
    else:
        parameters = (dsi_studio, 'src', filename, file_src, b_table)
        os.system(cmd_src_btable % parameters)

    patterns = [
        os.path.join(dir_src, filename_base + '*.src.gz.sz'),
        os.path.join(dir_src, filename_base + '*.sz'),
        os.path.join(dir_src, filename_base + '*.src.gz'),
    ]

    src_candidates = []
    for pattern in patterns:
        src_candidates.extend(sorted(glob.glob(pattern)))

    src_candidates = list(dict.fromkeys(src_candidates))

    if not src_candidates:
        raise FileNotFoundError(
            f"No DSI Studio source file found for base '{filename_base}' in {dir_src}"
        )

    file_src_real = src_candidates[0]

    # If unrealistic streamlines cross top of cortex are present due to an oversized mask, erode mask
    mask_erosion = 0
    if mask_erosion > 0:
        # Erode mask by 1 voxel
        dir_msk = erode_mask(dir_msk, os.path.join(dir_msk, 'eroded_mask.nii.gz'), n_voxels=mask_erosion)
        print(f'Eroded mask saved to {dir_msk}')

    # create fib files
    file_msk = dir_msk

    param_zero='1.25'
    # Select reconstruction parameters
    if vivo == "ex_vivo":
        param_zero='0.60'
        print(f'Using param0 value {param_zero} recommended for ex vivo data')
    elif vivo == "in_vivo":
        print(f'Using param0 value {param_zero} recommended for in vivo data')

    # get voxel size from nifti header to resample if make_isotropic == 'auto', use function get_min_voxel_size_mm
    min_vox_size_mm = get_min_voxel_size_mm(filename)
    if make_isotropic == "auto":
        iso_value = float(min_vox_size_mm)
    else:
        iso_value = float(make_isotropic)
    # Dev note: The parcellation image must be registered to the resampled (and flipped) diffusion image if this is the case
    #           We need to register the T2-weighted image to the resampled diffusion image to guide the registration of the parcellation image
    #           This should be done using the equivalent commands from registration_DTI.py

    additional_cmd = ''
    if iso_value > 0:
        additional_cmd = f'[Step T2][Edit][Resample]={iso_value}'
        print(f'Resampling to {iso_value} mm isotropic voxel size')
    
    use_eddy_correct = False
    use_dsi_topup = False
    rev_pe_image = '' # path to reverse phase encoding image if using dsi topup
    if use_eddy_correct and use_dsi_topup:
        additional_cmd = f'[Step T2][Corrections][TOPUP EDDY]={rev_pe_image}+{additional_cmd}'
    elif use_eddy_correct and not use_dsi_topup:
        additional_cmd = f'[Step T2][Corrections][EDDY]+{additional_cmd}'
    if use_dsi_topup:
        additional_cmd = f'[Step T2][Edit][TOPUP]={rev_pe_image}+{additional_cmd}'

    # default method value for DTI 
    method_rec=1
    if recon_method == "gqi":
        method_rec=4

    # Select template for DSI Studio
    # 1: C57BL6_mouse, 5: WHS_SD_rat 
    if template == "Mouse":
        template = 6
    elif template == "Rat":
        template = 5

    file_fib = os.path.join(dir_fib, filename_base + ext_fib)

    if additional_cmd:
        parameters = (dsi_studio, 'rec', file_src_real, file_msk, method_rec, file_fib, 0, '"'f'{additional_cmd}"')
        os.system(cmd_rec % parameters)
    else:
        cmd_rec_no_cmd = r'%s --action=%s --source=%s --mask=%s --method=%d --other_output=all --output=%s --check_btable=%d'
        parameters = (dsi_studio, 'rec', file_src_real, file_msk, method_rec, file_fib, 0)
        os.system(cmd_rec_no_cmd % parameters)

    # move fib to corresponding folders
    move_files(dir_src, dir_fib, f'/*{ext_fib}')

    # extracts maps: 2 ways:
    cmd_exp = r'%s --action=%s --source=%s --export=%s'
    file_fib = glob.glob(dir_fib+f'/*{ext_fib}')[0]
    if recon_method == "gqi":
        parameters = (dsi_studio, 'exp', file_fib, 'dti_fa')
    elif recon_method == "dti":
        parameters = (dsi_studio, 'exp', file_fib, 'fa')
    os.system(cmd_exp % parameters)

    # extracts maps: 2 ways:
    cmd_exp = r'%s --action=%s --source=%s --export=%s'
    file_fib = glob.glob(dir_fib + f'/*{ext_fib}')[0]
    parameters = (dsi_studio, 'exp', file_fib, 'md')
    os.system(cmd_exp % parameters)

    # extracts maps: 2 ways:
    cmd_exp = r'%s --action=%s --source=%s --export=%s'
    file_fib = glob.glob(dir_fib + f'/*{ext_fib}')[0]
    parameters = (dsi_studio, 'exp', file_fib, 'ad')
    os.system(cmd_exp % parameters)

    # extracts maps: 2 ways:
    cmd_exp = r'%s --action=%s --source=%s --export=%s'
    file_fib = glob.glob(dir_fib + f'/*{ext_fib}')[0]
    parameters = (dsi_studio, 'exp', file_fib, 'rd')
    os.system(cmd_exp % parameters)

    move_files(dir_fib, dir_qa, '/*qa.nii.gz')
    move_files(dir_fib, dir_qa, '/*fa.nii.gz')
    move_files(dir_fib, dir_qa, '/*md.nii.gz')
    move_files(dir_fib, dir_qa, '/*ad.nii.gz')
    move_files(dir_fib, dir_qa, '/*rd.nii.gz')

    # Generate PNG images for each NIfTI file using FSL's slicer tool
    # Collect all NIfTI files in dir_qa, including .nii and .nii.gz
    all_nifti_files = glob.glob(os.path.join(dir_qa, "*.nii")) + glob.glob(os.path.join(dir_qa, "*.nii.gz"))

    # Remove duplicates and ensure unique file paths
    all_nifti_files = list(set(all_nifti_files))

    nifti_files = all_nifti_files

    for nifti_path in nifti_files:
        base_name = os.path.splitext(os.path.basename(nifti_path))[0]
        # Remove .gz if present
        if base_name.endswith('.nii'):
            base_name = base_name[:-4]
        img = nii.load(nifti_path)
        png_slice_path = os.path.join(dir_qa, f"{base_name}.png")
        cmd = ["slicer", nifti_path, "-L -a", png_slice_path]
        subprocess.run(cmd, check=True)

    return min_vox_size_mm # , flip_image_y

def tracking(dsi_studio, dir_in, track_param='default', min_voxel_size_mm=0.1, thread_count=1, legacy=False):
    """
    Performs seed-based fiber-tracking.
    Default parameters are used unless a custom parameter is specified.
    """
    if not os.path.exists(dir_in):
        sys.exit("Input directory \"%s\" does not exist." % (dir_in,))

    if legacy:
        ext_fib = '.fib.gz'
    else:
        ext_fib = '.fz'

    # Define parameter sets
    param_sets = {
        'default':        ['0AD7A33C9A99193FE8D5123F0AD7233CCDCCCC3D9A99993EbF04240420FdcaCDCC4C3Ec'],
        'aida_optimized': [1000000, 0, '.01', '55', 0, '.02', '.1', '.3', '120.0'],
        'rat':            [1000000, 0, '.01', '60', 0, '.02', '.1', '.3', '20.0'],
        'mouse':          [1000000, 0, '.01', '45', 0, '.02', '.1', '.3', '15.0'],
    }

    if isinstance(track_param, str):
        params = param_sets.get(track_param.lower())
        if params is None:
            sys.exit(f'Unknown track_param set: {track_param}')
    elif isinstance(track_param, list):
        params = track_param
    else:
        sys.exit('track_param must be "default", "aida_optimized", "rat", "mouse", or a list of parameter values for fiber_count, interpolation, step_size, turning_angle, check_ending, fa_threshold, smoothing, min_length, max_length.')

    # change to input directory
    os.chdir(os.path.dirname(dir_in))

    filename = glob.glob(dir_in + f'/*{ext_fib}')[0]

    # Set tracking based on track_param:
    if track_param == 'default':
        print('Using DSI Studio default tracking parameters')
        # Use this tracking parameters in the form of parameter_id that you can get directly from the dsi_studio gui console. (this is here now the defualt mode)
        cmd_trk = r'%s --action=%s --source=%s --output=%s --parameter_id=%s'
        # Use this tracking parameters in the form of parameter_id that you can get directly from the dsi_studio gui console. (this is here now the defualt mode)
        parameters = (dsi_studio, 'trk', filename, os.path.join(dir_in, filename+'.trk.gz'), '0AD7A33C9A99193FE8D5123F0AD7233CCDCCCC3D9A99993EbF04240420FdcaCDCC4C3Ec')
    else:
        # Use this tracking parameters if you want to specify each tracking parameter separately.
        # cmd_trk = r'%s --action=%s --source=%s --output=%s --fiber_count=%d --interpolation=%d --step_size=%s --turning_angle=%s --check_ending=%d --fa_threshold=%s --smoothing=%s --min_length=%s --max_length=%s --thread_count=%s'
        cmd_trk = r'%s --action=%s --source=%s --output=%s --fiber_count=%d --interpolation=%d --step_size=%s --turning_angle=%s --check_ending=%d --fa_threshold=%s --smoothing=%s --min_length=%s --max_length=%s --thread_count=%s --export=tdi:color' # The tract-density image saved here may not be viewable in FSLeyes or ITK-SNAP, but is compatible with Mango viewer. 
        #parameters = (dsi_studio, 'trk', filename, os.path.join(dir_in, filename+'.trk.gz'), 1000000, 0, '.5', '55', 0, '.02', '.1', '.5', '12.0') #Our Old parameters
        #parameters = (dsi_studio, 'trk', filename, os.path.join(dir_in, filename+'.trk.gz'), 1000000, 0, '.01', '55', 0, '.02', '.1', '.3', '120.0') #Here are the optimized parameters (fatemeh)
        if track_param != 'aida_optimized':
            # step size = 1/2 (voxel size)
            params[2] = min_voxel_size_mm / 2
            # min streamline length = 2 * (voxel_size)
            params[7] = min_voxel_size_mm * 2
        parameters = (dsi_studio, 'trk', filename, os.path.join(dir_in, filename+'.trk.gz'), *params, thread_count)
    
    os.system(cmd_trk % parameters)

def merge_bval_bvec_to_btable(folder_path, preferred_stem=None):
    gradient_pair, error_message = find_matching_gradient_pair(folder_path, preferred_stem=preferred_stem)
    if gradient_pair is None:
        print(error_message)
        return False

    file_name = gradient_pair['stem']
    bval_file = gradient_pair['bval']
    bvec_file = gradient_pair['bvec']

    try:
        with open(bval_file, 'r') as bval_handle:
            bval_contents = bval_handle.read()
            # Split the content into a list of values (assuming it's space-separated)
            bval_values = bval_contents.strip().split()
            # Convert the list to a Pandas DataFrame and cast the 'bval' column to integers
            bval_table = pd.DataFrame({'bval': bval_values}).astype(float)

        with open(bvec_file, 'r') as bvec_handle:
            # Read lines and split each line into values
            bvec_lines = bvec_handle.readlines()
            bvec_values = [line.strip().split() for line in bvec_lines]

            # Create a Pandas DataFrame from the values
            bvec_table = pd.DataFrame(bvec_values, columns=[f'bvec_{i+1}' for i in range(len(bvec_values[0]))])
            # Transpose the bvec_table
            bvec_table = bvec_table.T

        # Merge bval_table and bvec_table
        merged_table = np.hstack((bval_table, bvec_table))
        # Convert the merged_table content to float
        merged_table = merged_table.astype(float)
        # Define the path for the final merged table
        final_path = os.path.join(folder_path, file_name + "_btable.txt")

        # Save the merged table to the final file
        np.savetxt(final_path, merged_table, fmt='%f', delimiter='\t')
        print(f"Merged table saved to {final_path}")
        return final_path
    except FileNotFoundError:
        print("One or both of the bval and bvec files were not found.")
        return False
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        return False

if __name__ == '__main__':
    pass
