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
import nibabel as nib
import numpy as np
import nipype.interfaces.fsl as fsl
import shutil
import subprocess
import pandas as pd

def scaleBy10(input_path, inv):
    data = nib.load(input_path)
    imgTemp = data.get_fdata()
    if inv is False:
        # create 4x4 scaling matrix and scale by 10 to match human like brain size
        scale = np.eye(4)
        scale[0, 0] = 10
        scale[1, 1] = 10
        scale[2, 2] = 10

        # Create new Nifti image with scaled affine
        scaled_affine = data.affine @ scale

        scaledNiiData = nib.Nifti1Image(imgTemp, scaled_affine)
        # overwrite old nifti
        fslPath = os.path.join(os.path.dirname(input_path), 'fslScaleTemp.nii.gz')
        nib.save(scaledNiiData, fslPath)
        return fslPath
    elif inv is True:
        #rescale nifti
        inv_scale = np.eye(4)
        inv_scale[0, 0] = 0.1
        inv_scale[1, 1] = 0.1
        inv_scale[2, 2] = 0.1

        unscaled_affine = data.affine @ inv_scale
        unscaledNiiData = nib.Nifti1Image(imgTemp, unscaled_affine)
        hdrOut = unscaledNiiData.header
        hdrOut.set_xyzt_units('mm')

        nib.save(unscaledNiiData, input_path)
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

def infer_slice_axis(nifti_path):
    img = nib.load(nifti_path)
    shape = img.shape[:3]
    zooms = img.header.get_zooms()[:3]

    # First, search for the smallest matrix dimension.
    # For typical 2D multislice DWI data, this is the slice axis.
    axis_index = min(range(3), key=lambda i: shape[i])
    axis_name = ["x", "y", "z"][axis_index]

    #print(f"Image shape: {shape}")
    #print(f"Voxel sizes: {zooms}")
    #print(f"Inferred slice axis for slice-wise MoCo: {axis_name}")

    return axis_name

def fsl_SeparateSliceMoCo(input_file, par_folder):
    # scale Nifti data by factor 10
    dataName = os.path.basename(input_file).split('.')[0]
    fslPath = scaleBy10(input_file, inv=False)

    aidamri_dir = os.getcwd()
    temp_dir = os.path.join(os.path.dirname(input_file), "temp")
    if not os.path.exists(temp_dir):
        os.mkdir(temp_dir)

    os.chdir(temp_dir)

    #find slice axis
    slice_axis = infer_slice_axis(fslPath)
    mySplit = fsl.Split(in_file=fslPath, dimension=slice_axis, out_base_name=dataName)
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
    myMerge = fsl.Merge(in_files=mcf_sliceFiles, dimension=slice_axis, merged_file=output_file)
    myMerge.run()

    for slc in mcf_sliceFiles: 
        os.remove(slc)

    # unscale result data by factor 10**(-1)
    output_file = scaleBy10(output_file, inv=True)
    
    os.chdir(aidamri_dir)
    if os.path.isdir(temp_dir):
        shutil.rmtree(temp_dir)

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
    nii_file = nib.load(nifti_path)
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


def reorient_nifti_to_lip(nifti_path):
    """
    Reorder a NIfTI image to LIP orientation in-place.
    """
    img = nib.load(nifti_path)
    current = nib.orientations.io_orientation(img.affine)
    target = nib.orientations.axcodes2ornt(("L", "I", "P"))

    if np.array_equal(current, target):
        return nifti_path

    transform = nib.orientations.ornt_transform(current, target)
    data = nib.orientations.apply_orientation(np.asanyarray(img.dataobj), transform)
    new_affine = img.affine @ nib.orientations.inv_ornt_aff(transform, img.shape)

    out_img = nib.Nifti1Image(data, new_affine, img.header)
    out_img.set_qform(new_affine, code=1)
    out_img.set_sform(new_affine, code=1)
    nib.save(out_img, nifti_path)
    return nifti_path


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

def srcgen(dsi_studio, dir_in, dir_msk, dir_out, b_table, recon_method='dti', vivo='in_vivo', make_isotropic=0, legacy=False, gradient_pair=None):
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

    input_dir = os.path.dirname(os.path.abspath(dir_in))

    if not os.path.isabs(dir_msk):
        dir_msk = os.path.join(input_dir, dir_msk)
    dir_msk = os.path.normpath(dir_msk)

    if not os.path.isfile(dir_msk):
        sys.exit("Mask file \"%s\" does not exist." % (dir_msk,))

    if not os.path.exists(dir_out):
        sys.exit("Output path \"%s\" does not exist." % (dir_out,))


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
    # method: 0:DSI, 1:DTI, 4:GQI 7:QSDR, param0: 1.25 (in vivo) diffusion sampling lenth ratio for GQI and QSDR reconstruction, 
    # check_btable: Set –check_btable=1 to test b-table orientation and apply automatic flippin, thread_count: number of multi-threads used to conduct reconstruction

    # create source files
    filename = os.path.basename(dir_in)
    filename_base = strip_nifti_suffix(filename)
    file_src = os.path.join(dir_src, filename_base + ext_src)
    if gradient_pair is not None:
        cmd = [
            dsi_studio,
            "--action=src",
            f"--source={filename}",
            f"--output={file_src}",
            f"--bval={gradient_pair['bval']}",
            f"--bvec={gradient_pair['bvec']}",
        ]
    else:
        cmd = [
            dsi_studio,
            "--action=src",
            f"--source={filename}",
            f"--output={file_src}",
            f"--b_table={b_table}",
        ]

    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)

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
        eroded_mask_path = os.path.join(
            os.path.dirname(dir_msk),
            strip_nifti_suffix(dir_msk) + "_eroded_mask.nii.gz"
        )

        dir_msk = erode_mask(dir_msk, eroded_mask_path, n_voxels=mask_erosion)
        print(f"Eroded mask saved to {dir_msk}")

    # create fib files
    file_msk = dir_msk

    #param0 still used in dsi_studio 2025?
    param_zero='1.25'
    # Select reconstruction parameters
    if vivo == "ex_vivo":
        param_zero='0.60'
        print(f'Using param0 value {param_zero} recommended for ex vivo data')
    elif vivo == "in_vivo":
        print(f'Using param0 value {param_zero} recommended for in vivo data')

    # get voxel size from nifti header to resample if make_isotropic == 'auto', use function get_min_voxel_size_mm
    min_vox_size_mm = get_min_voxel_size_mm(dir_in)

    if str(make_isotropic).lower() == "auto":
        iso_value = float(min_vox_size_mm)
    else:
        try:
            iso_value = float(make_isotropic)
        except ValueError:
            sys.exit(
                f"Invalid make_isotropic value: {make_isotropic}. "
                'Use 0, "auto", or a voxel size in mm, e.g. 0.2.'
            )
    # Dev note: The parcellation image must be registered to the resampled (and flipped) diffusion image if this is the case
    #           We need to register the T2-weighted image to the resampled diffusion image to guide the registration of the parcellation image
    #           This should be done using the equivalent commands from registration_DTI.py

    additional_cmd = ''
    #Atlas has to be in the same space resampled space if iso_value above 0
    if iso_value > 0:
        additional_cmd = f'[Step T2][Edit][Resample]={iso_value}'
        print(f'Resampling to {iso_value} mm isotropic voxel size')

    # Optional future DSI Studio corrections. Currently disabled because AIDAmri
    # performs slice-wise motion correction before source generation.
    use_eddy_correct = False
    use_dsi_topup = False
    rev_pe_image = '' # path to reverse phase encoding image if using dsi topup
    if use_eddy_correct and use_dsi_topup:
        additional_cmd = f'[Step T2][Corrections][TOPUP EDDY]={rev_pe_image}+{additional_cmd}'
    elif use_eddy_correct and not use_dsi_topup:
        additional_cmd = f'[Step T2][Corrections][EDDY]+{additional_cmd}'
    if use_dsi_topup:
        additional_cmd = f'[Step T2][Edit][TOPUP]={rev_pe_image}+{additional_cmd}'

    # method for reconstruction algorithm
    if recon_method == "dti":
        method_rec = 1
    elif recon_method == "gqi":
        method_rec = 4
    else:
        sys.exit(f"Unknown reconstruction method: {recon_method}")

    file_fib = os.path.join(dir_fib, filename_base + ext_fib)

    #Reconstruction command
    cmd = [
        dsi_studio,
        "--action=rec",
        f"--source={file_src_real}",
        f"--mask={file_msk}",
        f"--method={method_rec}",
        "--other_output=all",
        f"--output={file_fib}",
        "--check_btable=0",
    ]

    if additional_cmd:
        cmd.append(f"--cmd={additional_cmd}")

    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)

    # move fib to corresponding folders
    move_files(dir_src, dir_fib, f'/*{ext_fib}')

    fib_candidates = sorted(glob.glob(os.path.join(dir_fib, f"*{ext_fib}")))
    if not fib_candidates:
        raise FileNotFoundError(f"No reconstructed FIB/FZ file found in {dir_fib}")

    file_fib = fib_candidates[0]

    if recon_method == "gqi":
        exports = ["dti_fa", "md", "ad", "rd"]
    elif recon_method == "dti":
        exports = ["fa", "md", "ad", "rd"]
    else:
        sys.exit(f"Unknown reconstruction method: {recon_method}")

    for metric in exports:
        cmd = [
            dsi_studio,
            "--action=exp",
            f"--source={file_fib}",
            f"--export={metric}",
        ]
        print("Running:", " ".join(cmd))
        subprocess.run(cmd, check=True)

    for metric in exports:
        for metric_path in glob.glob(os.path.join(dir_fib, f"*.{metric}.nii.gz")):
            print(f"Reorienting {metric_path} to LIP")
            reorient_nifti_to_lip(metric_path)

    move_files(dir_fib, dir_qa, '/*fa.nii.gz')
    move_files(dir_fib, dir_qa, '/*md.nii.gz')
    move_files(dir_fib, dir_qa, '/*ad.nii.gz')
    move_files(dir_fib, dir_qa, '/*rd.nii.gz')

    #PNG generation creating QA images.
    all_nifti_files = glob.glob(os.path.join(dir_qa, "*.nii")) + glob.glob(os.path.join(dir_qa, "*.nii.gz"))
    #delete duplicates
    nifti_files = sorted(set(all_nifti_files))

    for nifti_path in nifti_files:
        base_name = os.path.basename(nifti_path)
        if base_name.endswith(".nii.gz"):
            base_name = base_name[:-7]
        elif base_name.endswith(".nii"):
            base_name = base_name[:-4]
        else:
            base_name = os.path.splitext(base_name)[0]

        png_slice_path = os.path.join(dir_qa, f"{base_name}.png")
        cmd = ["slicer", nifti_path, "-L", "-a", png_slice_path]
        print("Running:", " ".join(cmd))
        subprocess.run(cmd, check=True)

    return float(min_vox_size_mm)

def tracking(dsi_studio, dir_in, track_param='default', min_voxel_size_mm=0.1, thread_count=1, legacy=False):
    """
    Performs seed-based fiber-tracking.
    Default parameters are used unless a custom parameter is specified.
    """
    if not os.path.isdir(dir_in):
        sys.exit(f"Input directory does not exist: {dir_in}")

    if legacy:
        ext_fib = '.fib.gz'
    else:
        ext_fib = '.fz'

    # Define parameter sets
    param_sets = {
        'default':        ['0AD7A33C9A99193FE8D5123F0AD7233CCDCCCC3D9A99993EbF04240420FdcaCDCC4C3Ec'],
        'aida_optimized': [1000000, '.01', '55', 0, '.02', '.1', '.3', '120.0'],
        'rat':            [1000000, '.01', '60', 0, '.02', '.1', '.3', '20.0'],
        'mouse':          [1000000, '.01', '45', 0, '.02', '.1', '.3', '15.0'],
        'test':           [10000, '.01', '45', 0, '.02', '.1', '.3', '15.0'],
        #tract_count, step_size, turning_angle, check_ending, fa_threshold, smoothing, min_length, max_length
    }

    if isinstance(track_param, str):
        track_param_key = track_param.lower()
        params = param_sets.get(track_param_key)
        if params is None:
            sys.exit(f'Unknown track_param set: {track_param}')
        params = list(params)
    elif isinstance(track_param, list):
        track_param_key = "custom"
        params = list(track_param)
    else:
        sys.exit(
            'track_param must be "default", "aida_optimized", "rat", "mouse", '
            "or a list of values for tract_count, step_size, "
            "turning_angle, check_ending, fa_threshold, smoothing, min_length, and max_length."
        )

    fib_candidates = sorted(glob.glob(os.path.join(dir_in, f"*{ext_fib}")))
    if not fib_candidates:
        raise FileNotFoundError(f"No reconstructed file '*{ext_fib}' found in {dir_in}")
    if len(fib_candidates) > 1:
        print(f"Warning: multiple reconstructed files found. Using: {fib_candidates[0]}")

    filename = fib_candidates[0]
    track_file = filename + ".trk.gz"

    # Set tracking based on track_param:
    if track_param_key == "default":
        print('Using DSI Studio default tracking parameters')
        # Use this tracking parameters in the form of parameter_id that you can get directly from the dsi_studio gui console. (this is here now the defualt mode)
        cmd = [
            dsi_studio,
            "--action=trk",
            f"--source={filename}",
            f"--output={track_file}",
            f"--parameter_id={params[0]}",
        ]

    else:
        # Use this tracking parameters if you want to specify each tracking parameter separately.
        if track_param_key != "aida_optimized":
            params[1] = min_voxel_size_mm / 2
            params[6] = min_voxel_size_mm * 2

        # The tract-density image saved here may not be viewable in FSLeyes or ITK-SNAP, but is compatible with Mango viewer.
        cmd = [
            dsi_studio,
            "--action=trk",
            f"--source={filename}",
            f"--output={track_file}",
            f"--tract_count={int(params[0])}",
            f"--step_size={params[1]}",
            f"--turning_angle={params[2]}",
            f"--check_ending={int(params[3])}",
            f"--fa_threshold={params[4]}",
            f"--smoothing={params[5]}",
            f"--min_length={params[6]}",
            f"--max_length={params[7]}",
            f"--thread_count={thread_count}",
            "--export=tdi:color",
        ]

        # parameters = (dsi_studio, 'trk', filename, os.path.join(dir_in, filename+'.trk.gz'), 1000000, 0, '.5', '55', 0, '.02', '.1', '.5', '12.0') #Our Old parameters
        # parameters = (dsi_studio, 'trk', filename, os.path.join(dir_in, filename+'.trk.gz'), 1000000, 0, '.01', '55', 0, '.02', '.1', '.3', '120.0') #Here are the optimized parameters (fatemeh)

    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)

    tdi_color_file = track_file + ".tdi:color.nii.gz"
    tdi_color_renamed = track_file + ".tdi_color.nii.gz"

    if os.path.exists(tdi_color_file):
        os.replace(tdi_color_file, tdi_color_renamed)

def connectivity(dsi_studio, dir_in, dir_seeds, dir_out, dir_con, make_isotropic=0, legacy=False):
    """
Calculates connectivity data (types: pass and end).
"""
    if not os.path.isdir(dir_in):
        sys.exit(f"Input directory does not exist: {dir_in}")

    # dir_seeds is a seed/ROI/atlas file, not a directory.
    # If it is relative, resolve it relative to the DWI folder.
    if not os.path.isabs(dir_seeds):
        dir_seeds = os.path.join(os.path.dirname(dir_in), "DSI_studio", dir_seeds)

    dir_seeds = os.path.normpath(dir_seeds)

    if not os.path.isfile(dir_seeds):
        sys.exit(f"Seed/ROI file does not exist: {dir_seeds}")

    if not os.path.exists(dir_out):
        sys.exit(f"Output path does not exist: {dir_out}")

    dir_con = make_dir(dir_out, dir_con)

    ext_fib = ".fib.gz" if legacy else ".fz"

    fib_candidates = sorted(glob.glob(os.path.join(dir_in, f"*{ext_fib}")))
    if not fib_candidates:
        raise FileNotFoundError(f"No reconstructed file '*{ext_fib}' found in {dir_in}")
    if len(fib_candidates) > 1:
        print(f"WARNING: Multiple reconstructed files found. Using: {fib_candidates[0]}")

    filename = fib_candidates[0]

    trk_candidates = sorted(glob.glob(os.path.join(dir_in, "*trk.gz")))
    if not trk_candidates:
        raise FileNotFoundError(f"No tract file '*trk.gz' found in {dir_in}")
    if len(trk_candidates) > 1:
        print(f"WARNING: Multiple tract files found. Using: {trk_candidates[0]}")

    file_trk = trk_candidates[0]
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
        resampled_seeds_path = os.path.join(
            dir_con,
            strip_nifti_suffix(file_seeds) + "_resampled.nii.gz"
        )

        cmd_resample = [
            "flirt",
            "-in", file_seeds,
            "-ref", file_seeds,
            "-applyisoxfm", str(iso_value),
            "-nosearch",
            "-interp", "nearestneighbour",
            "-out", resampled_seeds_path,
        ]
        subprocess.run(cmd_resample, check=True)

        print(f'Resampling seeds image to {iso_value} mm isotropic voxel size')

        #PNG generation creating QA images.
        qc_orig_png = os.path.join(dir_con, "qc_seeds_orig.png")
        qc_resampled_png = os.path.join(dir_con, "qc_seeds_resampled.png")
        qc_combined_png = os.path.join(dir_con, "qc_resampled_seeds_combined.png")

        print("Creating quality-control image for original seed/ROI image")
        subprocess.run(["slicer", file_seeds, "-L", "-a", qc_orig_png], check=True)

        print("Creating quality-control image for resampled seed/ROI image")
        subprocess.run(["slicer", resampled_seeds_path, "-L", "-a", qc_resampled_png], check=True)

        print("Combining seed/ROI quality-control images")
        subprocess.run(
            ["pngappend", qc_orig_png, "-", qc_resampled_png, qc_combined_png],
            check=True
        )

        # Use resampled seed/ROI image for connectivity.
        file_seeds = resampled_seeds_path
        # # command needs to use the --t1t2 t2_rare.nii.gz to align the ROI files
        # cmd_ana = r'%s --action=%s --source=%s --tract=%s --connectivity=%s --connectivity_value=%s --connectivity_type=%s --t1t2=%s'
        # # Inverse scale the file_seeds by 10
        # # file_seeds = scaleBy10(file_seeds, inv=True)
    # Performs analysis on every connectivity value and type. DSI Studio reuses
    # the same default connectivity filename, so each result must be renamed
    # immediately after the corresponding command finishes.
    connect_vals = ['qa', 'count']
    connect_types = ['pass', 'end']
    tract_dir = os.path.dirname(file_trk)
    tract_base = os.path.basename(file_trk)
    roi_base = strip_nifti_suffix(os.path.basename(file_seeds))

    for value in connect_vals:
        for connectivity_type in connect_types:
            cmd = [
                dsi_studio,
                "--action=ana",
                f"--source={filename}",
                f"--tract={file_trk}",
                f"--connectivity={file_seeds}",
                f"--connectivity_value={value}",
                f"--connectivity_type={connectivity_type}",
            ]

            print("Running:", " ".join(cmd))
            subprocess.run(cmd, check=True)

            connectivity_outputs = sorted(glob.glob(os.path.join(
                tract_dir,
                f"{tract_base}.{roi_base}.connectivity.*"
            )))

            if not connectivity_outputs:
                print(
                    f"WARNING: No connectivity output files found for "
                    f"tract='{tract_base}', roi='{roi_base}', "
                    f"value='{value}', type='{connectivity_type}'."
                )

            for output_path in connectivity_outputs:
                output_dir = os.path.dirname(output_path)
                output_name = os.path.basename(output_path)

                if f".{value}.{connectivity_type}.connectivity." in output_name:
                    continue

                renamed_name = output_name.replace(
                    ".connectivity.",
                    f".{value}.{connectivity_type}.connectivity.",
                    1
                )

                os.replace(output_path, os.path.join(output_dir, renamed_name))

            move_files(tract_dir, dir_con, "/*.txt")
            move_files(tract_dir, dir_con, "/*.mat")

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
