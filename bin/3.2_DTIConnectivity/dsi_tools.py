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
    image_data = data.get_fdata()
    if inv is False:
        # Scale only the affine. The voxel data stay unchanged, but FSL sees a
        # brain size closer to human dimensions during motion correction.
        scale = np.eye(4)
        scale[0, 0] = 10
        scale[1, 1] = 10
        scale[2, 2] = 10

        scaled_affine = data.affine @ scale

        scaled_image = nib.Nifti1Image(image_data, scaled_affine)
        scaled_path = os.path.join(os.path.dirname(input_path), 'fslScaleTemp.nii.gz')
        nib.save(scaled_image, scaled_path)
        return scaled_path
    elif inv is True:
        # Undo the temporary affine scaling after FSL processing.
        inv_scale = np.eye(4)
        inv_scale[0, 0] = 0.1
        inv_scale[1, 1] = 0.1
        inv_scale[2, 2] = 0.1

        unscaled_affine = data.affine @ inv_scale
        unscaled_image = nib.Nifti1Image(image_data, unscaled_affine)
        output_header = unscaled_image.header
        output_header.set_xyzt_units('mm')

        nib.save(unscaled_image, input_path)
        return input_path
    else:
        sys.exit("Error: inv - parameter should be a boolean.")


def findSlicesData(path, pre):
    slice_files = []
    matching_files = glob.iglob(path + '/' + pre + '*.nii.gz', recursive=True)
    for filename in matching_files:
        slice_files.append(filename)
    slice_files.sort()
    return slice_files

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
    # Temporarily scale the affine before FSL slice-wise motion correction.
    input_stem = os.path.basename(input_file).split('.')[0]
    scaled_input_file = scaleBy10(input_file, inv=False)

    original_cwd = os.getcwd()
    temp_dir = os.path.join(os.path.dirname(input_file), "temp")
    if not os.path.exists(temp_dir):
        os.mkdir(temp_dir)

    os.chdir(temp_dir)

    slice_axis = infer_slice_axis(scaled_input_file)
    split_interface = fsl.Split(in_file=scaled_input_file, dimension=slice_axis, out_base_name=input_stem)
    split_interface.run()
    os.remove(scaled_input_file)

    # Split the volume along the inferred slice axis and correct every slice.
    slice_files = findSlicesData(os.getcwd(), input_stem)

    for slc in slice_files:
        output_file = os.path.join(par_folder, os.path.basename(slc))
        mcflirt = fsl.preprocess.MCFLIRT(in_file=slc, out_file=output_file, save_plots=True, terminal_output='none')
        mcflirt.run()
        os.remove(slc)

    # Merge corrected slices back to a single volume.
    corrected_slice_files = findSlicesData(par_folder, input_stem)
    output_file = os.path.join(os.path.dirname(input_file),
                               os.path.basename(input_file).split('.')[0]) + '_mcf.nii.gz'
    merge_interface = fsl.Merge(in_files=corrected_slice_files, dimension=slice_axis, merged_file=output_file)
    merge_interface.run()

    for slc in corrected_slice_files:
        os.remove(slc)

    # Restore the original affine scale.
    output_file = scaleBy10(output_file, inv=True)
    
    os.chdir(original_cwd)
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
    """
    Move files matching pattern from dir_in to dir_out.

    Existing call sites pass patterns such as '/*.nii.gz'. Strip the leading
    slash before joining so the source and destination paths are explicit.
    """
    time.sleep(1.0)
    file_list = glob.glob(os.path.join(dir_in, pattern.lstrip("/")))
    file_list.sort()

    time.sleep(1.0)
    for source_path in file_list:
        destination_path = os.path.join(dir_out, os.path.basename(source_path))
        shutil.copy(source_path, destination_path)
        os.remove(source_path)


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
    Create DSI Studio source/reconstruction files and export scalar metrics.
    """
    source_image_file = dir_in
    mask_file = dir_msk
    output_anchor = dir_out

    src_dir_name = r'src'
    fib_dir_name = r'fib_map'
    dsi_metrics_dir_name = r'DSI_studio'

    # Support for backwards compatibility with pre-2024 DSI Studio (AIDAmri <= v2.0)
    if legacy:
        ext_src = '.src.gz'
        ext_fib = '.fib.gz'
    else:
        ext_src = '.sz'
        ext_fib = '.fz'

    if not os.path.exists(source_image_file):
        sys.exit("Input file \"%s\" does not exist." % (source_image_file,))

    input_dir = os.path.dirname(os.path.abspath(source_image_file))

    #Resolve mask_file relative to the DWI input directory and normalize the result.
    if not os.path.isabs(mask_file):
        mask_file = os.path.join(input_dir, mask_file)
    mask_file = os.path.normpath(mask_file)

    if not os.path.isfile(mask_file):
        sys.exit("Mask file \"%s\" does not exist." % (mask_file,))

    if not os.path.exists(output_anchor):
        sys.exit("Output path \"%s\" does not exist." % (output_anchor,))

    # Validate the b-table or bval/bvec information that will be passed to DSI Studio.
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

    output_root = os.path.dirname(output_anchor)
    src_dir = make_dir(output_root, src_dir_name)
    fib_dir = make_dir(output_root, fib_dir_name)
    dsi_metrics_dir = make_dir(output_root, dsi_metrics_dir_name)

    # change to input directory
    os.chdir(input_dir)

    # Prefer explicit gradient files when they are available. This keeps the
    # source image coupled to the original gradients even after motion correction
    # changed the NIfTI filename to *_mcf.nii.gz.
    input_filename = os.path.basename(source_image_file)
    input_stem = strip_nifti_suffix(input_filename)
    expected_src_file = os.path.join(src_dir, input_stem + ext_src)
    if gradient_pair is not None:
        cmd = [
            dsi_studio,
            "--action=src",
            f"--source={input_filename}",
            f"--output={expected_src_file}",
            f"--bval={gradient_pair['bval']}",
            f"--bvec={gradient_pair['bvec']}",
        ]
    else:
        cmd = [
            dsi_studio,
            "--action=src",
            f"--source={input_filename}",
            f"--output={expected_src_file}",
            f"--b_table={b_table}",
        ]

    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)

    patterns = [
        os.path.join(src_dir, input_stem + '*.src.gz.sz'),
        os.path.join(src_dir, input_stem + '*.sz'),
        os.path.join(src_dir, input_stem + '*.src.gz'),
    ]

    src_candidates = []
    for pattern in patterns:
        src_candidates.extend(sorted(glob.glob(pattern)))

    src_candidates = list(dict.fromkeys(src_candidates))

    if not src_candidates:
        raise FileNotFoundError(
            f"No DSI Studio source file found for base '{input_stem}' in {src_dir}"
        )

    src_file = src_candidates[0]

    # If unrealistic streamlines cross top of cortex are present due to an oversized mask, erode mask
    mask_erosion = 0
    if mask_erosion > 0:
        eroded_mask_path = os.path.join(
            os.path.dirname(mask_file),
            strip_nifti_suffix(mask_file) + "_eroded_mask.nii.gz"
        )

        mask_file = erode_mask(mask_file, eroded_mask_path, n_voxels=mask_erosion)
        print(f"Eroded mask saved to {mask_file}")

    # param0 is currently reported for traceability only. The active DSI Studio
    # command below does not pass it unless the reconstruction command is extended.
    param_zero='1.25'
    if vivo == "ex_vivo":
        param_zero='0.60'
        print(f'Using param0 value {param_zero} recommended for ex vivo data')
    elif vivo == "in_vivo":
        print(f'Using param0 value {param_zero} recommended for in vivo data')

    min_vox_size_mm = get_min_voxel_size_mm(source_image_file)

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
    # If diffusion data are resampled here, all ROI/parcellation images must be
    # created in the same diffusion space before connectivity is calculated.

    additional_cmd = ''
    # Atlas/ROI images must match this resampled space when iso_value > 0.
    if iso_value > 0:
        additional_cmd = f'[Step T2][Edit][Resample]={iso_value}'
        print(f'Resampling to {iso_value} mm isotropic voxel size')

    # Optional future DSI Studio corrections. Currently disabled because AIDAmri
    # performs slice-wise motion correction before source generation.
    use_eddy_correct = False
    use_dsi_topup = False
    rev_pe_image = '' # path to reverse phase encoding image if DSI Studio topup is enabled
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

    fib_file = os.path.join(fib_dir, input_stem + ext_fib)

    # Reconstruction command. method: 1=DTI, 4=GQI.
    cmd = [
        dsi_studio,
        "--action=rec",
        f"--source={src_file}",
        f"--mask={mask_file}",
        f"--method={method_rec}",
        "--other_output=all", #diffusion metrics to compute. 'all' for every possible measure(fa,rd,rdi)
        f"--output={fib_file}",
        "--check_btable=0",#if 1 checks the gradient table and flips/swaps to fix gradient directions
    ]

    if additional_cmd:
        cmd.append(f"--cmd={additional_cmd}")

    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)

    # Older DSI Studio releases may write reconstructed files next to the source.
    move_files(src_dir, fib_dir, f'/*{ext_fib}')

    fib_candidates = sorted(glob.glob(os.path.join(fib_dir, f"*{ext_fib}")))
    if not fib_candidates:
        raise FileNotFoundError(f"No reconstructed FIB/FZ file found in {fib_dir}")

    fib_file = fib_candidates[0]

    if recon_method == "gqi":
        # DSI Studio exports the DTI-derived FA metric from GQI reconstruction as
        # dti_fa to distinguish it from GQI-specific anisotropy metrics.
        exports = ["dti_fa", "md", "ad", "rd"]
    elif recon_method == "dti":
        exports = ["fa", "md", "ad", "rd"]
    else:
        sys.exit(f"Unknown reconstruction method: {recon_method}")

    for metric in exports:
        cmd = [
            dsi_studio,
            "--action=exp",
            f"--source={fib_file}",
            f"--export={metric}",
        ]
        print("Running:", " ".join(cmd))
        subprocess.run(cmd, check=True)

    for metric in exports:
        for metric_path in glob.glob(os.path.join(fib_dir, f"*.{metric}.nii.gz")):
            print(f"Reorienting {metric_path} to LIP")
            reorient_nifti_to_lip(metric_path)

    move_files(fib_dir, dsi_metrics_dir, '/*fa.nii.gz')
    move_files(fib_dir, dsi_metrics_dir, '/*md.nii.gz')
    move_files(fib_dir, dsi_metrics_dir, '/*ad.nii.gz')
    move_files(fib_dir, dsi_metrics_dir, '/*rd.nii.gz')

    # Generate PNG quick-look images for all NIfTI outputs in DSI_studio.
    all_nifti_files = glob.glob(os.path.join(dsi_metrics_dir, "*.nii")) + glob.glob(os.path.join(dsi_metrics_dir, "*.nii.gz"))
    nifti_files = sorted(set(all_nifti_files))

    for nifti_path in nifti_files:
        base_name = os.path.basename(nifti_path)
        if base_name.endswith(".nii.gz"):
            base_name = base_name[:-7]
        elif base_name.endswith(".nii"):
            base_name = base_name[:-4]
        else:
            base_name = os.path.splitext(base_name)[0]

        png_slice_path = os.path.join(dsi_metrics_dir, f"{base_name}.png")
        cmd = ["slicer", nifti_path, "-L", "-a", png_slice_path]
        print("Running:", " ".join(cmd))
        subprocess.run(cmd, check=True)

    return float(min_vox_size_mm)

def tracking(dsi_studio, dir_in, track_param='default', min_voxel_size_mm=0.1, thread_count=1, legacy=False):
    """
    Performs seed-based fiber-tracking.
    Default parameters are used unless a custom parameter is specified.
    """
    fib_dir = dir_in
    if not os.path.isdir(fib_dir):
        sys.exit(f"Input directory does not exist: {fib_dir}")

    if legacy:
        ext_fib = '.fib.gz'
    else:
        ext_fib = '.fz'

    # Parameter sets:
    # tract_count, step_size, turning_angle, check_ending, fa_threshold,
    # smoothing, min_length, max_length
    param_sets = {
        'default':        ['0AD7A33C9A99193FE8D5123F0AD7233CCDCCCC3D9A99993EbF04240420FdcaCDCC4C3Ec'],
        'aida_optimized': [1000000, '.01', '55', 0, '.02', '.1', '.3', '120.0'],
        'rat':            [1000000, '.01', '60', 0, '.02', '.1', '.3', '20.0'],
        'mouse':          [1000000, '.01', '45', 0, '.02', '.1', '.3', '15.0'],
        'test':           [10000, '.01', '45', 0, '.02', '.1', '.3', '15.0'],
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

    fib_candidates = sorted(glob.glob(os.path.join(fib_dir, f"*{ext_fib}")))
    if not fib_candidates:
        raise FileNotFoundError(f"No reconstructed file '*{ext_fib}' found in {fib_dir}")
    if len(fib_candidates) > 1:
        sys.exit(
            f"Multiple reconstructed files found in {fib_dir}. "
            "Remove duplicates or run in a clean output folder."
        )

    fib_file = fib_candidates[0]
    track_file = fib_file + ".trk.gz"

    # Set tracking based on track_param:
    if track_param_key == "default":
        print('Using DSI Studio default tracking parameters')
        # DSI Studio can replay GUI tracking settings through a parameter_id.
        cmd = [
            dsi_studio,
            "--action=trk",
            f"--source={fib_file}",
            f"--output={track_file}",
            f"--parameter_id={params[0]}",
        ]

    else:
        # Custom/predefined sets pass each tracking parameter explicitly.
        if track_param_key != "aida_optimized":
            params[1] = min_voxel_size_mm / 2
            params[6] = min_voxel_size_mm * 2

        # Export a color TDI next to the tract file for visual QC.
        cmd = [
            dsi_studio,
            "--action=trk",
            f"--source={fib_file}",
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

    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)

    tdi_color_file = track_file + ".tdi:color.nii.gz"
    tdi_color_renamed = track_file + ".tdi_color.nii.gz"

    if os.path.exists(tdi_color_file):
        os.replace(tdi_color_file, tdi_color_renamed)

def connectivity(dsi_studio, dir_in, dir_seeds, dir_out, dir_con, make_isotropic=0, legacy=False):
    """
    Calculates connectivity data for pass-through and endpoint definitions.
    """
    fib_dir = dir_in
    seed_file = dir_seeds
    output_root = dir_out
    connectivity_dir_name = dir_con

    if not os.path.isdir(fib_dir):
        sys.exit(f"Input directory does not exist: {fib_dir}")

    # seed_file is a seed/ROI/atlas file, not a directory. Relative paths are
    # kept for backwards compatibility and resolved against DSI_studio.
    if not os.path.isabs(seed_file):
        seed_file = os.path.join(os.path.dirname(fib_dir), "DSI_studio", seed_file)

    seed_file = os.path.normpath(seed_file)

    if not os.path.isfile(seed_file):
        sys.exit(f"Seed/ROI file does not exist: {seed_file}")

    if not os.path.exists(output_root):
        sys.exit(f"Output path does not exist: {output_root}")

    connectivity_dir = make_dir(output_root, connectivity_dir_name)

    ext_fib = ".fib.gz" if legacy else ".fz"

    # Search for fib_files
    fib_candidates = sorted(glob.glob(os.path.join(fib_dir, f"*{ext_fib}")))
    if not fib_candidates:
        raise FileNotFoundError(f"No reconstructed file '*{ext_fib}' found in {fib_dir}")
    if len(fib_candidates) > 1:
        sys.exit(
            f"Multiple reconstructed files found in {fib_dir}. "
            "Remove duplicates or run in a clean output folder."
        )

    fib_file = fib_candidates[0]

    #Search for tracking file
    trk_candidates = sorted(glob.glob(os.path.join(fib_dir, "*trk.gz")))
    if not trk_candidates:
        raise FileNotFoundError(f"No tract file '*trk.gz' found in {fib_dir}")
    if len(trk_candidates) > 1:
        sys.exit(
            f"Multiple tract files found in {fib_dir}. "
            "Remove duplicates or run in a clean output folder."
        )

    tract_file = trk_candidates[0]

    # If srcgen resampled diffusion data, use a nearest-neighbor resampled copy
    # of the seed/ROI image for connectivity.
    iso_value = None
    if make_isotropic == "auto":
        iso_value = float(get_min_voxel_size_mm(seed_file))
    else:
        iso_value = float(make_isotropic)
    if iso_value is not None and iso_value > 0:
        resampled_seeds_path = os.path.join(
            connectivity_dir,
            strip_nifti_suffix(seed_file) + "_resampled.nii.gz"
        )

        cmd_resample = [
            "flirt",
            "-in", seed_file,
            "-ref", seed_file,
            "-applyisoxfm", str(iso_value),
            "-nosearch",
            "-interp", "nearestneighbour",
            "-out", resampled_seeds_path,
        ]
        subprocess.run(cmd_resample, check=True)

        print(f'Resampling seeds image to {iso_value} mm isotropic voxel size')

        # Generate PNG quick-look images for the original and resampled seed/ROI.
        qc_orig_png = os.path.join(connectivity_dir, "qc_seeds_orig.png")
        qc_resampled_png = os.path.join(connectivity_dir, "qc_seeds_resampled.png")
        qc_combined_png = os.path.join(connectivity_dir, "qc_resampled_seeds_combined.png")

        print("Creating quality-control image for original seed/ROI image")
        subprocess.run(["slicer", seed_file, "-L", "-a", qc_orig_png], check=True)

        print("Creating quality-control image for resampled seed/ROI image")
        subprocess.run(["slicer", resampled_seeds_path, "-L", "-a", qc_resampled_png], check=True)

        print("Combining seed/ROI quality-control images")
        subprocess.run(
            ["pngappend", qc_orig_png, "-", qc_resampled_png, qc_combined_png],
            check=True
        )

        seed_file = resampled_seeds_path
    # Performs analysis on every connectivity value and type. DSI Studio reuses
    # the same default connectivity filename, so each result must be renamed
    # immediately after the corresponding command finishes.
    connectivity_values = ['qa', 'count']
    connectivity_types = ['pass', 'end']
    tract_dir = os.path.dirname(tract_file)
    tract_base = os.path.basename(tract_file)
    roi_base = strip_nifti_suffix(os.path.basename(seed_file))

    for value in connectivity_values:
        for connectivity_type in connectivity_types:
            cmd = [
                dsi_studio,
                "--action=ana",
                f"--source={fib_file}",
                f"--tract={tract_file}",
                f"--connectivity={seed_file}",
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

            move_files(tract_dir, connectivity_dir, "/*.txt")
            move_files(tract_dir, connectivity_dir, "/*.mat")

if __name__ == '__main__':
    pass
