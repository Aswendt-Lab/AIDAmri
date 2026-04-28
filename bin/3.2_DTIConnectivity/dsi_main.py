#!/opt/env/bin/python
"""
Created on 10/08/2017

@author: Niklas Pallast
Neuroimaging & Neuroengineering
Department of Neurology
University Hospital Cologne
"""
from __future__ import print_function

import atexit
import argparse
import os
import glob
import dsi_tools
import shutil
import gzip
import subprocess
import sys


def enable_process_log(log_path):
    # Mirror stdout/stderr into a dataset-local log file so debugging works the
    # same whether the script is launched from the image or a mounted checkout.
    tee_proc = subprocess.Popen(["tee", "-a", log_path], stdin=subprocess.PIPE)
    os.dup2(tee_proc.stdin.fileno(), sys.stdout.fileno())
    os.dup2(tee_proc.stdin.fileno(), sys.stderr.fileno())
    sys.stdout = os.fdopen(sys.stdout.fileno(), "w", buffering=1, closefd=False)
    sys.stderr = os.fdopen(sys.stderr.fileno(), "w", buffering=1, closefd=False)

    def _cleanup():
        try:
            sys.stdout.flush()
            sys.stderr.flush()
        finally:
            try:
                tee_proc.stdin.close()
            except Exception:
                pass
            try:
                tee_proc.wait(timeout=5)
            except Exception:
                pass

    atexit.register(_cleanup)


def should_enable_process_log():
    # batchProc.py already redirects stdout/stderr into its own step log. In
    # that case dsi_main.py should not create an additional process.log or tee
    # the same lines twice. Allow an explicit override via environment variable
    # and otherwise fall back to a TTY check for direct interactive runs.
    if os.environ.get("AIDAMRI_DISABLE_PROCESS_LOG", "").lower() in {"1", "true", "yes"}:
        return False
    return os.isatty(sys.stdout.fileno()) and os.isatty(sys.stderr.fileno())

if __name__ == '__main__':
    # Resolve helper files relative to this script so the mounted repo can be
    # executed directly without depending on the shell's current directory.
    script_dir = os.path.dirname(os.path.abspath(__file__))
    # Keep the DSI Studio binary path in a sidecar file next to this script.
    dsi_path_file = os.path.join(script_dir, "dsi_studioPath.txt")
    with open(dsi_path_file, "r") as f:
        dsi_studio = f.read().splitlines()[0]

    # Use explicit b-table paths only when the user asks for them. In auto
    # mode, the pipeline now requires a real .bval/.bvec pair.
    b_table = None
    gradient_pair = None

    # default connectivity directory relative to input directory
    dir_con = r'connectivity'

    # Defining CLI flags
    parser = argparse.ArgumentParser(description='Get connectivity of DTI dataset')
    requiredNamed = parser.add_argument_group('Required named arguments')
    requiredNamed.add_argument('-i',
                               '--file_in',
                               help = 'path to the raw NIfTI DTI file (ends with *dwi.nii.gz)',
                               required=True
                               )
    parser.add_argument('-b',
                        '--b_table',
                        default='auto',  # Default to 'auto' for automatic selection
                        help='Specify the diffusion gradient source: "auto" requires matching .bval/.bvec files. Any other value is treated as a b-table path.'
                        )
    parser.add_argument('-r',
                        '--recon_method',
                        default='dti',
                        type=str.lower,
                        choices=['dti', 'gqi'],
                        help='Specify diffusion reconstruction method ("gqi" or default "dti").',
                        required=False
                       )
    parser.add_argument('-v',
                        '--vivo',
                        default='in_vivo',
                        type=str.lower,
                        choices=['in_vivo', 'ex_vivo'],
                        help='Specify in vivo or ex vivo data to adjust sampling length ratio (param0). "in_vivo" param0=1.25 (default), "ex_vivo" param0=0.60.',
                        required=False
                       )
    parser.add_argument('-m',
                        '--make_isotropic',
                        default='0',
                        help='Specify an isotropic voxel size in mm for resampling. Default 0 = no resampling. "auto" uses nibabel to read the NIFTI header for the minimum voxel size',
                        required=False
                       )
    parser.add_argument('-t',
                        '--track_params',
                        default='default',
                        help='Specify tracking parameters from a pre-defined set ("aida_optimized", "rat", or "mouse") or as a list of values for fiber_count, interpolation, step_size, turning_angle, check_ending, fa_threshold, smoothing, min_length, and max_length.',
                        required=False
                       )
    parser.add_argument('-template',
                        '--template',
                        default='mouse',
                        choices = ['mouse', 'rat'],
                        type=str.lower,
                        help='Specify the template to use for the reconstruction. Default is mouse. Other options is "rat"',
                        required=False
                       )
    parser.add_argument('-thread_count',
                        '--thread_count',
                        type=int,
                        default=1,
                        help='Specify the number of threads to use for fiber tracking. Default is 1.',
                        required=False
                        )
    parser.add_argument('-l',
                        '--legacy',
                        help='Legacy file types for DSI-Studio releases before 2024. Default is False (uses new more storage-efficient ".sz" and ".fz" file types)',
                        action = 'store_true'
                        )
    parser.add_argument('-nomcf',
                        '--no_motion_correction',
                        action='store_true',
                        help='Specify whether to skip motion correction. Default is False (perform motion correction). Set to "true" to skip motion correction.',
                        required=False
                        )
    parser.add_argument('-o',
                        '--optional',
                        nargs = '*',
                        help = 'Optional arguments.\n\t"fa0": Renames the FA metric data to former DSI naming convention.\n\t"nii_gz": Converts ROI labeling relating files from .nii to .nii.gz format to match former data structures.'
                        )
    args = parser.parse_args()

    file_cur = os.path.dirname(args.file_in)
    process_log = os.path.join(file_cur, "process.log")
    if should_enable_process_log():
        enable_process_log(process_log)
        print(f"Writing process log to {process_log}")
        
    # Determine the gradient source based on the -b option.
    if str(args.b_table).lower() == 'auto':
        input_dir = os.path.dirname(os.path.abspath(args.file_in))
        input_stem = dsi_tools.strip_nifti_suffix(args.file_in)
        gradient_pair, gradient_error = dsi_tools.find_matching_gradient_pair(
            input_dir,
            preferred_stem=input_stem,
        )
        if gradient_pair is None:
            sys.exit(f"ERROR: {gradient_error}")
        else:
            print(f"Using gradient files: {gradient_pair['bval']} and {gradient_pair['bvec']}")
    else:
        b_table = args.b_table
        print(f"Using explicit b-table: {b_table}")

    # Preparing directories
    dsi_path = os.path.join(file_cur, 'DSI_studio')
    mcf_path = os.path.join(file_cur, 'mcf_Folder')
    dir_mask = sorted(glob.glob(os.path.join(dsi_path, '*BetMask_scaled.nii')))
    if not dir_mask:
        dir_mask = sorted(glob.glob(os.path.join(dsi_path, '*BetMask_scaled.nii.gz'))) # check for ending (either .nii or .nii.gz)
        if not dir_mask:
            # check for mask without scaled in name
            dir_mask = sorted(glob.glob(os.path.join(dsi_path, '*BetMask.nii.gz')))
    if not dir_mask:
        raise FileNotFoundError("No BET mask found in DSI_studio folder.")

    dir_mask = dir_mask[0]

    dir_out = args.file_in

    if str(args.make_isotropic).lower() == 'auto':
        make_isotropic = 'auto'
    else:
        make_isotropic = float(args.make_isotropic)
    flip_image_y = False

    if args.template.lower() == 'rat' or str(args.template) == '5':
        template = 5
    elif args.template.lower() == 'mouse' or str(args.template) == '1':
        template = 1
    else:
        sys.exit(
            f"Error: Invalid template value: {args.template}. "
            "Allowed values are: 'mouse' or '1', 'rat' or '5'."
        )

    # if it exists, find the denoised dwi data and use it as file_in
    file_in = args.file_in
    if os.path.exists(file_cur):
        denoised = sorted(glob.glob(os.path.join(file_cur, '*Denoised.nii*')))
        if denoised:
            file_in = denoised[0]

    if os.path.exists(mcf_path):
        shutil.rmtree(mcf_path)
   
    if args.no_motion_correction:
        print("Skipping motion correction")
    else:
        print("Performing slice-wise motion correction")
        os.mkdir(mcf_path)
        # Use FSL's SeparateSliceMoCo to perform motion correction
        if not os.path.exists(dsi_path):
            os.makedirs(dsi_path)
        file_in = dsi_tools.fsl_SeparateSliceMoCo(file_in, mcf_path)

    voxel_size = dsi_tools.srcgen(
        dsi_studio,
        file_in,
        dir_mask,
        dir_out,
        b_table,
        args.recon_method,
        args.vivo,
        make_isotropic,
        flip_image_y,
        template,
        args.legacy,
        gradient_pair=gradient_pair,
    )
    file_in = os.path.join(file_cur,'fib_map')

    track_param = args.track_params

    # Fiber tracking
    dir_out = os.path.dirname(args.file_in)
    dsi_tools.tracking(dsi_studio, file_in, track_param, voxel_size, args.thread_count, args.legacy)

    # Calculating connectivity
    suffixes = ['*StrokeMask_scaled.nii', '*parental_Mask_scaled.nii', '*Anno_scaled.nii', '*AnnoSplit_parental_scaled.nii']
    # if bet4animal is True:
    #     suffixes = ['*StrokeMask.nii', '*parental_Mask.nii', '*Anno.nii', '*AnnoSplit_parental.nii']
    for f in suffixes:
        dir_seeds = sorted(glob.glob(os.path.join(file_cur, 'DSI_studio', f)))
        if not dir_seeds:
            dir_seeds = sorted(glob.glob(os.path.join(file_cur, 'DSI_studio', f + '.gz'))) # check for ending (either .nii or .nii.gz)
        if not dir_seeds:
            continue
        dir_seeds = dir_seeds[0]
        dsi_tools.connectivity(dsi_studio, file_in, dir_seeds, dir_out, dir_con, make_isotropic, flip_image_y, args.legacy)

    # rename files to reduce path length
    confiles = os.path.join(file_cur, dir_con)
    if not os.path.isdir(confiles):
        raise FileNotFoundError(f"Connectivity folder not found: {confiles}")
    data_list = os.listdir(confiles)
    for filename in data_list:
        if args.recon_method == "dti":
            if args.legacy:
                splittedName = filename.split('.src.gz.dti.fib.gz.trk.gz.')
            else:
                splittedName = filename.split('.sz.dti.fz.trk.gz.')
        elif args.recon_method == "gqi":
            if args.legacy:
                splittedName = filename.split('.src.gz.gqi.fib.gz.trk.gz.')
            else:
                splittedName = filename.split('.sz.gqi.fz.trk.gz.')
        if len(splittedName)>1:
            newName = splittedName[1]
            newName = os.path.join(confiles,newName)
            if os.path.isfile(newName):
                os.remove(newName)
            oldName = os.path.join(confiles,filename)
            os.rename(oldName,newName)

    # Including optional arguments regarding deprecated terminology
    if args.optional is not None:
        opts = [s.lower() for s in args.optional]

        file_list = os.listdir(dsi_path)
        for f in file_list:

            # fa0 was a former term used in earlier DSI-studio versions; the '0' in fa0 referred to the first fiber track. However, DTI can only result in one track, therefore only one fractional anisotropy value per voxel is given, thus the collective values are referred to as fa. With the 'fa0' flag toggled on, the 'fa' data file is renamed to the former naming convention (fa0).
            if 'fa0' in opts and f.endswith('fa.nii.gz'):
                newName = f.split('fa.nii.gz')[0] + 'fa0.nii.gz'
                newName = os.path.join(dsi_path, newName)
                oldName = os.path.join(dsi_path, f)
                if os.path.isfile(newName):
                    os.remove(newName)
                os.rename(oldName, newName)
            
            # Due to changes in ROI annotations the corresponding files are saved as '.nii' files as opposed to '.nii.gz' files in earlier versions of DSI studio. With the 'nii_gz' flag toggled on, the '.nii' files are renamed to '.nii.gz'.
            if 'nii_gz' in opts and f.endswith('.nii'):
                oldName = os.path.join(dsi_path, f)
                newName = os.path.join(dsi_path, f + '.gz')

                if os.path.isfile(newName):
                    os.remove(newName)

                with open(oldName, 'rb') as f_in:
                    with gzip.open(newName, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)

                os.remove(oldName)
