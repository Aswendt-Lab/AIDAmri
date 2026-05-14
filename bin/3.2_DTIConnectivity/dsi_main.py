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


def normalize_track_params(values):
    if len(values) == 1:
        return values[0]
    if len(values) != 8:
        raise argparse.ArgumentTypeError(
            "-t/--track_params expects one preset name or exactly 8 custom values: "
            "tract_count step_size turning_angle check_ending fa_threshold smoothing min_length max_length"
        )
    return values


def positive_int(value):
    try:
        parsed = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{value} is not an integer")
    if parsed < 1:
        raise argparse.ArgumentTypeError("value must be >= 1")
    return parsed


if __name__ == '__main__':
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
                        nargs='+',
                        default=['default'],
                        help='Specify tracking parameters from a pre-defined set ("aida_optimized", "rat", or "mouse") or as a list of values for tract_count, step_size, turning_angle, check_ending, fa_threshold, smoothing, min_length, and max_length.',
                        required=False
                       )
    parser.add_argument('--thread_count',
                        type=positive_int,
                        default=1,
                        help='Specify the number of threads to use for fiber tracking. Default is 1.',
                        required=False
                        )
    parser.add_argument('-l',
                        '--legacy',
                        help='Legacy file types for DSI-Studio releases before 2024. Default is False (uses new more storage-efficient ".sz" and ".fz" file types)',
                        action = 'store_true'
                        )
    parser.add_argument('--skip_motion_correction',
                        action='store_true',
                        help='Specify whether to skip motion correction. Default is False (perform motion correction). Set to "true" to skip motion correction.',
                        required=False
                        )
    parser.add_argument('-o',
                        '--optional',
                        nargs = '*',
                        choices=['fa0', 'nii_gz'],
                        help = 'Optional arguments.\n\t"fa0": Renames the FA metric data to the former DSI naming convention.\n\t"nii_gz": Compresses legacy .nii files in DSI_studio to .nii.gz.'
                        )
    args = parser.parse_args()

    try:
        args.track_params = normalize_track_params(args.track_params)
    except argparse.ArgumentTypeError as e:
        parser.error(str(e))

    if str(args.make_isotropic).lower() == 'auto':
        make_isotropic = 'auto'
    else:
        try:
            make_isotropic = float(args.make_isotropic)
        except ValueError:
            parser.error(
                f'Invalid --make_isotropic value "{args.make_isotropic}". '
                'Use 0, "auto", or a voxel size in mm, e.g. 0.2.'
            )

    import dsi_tools

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

    # Connectivity outputs are written relative to the DWI directory.
    connectivity_dir_name = r'connectivity'

    dwi_dir = os.path.dirname(args.file_in)
    process_log = os.path.join(dwi_dir, "process.log")
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

    # Resolve the dataset-local folders and the brain mask used by DSI Studio.
    dsi_metrics_dir = os.path.join(dwi_dir, 'DSI_studio')
    motion_correction_dir = os.path.join(dwi_dir, 'mcf_Folder')
    mask_candidates = sorted(glob.glob(os.path.join(dwi_dir, '*Bet_mask.nii.gz')))
    if not mask_candidates:
        mask_candidates = sorted(glob.glob(os.path.join(dwi_dir, '*Bet_mask.nii')))
    if not mask_candidates:
        raise FileNotFoundError("No BET mask found in DWI folder.")

    brain_mask_file = mask_candidates[0]
    dwi_output_anchor = args.file_in

    # Prefer denoised DWI data when the preprocessing step created it.
    dwi_input_file = args.file_in
    if os.path.exists(dwi_dir):
        denoised = sorted(glob.glob(os.path.join(dwi_dir, '*Denoised.nii*')))
        if denoised:
            dwi_input_file = denoised[0]

    if os.path.exists(motion_correction_dir):
        shutil.rmtree(motion_correction_dir)
   
    if args.skip_motion_correction:
        print("Skipping motion correction")
    else:
        print("Performing slice-wise motion correction")
        os.mkdir(motion_correction_dir)
        dwi_input_file = dsi_tools.fsl_SeparateSliceMoCo(dwi_input_file, motion_correction_dir)

    voxel_size = dsi_tools.srcgen(
        dsi_studio,
        dwi_input_file,
        brain_mask_file,
        dwi_output_anchor,
        b_table,
        args.recon_method,
        args.vivo,
        make_isotropic,
        args.legacy,
        gradient_pair=gradient_pair,
    )
    fib_dir = os.path.join(dwi_dir, 'fib_map')

    tracking_params = args.track_params

    # Fiber tracking
    connectivity_output_root = os.path.dirname(args.file_in)
    dsi_tools.tracking(
        dsi_studio,
        fib_dir,
        tracking_params,
        voxel_size,
        args.thread_count,
        args.legacy
    )

    # Calculating connectivity
    seed_patterns = [
        '*Stroke_mask_anno.nii.gz',
        '*_AnnoSplit.nii.gz',
        '*_AnnoSplit_parental.nii.gz',
    ]

    for seed_pattern in seed_patterns:
        seed_candidates = sorted(glob.glob(os.path.join(dwi_dir, seed_pattern)))
        if not seed_candidates and seed_pattern.endswith('.nii.gz'):
            seed_candidates = sorted(glob.glob(os.path.join(dwi_dir, seed_pattern[:-3])))
        if not seed_candidates:
            print(
                f"WARNING: No connectivity seed/ROI file found for pattern "
                f"{seed_pattern} in {dwi_dir}. Skipping."
            )
            continue
        seed_file = seed_candidates[0]

        # Connectivity
        dsi_tools.connectivity(
            dsi_studio,
            fib_dir,
            seed_file,
            connectivity_output_root,
            connectivity_dir_name,
            make_isotropic,
            args.legacy
        )

    # Rename connectivity files to reduce path length.
    connectivity_dir = os.path.join(dwi_dir, connectivity_dir_name)
    if not os.path.isdir(connectivity_dir):
        raise FileNotFoundError(f"Connectivity folder not found: {connectivity_dir}")

    if args.recon_method == "dti":
        split_token = ".src.gz.dti.fib.gz.trk.gz." if args.legacy else ".sz.dti.fz.trk.gz."
    elif args.recon_method == "gqi":
        split_token = ".src.gz.gqi.fib.gz.trk.gz." if args.legacy else ".sz.gqi.fz.trk.gz."
    else:
        raise ValueError(f"Unknown reconstruction method: {args.recon_method}")

    for filename in os.listdir(connectivity_dir):
        if split_token not in filename:
            continue

        new_filename = filename.split(split_token, 1)[1]
        old_path = os.path.join(connectivity_dir, filename)
        new_path = os.path.join(connectivity_dir, new_filename)

        if os.path.isfile(new_path):
            os.remove(new_path)

        os.rename(old_path, new_path)

    # Optional compatibility renames for older AIDAmri/DSI Studio outputs.
    if args.optional is not None:
        opts = [s.lower() for s in args.optional]

        file_list = os.listdir(dsi_metrics_dir)
        for output_name in file_list:

            # Older DSI Studio versions used fa0 for the first-fiber FA metric.
            # DTI exports only one FA value per voxel, so current outputs use fa.
            if 'fa0' in opts and output_name.endswith('.fa.nii.gz'):
                new_metric_path = output_name.split('fa.nii.gz')[0] + 'fa0.nii.gz'
                new_metric_path = os.path.join(dsi_metrics_dir, new_metric_path)
                old_metric_path = os.path.join(dsi_metrics_dir, output_name)
                if os.path.isfile(new_metric_path):
                    os.remove(new_metric_path)
                os.rename(old_metric_path, new_metric_path)
            
            # Compress legacy .nii files in DSI_studio when the old output
            # layout is explicitly requested.
            if 'nii_gz' in opts and output_name.endswith('.nii'):
                old_nifti_path = os.path.join(dsi_metrics_dir, output_name)
                new_nifti_path = os.path.join(dsi_metrics_dir, output_name + '.gz')

                if os.path.isfile(new_nifti_path):
                    os.remove(new_nifti_path)

                with open(old_nifti_path, 'rb') as f_in:
                    with gzip.open(new_nifti_path, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)

                os.remove(old_nifti_path)
