"""
Created on 18/11/2020

@author: Marc Schneider
AG Neuroimaging and Neuroengineering of Experimental Stroke
Department of Neurology, University Hospital Cologne

This script runs every needed script for all (pre-)processing and registration
steps. The data needs to be ordered like after Bruker2NIfTI conversion:
project_folder/days/groups/subjects/.
For the script to work, it needs to be placed within the /bin folder of AIDAmri.

Example:
python batchProc.py -i /Volumes/Desktop/MRI/proc_data -t anat dwi func t2map
"""

import glob
import os
import fnmatch
import shutil
from pathlib import Path
import nibabel as nii
import concurrent.futures
import functools
import subprocess
from tqdm import tqdm
import multiprocessing
import logging
import shlex
import time


def findData(projectPath, sessions, data_types):
    if not data_types:
        data_types = ["anat", "dwi", "func", "t2map"]
    # This function screens all existing paths. Within these paths, this function collects all subject
    # folders, which are all folders that are not named 'Physio'.
    full_path_list = os.listdir(projectPath)
    all_wanted_paths, anat_files, dwi_files, func_files, t2map_files = [], [], [], [], []

    # collect ses paths
    for path in full_path_list:
        if path.startswith("sub-"):
            sub_root = os.path.join(projectPath, path)
            wanted_paths = os.listdir(sub_root)
            wanted_paths = [
                os.path.join(sub_root, wp)
                for wp in wanted_paths
                if wp.startswith("ses-")
            ]
            all_wanted_paths.extend(wanted_paths)

    # filter sessions (exact match on path components)
    if sessions:
        wanted = {f"ses-{s}" for s in sessions}
        matching_paths = []
        for p in all_wanted_paths:
            parts = os.path.normpath(p).split(os.sep)
            if any(part in wanted for part in parts):
                matching_paths.append(p)
        all_wanted_paths = matching_paths

    # collect datatype folders
    for path in all_wanted_paths:
        for sub_dir in os.listdir(path):
            if sub_dir == "anat" and "anat" in data_types:
                anat_files.append(os.path.join(path, sub_dir))

            elif sub_dir == "dwi" and "dwi" in data_types:
                dwi_files.append(os.path.join(path, sub_dir))

            elif sub_dir == "func" and "func" in data_types:
                func_files.append(os.path.join(path, sub_dir))

            elif sub_dir == "t2map" and "t2map" in data_types:
                t2map_files.append(os.path.join(path, sub_dir))

    return {"anat": anat_files, "dwi": dwi_files, "func": func_files, "t2map": t2map_files}

def _get_arg_after(flags, argv):
    for f in flags:
        if f in argv:
            i = argv.index(f)
            if i + 1 < len(argv):
                return argv[i + 1]
    return None

def _log_base_from_input(input_path: str) -> str:
    # If input is a dir -> log in that dir
    # If input is a file -> log in its parent dir
    return input_path if os.path.isdir(input_path) else os.path.dirname(input_path)

def run_subprocess(command, datatype, step, anat_process=False):
    timeout = 3600
    command_args = shlex.split(command)

    inp = _get_arg_after(["-i", "--input", "--input_file"], command_args)
    if inp is None:
        inp = next((a for a in reversed(command_args)
                    if a.endswith(".nii") or a.endswith(".nii.gz")), command_args[-1])

    base = _log_base_from_input(inp)

    # default location
    log_file = os.path.join(base, f"{step}.log")

    # special case: anat/process wants different filenames
    if datatype == "anat" and step == "process":
        log_name = f"{step}.log" if anat_process else f"{step}_par.log"
        log_file = os.path.join(base, log_name)

    #Determine sub / ses
    normalized_path = os.path.normpath(inp)
    directories = normalized_path.split(os.path.sep)
    sub = next((d for d in directories if d.startswith("sub-")), "sub-UNKNOWN")
    ses = next((d for d in directories if d.startswith("ses-")), "ses-UNKNOWN")

    try:
        logging.info(f"Running command: {command}.\nCheck {log_file} for further information.")
        if os.path.exists(log_file):
            os.remove(log_file)    
        with open(log_file, 'w') as outfile:
            time.sleep(2) # make sure logging file is created before starting the subprocess
            result = subprocess.run(command_args, stdout=outfile, stderr=outfile, text=True, timeout=timeout)
            if result.returncode != 0:
                return sub,ses,datatype,step
            else:
                return 0
    except subprocess.TimeoutExpired:
        logging.error(f'Timeout expired for command: {command_args}')
        return sub,ses,datatype,step
    except Exception as e:
        logging.error(f'Error while executing the command: {command_args} Errorcode: {str(e)}')
        raise
    

def executeScripts(currentPath_wData, dataFormat, step, cfg, stc=False):
    # For every datatype (T2w, fMRI, DTI), go in all days/group/subjects folders
    # and execute the respective (pre-)processing/registration-scripts.
    # If a certain file does not exist, a note will be created in the errorList.
    # cwd should contain the path of the /bin folder (the user needs to navigate to the /bin folder before executing this script)
    #KEEP IN MIND DUE TO PARALLEL COMPUTING NO ERRORS IN THIS FUNCTION WILL BE PRINTED OUT => GREY ZONE
    errorList = [];
    message = '';
    cwd = str(Path(__file__).resolve().parent)
    currentPath_wData = Path(currentPath_wData)
    # currentPath_wData = projectfolder/sub/ses/dataFormat (e.g. anat, func, dwi)
    if os.path.isdir(currentPath_wData):
        if dataFormat == 'anat':
            if step == "preprocess":
                os.chdir(os.path.join(cwd, '2.1_T2PreProcessing'))
                currentFile = list(currentPath_wData.glob("*T2w.nii.gz"))
                if len(currentFile) > 0:
                    command = f'python preProcessing_T2.py -i {currentFile[0]}'
                    # Bias: skip + method (T2 has both)
                    if cfg.get("t2_bias_skip") is not None:
                        command += f' -b {cfg["t2_bias_skip"]}'

                    # Only set method if specified
                    if cfg.get("t2_bias_method"):
                        command += f' --bias_method {cfg["t2_bias_method"]}'
                    # fallback to old biasfieldcorr
                    elif cfg.get("biasfieldcorr"):
                        command += f' --bias_method {cfg["biasfieldcorr"]}'

                    # BET-Parameter
                    if cfg.get("t2_frac") is not None:
                        command += f' -f {cfg["t2_frac"]}'
                    if cfg.get("t2_radius") is not None:
                        command += f' -r {cfg["t2_radius"]}'
                    if cfg.get("t2_vertical_gradient") is not None:
                        command += f' -g {cfg["t2_vertical_gradient"]}'
                    if cfg.get("t2_center") is not None:
                        cx, cy, cz = cfg["t2_center"]
                        command += f' -c {cx} {cy} {cz}'

                    # bet4animal
                    if cfg.get("bet4animal"):
                        command += ' --use_bet4animal'

                    result = run_subprocess(command, dataFormat, step)
                    if result != 0:
                        errorList.append(result)
                else:
                    message = f'Could not find *T2w.nii.gz in {str(currentPath_wData)}'
                    logging.error(message)
                    errorList.append(message)
                os.chdir(cwd)

            elif step == "registration":
                os.chdir(os.path.join(cwd, '2.1_T2PreProcessing'))
                currentFile = list(currentPath_wData.glob("*Bet.nii.gz"))
                if len(currentFile) > 0:
                    command = f'python registration_T2.py -i {currentFile[0]}'
                    r1 = run_subprocess(f'python registration_T2.py -i {currentFile[0]}', dataFormat, step)
                    if r1 != 0:
                        errorList.append(r1)
                    command = f'python t2_value_extraction.py -i {currentFile[0]}'
                    r2 = run_subprocess(f'python t2_value_extraction.py -i {currentFile[0]}', dataFormat, step)
                    if r2 != 0:
                        errorList.append(r2)
                else:
                    message = f'Could not find *Bet.nii.gz in {str(currentPath_wData)}'
                    logging.error(message)
                    errorList.append(message)
                os.chdir(cwd)

            elif step == "process":
                has_stroke_mask = any(currentPath_wData.glob("**/*Stroke_mask.nii.gz"))
                if not has_stroke_mask:
                    message = f"No stroke mask found for {str(currentPath_wData)}, proceeding without mask."
                    logging.info(message)  #write in log-file
                    #print(message, flush=True)
                    return 0
                os.chdir(os.path.join(cwd, '3.1_T2Processing'))
                command = f'python getIncidenceSize_par.py -i {str(currentPath_wData)}'
                result = run_subprocess(command, dataFormat, step)
                if isinstance(result, tuple) and len(result) == 4:
                    os.chdir(os.path.join(cwd, '3.1_T2Processing'))

                    r_par = run_subprocess(f'python getIncidenceSize_par.py -i {str(currentPath_wData)}',
                                           dataFormat, step)

                    if r_par != 0:
                        errorList.append(r_par)
                        r_ser = run_subprocess(f'python getIncidenceSize.py -i {str(currentPath_wData)}',
                                               dataFormat, step, anat_process=True)
                        if r_ser != 0:
                            errorList.append(r_ser)

                    os.chdir(cwd)


        elif dataFormat == 'func':
            if step == "preprocess":
                os.chdir(os.path.join(cwd, '2.3_fMRIPreProcessing'))
                currentFile = list(currentPath_wData.glob("*EPI.nii.gz"))
                if len(currentFile)>0:
                    command = f'python preProcessing_fMRI.py -i {currentFile[0]}'
                    result = run_subprocess(command,dataFormat,step)
                    if result != 0:
                        errorList.append(result)
                else:
                    message = f'Could not find *EPI.nii.gz in {str(currentPath_wData)}';
                    logging.error(message)
                    errorList.append(message)
                os.chdir(cwd)
            elif step == "registration":
                os.chdir(os.path.join(cwd, '2.3_fMRIPreProcessing'))
                currentFile = list(currentPath_wData.glob("*SmoothBet.nii.gz"))
                if len(currentFile)>0:
                    command = f'python registration_rsfMRI.py -i {currentFile[0]}'
                    result = run_subprocess(command,dataFormat,step)
                    if result != 0:
                        errorList.append(result)
                else:
                    message = f'Could not find *SmoothBet.nii.gz in {str(currentPath_wData)}';
                    logging.error(message)
                    errorList.append(message)
                os.chdir(cwd)
            elif step == "process":
                currentFile = list(currentPath_wData.glob("*EPI.nii.gz"))
                if len(currentFile)>0:
                    os.chdir(os.path.join(cwd, '3.3_fMRIActivity'))
                    command = f'python process_fMRI.py -i {currentFile[0]} -stc {stc}'
                    result = run_subprocess(command,dataFormat,step)
                    if result != 0:
                        errorList.append(result)
                    os.chdir(cwd)
        elif dataFormat == 't2map':
            if step == "preprocess":
                os.chdir(os.path.join(cwd, '4.1_T2mapPreProcessing'))
                currentFile = list(currentPath_wData.glob("*MEMS.nii.gz"))
                if len(currentFile)>0:
                    command = f'python preProcessing_T2MAP.py -i {currentFile[0]}'
                    result = run_subprocess(command,dataFormat,step)
                    if result != 0:
                        errorList.append(result)
                else:
                    message = f'Could not find *MEMS.nii.gz in {str(currentPath_wData)}';
                    logging.error(message)
                    errorList.append(message)
                os.chdir(cwd)
            elif step == "registration":
                os.chdir(os.path.join(cwd, '4.1_T2mapPreProcessing'))
                currentFile = list(currentPath_wData.glob("*SmoothMicoBet.nii.gz"))
                if len(currentFile)>0:
                    command = f'python registration_T2MAP.py -i {currentFile[0]}'
                    result = run_subprocess(command,dataFormat,step)
                    if result != 0:
                        errorList.append(result)
                else:
                    message = f'Could not find *SmoothMicoBet.nii.gz in {str(currentPath_wData)}';
                    print(message)
                    errorList.append(message)
                os.chdir(cwd)
            elif step == "process":
                currentFile = list(currentPath_wData.glob("*T2w_MAP.nii.gz"))
                if len(currentFile)>0:
                    command = f'python t2map_data_extract.py -i {currentFile[0]}'
                    result = run_subprocess(command,dataFormat,step)
                    if result != 0:
                        errorList.append(result)
                else:
                    message = f'Could not find *T2w_MAP.nii.gz in {str(currentPath_wData)}';
                    logging.error(message)
                    errorList.append(message)
                os.chdir(cwd)
        elif dataFormat == 'dwi':
            if step == "preprocess":
                os.chdir(os.path.join(cwd, '2.2_DTIPreProcessing'))
                currentFile = list(currentPath_wData.glob("*dwi.nii.gz"))
                if len(currentFile) > 0:
                    command = f'python preProcessing_DTI.py -i {currentFile[0]}'

                    # DWI BET parameter (only append if set, otherwise use script defaults)
                    if cfg.get("dwi_frac") is not None:
                        command += f' -f {cfg["dwi_frac"]}'
                    if cfg.get("dwi_radius") is not None:
                        command += f' -r {cfg["dwi_radius"]}'
                    if cfg.get("dwi_vertical_gradient") is not None:
                        command += f' -g {cfg["dwi_vertical_gradient"]}'

                    # Bias field (DWI: none/micro/ANTS)
                    # dwi_bias_method with choices [“none”,‘mico’,“ants”], default="none"
                    if cfg.get("dwi_bias_method", "none") != "none":
                        command += f' -b {cfg["dwi_bias_method"]}'

                    # Denoiser
                    if cfg.get("dwi_denoiser"):
                        command += f' --denoiser {cfg["dwi_denoiser"]}'

                    # Flags
                    if cfg.get("bet4animal"):
                        command += ' --use_bet4animal'

                    if cfg.get("dwi_average_b0"):
                        command += ' --average_b0'

                    if cfg.get("dwi_skip_min"):
                        command += ' --skip_min'

                    if cfg.get("dwi_deoblique"):
                        command += " --deoblique"

                    result = run_subprocess(command, dataFormat, step)
                    if result != 0:
                        errorList.append(result)
                else:
                    message = f'Could not find *dwi.nii.gz in {str(currentPath_wData)}';
                    logging.error(message)
                    errorList.append(message)
                os.chdir(cwd)
            elif step == "registration":
                os.chdir(os.path.join(cwd, '2.2_DTIPreProcessing'))
                currentFile = list(currentPath_wData.glob("*Smooth*Bet.nii.gz"))
                if len(currentFile)>0:
                    command = f'python registration_DTI.py -i {currentFile[0]}'
                    result = run_subprocess(command,dataFormat,step)
                    if result != 0:
                        errorList.append(result)
                else:
                    message = f'Could not find *Smooth*Bet.nii.gz in {currentPath_wData}';
                    logging.error(message)
                    errorList.append(message)
                os.chdir(cwd)
            elif step == "process":
                currentFile = list(currentPath_wData.glob("*dwi.nii.gz"))
                if cfg.get("dwi_denoiser") == "patch2self" or cfg.get("denoiser") == "patch2self":
                    currentFile = list(currentPath_wData.glob("*Patch2SelfDenoised.nii.gz"))
                # Appends optional (fa0, nii_gz) flags to DTI main process if passed
                if len(currentFile)>0:
                    # Pull values from cfg (with defaults)
                    track_param = cfg.get("dsi_track_param", cfg.get("track_param", "default"))
                    recon_method = cfg.get("dsi_recon_method", cfg.get("recon_method", "dti"))
                    vivo = cfg.get("dsi_vivo", cfg.get("vivo", "in_vivo"))
                    make_iso = cfg.get("dsi_make_isotropic", cfg.get("make_isotropic", 0))
                    flip_y = 1 if cfg.get("dsi_flip_image_y", cfg.get("flip_image_y", False)) else 0
                    template_val = str(cfg.get("dsi_template", cfg.get("template", "mouse")))
                    thread_count = cfg.get("num_processes", 1)
                    legacy = bool(cfg.get("dsi_legacy", cfg.get("legacy", False)))
                    no_mcf = bool(cfg.get("dsi_no_mcf", cfg.get("no_mcf", False)))

                    cli_str = (
                        f'dsi_main.py -i {currentFile[0]} '
                        f'-t {track_param} -r {recon_method} -v {vivo} -m {make_iso} '
                        f'-y {flip_y} -template {template_val} -thread_count {thread_count}'
                    )
                    if legacy:
                        cli_str += ' -l'
                    if no_mcf:
                        cli_str += ' -nomcf'

                    os.chdir(cwd + '/3.2_DTIConnectivity')
                    command = f'python {cli_str}'
                    result = run_subprocess(command,dataFormat,step)
                    if result != 0:
                        errorList.append(result)
                os.chdir(cwd)
        else:
            message = 'The data folders'' names do not match anat, dwi, func or t2map';
            logging.error(message);
            errorList.append(message)
    else:
        message = 'The folder ' + dataFormat + ' does not exist in ' + str(currentPath_wData)
        logging.error(message)
        errorList.append(message)
    
    if errorList:
        return errorList
    else:
        return 0
    
 
def find(pattern, path):
    # This function finds all files with a specified fragment within
    # the given path
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "Batch processing of all data (AIDAmri). "
            "Runs preprocessing, registration and processing steps for T2, DWI, fMRI and T2map.\n\n"
            "Example:\n"
            "python batchProc.py -i /path/to/proc_data -t anat dwi "
            "--t2-frac 0.15 --t2-bias-method mico "
            "--dwi-denoiser patch2self "
            "--dsi-template mouse"
        ),
        formatter_class=argparse.RawTextHelpFormatter
    )

    # ============================================================
    # REQUIRED
    # ============================================================
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-i", "--input",
        required=True,
        help="Path to the parent project folder (e.g. proc_data)"
    )

    # ============================================================
    # GLOBAL / BATCH CONTROL
    # ============================================================
    batch = parser.add_argument_group("batch control")
    batch.add_argument(
        "-s", "--sessions",
        nargs="+",
        help="Process only selected sessions (e.g. Baseline P7 P14)"
    )
    batch.add_argument(
        "-t", "--data-types",
        nargs="+",
        help="Data types to process (anat, dwi, func, t2map). Default: all"
    )
    batch.add_argument(
        "-d", "--debug-steps",
        dest="debug_steps",
        nargs="+",
        help="Processing steps to run (preprocess registration process). Default: all"
    )
    batch.add_argument(
        "--slice-time-correction",
        action="store_true",
        help="Enable slice time correction for fMRI"
    )
    # ============================================================
    # CPU / PARALLELIZATION
    # ============================================================
    cpu = parser.add_argument_group("cpu / parallelization")
    cpu.add_argument(
        "-c", "--cpu-cores",
        default="Half",
        choices=["Min", "Half", "Max"],
        help="CPU usage preset (Min, Half, Max)"
    )
    cpu.add_argument(
        "-e", "--expert-cpu",
        type=int,
        help="Explicit number of parallel processes"
    )

    # ============================================================
    # T2 PREPROCESSING (preProcessing_T2.py)
    # ============================================================
    t2 = parser.add_argument_group("T2 preprocessing (preProcessing_T2.py)")
    t2.add_argument(
        "--t2-bias-method",
        choices=["mico", "ants"],
        help="Bias field correction method for T2 (mico or ants)"
    )
    t2.add_argument(
        "--t2-bias-skip",
        type=float,
        help="Skip T2 bias field correction (1 = skip, 0 = apply)"
    )
    t2.add_argument(
        "--t2-frac",
        type=float,
        help="BET fractional intensity threshold (default in script: 0.15)"
    )
    t2.add_argument(
        "--t2-radius",
        type=int,
        help="BET head radius in mm (default in script: 45)"
    )
    t2.add_argument(
        "--t2-vertical-gradient",
        type=float,
        help="BET vertical gradient (default in script: 0.0)"
    )
    t2.add_argument(
        "--t2-center",
        nargs=3,
        type=int,
        metavar=("X", "Y", "Z"),
        help="BET center in voxel coordinates"
    )

    # ============================================================
    # DWI PREPROCESSING (preProcessing_DTI.py)
    # ============================================================
    dwi = parser.add_argument_group("DWI preprocessing (preProcessing_DTI.py)")
    dwi.add_argument(
        "--dwi-denoiser",
        choices=["patch2self"],
        help="DWI denoising method"
    )
    dwi.add_argument(
        "--dwi-average-b0",
        action="store_true",
        help="Average b0 volumes before DWI processing"
    )
    dwi.add_argument(
        "--dwi-skip-min",
        action="store_true",
        help="Skip minimum intensity projection step"
    )
    dwi.add_argument(
        "--dwi-frac",
        type=float,
        help="BET fractional intensity threshold for DWI"
    )
    dwi.add_argument(
        "--dwi-radius",
        type=int,
        help="BET head radius (mm) for DWI"
    )
    dwi.add_argument(
        "--dwi-vertical-gradient",
        type=float,
        help="BET vertical gradient for DWI"
    )
    dwi.add_argument(
        "--dwi-bias-method",
        choices=["none", "mico", "ants"],
        default="none",
        help="Bias field correction for DWI: none|mico|ants (maps to preProcessing_DTI.py --biasfieldcorr)"
    )
    dwi.add_argument(
        '--dwi-deoblique',
        help='Deoblique input using AFNI 3dWarp -deoblique',
        action='store_true'
    )
    # ============================================================
    # BET / ANIMAL-SPECIFIC
    # ============================================================
    bet = parser.add_argument_group("BET / animal settings")
    bet.add_argument(
        "--bet4animal",
        action="store_true",
        help="Use BET tuned for animal brains (bet4animal)"
    )

    # ============================================================
    # DSI STUDIO / TRACTOGRAPHY (dsi_main.py)
    # ============================================================
    dsi = parser.add_argument_group("DSI Studio / tractography (dsi_main.py)")
    dsi.add_argument(
        "--dsi-recon-method",
        default="dti",
        choices=["dti", "gqi"],
        help="DSI reconstruction method"
    )
    dsi.add_argument(
        "--dsi-vivo",
        default="in_vivo",
        choices=["in_vivo", "ex_vivo"],
        help="In vivo or ex vivo data (controls sampling length)"
    )
    dsi.add_argument(
        "--dsi-make-isotropic",
        type=float,
        default=0.0,
        help="Voxel size (mm) for isotropic resampling (0 = off, auto = header)"
    )
    dsi.add_argument(
        "--dsi-flip-image-y",
        action="store_true",
        help="Flip image in Y direction before DSI processing"
    )
    dsi.add_argument(
        "--dsi-template",
        default="mouse",
        help="DSI template (mouse, rat or numeric ID)"
    )
    dsi.add_argument(
        "--dsi-track-param",
        default="default",
        help="Tracking parameter preset (default, mouse, rat, aida_optimized)"
    )
    dsi.add_argument(
        "--dsi-no-mcf",
        action="store_true",
        help="Skip slice-wise MCFLIRT motion correction"
    )
    dsi.add_argument(
        "--dsi-legacy",
        action="store_true",
        help="Enable legacy .fib.gz / .src.gz support"
    )

    args = parser.parse_args()

    pathToData = args.input
    sessions = args.sessions
    
    #configurate the logging module
    log_file_path = os.path.join(pathToData, "batchproc_log.txt")
    logging.basicConfig(filename=log_file_path, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', force =True)

    stc = args.slice_time_correction

    if args.data_types is None:
        data_types = ["anat", "dwi", "func", "t2map"]
    else:
        data_types = args.data_types

    if args.debug_steps is None:
        steps = ["preprocess","registration","process"]
    else:
        steps = args.debug_steps
    
    print('Entered information:')
    print(pathToData)
    print('data_types %s' % data_types)
    print('Slice time correction [%s]' % stc)
    print('Steps %s' % steps)
    print()

    all_files = findData(pathToData, sessions, data_types)

    num_processes = 1

    if args.cpu_cores.upper() == "MIN":
        num_processes = 1
    elif args.cpu_cores.upper() == "HALF":
        num_processes = int(multiprocessing.cpu_count() / 2)
    elif args.cpu_cores.upper() == "MAX":
        num_processes = multiprocessing.cpu_count()

    print(args)
    
    if args.expert_cpu:
        num_processes = int(args.expert_cpu)
    
    print(f"Running with {num_processes} parallel processes!")
    logging.info(f"Using {num_processes} CPUs for the parallelization")
    logging.info(f"Processing following datasets:\n{all_files}")
    # turns argparse.Namespace into a dict
    cfg = vars(args)
    cfg["num_processes"] = num_processes
    logging.info(
        "DSI settings: recon=%s vivo=%s make_isotropic=%s flip_y=%s template=%s track_param=%s no_mcf=%s legacy=%s",
        args.dsi_recon_method,
        args.dsi_vivo,
        args.dsi_make_isotropic,
        args.dsi_flip_image_y,
        args.dsi_template,
        args.dsi_track_param,
        args.dsi_no_mcf,
        args.dsi_legacy)

    for key, value in all_files.items():
        if value:
            error_list_all = []
            print()
            print(f"Entered {key} data: \n{value}")
            print()
            print(f"\n{key} processing \33[5m...\33[0m (wait!)") 
            print()
            for step in steps: 
                error_list_step = []
                progress_bar = tqdm(total=len(value), desc=f"{step} {key} data")
                with concurrent.futures.ProcessPoolExecutor(max_workers=num_processes) as executor:
                    futures = [executor.submit(executeScripts, path, key, step, cfg, stc) for path in value]

                    # --- collect errors robustly ---
                    flat_errors_step = []

                    for future in concurrent.futures.as_completed(futures):
                        progress_bar.update(1)

                        res = future.result()

                        # normalize result into a flat list
                        if res == 0 or res is None:
                            continue

                        if isinstance(res, list):
                            flat_errors_step.extend(res)
                        else:
                            flat_errors_step.append(res)

                    concurrent.futures.wait(futures)
                    progress_bar.close()

                    # keep a per-step and per-datatype summary
                    if not flat_errors_step:
                        print(f"{key} {step}  \033[0;30;42m COMPLETED \33[0m")
                    else:
                        print(f"{key} {step}  \033[0;30;41m INCOMPLETE \33[0m")
                        error_list_all.extend(flat_errors_step)

                    logging.info(f"{key} {step} processing completed")

            logging.error(f"Following errors were occurring: {error_list_all}")
            logging.info(f"{key} processing completed")

            if not error_list_all:
                print(f"\n{key} processing \033[0;30;42m COMPLETED \33[0m")
            else:
                print(f"\n{key} processing \033[0;30;41m INCOMPLETE \33[0m")
                print()
                for err in error_list_all:
                    if isinstance(err, tuple) and len(err) == 4:
                        sub, ses, dtype, stepname = err
                        print(
                            f"Error in sub: {sub} in session: {ses} in datatype: {dtype} and step: {stepname}. Check log.")
                    else:
                        # strings or unexpected types
                        print(f"Error: {err}")

            
                

 
