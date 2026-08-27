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

import argparse
import ast
import os
import fnmatch
import csv
from calendar import month_name
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo
import concurrent.futures
import subprocess
from tqdm import tqdm
import multiprocessing
import logging
import shlex
import time
import sys
import shutil

from common.artifact_manifest import OutputTracker

FATAL_LIP_HEADER_EXIT_CODE = 86
REPORT_TIMEZONE = ZoneInfo("Europe/Berlin")
AIDAMRI_GIT_INFO_SOURCES = [
    "/aida/build/AIDAmri_git_information.txt",
    "/aida/DATA/AIDAmri_git_information.txt",
]
AIDAMRI_GIT_INFO_FILENAME = "AIDAmri_git_information.txt"


class BatchProcFormatter(logging.Formatter):
    def formatTime(self, record, datefmt=None):
        timestamp = datetime.fromtimestamp(record.created, REPORT_TIMEZONE)
        return (
            f"{timestamp.day:02d} {month_name[timestamp.month]} {timestamp.year} "
            f"{timestamp:%H:%M:%S} {timestamp.tzname()}"
        )


def configure_logging(log_file_path):
    handler = logging.FileHandler(log_file_path)
    handler.setFormatter(
        BatchProcFormatter("%(asctime)s - %(levelname)s - %(message)s")
    )
    logging.basicConfig(level=logging.INFO, handlers=[handler], force=True)


def copy_aidamri_git_information_to_proc(proc_dir):
    """Best-effort copy of build provenance into the processed data folder."""
    target = os.path.join(proc_dir, AIDAMRI_GIT_INFO_FILENAME)
    for source in AIDAMRI_GIT_INFO_SOURCES:
        if not os.path.isfile(source):
            continue
        try:
            if os.path.abspath(source) != os.path.abspath(target):
                shutil.copyfile(source, target)
        except OSError:
            pass
        return


def findData(projectPath, sessions, data_types):
    if not data_types:
        data_types = ["anat", "dwi", "func", "t2map"]
    # This function screens all existing paths. Within these paths, this function collects all subject
    # folders, which are all folders that are not named 'Physio'.
    full_path_list = sorted(os.listdir(projectPath))
    all_wanted_paths, anat_files, dwi_files, func_files, t2map_files = [], [], [], [], []

    # collect ses paths
    for path in full_path_list:
        if path.startswith("sub-"):
            sub_root = os.path.join(projectPath, path)
            wanted_paths = sorted(os.listdir(sub_root))
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
        for sub_dir in sorted(os.listdir(path)):
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

def _quote(value) -> str:
    '''
    Safely quote a value for shell command usage, handling spaces and special characters.
    '''
    return shlex.quote(str(value))

def run_subprocess(command, datatype, step, anat_process=False):
    timeout = 5400 #timeout (sec) for subprocess
    command_args = shlex.split(command)

    inp = _get_arg_after(["-i", "--input", "--input-file"], command_args)
    if inp is None:
        inp = next((a for a in reversed(command_args)
                    if a.endswith(".nii") or a.endswith(".nii.gz")), command_args[-1])

    base = _log_base_from_input(inp)
    #starting aidamri output tracker
    output_tracker = None
    if datatype in {"anat", "dwi", "func", "t2map"}:
        output_tracker = OutputTracker.start(base, datatype, step)

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
        with open(log_file, 'w') as outfile:
            time.sleep(2) # make sure logging file is created before starting the subprocess
            child_env = os.environ.copy()
            child_env["AIDAMRI_DISABLE_SCRIPT_LOG"] = "1"
            # dsi_main.py can create its own process.log during interactive
            # runs. Disable that side log here because batchProc.py already
            # captures stdout/stderr into the step-specific batch log.
            if any(arg.endswith("dsi_main.py") for arg in command_args):
                child_env["AIDAMRI_DISABLE_PROCESS_LOG"] = "1"
            child_env["AIDAMRI_DISABLE_SPINNER"] = "1"
            result = subprocess.run(
                command_args,
                stdout=outfile,
                stderr=outfile,
                text=True,
                timeout=timeout,
                env=child_env,
            )
            if result.returncode != 0:
                if (
                    result.returncode == FATAL_LIP_HEADER_EXIT_CODE
                    and datatype in {"anat", "dwi"}
                    and step == "preprocess"
                ):
                    raise RuntimeError(
                        f"Fatal header check failure in {inp}. Expected LIP orientation."
                    )
                return sub,ses,datatype,step
            else:
                return 0
    except subprocess.TimeoutExpired:
        logging.error(f'Timeout expired for command: {command_args}')
        return sub,ses,datatype,step
    except Exception as e:
        logging.error(f'Error while executing the command: {command_args} Errorcode: {str(e)}')
        raise
    finally:
        if output_tracker is not None:
            output_tracker.finalize()
    

def executeScripts(currentPath_wData, dataFormat, step, cfg):
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
                currentFile = sorted(currentPath_wData.glob("*T2w.nii.gz"))
                if len(currentFile) > 0:
                    command = f'python preProcessing_T2.py -i {_quote(currentFile[0])}'

                    # Bias field correction for T2: none | mico | ants
                    if cfg.get("t2_bias_method") is not None:
                        command += f' -b {cfg["t2_bias_method"]}'

                    if cfg.get("t2_bet") is not None:
                        command += f' --bet {cfg["t2_bet"]}'

                    # BET-Parameter
                    if cfg.get("t2_frac") is not None:
                        command += f' -f {cfg["t2_frac"]}'
                    if cfg.get("t2_radius") is not None:
                        command += f' -r {cfg["t2_radius"]}'
                    if cfg.get("t2_gradient") is not None:
                        command += f' -g {cfg["t2_gradient"]}'
                    if cfg.get("t2_center") is not None:
                        cx, cy, cz = cfg["t2_center"]
                        command += f' -c {cx} {cy} {cz}'

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
                currentFile = sorted(currentPath_wData.glob("*Bet.nii.gz"))
                if len(currentFile) > 0:
                    r1 = run_subprocess(f'python registration_T2.py -i {_quote(currentFile[0])}', dataFormat, step)
                    if r1 != 0:
                        errorList.append(r1)
                    r2 = run_subprocess(f'python t2_value_extraction.py -i {_quote(currentFile[0])}', dataFormat, step)
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
                if cfg.get("t2_incidence_script") == "detailed":
                    command = f'python getIncidenceSize.py -i {_quote(currentPath_wData)}'
                    result = run_subprocess(command, dataFormat, step, anat_process=True)
                else:
                    command = f'python getIncidenceSize_par.py -i {_quote(currentPath_wData)}'
                    result = run_subprocess(command, dataFormat, step)

                if result != 0:
                    errorList.append(result)

                os.chdir(cwd)


        elif dataFormat == 'func':
            if step == "preprocess":
                os.chdir(os.path.join(cwd, '2.3_fMRIPreProcessing'))
                currentFile = sorted(currentPath_wData.glob("*EPI.nii.gz"))
                if len(currentFile)>0:
                    command = f'python preProcessing_fMRI.py -i {_quote(currentFile[0])}'
                    if cfg.get("func_bias_method") is not None:
                        command += f' -b {cfg["func_bias_method"]}'
                    if cfg.get("func_skip_smoothing") is True:
                        command += ' --skip-smoothing'
                    if cfg.get("func_bet") is not None:
                        command += f' --bet {cfg["func_bet"]}'
                    if cfg.get("func_frac") is not None:
                        command += f' -f {cfg["func_frac"]}'
                    if cfg.get("func_radius") is not None:
                        command += f' -r {cfg["func_radius"]}'
                    if cfg.get("func_gradient") is not None:
                        command += f' -g {cfg["func_gradient"]}'
                    if cfg.get("func_center") is not None:
                        cx, cy, cz = cfg["func_center"]
                        command += f' -c {cx} {cy} {cz}'
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
                currentFile = sorted(currentPath_wData.glob("*Bet.nii.gz"))
                if len(currentFile)>0:
                    command = f'python registration_rsfMRI.py -i {_quote(currentFile[0])}'
                    if cfg.get("func_atlas_mask_t2") is True:
                        command += " --atlas-mask-t2"
                    result = run_subprocess(command,dataFormat,step)
                    if result != 0:
                        errorList.append(result)
                else:
                    message = f'Could not find *Bet.nii.gz in {str(currentPath_wData)}';
                    logging.error(message)
                    errorList.append(message)
                os.chdir(cwd)
            elif step == "process":
                currentFile = sorted(currentPath_wData.glob("*EPI.nii.gz"))
                if len(currentFile)>0:
                    os.chdir(os.path.join(cwd, '3.3_fMRIActivity'))
                    command = f'python process_fMRI.py -i {_quote(currentFile[0])} --bet {cfg["func_bet"]}'
                    if cfg.get("func_stc") is True:
                        command += ' -stc'
                    if cfg.get("func_frac") is not None:
                        command += f' --bet-frac {cfg["func_frac"]}'
                    if cfg.get("func_radius") is not None:
                        command += f' --bet-radius {cfg["func_radius"]}'
                    if cfg.get("func_gradient") is not None:
                        command += f' --bet-gradient {cfg["func_gradient"]}'
                    if cfg.get("func_center") is not None:
                        cx, cy, cz = cfg["func_center"]
                        command += f' -ctr {cx} {cy} {cz}'
                    result = run_subprocess(command,dataFormat,step)
                    if result != 0:
                        errorList.append(result)
                    os.chdir(cwd)
                else:
                    message = f'Could not find *EPI.nii.gz in {str(currentPath_wData)}';
                    logging.error(message)
                    errorList.append(message)
        elif dataFormat == 't2map':
            if step == "preprocess":
                os.chdir(os.path.join(cwd, '4.1_T2mapPreProcessing'))
                currentFile = sorted(currentPath_wData.glob("*MEMS.nii.gz"))
                if len(currentFile)>0:
                    command = f'python preProcessing_T2MAP.py -i {_quote(currentFile[0])}'
                    if cfg.get("t2map_bias_method") is not None:
                        command += f' -b {cfg["t2map_bias_method"]}'
                    if cfg.get("t2map_bet") is not None:
                        command += f' --bet {cfg["t2map_bet"]}'
                    if cfg.get("t2map_frac") is not None:
                        command += f' -f {cfg["t2map_frac"]}'
                    if cfg.get("t2map_radius") is not None:
                        command += f' -r {cfg["t2map_radius"]}'
                    if cfg.get("t2map_gradient") is not None:
                        command += f' -g {cfg["t2map_gradient"]}'
                    if cfg.get("t2map_center") is not None:
                        cx, cy, cz = cfg["t2map_center"]
                        command += f' -c {cx} {cy} {cz}'
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
                currentFile = sorted(
                    currentPath_wData.glob("*Bet.nii.gz"),
                    key=lambda path: path.stat().st_mtime,
                    reverse=True,
                )
                if len(currentFile)>0:
                    command = f'python registration_T2MAP.py -i {_quote(currentFile[0])}'
                    result = run_subprocess(command,dataFormat,step)
                    if result != 0:
                        errorList.append(result)
                else:
                    message = f'Could not find *Bet.nii.gz in {str(currentPath_wData)}';
                    print(message)
                    errorList.append(message)
                os.chdir(cwd)
            elif step == "process":
                os.chdir(os.path.join(cwd, '4.1_T2mapPreProcessing'))
                currentFile = sorted(currentPath_wData.glob("*T2w_MAP.nii.gz"))
                if len(currentFile)>0:
                    command = f'python t2map_data_extract.py -i {_quote(currentFile[0])}'
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
                currentFile = sorted(currentPath_wData.glob("*dwi.nii.gz"))
                if len(currentFile) > 0:
                    command = f'python preProcessing_DTI.py -i {_quote(currentFile[0])}'

                    # DWI BET parameter (only append if set, otherwise use script defaults)
                    if cfg.get("dwi_frac") is not None:
                        command += f' -f {cfg["dwi_frac"]}'
                    if cfg.get("dwi_radius") is not None:
                        command += f' -r {cfg["dwi_radius"]}'
                    if cfg.get("dwi_gradient") is not None:
                        command += f' -g {cfg["dwi_gradient"]}'

                    # Bias field
                    if cfg.get("dwi_bias_method") is not None:
                        command += f' -b {cfg["dwi_bias_method"]}'

                    # Denoiser
                    if cfg.get("dwi_denoiser") is not None:
                        command += f' --denoiser {cfg["dwi_denoiser"]}'

                    if cfg.get("dwi_bet") is not None:
                        command += f' --bet {cfg["dwi_bet"]}'

                    if cfg.get("dwi_average_b0") is True:
                        command += ' --average-b0'

                    if cfg.get("dwi_skip_smoothing") is True:
                        command += ' --skip-smoothing'

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
                currentFile = sorted(currentPath_wData.glob("*Bet.nii.gz"))
                if len(currentFile)>0:
                    command = f'python registration_DTI.py -i {_quote(currentFile[0])}'
                    result = run_subprocess(command,dataFormat,step)
                    if result != 0:
                        errorList.append(result)
                else:
                    message = f'Could not find *Bet.nii.gz in {currentPath_wData}';
                    logging.error(message)
                    errorList.append(message)
                os.chdir(cwd)
            elif step == "process":
                currentFile = sorted(currentPath_wData.glob("*dwi.nii.gz"))
                if cfg.get("dwi_denoiser") == "patch2self":
                    currentFile = sorted(currentPath_wData.glob("*Patch2SelfDenoised.nii.gz"))
                # Appends optional (fa0, nii_gz) flags to DTI main process if passed
                if len(currentFile)>0:
                    # Pull values from cfg (with defaults)
                    optional = cfg.get("dsi_optional")
                    thread_count = cfg.get("num_processes", 1)

                    cli_str = (
                        f'dsi_main.py -i {_quote(currentFile[0])} '
                        f'--thread-count {thread_count}'
                    )
                    if cfg.get("dsi_b_table") is not None:
                        cli_str += f' -b {_quote(cfg["dsi_b_table"])}'
                    if cfg.get("dsi_track_param") is not None:
                        track_param = cfg["dsi_track_param"]
                        if isinstance(track_param, (list, tuple)):
                            track_param_args = ' '.join(_quote(item) for item in track_param)
                        else:
                            track_param_args = _quote(track_param)
                        cli_str += f' -t {track_param_args}'
                    if cfg.get("dsi_recon_method") is not None:
                        cli_str += f' -r {_quote(cfg["dsi_recon_method"])}'
                    if cfg.get("dsi_vivo") is not None:
                        cli_str += f' -v {_quote(cfg["dsi_vivo"])}'
                    if cfg.get("dsi_make_isotropic") is not None:
                        cli_str += f' -m {_quote(cfg["dsi_make_isotropic"])}'
                    if cfg.get("dsi_legacy") is True:
                        cli_str += ' -l'
                    if cfg.get("dsi_skip_motion_correction") is True:
                        cli_str += ' --skip-motion-correction'
                    if optional is not None:
                        cli_str += ' -o'
                        if len(optional) > 0:
                            cli_str += ' ' + ' '.join(_quote(item) for item in optional)

                    os.chdir(cwd + '/3.2_DTIConnectivity')
                    command = f'python {cli_str}'
                    result = run_subprocess(command,dataFormat,step)
                    if result != 0:
                        errorList.append(result)
                else:
                    message = f'Could not find DWI input for DSI processing in {str(currentPath_wData)}';
                    logging.error(message)
                    errorList.append(message)
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


def create_qc_reports(project_path, steps, custom_parameters=None):
    requested_steps = set(steps)

    try:
        from helper_tools.batch_qc_reports import (
            build_bet_qc_report,
            build_registration_qc_report,
        )
    except Exception as exc:
        logging.warning("Could not import batch report tools: %s", exc)
        print(f"Report generation skipped: {exc}")
        return

    if "preprocess" in requested_steps:
        try:
            print("Creating BET report...")
            logging.info("Creating BET report...")
            html_path, count = build_bet_qc_report(
                project_path,
                n_slices=7,
                custom_parameters=custom_parameters,
            )
            if html_path:
                print(f"BET report written to {html_path} ({count} image(s))")
                logging.info("BET report written to %s (%s images)", html_path, count)
            else:
                print("BET report skipped: no BET files found.")
                logging.info("BET report skipped: no BET files found.")
        except Exception as exc:
            logging.warning("BET report generation failed: %s", exc)
            print(f"BET report generation failed: {exc}")

    if "registration" in requested_steps:
        try:
            print("Creating Registration report...")
            logging.info("Creating Registration report...")
            html_path, count = build_registration_qc_report(
                project_path,
                n_slices=7,
                custom_parameters=custom_parameters,
            )
            if html_path:
                print(f"Registration report written to {html_path} ({count} image(s))")
                logging.info("Registration report written to %s (%s images)", html_path, count)
            else:
                print("Registration report skipped: no BET/AnnoSplit_parental pairs found.")
                logging.info("Registration report skipped: no BET/AnnoSplit_parental pairs found.")
        except Exception as exc:
            logging.warning("Registration report generation failed: %s", exc)
            print(f"Registration report generation failed: {exc}")


def format_step_label(step):
    labels = {
        "preprocess": "Preprocessing",
        "registration": "Registration",
        "process": "Processing",
    }
    return labels.get(step, step.capitalize())


CLI_DEFAULT_DESCRIPTIONS = {
    "sessions": "all sessions",
    "data_types": "anat, dwi, func, t2map",
    "debug_steps": "preprocess, registration, process",
    "exemptionlist": "no exemptionlist",
    "cpu_percent": "use cpu_cores",
    "func_atlas_mask_t2": "disabled",
}

CLI_DEFAULT_SOURCES = {
    "t2_bias_method": [("2.1_T2PreProcessing/preProcessing_T2.py", "--bias-method")],
    "t2_bet": [("2.1_T2PreProcessing/preProcessing_T2.py", "--bet")],
    "t2_frac": [("2.1_T2PreProcessing/preProcessing_T2.py", "--frac")],
    "t2_radius": [("2.1_T2PreProcessing/preProcessing_T2.py", "--radius")],
    "t2_gradient": [("2.1_T2PreProcessing/preProcessing_T2.py", "--horizontal-gradient")],
    "t2_center": [("2.1_T2PreProcessing/preProcessing_T2.py", "--center")],
    "dwi_denoiser": [("2.2_DTIPreProcessing/preProcessing_DTI.py", "--denoiser")],
    "dwi_average_b0": [("2.2_DTIPreProcessing/preProcessing_DTI.py", "--average-b0")],
    "dwi_bet": [("2.2_DTIPreProcessing/preProcessing_DTI.py", "--bet")],
    "dwi_skip_smoothing": [("2.2_DTIPreProcessing/preProcessing_DTI.py", "--skip-smoothing")],
    "dwi_frac": [("2.2_DTIPreProcessing/preProcessing_DTI.py", "--frac")],
    "dwi_radius": [("2.2_DTIPreProcessing/preProcessing_DTI.py", "--radius")],
    "dwi_gradient": [("2.2_DTIPreProcessing/preProcessing_DTI.py", "--horizontal-gradient")],
    "dwi_bias_method": [("2.2_DTIPreProcessing/preProcessing_DTI.py", "--bias-method")],
    "func_bias_method": [("2.3_fMRIPreProcessing/preProcessing_fMRI.py", "--bias-method")],
    "func_skip_smoothing": [("2.3_fMRIPreProcessing/preProcessing_fMRI.py", "--skip-smoothing")],
    "func_bet": [
        ("2.3_fMRIPreProcessing/preProcessing_fMRI.py", "--bet"),
        ("3.3_fMRIActivity/process_fMRI.py", "--bet"),
    ],
    "func_frac": [
        ("2.3_fMRIPreProcessing/preProcessing_fMRI.py", "--frac"),
        ("3.3_fMRIActivity/process_fMRI.py", "--bet-frac"),
    ],
    "func_radius": [
        ("2.3_fMRIPreProcessing/preProcessing_fMRI.py", "--radius"),
        ("3.3_fMRIActivity/process_fMRI.py", "--bet-radius"),
    ],
    "func_gradient": [
        ("2.3_fMRIPreProcessing/preProcessing_fMRI.py", "--horizontal-gradient"),
        ("3.3_fMRIActivity/process_fMRI.py", "--bet-gradient"),
    ],
    "func_center": [
        ("2.3_fMRIPreProcessing/preProcessing_fMRI.py", "--center"),
        ("3.3_fMRIActivity/process_fMRI.py", "--center"),
    ],
    "func_stc": [("3.3_fMRIActivity/process_fMRI.py", "--slicetimecorrection")],
    "t2map_bet": [("4.1_T2mapPreProcessing/preProcessing_T2MAP.py", "--bet")],
    "t2map_bias_method": [("4.1_T2mapPreProcessing/preProcessing_T2MAP.py", "--bias-method")],
    "t2map_frac": [("4.1_T2mapPreProcessing/preProcessing_T2MAP.py", "--frac")],
    "t2map_radius": [("4.1_T2mapPreProcessing/preProcessing_T2MAP.py", "--radius")],
    "t2map_gradient": [("4.1_T2mapPreProcessing/preProcessing_T2MAP.py", "--horizontal-gradient")],
    "t2map_center": [("4.1_T2mapPreProcessing/preProcessing_T2MAP.py", "--center")],
    "dsi_b_table": [("3.2_DTIConnectivity/dsi_main.py", "--b-table")],
    "dsi_recon_method": [("3.2_DTIConnectivity/dsi_main.py", "--recon-method")],
    "dsi_vivo": [("3.2_DTIConnectivity/dsi_main.py", "--vivo")],
    "dsi_make_isotropic": [("3.2_DTIConnectivity/dsi_main.py", "--make-isotropic")],
    "dsi_track_param": [("3.2_DTIConnectivity/dsi_main.py", "--track-params")],
    "dsi_skip_motion_correction": [("3.2_DTIConnectivity/dsi_main.py", "--skip-motion-correction")],
    "dsi_legacy": [("3.2_DTIConnectivity/dsi_main.py", "--legacy")],
    "dsi_optional": [("3.2_DTIConnectivity/dsi_main.py", "--optional")],
}

SCRIPT_ARGUMENT_DEFAULT_CACHE = {}


def literal_ast_value(node):
    try:
        return ast.literal_eval(node)
    except (ValueError, TypeError):
        return None


def get_script_argument_defaults(script_relative_path):
    if script_relative_path in SCRIPT_ARGUMENT_DEFAULT_CACHE:
        return SCRIPT_ARGUMENT_DEFAULT_CACHE[script_relative_path]

    script_path = Path(__file__).resolve().parent / script_relative_path
    defaults = {}

    try:
        tree = ast.parse(script_path.read_text())
    except OSError:
        SCRIPT_ARGUMENT_DEFAULT_CACHE[script_relative_path] = defaults
        return defaults

    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        if not isinstance(node.func, ast.Attribute) or node.func.attr != "add_argument":
            continue

        option_strings = []
        for arg in node.args:
            value = literal_ast_value(arg)
            if isinstance(value, str):
                option_strings.append(value)
        if not option_strings:
            continue

        keyword_values = {keyword.arg: keyword.value for keyword in node.keywords}
        if "default" in keyword_values:
            default = literal_ast_value(keyword_values["default"])
        else:
            action = literal_ast_value(keyword_values["action"]) if "action" in keyword_values else None
            if action == "store_true":
                default = False
            elif action == "store_false":
                default = True
            else:
                default = None

        for option in option_strings:
            defaults[option] = default

    SCRIPT_ARGUMENT_DEFAULT_CACHE[script_relative_path] = defaults
    return defaults


def format_cli_value(value):
    if value is None:
        return "not set"
    if isinstance(value, bool):
        return "enabled" if value else "disabled"
    if isinstance(value, (list, tuple)):
        if not value:
            return "[]"
        return ", ".join(str(item) for item in value)
    return str(value)


def format_cli_default(dest):
    if dest in CLI_DEFAULT_DESCRIPTIONS:
        return f"default: {CLI_DEFAULT_DESCRIPTIONS[dest]}"

    sources = CLI_DEFAULT_SOURCES.get(dest)
    if not sources:
        return "default"

    parts = []
    for script_relative_path, option in sources:
        defaults = get_script_argument_defaults(script_relative_path)
        default = defaults.get(option)
        parts.append(f"{Path(script_relative_path).name}={format_cli_value(default)}")

    return "default: " + "; ".join(parts)


def print_cli_parameters(args):
    values = vars(args)
    parameters = []
    for label, value in values.items():
        display_value = format_cli_default(label) if value is None else format_cli_value(value)
        parameters.append(f"{label}={display_value}")

    message = "Command line (global) parameters : " + ", ".join(parameters)
    print(message)
    logging.info(message)


def get_explicit_cli_parameters(parser, args, argv):
    """Return parsed values for options explicitly supplied on the command line."""
    supplied_options = {value.split("=", 1)[0] for value in argv if value.startswith("-")}
    parameters = []

    for action in parser._actions:
        if action.dest in {"help", "input"}:
            continue
        if not supplied_options.intersection(action.option_strings):
            continue

        long_option = next(
            (option for option in action.option_strings if option.startswith("--")),
            action.option_strings[0],
        )
        parameters.append((long_option, getattr(args, action.dest)))

    return parameters


NON_SAMPLE_OVERRIDE_DESTS = {
    "help",
    "input",
    "sessions",
    "data_types",
    "debug_steps",
    "cpu_cores",
    "cpu_percent",
    "exemptionlist",
}


def get_study_id_from_path(path):
    for part in Path(path).parts:
        if part.startswith("sub-"):
            return part
    return None


def get_session_from_path(path):
    for part in Path(path).parts:
        if part.startswith("ses-"):
            return part
    return None


def normalize_exemption_column(column):
    normalized = column.strip()
    if normalized.startswith("--"):
        normalized = normalized[2:]

    if normalized == "StudyID":
        return "StudyID"
    if normalized == "session":
        return "session"
    return normalized


def get_sample_override_actions(parser):
    actions = {}
    for action in parser._actions:
        if action.dest in NON_SAMPLE_OVERRIDE_DESTS:
            continue
        long_option = next(
            (option for option in action.option_strings if option.startswith("--")),
            None,
        )
        if long_option:
            actions[normalize_exemption_column(long_option)] = action
    return actions


def split_exemption_values(value):
    if "," in value or ";" in value:
        return [item.strip() for item in value.replace(";", ",").split(",") if item.strip()]
    return value.split()


def csv_cell_has_value(value):
    if isinstance(value, list):
        return any((item or "").strip() for item in value)
    return bool((value or "").strip())


def convert_exemption_item(action, value, column):
    converter = action.type or str
    try:
        converted = converter(value)
    except Exception as exc:
        raise ValueError(
            f"Invalid value {value!r} in exemptionlist column {column!r}: {exc}"
        ) from exc

    if action.choices is not None and converted not in action.choices:
        raise ValueError(
            f"Invalid value {value!r} in exemptionlist column {column!r}. "
            f"Allowed values: {', '.join(map(str, action.choices))}"
        )
    return converted


def convert_exemption_value(action, raw_value, column):
    value = raw_value.strip()
    if value == "":
        return None

    if action.nargs == 0:
        lowered = value.lower()
        if lowered in {"1", "true", "yes", "y", "on"}:
            return True
        if lowered in {"0", "false", "no", "n", "off"}:
            return False
        raise ValueError(
            f"Invalid boolean value {value!r} in exemptionlist column {column!r}. "
            "Use true/false, yes/no or 1/0."
        )

    if action.nargs is None:
        return convert_exemption_item(action, value, column)

    values = split_exemption_values(value)
    if isinstance(action.nargs, int) and len(values) != action.nargs:
        raise ValueError(
            f"Column {column!r} expects {action.nargs} value(s), got {len(values)}: {value!r}"
        )
    if action.nargs == "+" and not values:
        raise ValueError(f"Column {column!r} expects at least one value")

    return [convert_exemption_item(action, item, column) for item in values]


def load_exemptionlist(path, parser):
    exemption_path = Path(path)
    if not exemption_path.is_file():
        raise ValueError(f"Exemptionlist does not exist: {path}")

    override_actions = get_sample_override_actions(parser)
    exemptions = {}

    with exemption_path.open(newline="") as csv_file:
        reader = csv.DictReader(csv_file)
        if not reader.fieldnames:
            raise ValueError("Exemptionlist is empty or has no header")

        normalized_columns = {
            column: normalize_exemption_column(column)
            for column in reader.fieldnames
            if column is not None
        }
        study_id_columns = [
            column for column, normalized in normalized_columns.items()
            if normalized == "StudyID"
        ]
        if len(study_id_columns) != 1:
            raise ValueError("Exemptionlist must contain exactly one StudyID column")
        study_id_column = study_id_columns[0]

        session_columns = [
            column for column, normalized in normalized_columns.items()
            if normalized == "session"
        ]
        if len(session_columns) != 1:
            raise ValueError("Exemptionlist must contain exactly one session column")
        session_column = session_columns[0]

        parameter_columns = []
        for column, normalized in normalized_columns.items():
            if column in {study_id_column, session_column}:
                continue
            if normalized not in override_actions:
                raise ValueError(
                    f"Unknown exemptionlist column {column!r}. "
                    "Use CLI parameter names without leading '--', e.g. func-frac."
                )
            parameter_columns.append((column, normalized, override_actions[normalized]))

        for line_number, row in enumerate(reader, start=2):
            if row.get(None) and csv_cell_has_value(row.get(None)):
                raise ValueError(f"Too many columns in exemptionlist line {line_number}")

            if not any(csv_cell_has_value(value) for value in row.values()):
                continue

            study_id = (row.get(study_id_column) or "").strip()
            if not study_id:
                raise ValueError(f"Missing StudyID in exemptionlist line {line_number}")
            session = (row.get(session_column) or "").strip()
            if not session:
                raise ValueError(f"Missing session in exemptionlist line {line_number}")

            exemption_key = (study_id, session)
            if exemption_key in exemptions:
                raise ValueError(
                    f"Duplicate StudyID/session pair {study_id!r}, {session!r} in exemptionlist"
                )

            overrides = {}
            for column, normalized, action in parameter_columns:
                converted = convert_exemption_value(action, row.get(column) or "", column)
                if converted is not None:
                    overrides[action.dest] = converted
            exemptions[exemption_key] = overrides

    return exemptions


def get_batch_exemption_keys(all_files):
    keys = set()
    for paths in all_files.values():
        for path in paths:
            study_id = get_study_id_from_path(path)
            session = get_session_from_path(path)
            if study_id and session:
                keys.add((study_id, session))
    return keys


def validate_exemption_samples(exemptions, all_files):
    batch_keys = get_batch_exemption_keys(all_files)
    missing = sorted(set(exemptions) - batch_keys)
    if missing:
        missing_labels = [f"{study_id}/{session}" for study_id, session in missing]
        raise ValueError(
            "Exemptionlist contains StudyID/session pair(s) that are not part of the current batch: "
            + ", ".join(missing_labels)
        )


def get_sample_cfg(base_cfg, path, exemptions):
    study_id = get_study_id_from_path(path)
    session = get_session_from_path(path)
    exemption_key = (study_id, session)
    if exemption_key not in exemptions:
        return base_cfg

    cfg = dict(base_cfg)
    cfg.update(exemptions[exemption_key])
    return cfg


TQDM_BAR_FORMAT = (
    "{desc}: {percentage:3.0f}%|{bar}| "
    "{n_fmt}/{total_fmt} done | elapsed {elapsed} | remaining {remaining}"
)

if __name__ == "__main__":
    def parse_cpu_percent(value):
        cpu_count = multiprocessing.cpu_count()
        value = str(value).strip()

        if value.endswith("%"):
            value = value[:-1]

        try:
            percent = float(value)
        except ValueError:
            raise argparse.ArgumentTypeError(
                "--cpu-percent must be a percentage from 1 to 100, e.g. 50 or 50%"
            )
        if percent <= 0 or percent > 100:
            raise argparse.ArgumentTypeError(
                "--cpu-percent must be greater than 0 and at most 100"
            )
        return max(1, int(cpu_count * percent / 100 + 0.5))

    def parse_cpu_cores(value):
        value = str(value).strip().lower()
        if value in {"min", "half", "max"}:
            return value
        try:
            cores = int(value)
        except ValueError:
            raise argparse.ArgumentTypeError(
                "--cpu-cores must be min, half, max or a positive integer"
            )
        if cores < 1:
            raise argparse.ArgumentTypeError("--cpu-cores must be at least 1")
        cpu_count = multiprocessing.cpu_count()
        if cores > cpu_count:
            raise argparse.ArgumentTypeError(
                f"--cpu-cores must not exceed the available CPU cores ({cpu_count})"
            )
        return cores

    parser = argparse.ArgumentParser(
        description=(
            "Batch processing of all data (AIDAmri). "
            "Runs preprocessing, registration and processing steps for T2, DWI, fMRI and T2map.\n\n"
            "Example:\n"
            "python batchProc.py -i /path/to/proc_data -t anat dwi "
            "--t2-frac 0.1 --t2-bias-method mico "
            "--dwi-denoiser patch2self"
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
        choices=["anat", "dwi", "func", "t2map"],
        help="Data types to process (anat, dwi, func, t2map). Default: all"
    )
    batch.add_argument(
        "-d", "--debug-steps",
        dest="debug_steps",
        nargs="+",
        choices=["preprocess", "registration", "process"],
        help="Processing steps to run (preprocess registration process). Default: all"
    )
    batch.add_argument(
        "--exemptionlist",
        help=(
            "CSV file with per-StudyID parameter overrides. "
            "Empty cells keep the normal CLI/default value. For boolean flags, use false, no, 0 or off to disable a CLI/default value. "
        )
    )
    # ============================================================
    # CPU / PARALLELIZATION
    # ============================================================
    cpu = parser.add_argument_group("cpu / parallelization")
    cpu.add_argument(
        "-c", "--cpu-cores",
        default="half",
        type=parse_cpu_cores,
        help="CPU usage preset (min, half, max) or explicit number of parallel processes"
    )
    cpu.add_argument(
        "-p", "--cpu-percent",
        dest="cpu_percent",
        type=parse_cpu_percent,
        help="CPU percentage for parallel processes, e.g. 50 or 50%%"
    )

    # ============================================================
    # T2 OPTIONS
    # ============================================================
    t2 = parser.add_argument_group("T2 options")
    t2.add_argument(
        "--t2-bias-method",
        choices=["skip", "mico", "ants"],
        type=str.lower,
        help="Bias field correction method for T2; omitted uses preProcessing_T2.py default"
    )
    t2.add_argument(
        "--t2-bet",
        choices=["skip", "bet", "bet4animal"],
        type=str.lower,
        help="Brain extraction method for T2; omitted uses preProcessing_T2.py default"
    )

    t2.add_argument(
        "--t2-frac",
        type=float,
        help="BET fractional intensity threshold"
    )
    t2.add_argument(
        "--t2-radius",
        type=int,
        help="BET head radius in mm"
    )
    t2.add_argument(
        "--t2-gradient",
        type=float,
        help="BET horizontal gradient"
    )
    t2.add_argument(
        "--t2-center",
        nargs=3,
        type=int,
        metavar=("X", "Y", "Z"),
        help="BET center in voxel coordinates"
    )
    t2.add_argument(
        "--t2-incidence-script",
        choices=["par", "detailed"],
        default="par",
        help="T2 incidence script for anat process: par (parental Atlas) runs getIncidenceSize_par.py, detailed runs getIncidenceSize.py. Default: par"
    )

    # ============================================================
    # DWI OPTIONS
    # ============================================================
    dwi = parser.add_argument_group("DWI options")
    dwi.add_argument(
        "--dwi-denoiser",
        choices=["patch2self"],
        type=str.lower,
        help="DWI denoising method"
    )
    dwi.add_argument(
        "--dwi-average-b0",
        action="store_true",
        default=None,
        help="Average b0 volumes before DWI processing"
    )
    dwi.add_argument(
        "--dwi-bet",
        choices=["skip", "bet", "bet4animal"],
        type=str.lower,
        help="Brain extraction method for DWI; omitted uses preProcessing_DTI.py default"
    )
    dwi.add_argument(
        "--dwi-skip-smoothing",
        action="store_true",
        default=None,
        help="Skip spatial median smoothing in dwi preprocessing; the 3D median reference image is still created"
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
        "--dwi-gradient",
        type=float,
        help="BET horizontal gradient for DWI"
    )
    dwi.add_argument(
        "--dwi-bias-method",
        choices=["skip", "mico", "ants"],
        type=str.lower,
        help="Bias field correction for DWI; omitted uses preProcessing_DTI.py default"
    )

    # ============================================================
    # fMRI OPTIONS
    # ============================================================
    func = parser.add_argument_group("fMRI options")
    func.add_argument(
        "--func-bias-method",
        choices=["skip", "ants"],
        type=str.lower,
        help="Bias field correction for fMRI; omitted uses preProcessing_fMRI.py default"
    )
    func.add_argument(
        "--func-skip-smoothing",
        action="store_true",
        default=None,
        help="Skip spatial median smoothing in fMRI preprocessing; the 3D median reference image is still created"
    )
    func.add_argument(
        "--func-bet",
        choices=["skip", "bet", "bet4animal"],
        type=str.lower,
        help="Brain extraction method for fMRI preprocess/process; omitted uses preProcessing_fMRI.py default"
    )
    func.add_argument(
        "--func-frac",
        type=float,
        help="BET fractional intensity threshold for fMRI"
    )
    func.add_argument(
        "--func-radius",
        type=int,
        help="BET head radius (mm) for fMRI"
    )
    func.add_argument(
        "--func-gradient",
        type=float,
        help="BET horizontal gradient for fMRI"
    )
    func.add_argument(
        "--func-center",
        nargs=3,
        type=float,
        metavar=("X", "Y", "Z"),
        help="BET center in voxel coordinates for fMRI"
    )
    func.add_argument(
        "--func-atlas-mask-t2",
        action="store_true",
        default=None,
        help="Mask the T2 BET with the T2 registered atlas annotation before fMRI registration"
    )
    func.add_argument(
        "--func-stc",
        action="store_true",
        default=None,
        help="Enable slice time correction for fMRI processing"
    )

    # ============================================================
    # T2MAP OPTIONS
    # ============================================================
    t2map = parser.add_argument_group("T2map options")
    t2map.add_argument(
        "--t2map-bet",
        choices=["skip", "bet", "bet4animal"],
        type=str.lower,
        help="Brain extraction method for T2map; omitted uses preProcessing_T2MAP.py default"
    )
    t2map.add_argument(
        "--t2map-bias-method",
        choices=["skip", "mico"],
        type=str.lower,
        help='Biasfield correction method for T2map; omitted uses preProcessing_T2MAP.py default'
    )
    t2map.add_argument(
        "--t2map-frac",
        type=float,
        help="BET fractional intensity threshold for T2map"
    )
    t2map.add_argument(
        "--t2map-radius",
        type=int,
        help="BET head radius (mm) for T2map"
    )
    t2map.add_argument(
        "--t2map-gradient",
        type=float,
        help="BET horizontal gradient for T2map"
    )
    t2map.add_argument(
        "--t2map-center",
        nargs=3,
        type=float,
        metavar=("X", "Y", "Z"),
        help="BET center in voxel coordinates for T2map"
    )

    # ============================================================
    # DSI STUDIO / TRACTOGRAPHY (dsi_main.py)
    # ============================================================
    dsi = parser.add_argument_group("DSI Studio / tractography (dsi_main.py)")
    dsi.add_argument(
        "--dsi-b-table",
        help='Diffusion gradient source; omitted uses dsi_main.py default'
    )
    dsi.add_argument(
        "--dsi-recon-method",
        type=str.lower,
        choices=["dti", "gqi"],
        help="DSI reconstruction method; omitted uses dsi_main.py default"
    )
    dsi.add_argument(
        "--dsi-vivo",
        type=str.lower,
        choices=["in_vivo", "ex_vivo"],
        help="In vivo or ex vivo data; omitted uses dsi_main.py default"
    )
    dsi.add_argument(
        "--dsi-make-isotropic",
        help="Voxel size (mm) for isotropic resampling; omitted uses dsi_main.py default"
    )
    dsi.add_argument(
        "--dsi-track-param",
        nargs="+",
        help="Tracking parameter preset or 8 custom values; omitted uses dsi_main.py default"
    )
    dsi.add_argument(
        "--dsi-skip-motion-correction",
        dest="dsi_skip_motion_correction",
        action="store_true",
        default=None,
        help="Skip slice-wise motion correction"
    )
    dsi.add_argument(
        "--dsi-legacy",
        action="store_true",
        default=None,
        help="Enable legacy .fib.gz / .src.gz support"
    )
    dsi.add_argument(
        "--dsi-optional",
        nargs="*",
        choices=["fa0", "nii_gz"],
        help="Optional dsi_main compatibility outputs (fa0, nii_gz)"
    )

    args = parser.parse_args()

    custom_parameters = get_explicit_cli_parameters(parser, args, sys.argv[1:])

    pathToData = args.input
    sessions = args.sessions
    
    #configurate the logging module
    log_file_path = os.path.join(pathToData, "batchproc_log.txt")
    configure_logging(log_file_path)
    copy_aidamri_git_information_to_proc(pathToData)

    if args.data_types is None:
        data_types = ["anat", "dwi", "func", "t2map"]
    else:
        data_types = args.data_types

    if args.debug_steps is None:
        steps = ["preprocess","registration","process"]
    else:
        steps = args.debug_steps
    
    all_files = findData(pathToData, sessions, data_types)
    exemptions = {}
    if args.exemptionlist:
        try:
            exemptions = load_exemptionlist(args.exemptionlist, parser)
            validate_exemption_samples(exemptions, all_files)
        except ValueError as exc:
            logging.error("Invalid exemptionlist: %s", exc)
            parser.error(str(exc))

        print(f"Loaded exemptionlist for {len(exemptions)} StudyID/session pair(s): {args.exemptionlist}")
        logging.info(
            "Loaded exemptionlist for %s StudyID/session pair(s): %s",
            len(exemptions),
            args.exemptionlist,
        )
        for (study_id, session), overrides in sorted(exemptions.items()):
            logging.info("Exemption parameters for %s/%s: %s", study_id, session, overrides)

    num_processes = 1

    if isinstance(args.cpu_cores, int):
        num_processes = args.cpu_cores
    elif args.cpu_cores == "min":
        num_processes = 1
    elif args.cpu_cores == "half":
        num_processes = int(multiprocessing.cpu_count() / 2)
    elif args.cpu_cores == "max":
        num_processes = multiprocessing.cpu_count()

    print_cli_parameters(args)
    
    if args.cpu_percent is not None:
        num_processes = args.cpu_percent
    
    print(f"Running with {num_processes} CPUs for the parallelization!")
    logging.info(f"Using {num_processes} CPUs for the parallelization")
    logging.info(f"Processing following datasets:\n{all_files}")
    # turns argparse.Namespace into a dict
    cfg = vars(args)
    cfg["num_processes"] = num_processes

    for key, value in all_files.items():
        if value:
            error_list_all = []
            print()
            print(f"Entered {key} data: \n{value}")
            print()
            print(f"\nStarting {key} pipeline \33[5m...\33[0m")
            print()
            for step in steps: 
                error_list_step = []
                step_label = format_step_label(step)
                step_display = f"{step_label} {key} data"
                progress_bar = tqdm(
                    total=len(value),
                    desc=step_display,
                    unit="dataset",
                    bar_format=TQDM_BAR_FORMAT,
                )
                with concurrent.futures.ProcessPoolExecutor(max_workers=num_processes) as executor:
                    futures = [
                        executor.submit(
                            executeScripts,
                            path,
                            key,
                            step,
                            get_sample_cfg(cfg, path, exemptions),
                        )
                        for path in value
                    ]

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
                        print(f"{step_display}  \033[0;30;42m Complete \33[0m")
                    else:
                        print(f"{step_display}  \033[0;30;41m Incomplete \33[0m")
                        error_list_all.extend(flat_errors_step)

                    logging.info(f"{key} {step} processing completed")

            logging.error(f"Following errors were occurring: {error_list_all}")
            logging.info(f"{key} pipeline completed")

            if not error_list_all:
                print(f"\n{key} pipeline \033[0;30;42m COMPLETED \33[0m")
            else:
                print(f"\n{key} pipeline \033[0;30;41m INCOMPLETE \33[0m")
                print()
                for err in error_list_all:
                    if isinstance(err, tuple) and len(err) == 4:
                        sub, ses, dtype, stepname = err
                        print(
                            f"Error in sub: {sub} in session: {ses} in datatype: {dtype} and step: {stepname}. Check log.")
                    else:
                        # strings or unexpected types
                        print(f"Error: {err}")

    create_qc_reports(pathToData, steps, custom_parameters=custom_parameters)

 
