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


def findData(projectPath, sessions, dataTypes):
    # This function screens all existing paths. Within these paths, this function collects all subject
    # folders, which are all folders that are not named 'Physio'.
    full_path_list = os.listdir(projectPath)
    all_wanted_paths, anat_files, dwi_files, func_files, t2map_files = [], [], [], [], []

    for path in full_path_list:  
        if "sub" in path and not ".DS_Store" in path:
            wanted_paths = os.listdir(os.path.join(projectPath, path))
            wanted_paths = [os.path.join(projectPath, path, wanted_path) for wanted_path in wanted_paths if "ses" in wanted_path]
            all_wanted_paths.extend(wanted_paths)
            
    if sessions:
        ses_path = []
        matching_paths = []
        for ses in sessions:
            ses_path.append("ses-" + ses)
        for path in all_wanted_paths:
            if any(ses in path for ses in ses_path):
                matching_paths.append(path)
        
        all_wanted_paths = matching_paths

    for path in all_wanted_paths:
        for sub_dir in os.listdir(path):
            if sub_dir == "anat" and "anat" in dataTypes:
                anat_files.append(os.path.join(path, sub_dir))

            elif sub_dir == "dwi" and "dwi" in dataTypes:
                dwi_files.append(os.path.join(path, sub_dir))

            elif sub_dir == "func" and "func" in dataTypes:
                func_files.append(os.path.join(path, sub_dir))

            elif sub_dir == "t2map" and "t2map" in dataTypes:
                t2map_files.append(os.path.join(path, sub_dir))

    all_files = {}
    all_files["anat"] = anat_files
    all_files["dwi"] = dwi_files
    all_files["func"] = func_files
    all_files["t2map"] = t2map_files
        

    return all_files
    
def run_subprocess(command,datatype,step,anat_process=False):
    timeout = 3600 # set maximum time in seconds after which the subprocess will be terminated
    command_args = shlex.split(command)
    file = command_args[-1]
    if datatype == "func" and step =="process":
        file = command_args[-3]
    log_file = os.path.join(os.path.dirname(file), step + ".log")
    if datatype == "anat" and step == "process":
        log_file = os.path.join(os.path.dirname(file), datatype, step + ".log")
        if anat_process == False:
            log_file = os.path.join(os.path.dirname(file), datatype, step + "_par" + ".log")

    # find current sub name
    normalized_path = os.path.normpath(file)
    directories = normalized_path.split(os.path.sep)
    sub = [directory for directory in directories if "sub-" in directory][0]
    ses = [directory for directory in directories if "ses-" in directory][0]
    
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
    

def executeScripts(currentPath_wData, dataFormat, step, stc=False, *optargs):
    # For every datatype (T2w, fMRI, DTI), go in all days/group/subjects folders
    # and execute the respective (pre-)processing/registration-scripts.
    # If a certain file does not exist, a note will be created in the errorList.
    # cwd should contain the path of the /bin folder (the user needs to navigate to the /bin folder before executing this script)
    #KEEP IN MIND DUE TO PARALLEL COMPUTING NO ERRORS IN THIS FUNCTION WILL BE PRINTED OUT => GREY ZONE
    errorList = [];
    message = '';
    cwd = str(os.getcwd())
    currentPath_wData = Path(currentPath_wData)
    # currentPath_wData = projectfolder/sub/ses/dataFormat (e.g. anat, func, dwi)
    #Find logging file
    root_path = Path(currentPath_wData).parents[2]
    log_file_path = os.path.join(root_path, "batchproc_log.txt")
    #Initialize logging only if no handler is active
    if not logging.getLogger().hasHandlers():
        logging.basicConfig(
            filename=log_file_path,
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
    if os.path.isdir(currentPath_wData):
        if dataFormat == 'anat':
            os.chdir(cwd + '/2.1_T2PreProcessing')
            if step == "preprocess":
                currentFile = list(currentPath_wData.glob("*T2w.nii.gz"))
                if len(currentFile)>0:
                    command = f'python preProcessing_T2.py -i {currentFile[0]}'
                    result = run_subprocess(command,dataFormat,step)
                    if result != 0:
                        errorList.append(result)
                else:
                    message = 'Could not find *T2w.nii.gz in {currentPath_wData}';
                    logging.error(message)
                    errorList.append(message)
                os.chdir(cwd)
            elif step == "registration":
                currentFile = list(currentPath_wData.glob("*Bet.nii.gz"))
                if len(currentFile)>0:
                    command = f'python registration_T2.py -i {currentFile[0]}'
                    result = run_subprocess(command,dataFormat,step)
                    command = f'python t2_value_extraction.py -i {currentFile[0]}'
                    result = run_subprocess(command,dataFormat,step)
                    if result != 0:
                        errorList.append(result)
                else:
                    message = 'Could not find *BiasBet.nii.gz in {currentPath_wData}';
                    logging.error(message)
                    errorList.append(message)
                os.chdir(cwd)
            elif step == "process":
                has_stroke_mask = any(currentPath_wData.glob("**/*Stroke_mask.nii.gz")) #search for stroke mask
                if not has_stroke_mask:
                    message = f"No stroke mask found for {currentPath_wData}, proceeding without mask."
                    logging.info(message)  #write in log-file
                    #print(message, flush=True)
                    return 0
                os.chdir(cwd + '/3.1_T2Processing')
                command = f'python getIncidenceSize_par.py -i {currentPath_wData}'
                result = run_subprocess(command,dataFormat,step)
                if isinstance(result, tuple) and len(result) == 4:
                    errorList.append(result)
                command = f'python getIncidenceSize.py -i {currentPath_wData}'
                result = run_subprocess(command,dataFormat,step,anat_process=True)
                if isinstance(result, tuple) and len(result) == 4:
                    errorList.append(result)
                os.chdir(cwd)
        elif dataFormat == 'func':
            os.chdir(cwd + '/2.3_fMRIPreProcessing')
            if step == "preprocess":
                currentFile = list(currentPath_wData.glob("*EPI.nii.gz"))
                if len(currentFile)>0:
                    command = f'python preProcessing_fMRI.py -i {currentFile[0]}'
                    result = run_subprocess(command,dataFormat,step)
                    if result != 0:
                        errorList.append(result)
                else:
                    message = 'Could not find *EPI.nii.gz in {currentPath_wData}';
                    logging.error(message)
                    errorList.append(message)
                os.chdir(cwd)
            elif step == "registration":
                currentFile = list(currentPath_wData.glob("*SmoothBet.nii.gz"))
                if len(currentFile)>0:
                    command = f'python registration_rsfMRI.py -i {currentFile[0]}'
                    result = run_subprocess(command,dataFormat,step)
                    if result != 0:
                        errorList.append(result)
                else:
                    message = 'Could not find *SmoothBet.nii.gz in {currentPath_wData}';
                    logging.error(message)
                    errorList.append(message)
                os.chdir(cwd)
            elif step == "process":
                currentFile = list(currentPath_wData.glob("*EPI.nii.gz"))
                if len(currentFile)>0:
                    os.chdir(cwd + '/3.3_fMRIActivity')
                    command = f'python process_fMRI.py -i {currentFile[0]} -stc {stc}'
                    result = run_subprocess(command,dataFormat,step)
                    if result != 0:
                        errorList.append(result)
                    os.chdir(cwd)
        elif dataFormat == 't2map':
            os.chdir(cwd + '/4.1_T2mapPreProcessing')
            if step == "preprocess":
                currentFile = list(currentPath_wData.glob("*MEMS.nii.gz"))
                if len(currentFile)>0:
                    command = f'python preProcessing_T2MAP.py -i {currentFile[0]}'
                    result = run_subprocess(command,dataFormat,step)
                    if result != 0:
                        errorList.append(result)
                else:
                    message = f'Could not find *MEMS.nii.gz in {currentPath_wData}';
                    logging.error(message)
                    errorList.append(message)
                os.chdir(cwd)
            elif step == "registration":
                currentFile = list(currentPath_wData.glob("*SmoothMicoBet.nii.gz"))
                if len(currentFile)>0:
                    command = f'python registration_T2MAP.py -i {currentFile[0]}'
                    result = run_subprocess(command,dataFormat,step)
                    if result != 0:
                        errorList.append(result)
                else:
                    message = f'Could not find *SmoothMicoBet.nii.gz in {currentPath_wData}';
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
                    message = f'Could not find *T2w_MAP.nii.gz in {currentPath_wData}';
                    logging.error(message)
                    errorList.append(message)
                os.chdir(cwd)
        elif dataFormat == 'dwi':
            os.chdir(cwd + '/2.2_DTIPreProcessing')
            if step == "preprocess":
                currentFile = list(currentPath_wData.glob("*dwi.nii.gz"))
                if len(currentFile)>0:
                    command = f'python preProcessing_DTI.py -i {currentFile[0]}'
                    result = run_subprocess(command,dataFormat,step)
                    if result != 0:
                        errorList.append(result)
                else:
                    message = f'Could not find *dwi.nii.gz in {currentPath_wData}';
                    logging.error(message)
                    errorList.append(message)
                os.chdir(cwd)
            elif step == "registration":
                currentFile = list(currentPath_wData.glob("*SmoothMicoBet.nii.gz"))
                if len(currentFile)>0:
                    command = f'python registration_DTI.py -i {currentFile[0]}'
                    result = run_subprocess(command,dataFormat,step)
                    if result != 0:
                        errorList.append(result)
                else:
                    message = f'Could not find *SmoothMicoBet.nii.gz in {currentPath_wData}';
                    logging.error(message)
                    errorList.append(message)
                os.chdir(cwd)
            elif step == "process":
                currentFile = list(currentPath_wData.glob("*dwi.nii.gz"))
                # Appends optional (fa0, nii_gz) flags to DTI main process if passed
                if len(currentFile)>0:
                    cli_str = f'dsi_main.py -i {currentFile[0]}'
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
        message = 'The folder '+dataFormat+' does not exist in '+currentPath
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
    parser = argparse.ArgumentParser(description='Batch processing of all data. This script runs every needed script for all registration and processing steps. The data needs to be ordered like after Bruker2NIfTI conversion: project_folder/days/groups/subjects/. For the script to work, it needs to be placed within the /bin folder of AIDAmri. Example: python batchProc.py -f /Volumes/Desktop/MRI/proc_data -g Treatment_C3a Treatment_PBS -d Baseline P7 P14 P28 P42 P56 -t T2w fMRI DTI')
    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument('-i', '--input', required=True,
                        help='Path to the parent project folder of the dataset, e.g. proc_data')


    optionalNamed = parser.add_argument_group('optional arguments')
    optionalNamed.add_argument('-s', '--sessions', required=False,
                        help='Select which sessions of your data should be processed, if no days are given all data will be used.', nargs='+')
    optionalNamed.add_argument('-stc', '--slicetimecorrection', default = "False", type=str,
                               help='Set True or False if a slice time correction should be performed. Only set true if you converted raw bruker data with conv2nifti.py from aidamri beforehand. Otherwise choose False')
    optionalNamed.add_argument('-t', '--dataTypes', required=False, nargs='+', help='Data types to be processed e.g. anat, dwi and/or func. Multiple specifications are possible.')
    optionalNamed.add_argument('-ds', '--debug_steps', required=False, nargs='+', help='Define which steps of the processing should be done. Default = [preprocess, registration, process]')
    optionalNamed.add_argument('-cpu', '--cpu_cores', required=False, default = "Half", help='Define how many parallel processes should be use to process your data. CAUTION: Too many processes will slow down your computer noticeably. Select between: ["Min", "Half", "Max"]')
    optionalNamed.add_argument('-e_cpu', '--expert_cpu', required=False, help='Define precisely how many parallel processes should be used. Enter a number.')
    

    args = parser.parse_args()
    pathToData = args.input
    sessions = args.sessions
    
    #Konfiguriere das Logging-Modul
    log_file_path = os.path.join(pathToData, "batchproc_log.txt")
    logging.basicConfig(filename=log_file_path, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    if args.slicetimecorrection is None:
        stc = False
    else:
        stc = args.slicetimecorrection
    if args.dataTypes is None:
        dataTypes = ["anat", "dwi", "func", "t2map"]
    else:
        dataTypes = args.dataTypes

    if args.debug_steps is None:
        steps = ["preprocess","registration","process"]
    else:
        steps = args.debug_steps
    
    print('Entered information:')
    print(pathToData)
    print('dataTypes %s' % dataTypes)
    print('Slice time correction [%s]' % stc)
    print('Steps %s' % steps)
    print()

    all_files = findData(pathToData, sessions, dataTypes)

    if args.cpu_cores.upper() == "MIN":
        num_processes = 1
    elif args.cpu_cores.upper() == "HALF":
        num_processes = int(multiprocessing.cpu_count() / 2)
    elif args.cpu_cores.upper() == "MAX":
        num_processes = multiprocessing.cpu_count()

    if args.expert_cpu:
        num_processes = int(args.expert_cpu)
    
    print(f"Running with {num_processes} parallel processes!")

    logging.info(f"Entered information:\n{pathToData}\n dataTypes {dataTypes}\n Slice time correction [{stc}]")
    logging.info(f"Using {num_processes} CPUs for the parallelization")
    logging.info(f"Processing following datasets:\n{all_files}")

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
                    futures = [executor.submit(executeScripts, path, key, step, stc) for path in value]

                    for future in concurrent.futures.as_completed(futures):
                        progress_bar.update(1)

                        errorList = future.result()
                        if errorList != 0:
                            if isinstance(errorList, list):
                                error_list_step.extend(errorList)
                            else:
                                error_list_step.append(errorList)
                        
                    concurrent.futures.wait(futures)
                progress_bar.close()
                error_list_all.extend(error_list_step)
                if not error_list_step:
                    print(f"{key} {step}  \033[0;30;42m COMPLETED \33[0m")
                else:
                    print(f"{key} {step}  \033[0;30;41m INCOMPLETE \33[0m")
                logging.info(f"{key} {step} processing completed")
                
                
            logging.error(f"Following errors were occuring {error_list_all}")
            logging.info(f"{key} processing completed")
            if not error_list_all:
                print(f"\n{key} processing \033[0;30;42m COMPLETED \33[0m")
            else:
                print(f"\n{key} processing \033[0;30;41m INCOMPLETE \33[0m")
            if error_list_all:
                print()
                for error in error_list_all:
                    if isinstance(error, tuple) and len(error) == 4:
                        sub, ses, datatype, step = error
                        print(
                            f"Error in sub: {sub} in session: {ses} in datatype: {datatype} and step: {step}. Check logging file for further information")
                    else:
                        print(f"Unrecognized error format: {error}")
                

 