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
    if datatype == "func" and step == "process":
        file = command_args[-3]
    elif datatype == "dwi" and step == "process":
        file = command_args[3]
        # print(file) # debug
    elif datatype == "dwi" and step == "preprocess":
        file = command_args[3]
    elif datatype == "anat" and step == "preprocess":
        file = command_args[3]
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
    cwd = str(Path(__file__).resolve().parent)
    currentPath_wData = Path(currentPath_wData)
    # currentPath_wData = projectfolder/sub/ses/dataFormat (e.g. anat, func, dwi)
    if os.path.isdir(currentPath_wData):
        if dataFormat == 'anat':
            if step == "preprocess":
                os.chdir(os.path.join(cwd, '2.1_T2PreProcessing'))
                currentFile = list(currentPath_wData.glob("*T2w.nii.gz"))
                if len(currentFile) > 0:
                    command = f'python preProcessing_T2.py -i {currentFile[0]} --bias_method {args.biasfieldcorr}'
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
                    result = run_subprocess(command, dataFormat, step)
                    command = f'python t2_value_extraction.py -i {currentFile[0]}'
                    result = run_subprocess(command, dataFormat, step)
                    if result != 0:
                        errorList.append(result)
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
                    os.chdir(cwd + '/3.1_T2Processing')
                    command = f'python getIncidenceSize_par.py -i {currentPath_wData}'
                    result = run_subprocess(command,dataFormat,step)
                if result != 0:
                    errorList.append(result)
                    command = f'python getIncidenceSize.py -i {str(currentPath_wData)}'
                    result = run_subprocess(command, dataFormat, step, anat_process=True)
                if isinstance(result, tuple) and len(result) == 4:
                    command = f'python getIncidenceSize.py -i {currentPath_wData}'
                    result = run_subprocess(command,dataFormat,step,anat_process=True)
                if result != 0:
                    errorList.append(result)
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
                if len(currentFile)>0:
                    command = (f'python preProcessing_DTI.py -i {currentFile[0]} -f 0.5 -b {args.biasfieldcorr} --denoiser {args.denoiser} --use_bet4animal {args.bet4animal} --average_b0 {args.average_b0} --skip_min {args.skip_min}')
                    result = run_subprocess(command,dataFormat,step)
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
                if args.denoiser == 'patch2self':
                    currentFile = list(currentPath_wData.glob("*Patch2SelfDenoised.nii.gz"))
                # Appends optional (fa0, nii_gz) flags to DTI main process if passed
                if len(currentFile)>0:
                    cli_str = f'dsi_main.py -i {currentFile[0]} -t {track_param} -r {recon_method} -v {vivo} -m {make_isotropic} -y {flip_image_y} -template {template} -thread_count {num_processes} -l {legacy} -nomcf {no_mcf}'
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
    optionalNamed.add_argument('-denoise', '--denoiser', required=False, default=None, help='Specify the denoising method to use. Options: "patch2self" for Patch2Self denoising.')
    optionalNamed.add_argument('-bet4animal', '--bet4animal', required=False, default=False, help='Use FSL BET tuned for animal data. Default is False. Set to True to use FSL BET tuned for animal data.')
    optionalNamed.add_argument('-average_b0', '--average_b0', required=False, default=False, help='Average b0 volumes in DTI data. Default is False. Set to True to average b0 volumes.')
    optionalNamed.add_argument('-skip_min', '--skip_min', required=False, default=False, help='Skip the minimum intensity projection step in DTI preprocessing. Default is False. Set to True to skip this step.')
    optionalNamed.add_argument('-b', '--biasfieldcorr', help='Biasfield correction method - default=None, other options are "mico" or "ants"', nargs='?', type=str,default=None)
    optionalNamed.add_argument('-no_mcf', '--no_mcf', required=False, default=False, help='Skip the slice-wise MCFLIRT motion and correction step in DTI processing. Default is False. Set to True to skip this step.')
    optionalNamed.add_argument('-r', '--recon_method', required=False, default='dti', help='Specify diffusion reconstruction for DSI Studio (Default="dti", "gqi").')
    optionalNamed.add_argument('-v', '--vivo', required=False, default='in_vivo', help='Specify in vivo or ex vivo data for diffusion sampling length param0 for DSI Studio (Default="in_vivo" : param0=1.25, "ex_vivo" : param0=0.60).')
    optionalNamed.add_argument('-m', '--make_isotropic', required=False, default=0, help='Provide voxel size (mm) for isotropic resampling of diffusion data in DSI Studio (Default=0 : no resampling, "auto" uses the NIFTI header to find the voxel size for resampling).')
    optionalNamed.add_argument('-f', '--flip_image_y', required=False, default=False, help='Specify whether to flip the image in the y-direction. Default is None (no flip). Set to "true" to flip the image.')
    optionalNamed.add_argument('-template', '--template', required=False, default=1, help='Specify the template to use for the reconstruction step T2 in DSI Studio. Default is 1 (mouse). Other options are "Rat" (5) or "Mouse" (1).')
    optionalNamed.add_argument('-track_param', '--track_param', required=False, default='default', help='Provide custom tracking parameter values for DSI Studio. Options: "default", "aida_optimized", "mouse", "rat", or a list of values for: --fiber_count --interpolation --step_size --turning_angle --check_ending --fa_threshold --smoothing --min_length --max_length')
    optionalNamed.add_argument('-l', '--legacy', required=False, default=False, help='Support for legacy file types in DSI-Studio. Default is False. Set to True to use with ".fib.gz" and ".src.gz" files.')
    

    args = parser.parse_args()
    pathToData = args.input
    sessions = args.sessions
    
    #configurate the logging module
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

    num_processes = 1

    if args.cpu_cores.upper() == "MIN":
        num_processes = 1
    elif args.cpu_cores.upper() == "HALF":
        num_processes = int(multiprocessing.cpu_count() / 2)
    elif args.cpu_cores.upper() == "MAX":
        num_processes = multiprocessing.cpu_count()

    print(args)

    no_mcf = False
    if args.no_mcf is True or str(args.no_mcf).lower() == 'true':
        no_mcf = True
        print(f"Skipping slice-wise MCFLIRT motion and correction step.")
        logging.info(f"Skipping slice-wise MCFLIRT motion and correction step.")

    if args.recon_method:
        recon_method = args.recon_method
    else:
        recon_method = 'dti'
        
    logging.info(f"Using DSI Studio option for reconstruction: {recon_method}")

    if args.vivo:
        vivo = args.vivo
    else:
        vivo = 'in_vivo'
    
    logging.info(f"Using DSI Studio option param0 = {recon_method}")


    if args.make_isotropic != 0:
        make_isotropic = args.make_isotropic
        logging.info(f"Using DSI Studio option for reconstruction: isotropic voxel size resampling {make_isotropic}")
    
    if args.expert_cpu:
        num_processes = int(args.expert_cpu)

    if args.track_param:
        track_param = args.track_param
    else:
        track_param = 'default'
    
    flip_image_y = False
    if args.flip_image_y is None:
        flip_image_y = False
    elif str(args.flip_image_y).lower() == 'true':
        flip_image_y = True

    template = 6 # new default for mouse
    if args.template.lower() == 'rat':
        template = 5
    elif args.template.lower() == 'mouse':
        template = 6 # new default for mouse, pre-2024 DSI Studio used 1 for mouse template
    else:
        try:
            template = int(args.template)
        except ValueError:
            print(f"Invalid template value: {args.template}. Using default template 6 (mouse).")
            logging.info(f"Using template: {template}")
            template = 6

    legacy = False
    if args.legacy is True:
        legacy = True
        print(f"Using legacy file types .fib.gz and .src.gz for DSI Studio")
        logging.info(f"Using legacy file types .fib.gz and .src.gz for DSI Studio")
    
    print(f"Running with {num_processes} parallel processes!")

    logging.info(f"Entered information:\n{pathToData}\n dataTypes {dataTypes}\n Slice time correction [{stc}]")
    logging.info(f"Using DSI Studio options reconstruction: {recon_method} for {vivo} data")
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
                        error = error[0]
                        print(f"Error in sub: {error[0]} in session: {error[1]} in datatype: {error[2]} and step: {error[3]}. Check logging file for further information")
            
                    print()
                    for error in error_list_all:
                        error = error[0]
                        print(f"Error in sub: {error[0]} in session: {error[1]} in datatype: {error[2]} and step: {error[3]}. Check logging file for further information")
            
                

 
