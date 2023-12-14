"""
Created on 18/11/2020

@author: Leon Scharwächter
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
    

def executeScripts(currentPath_wData, dataFormat, stc=False, *optargs):
    # For every datatype (T2w, fMRI, DTI), go in all days/group/subjects folders
    # and execute the respective (pre-)processing/registration-scripts.
    # If a certain file does not exist, a note will be created in the errorList.
    # cwd should contain the path of the /bin folder (the user needs to navigate to the /bin folder before executing this script)
    errorList = [];
    message = '';
    cwd = str(os.getcwd())
    # currentPath_wData = projectfolder/sub/ses/dataFormat (e.g. anat, func, dwi)
    if os.path.isdir(currentPath_wData):
        if dataFormat == 'anat':
            os.chdir(cwd + '/2.1_T2PreProcessing')
            currentFile = find("*T2w.nii.gz", currentPath_wData)
            if len(currentFile)>0:
                #print('Run python 2.1_T2PreProcessing/preProcessing_T2.py -i '+currentFile[0])
                os.system('python preProcessing_T2.py -i '+currentFile[0])
            else:
                message = 'Could not find *T2w.nii.gz in '+currentPath_wData;
                print(message)
                errorList.append(message)
            currentFile = find("*Bet.nii.gz", currentPath_wData)
            if len(currentFile)>0:
                #print('Run python 2.1_T2PreProcessing/registration_T2.py -i '+currentFile[0])
                os.system('python registration_T2.py -i '+currentFile[0])
            else:
                message = 'Could not find *BiasBet.nii.gz in '+currentPath_wData;
                print(message)
                errorList.append(message)
            #print('Run python 3.1_T2Processing/getIncidenceSize_par.py -i '+currentPath_wData)
            os.chdir(cwd + '/3.1_T2Processing')
            os.system('python getIncidenceSize_par.py -i '+currentPath_wData)
            #print('Run python 3.1_T2Processing/getIncidenceSize.py -i '+currentPath_wData)
            os.system('python getIncidenceSize.py -i '+currentPath_wData)
            os.chdir(cwd)
        elif dataFormat == 'func':
            os.chdir(cwd + '/2.3_fMRIPreProcessing')
            currentFile = find("*EPI.nii.gz", currentPath_wData)
            if len(currentFile)>0:
                print('Run python 2.3_fMRIPreProcessing/preProcessing_fMRI.py -i '+currentFile[0])
                os.system('python preProcessing_fMRI.py -i '+currentFile[0])
            else:
                message = 'Could not find *EPI.nii.gz in '+currentPath_wData;
                print(message)
                errorList.append(message)
            currentFile = find("*SmoothBet.nii.gz", currentPath_wData)
            if len(currentFile)>0:
                print('Run python 2.3_fMRIPreProcessing/registration_rsfMRI.py -i '+currentFile[0])
                os.system('python registration_rsfMRI.py -i '+currentFile[0])
            else:
                message = 'Could not find *SmoothBet.nii.gz in '+currentPath_wData;
                print(message)
                errorList.append(message)
            currentFile = find("*EPI.nii.gz", currentPath_wData)
            if len(currentFile)>0:
                print('Run python 3.3_fMRIActivity/process_fMRI.py -i '+ currentFile[0] + ' -stc ' + str(stc))
                os.chdir(cwd + '/3.3_fMRIActivity')
                os.system('python process_fMRI.py -i '+ currentFile[0] + ' -stc ' + str(stc))
            os.chdir(cwd)
        elif dataFormat == 't2map':
            os.chdir(cwd + '/4.1_T2mapPreProcessing')
            currentFile = find("*MEMS.nii.gz", currentPath_wData)
            if len(currentFile)>0:
                print('Run python 4.1_T2mapPreProcessing/preProcessing_T2MAP.py -i '+currentFile[0])
                os.system('python preProcessing_T2MAP.py -i '+currentFile[0])
            else:
                message = 'Could not find *MEMS.nii.gz in '+currentPath_wData;
                print(message)
                errorList.append(message)
            currentFile = find("*SmoothMicoBet.nii.gz", currentPath_wData)
            if len(currentFile)>0:
                print('Run python 4.1_T2mapPreProcessing/registration_T2MAP.py -i '+currentFile[0])
                os.system('python registration_T2MAP.py -i '+currentFile[0])
            else:
                message = 'Could not find *SmoothMicoBet.nii.gz in '+currentPath_wData;
                print(message)
                errorList.append(message)
            currentFile = find("*T2w_MAP.nii.gz", currentPath_wData)
            rois_file = find("*AnnoSplit_t2map.nii.gz", currentPath_wData)
            if len(currentFile)>0 and len(rois_file)>0:
                print('Run python 4.1_T2mapPreProcessing/t2map_data_extract.py -i '+currentFile[0] + ' -r ' + rois_file[0])
                os.system('python t2map_data_extract.py -i '+currentFile[0] + ' -r ' + rois_file[0])
            else:
                message = 'Could not find *T2w_MAP.nii.gz in '+currentPath_wData;
                print(message)
                errorList.append(message)
        elif dataFormat == 'dwi':
            os.chdir(cwd + '/2.2_DTIPreProcessing')
            currentFile = find("*dwi.nii.gz", currentPath_wData)
            if len(currentFile)>0:
                print('Run python 2.2_DTIPreProcessing/preProcessing_DTI.py -i '+currentFile[0])
                os.system('python preProcessing_DTI.py -i '+currentFile[0])
            else:
                message = 'Could not find *dwi.nii.gz in '+currentPath_wData;
                print(message)
                errorList.append(message)
            currentFile = find("*SmoothMicoBet.nii.gz", currentPath_wData)
            if len(currentFile)>0:
                print('Run python 2.2_DTIPreProcessing/registration_DTI.py -i '+currentFile[0])
                os.system('python registration_DTI.py -i '+currentFile[0])
            else:
                message = 'Could not find *SmoothMicoBet.nii.gz in '+currentPath_wData;
                print(message)
                errorList.append(message)
            currentFile = find("*dwi.nii.gz", currentPath_wData)

            # Appends optional (fa0, nii_gz) flags to DTI main process if passed
            if len(currentFile)>0:
                cli_str = r'dsi_main.py -i %s' % currentFile[0] 
                print('Run python 3.2_DTIConnectivity/%s' % cli_str)
                os.chdir(cwd + '/3.2_DTIConnectivity')
                os.system('python %s' % cli_str)
            os.chdir(cwd)
        else:
            message = 'The data folders'' names do not match T2w, fMRI or DTI';
            print(message);
            errorList.append(message)
    else:
        message = 'The folder '+dataFormat+' does not exist in '+currentPath
        print(message)
        errorList.append(message)
    print('')
    print('Errors:')
    print(errorList)
 
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

    args = parser.parse_args()
    pathToData = args.input
    sessions = args.sessions

    if args.slicetimecorrection is None:
        stc = False
    else:
        stc = args.slicetimecorrection
    if args.dataTypes is None:
        dataTypes = ["anat", "dwi", "func", "t2map"]
    else:
        dataTypes = args.dataTypes
    
    print('Entered information:')
    print(pathToData)
    print('dataTypes %s' % dataTypes)
    print('Slice time correction [%s]' % stc)
    print()

    all_files = findData(pathToData, sessions, dataTypes)

    num_processes = multiprocessing.cpu_count()

    for key, value in all_files.items():
        if value:
            print(f"Entered {key} data: \n{value}")
            print()
            with concurrent.futures.ProcessPoolExecutor(max_workers=num_processes) as executor:

                progress_bar = tqdm(total=len(value), desc=f"Processing {key} data")

                print(f"\n{key} processing \33[5m...\33[0m (wait!)", end="\r")
                
                futures = [executor.submit(executeScripts, path, key, stc) for path in value]

                for future in futures:
                    future.result()
                    progress_bar.update(1)

                concurrent.futures.wait(futures)

                progress_bar.close()

                print(f"{key} processing  \033[0;30;42m COMPLETED \33[0m")

 