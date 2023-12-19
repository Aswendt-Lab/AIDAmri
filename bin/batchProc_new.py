"""
Created on 11/18/2020
Updated on 12/19/2023

@author: Leon Scharwächter, Marc Schneider
AG Neuroimaging and Neuroengineering of Experimental Stroke
Department of Neurology, University Hospital Cologne

This script runs every needed script for all (pre-)processing and registration
steps. The data needs to be ordered like after Bruker2NIfTI conversion:
project_folder/days/groups/subjects/.
For the script to work, it needs to be placed within the /bin folder of AIDAmri.

Example:
python batchProc.py -i /Volumes/Desktop/MRI/proc_data -t anat dwi func t2map
"""

import os
import glob
import fnmatch
import shutil
import subprocess
from pathlib import Path
import concurrent.futures
from tqdm import tqdm
import multiprocessing

def findData(projectPath, sessions, dataTypes):
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
    errorList = []
    message = ''
    cwd = str(Path.cwd())
    currentPath_wData = Path(currentPath_wData)
    
    if os.path.isdir(currentPath_wData):
        if dataFormat == 'anat':
            os.chdir(cwd + '/2.1_T2PreProcessing')
            currentFile = list(currentPath_wData.glob("*T2w.nii.gz"))
            print(f"currentPath_wData: {currentPath_wData}")
            print(f"currentFile: {currentFile}")

            if len(currentFile) > 0:
                print(f"Found *T2w.nii.gz in {currentPath_wData}")
                os.system(f'python preProcessing_T2.py -i {currentFile[0]}')
            else:
                message = f'Could not find *T2w.nii.gz in {currentPath_wData}'
                print(message)
                errorList.append(message)

            os.chdir(cwd + '/3.1_T2Processing')
            os.system('python getIncidenceSize_par.py -i '+currentPath_wData)
            os.system('python getIncidenceSize.py -i '+currentPath_wData)
            os.chdir(cwd)
        elif dataFormat == 'func':
            os.chdir(cwd + '/2.3_fMRIPreProcessing')
            currentFile = glob.glob(os.path.join(currentPath_wData,"*EPI.nii.gz"))
            if len(currentFile) > 0:
                os.system('python preProcessing_fMRI.py -i '+currentFile[0])
            else:
                message = 'Could not find *EPI.nii.gz in '+currentPath_wData;
                print(message)
                errorList.append(message)

            currentFile = glob.glob(os.path.join(currentPath_wData,"*SmoothBet.nii.gz"))
            if len(currentFile) > 0:
                os.system('python registration_rsfMRI.py -i '+currentFile[0])
            else:
                message = 'Could not find *SmoothBet.nii.gz in '+currentPath_wData;
                print(message)
                errorList.append(message)

            currentFile = glob.glob(os.path.join(currentPath_wData,"*EPI.nii.gz"))
            if len(currentFile) > 0:
                os.chdir(cwd + '/3.3_fMRIActivity')
                os.system('python process_fMRI.py -i '+ currentFile[0] + ' -stc ' + str(stc))
            os.chdir(cwd)
        elif dataFormat == 't2map':
            os.chdir(cwd + '/4.1_T2mapPreProcessing')
            currentFile = glob.glob(os.path.join(currentPath_wData,"*MEMS.nii.gz"))
            if len(currentFile) > 0:
                os.system('python preProcessing_T2MAP.py -i '+currentFile[0])
            else:
                message = 'Could not find *MEMS.nii.gz in '+currentPath_wData;
                print(message)
                errorList.append(message)

            currentFile = glob.glob(os.path.join(currentPath_wData,"*SmoothMicoBet.nii.gz"))
            if len(currentFile) > 0:
                os.system('python registration_T2MAP.py -i '+currentFile[0])
            else:
                message = 'Could not find *SmoothMicoBet.nii.gz in '+currentPath_wData;
                print(message)
                errorList.append(message)

            currentFile = glob.glob(os.path.join(currentPath_wData,"*T2w_MAP.nii.gz"))
            rois_file = find("*AnnoSplit_t2map.nii.gz", currentPath_wData)
            if len(currentFile) > 0 and len(rois_file) > 0:
                print('Run python 4.1_T2mapPreProcessing/t2map_data_extract.py -i '+currentFile[0] + ' -r ' + rois_file[0])
                os.system('python t2map_data_extract.py -i '+currentFile[0] + ' -r ' + rois_file[0])
            else:
                message = 'Could not find *T2w_MAP.nii.gz in '+currentPath_wData;
                print(message)
                errorList.append(message)
        elif dataFormat == 'dwi':
            os.chdir(cwd + '/2.2_DTIPreProcessing')
            currentFile = glob.glob(os.path.join(currentPath_wData,"*dwi.nii.gz"))
            if len(currentFile) > 0:
                os.system('python preProcessing_DTI.py -i '+currentFile[0])
            else:
                message = 'Could not find *dwi.nii.gz in '+currentPath_wData;
                print(message)
                errorList.append(message)

            currentFile = glob.glob(os.path.join(currentPath_wData, "*SmoothMicoBet.nii.gz"))
            if len(currentFile) > 0:
                os.system('python registration_DTI.py -i '+currentFile[0])
            else:
                message = 'Could not find *SmoothMicoBet.nii.gz in '+currentPath_wData;
                print(message)
                errorList.append(message)

            currentFile = glob.glob(os.path.join(currentPath_wData,"*dwi.nii.gz"))

            if len(currentFile) > 0:
                cli_str = r'dsi_main.py -i %s' % currentFile[0]
                print('Run python 3.2_DTIConnectivity/%s' % cli_str)
                os.chdir(cwd + '/3.2_DTIConnectivity')
                os.system('python %s' % cli_str)
            os.chdir(cwd)
        else:
            message = 'The data folders'' names do not match anat, dwi, func or t2map';
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
    optionalNamed.add_argument('-stc', '--slicetimecorrection', default="False", type=str,
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
