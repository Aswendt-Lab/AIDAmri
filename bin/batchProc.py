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
python batchProc.py -f /Volumes/Desktop/MRI/proc_data -g Treatment_C3a Treatment_PBS -d Baseline P7 P14 P28 P42 P56 -t fMRI DTI
"""

import glob
import os
import fnmatch

def findData(projectPath, days, groups):
    # This function screens all existing paths based on the specified days
    # and groups. Within these paths, this function collects all subject
    # folders, which are all folders that are not named 'Physio'.
    fullPath_list = []
    path_list = [os.path.join(projectPath, days[d], groups[g]) for d in range(len(days)) for g in range(len(groups))]
    for checkPath in path_list:
        list_subfolders = [f.name for f in os.scandir(checkPath) if f.is_dir() and f.name.lower() != 'physio']
        for subject in list_subfolders:
            fullPath_list.append(os.path.join(checkPath, subject))
    return(fullPath_list)

def executeScripts(fullPath, dataTypeInput, *optargs):
    # For every datatype (T2w, fMRI, DTI), go in all days/group/subjects folders
    # and execute the respective (pre-)processing/registration-scripts.
    # If a certain file does not exist, a note will be created in the errorList.
    # cwd should contain the path of the /bin folder (the user needs to navigate to the /bin folder before executing this script)
    errorList = [];
    message = '';
    cwd = str(os.getcwd())
    for dataFormat in dataTypeInput:
        for currentPath in fullPath:
            currentPath_wData = os.path.join(currentPath,dataFormat)
            # currentPath_wData = projectfolder/day/group/subject/dataFormat (e.g. T2w, fMRI, DTI)
            if os.path.isdir(currentPath_wData):
                if dataFormat == 'T2w':
                    os.chdir(cwd + '/2.1_T2PreProcessing')
                    currentFile = find("*1.nii.gz", currentPath_wData)
                    if len(currentFile)>0:
                        print('Run python 2.1_T2PreProcessing/preProcessing_T2.py -i '+currentFile[0])
                        os.system('python preProcessing_T2.py -i '+currentFile[0])
                    else:
                        message = 'Could not find *1.nii.gz in '+currentPath_wData;
                        print(message)
                        errorList.append(message)
                    currentFile = find("*BiasBet.nii.gz", currentPath_wData)
                    if len(currentFile)>0:
                        print('Run python 2.1_T2PreProcessing/registration_T2.py -i '+currentFile[0])
                        os.system('python registration_T2.py -i '+currentFile[0])
                    else:
                        message = 'Could not find *BiasBet.nii.gz in '+currentPath_wData;
                        print(message)
                        errorList.append(message)
                    print('Run python 3.1_T2Processing/getIncidenceSize_par.py -i '+currentPath_wData)
                    os.chdir(cwd + '/3.1_T2Processing')
                    os.system('python getIncidenceSize_par.py -i '+currentPath_wData)
                    print('Run python 3.1_T2Processing/getIncidenceSize.py -i '+currentPath_wData)
                    os.system('python getIncidenceSize.py -i '+currentPath_wData)
                    os.chdir(cwd)
                elif dataFormat == 'fMRI':
                    os.chdir(cwd + '/2.3_fMRIPreProcessing')
                    currentFile = find("*1.nii.gz", currentPath_wData)
                    if len(currentFile)>0:
                        print('Run python 2.3_fMRIPreProcessing/preProcessing_fMRI.py -i '+currentFile[0])
                        os.system('python preProcessing_fMRI.py -i '+currentFile[0])
                    else:
                        message = 'Could not find *1.nii.gz in '+currentPath_wData;
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
                    currentFile = find("*1.nii.gz", currentPath_wData)
                    if len(currentFile)>0:
                        print('Run python 3.3_fMRIActivity/process_fMRI.py -i '+ currentFile[0])
                        os.chdir(cwd + '/3.3_fMRIActivity')
                        os.system('python process_fMRI.py -i '+ currentFile[0])
                    os.chdir(cwd)
                elif dataFormat == 'DTI':
                    os.chdir(cwd + '/2.2_DTIPreProcessing')
                    currentFile = find("*1.nii.gz", currentPath_wData)
                    if len(currentFile)>0:
                        print('Run python 2.2_DTIPreProcessing/preProcessing_DTI.py -i '+currentFile[0])
                        os.system('python preProcessing_DTI.py -i '+currentFile[0])
                    else:
                        message = 'Could not find *1.nii.gz in '+currentPath_wData;
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
                    currentFile = find("*1.nii.gz", currentPath_wData)

                    # Appends optional (fa0, nii_gz) flags to DTI main process if passed
                    if len(currentFile)>0:
                        cli_str = r'dsi_main.py -i %s' % currentFile[0]
                        if len(optargs) > 0:
                            cli_str += ' -o'
                            for arg in optargs:
                                cli_str += ' %s' % arg
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
    requiredNamed.add_argument('-f', '--folder', required=True,
                        help='Path to the parent project folder of the dataset, e.g. proc_data')
    requiredNamed.add_argument('-g', '--groups', required=True, nargs='+', help='Group names as in the Bruker2NIfTI processed project folder')
    requiredNamed.add_argument('-d', '--days', required=True, nargs='+', help='Day names as in the Bruker2NIfTI processed project folder')
    requiredNamed.add_argument('-t', '--dataTypes', required=True, nargs='+', help='Data types to be processed e.g. T2w, DTI and/or fMRI. Multiple specifications are possible.')

    optionalNamed = parser.add_argument_group('optional arguments')
    optionalNamed.add_argument('-o', '--optional', nargs = '*', help = 'Optional arguments (e.g. deprecated terminology)')

    args = parser.parse_args()
    pathToData = args.folder
    groupNames = args.groups
    dayNames = args.days
    dataTypes = args.dataTypes
    if args.optional is not None:
        optionals = ' '.join(args.optional)
    else:
        optionals = args.optional
    
    print('Entered information:')
    print(groupNames)
    print(dayNames)
    print(dataTypes)
    print('Optional arguments: [%s]' % optionals)

    listMr = findData(pathToData, dayNames, groupNames)
    if optionals is not None:
        executeScripts(listMr, dataTypes, optionals)
    else:
        executeScripts(listMr, dataTypes)
