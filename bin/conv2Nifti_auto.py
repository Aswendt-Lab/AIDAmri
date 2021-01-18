"""
Created on 27/10/2020

@author: Leon Scharwächter
AG Neuroimaging and Neuroengineering of Experimental Stroke
Department of Neurology, University Hospital Cologne

This script automates the conversion from the raw bruker data format to the NIfTI
format for the whole dataset using 1_PV2NIfTiConverter/pv_conv2Nifti.py. The raw
data needs to be in the following structure: projectfolder/days/subjects/data/.
For this script to work, the groupMapping.csv needs to be adjusted, where the group
name of every subject's folder in the raw data structure needs to be specified.
This script computes the conversion either for all data in the raw project folder
or for certain days and/or groups specified through the optional
arguments -d and -g . During the processing a new folder called proc_data is being
created in the same directory where the raw data folder is located.

Example:
python conv2Nifti_auto.py -f /Volumes/Desktop/MRI/raw_data -d Baseline P1 P7 P14
"""

import os
import csv
import shutil

def findData(projectPath, days):
# This function screens all existing paths based on the specified days
# and groups. Within these paths, this function collects all subject
# folders, which are all folders that are not named 'Physio'.
    fullPath_list = []
    path_list = [os.path.join(projectPath, days[d]) for d in range(len(days))]
    for checkPath in path_list:
        list_subfolders = [f.name for f in os.scandir(checkPath) if f.is_dir() and f.name.lower() != 'physio']
        for subject in list_subfolders:
            fullPath_list.append(os.path.join(checkPath, subject))
    return(fullPath_list)
    
def convertToNifti(subjectFolder):
    print('Run python 1_PV2NIfTiConverter/pv_conv2Nifti.py -i '+subjectFolder)
    os.system('python 1_PV2NIfTiConverter/pv_conv2Nifti.py -i '+subjectFolder)
    print('Done')
    
def moveFolder(newSubjectFolder, destination):
    shutil.move(newSubjectFolder, destination)
    
def getGroupName(subjectName):
# This function finds the index of the subject from the csv-table and
# returns the group-value at that position
    return csv_listGroups[csv_listSubjects.index(subjectName)]
    
def getSubjectAndDay(subjectPath):
# This function finds the subject name and the respective day from a
# specific subject path. Necessary structure: rawprojectfolder/day/subject
    subjectName = os.path.basename(subjectPath)
    dayPath = os.path.dirname(subjectPath)
    dayName = os.path.basename(dayPath)
    return subjectName, dayName

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='This script automates the conversion from the raw bruker data format to the NIfTI format using 1_PV2NIfTiConverter/pv_conv2Nifti.py. The raw data needs to be in the following structure: projectfolder/days/subjects/data/. For this script to work, the groupMapping.csv needs to be adjusted, where the group name of every subject''s folder in the raw data structure needs to be specified. This script computes the converison either for all data in the raw project folder or for certain days and/or groups specified through the optional arguments -d and -g. During the processing a new folder called proc_data is being created in the same directory where the raw data folder is located. Example: python conv2Nifti_auto.py -f /Volumes/Desktop/MRI/raw_data -d Baseline P1 P7 P14 P28')
    parser.add_argument('-f', '--folder', required=True,
                        help='Path to the parent project folder of the dataset, e.g. raw_data')
    parser.add_argument('-g', '--groups', nargs='+', type=str, required=False, help='Group names as in the bruker raw project folder')
    parser.add_argument('-d', '--days', nargs='+', type=str, required=False, help='Day names as in the bruker raw project folder')

    args = parser.parse_args()
    pathToRawData = args.folder
    groupNames = args.groups
    dayNames = args.days

print('Entered information:')

if dayNames == None:
    # if no days were specified, collect all subfolders (day folders) within the
    # raw project folder
    dayNames = [f.name for f in os.scandir(pathToRawData) if f.is_dir()]
print('Days to process: ', dayNames)

# Open the *.csv-file, read the subject/group pairs and save both columns in separate lists
csv_listSubjects = []
csv_listGroups = []
with open('groupMapping.csv', newline='') as csv_file:
    reader = csv.reader(csv_file, delimiter=';')
    next(reader, None) # Skip the header
    for subjects, groups in reader:
        csv_listSubjects.append(subjects)
        csv_listGroups.append(groups)
        
if groupNames == None:
    # If no groups were specified, get all unique names of the .csv group column
    # using set()
    groupNames = list(set(csv_listGroups))
print('Groups to process: ', groupNames)

# Get all subject paths (projectfolder/days/subjects/)
rawData_subfolders = findData(pathToRawData,dayNames)

# Get the number of all subject paths
numberOfSubjectFolders = len(rawData_subfolders)

#countSubjects = 0
#countSubjects = [csv_listGroups.count(group) for group in groupNames]
#print('In Total: '+str(sum(countSubjects))+' subject folders to process...')

# Create the proc_data folder in the same hierarchical level where the raw data is located
rawDataDirectory = os.path.dirname(os.path.join(pathToRawData))
procDataFolder = os.path.join(rawDataDirectory, 'proc_data')
if not os.path.isdir(procDataFolder):
    os.mkdir(procDataFolder)

# Now comes the serious part: Convert every subject folder to the NIfTI format,
# if the subject's group was specified or no group specification was made.
# Move the new generated, processed subject folder to the new proc_data folder.
# Generate the corresponding path beforehand, if necessary.
for i in range(numberOfSubjectFolders):
    current_subfolder = rawData_subfolders[i]
    subName, subDay = getSubjectAndDay(current_subfolder)
    subGroup = getGroupName(subName)
    if subGroup in groupNames:
        convertToNifti(current_subfolder)
        newDest_Day = os.path.join(procDataFolder,subDay)
        newDest_DayGroup = os.path.join(procDataFolder,subDay,subGroup)
        if not os.path.isdir(newDest_Day):
            os.mkdir(newDest_Day)
        if not os.path.isdir(newDest_DayGroup):
            os.mkdir(newDest_DayGroup)
        currentDest = os.path.join(current_subfolder,subName)
        moveFolder(currentDest,newDest_DayGroup)
