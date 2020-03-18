"""
Created on 10/08/2017

@author: Niklas Pallast
Neuroimaging & Neuoengineering, Cologne

"""


import glob
import os


def findData(path, data, post):
    regMR_list = []
    if post is not None and data is not None:
        fullString = os.path.join(path, data, post)
        fileALL = glob.iglob(fullString, recursive=True)
    elif data is not None:
        fullString = os.path.join(path, data)
        fileALL = glob.iglob(fullString, recursive=True)
    else:
        fullString = os.path.join(path)
        fileALL = glob.iglob(fullString, recursive=False)

    for filename in fileALL:
        regMR_list.append(filename)

    return regMR_list

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Batch proscessing of all data')
    parser.add_argument('-c', '--command',
                        help='full path to the python command e.g. /Users/aswendtm/AIDA/bin/2.1_T2PreProcessing/registration_T2.py')
    parser.add_argument('-f', '--folder',
                        help='parent folder of dataset e.g. /Volumes/AG_Aswendt_Projects/TVA_GFAP_Vimentin_Goeteborg/MRI/proc_data/P1/GFAP_Vim/GV*')
    parser.add_argument('-d', '--dataType', help='folder of data e.g. DTI, T2w, or fMRI')
    parser.add_argument('-p', '--postfix', help='postfix of filename e.g. *Bet.nii.gz or *1.nii.gz', default=None)
    args = parser.parse_args()

    pathData = args.folder
    dataType = args.dataType
    postfix = args.postfix
    command = os.path.basename(args.command)
    cmd_path = os.path.dirname(args.command)

    listMr = findData(pathData, dataType, postfix)
    for x in range(len(listMr)):
        print(listMr[x])
    os.chdir(cmd_path)
    for i in listMr:
        print('python ' + command + ' -i ' + i)  # + ' -d')

        os.system('python ' + command + ' -i ' + i)  # + ' -d')

    os.chdir(os.path.dirname(cmd_path))
