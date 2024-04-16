'''
Created on 08.04.2019
Updated: 26.09.2020
@author: Niklas Pallast and Markus Aswendt

process all DTI data
'''


import glob
import os
import numpy as np

def findData(path):


    regAtlas_list = []
    fileALL = glob.iglob(path + '/GV*/DTI/DSI_studio/*_rsfMRISplit_scaled.nii.gz', recursive=True)
    for filename in fileALL:
        regAtlas_list.append(filename)



    return regAtlas_list

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Find all related DTI data')
    parser.add_argument('-p','--pathData', help='Path to study')

    args = parser.parse_args()

    pathData = args.pathData

    listAtlas = findData(pathData)
    print(listAtlas)
    for i in range(np.size(listAtlas)):
        print(listAtlas[i])
        curPath = os.path.dirname(listAtlas[i])
        dti = glob.glob(curPath + '/*.md.nii.gz')
        if dti:
            print('python DTIdata_extract.py ' + dti[0] + ' ' + listAtlas[i] + ' -t ./acronyms_ARA.txt')
            os.system('python DTIdata_extract.py ' + dti[0] + ' ' + listAtlas[i] + ' -t ./acronyms_splitted_ARA.txt')
    for i in range(np.size(listAtlas)):
        print(listAtlas[i])
        curPath = os.path.dirname(listAtlas[i])
        dti = glob.glob(curPath + '/*.fa0.nii.gz')
        if dti:
            print('python DTIdata_extract.py ' + dti[0] + ' ' + listAtlas[i] + ' -t ./acronyms_ARA.txt')
            os.system('python DTIdata_extract.py ' + dti[0] + ' ' + listAtlas[i] + ' -t ./acronyms_splitted_ARA.txt')
    for i in range(np.size(listAtlas)):
        print(listAtlas[i])
        curPath = os.path.dirname(listAtlas[i])
        dti = glob.glob(curPath + '/*.rd.nii.gz')
        if dti:
            print('python DTIdata_extract.py ' + dti[0] + ' ' + listAtlas[i] + ' -t ./acronyms_ARA.txt')
            os.system('python DTIdata_extract.py ' + dti[0] + ' ' + listAtlas[i] + ' -t ./acronyms_splitted_ARA.txt')
    for i in range(np.size(listAtlas)):
        print(listAtlas[i])
        curPath = os.path.dirname(listAtlas[i])
        dti = glob.glob(curPath + '/*.ad.nii.gz')
        if dti:
            print('python DTIdata_extract.py ' + dti[0] + ' ' + listAtlas[i] + ' -t ./acronyms_ARA.txt')
            os.system('python DTIdata_extract.py ' + dti[0] + ' ' + listAtlas[i] + ' -t ./acronyms_splitted_ARA.txt')
