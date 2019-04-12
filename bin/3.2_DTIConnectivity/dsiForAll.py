"""
Created on 10/08/2017

@author: Niklas Pallast
Neuroimaging & Neuoengineering, Cologne

"""


import glob
import os


def findRegisteredData(path):
    regMR_list = []
    fileALL = glob.iglob(path + '/P27/SP_V1_2_2*/DTI/*.1.nii.gz', recursive=True)
    for filename in fileALL:
        regMR_list.append(filename)

    return regMR_list

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Exponential curve fitting function  ')
    parser.add_argument('-p','--pathData', help='file name:Brain extracted input data')

    args = parser.parse_args()

    pathData = args.pathData

    listMr = findRegisteredData(pathData)
    print(listMr)
    for i in listMr:
        print('python dsi_main.py -i '+i)
        os.system('python dsi_main.py -i '+i)