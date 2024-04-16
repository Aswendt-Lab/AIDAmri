"""
Created on 10/08/2017

@author: Niklas Pallast
Neuroimaging & Neuroengineering
Department of Neurology
University Hospital Cologne

"""

import sys,os
import numpy as np
import glob
import parReader
import i32Reader

def findData(path,addon):
    reg_list = []
    fileALL = glob.iglob(path+'/'+addon, recursive=True)
    for filename in fileALL:
        reg_list.append(filename)
    return reg_list

def getRegrTable(file_name):
    # proof par Folder
    par_Path = os.path.join(file_name,'rs-fMRI_mcf')
    if not os.path.exists(par_Path):
        sys.exit("Error: %s is not an existing directory or file." % (par_Path,))

    # proof i32 Folder
    i32_Path =  os.path.join(file_name,'rawMonData')
    if not os.path.exists(i32_Path):
        sys.exit("Error: %s is not an existing directory or file." % (i32_Path,))

    # generate target Folder
    target_folder = os.path.join(file_name, 'txtRegrPython')
    if not os.path.exists(target_folder):
        os.mkdir(target_folder)

    # get all file entries and compare the length
    listofPar_names = findData(par_Path,'*mcf.mat')
    listofI32_names = findData(i32_Path,'*I32')
    #if not  len(listofPar_names) == len(listofI32_names):
    #   print('\x1b[00;37;43m' + 'Some Data of I32 have no corresponding par data!' + '\x1b[0m')

    headlineStr = ['#Resp. BLC(1)','Resp. Deriv.(2)','Card. BLC(3)',
                  'Card. Deriv.(4)','RotX(5)','RotY(6)','RotZ(7)',
                  'dX(8)','dY(9)','dZ(10)','1st Order Drift(11)',
                  '2nd Order Drift(12)','3rd Order Drift(13)']

    for i in range(len(listofPar_names)):
        #  get par folder info
        par_folder_path = os.path.join(listofPar_names[i])
        fullpar_folder_name = os.path.basename(par_folder_path).split('.')
        par_folder_name2comp = fullpar_folder_name[1]+'.PHYSWAV.'+fullpar_folder_name[2]


        # get par-Folder content
        cur_contentOfPar = findData(par_folder_path, '*.par')
        numberOfSlices = len(cur_contentOfPar)

        # find corresponding I32 dataset
        str_indexI32 = [s for s in listofI32_names if par_folder_name2comp in s]
        if len(str_indexI32) < 1:
            print('\x1b[00;37;43m' + 'No corresponding .par - .i32 file for %s'  % (par_folder_name2comp,) + '\x1b[0m')
            continue
        else:
            i32_folder_path = str_indexI32[0]

        # read the first par Table to get the real number of Repition
        cur_par_file_path = cur_contentOfPar[0]
        parTestTable = parReader.getPar(cur_par_file_path)
        # delete the first five measurements
        numberOfAllRepitionsParTable = len(parTestTable) - 5


        # get i32 - Data
        trigger,i32Table = i32Reader.getI32(i32_folder_path,numberOfSlices,numberOfAllRepitionsParTable)
        numberOfAllRepitionsI32 = len(trigger)/numberOfSlices

        # generate drifts
        driftTable = np.zeros([numberOfAllRepitionsParTable,3])
        x = np.linspace(-1,1,numberOfAllRepitionsParTable)
        driftTable[:,0] = x
        driftTable[:,1] = x**2
        driftTable[:,2] = x**3

        for j in range(len(cur_contentOfPar)):
            tempRgrName = os.path.basename(cur_contentOfPar[j]).split('_slc')
            cur_slc = int(tempRgrName[1][0:4])
            rgr_folder_name = tempRgrName[0] + '_mcf_slice_' + tempRgrName[1][0:4] + '.txt'
            rgr_folder_path = os.path.join(target_folder,rgr_folder_name)

            # get par - Data
            cur_par_file_path = cur_contentOfPar[j]
            parTable = parReader.getPar(cur_par_file_path)

            # get I32 entries for cur_slc
            cur_slc_i32entries = i32Table[trigger[cur_slc::numberOfSlices]]

            # merge i32Table, parTable and driftTable
            fid = open(rgr_folder_path,'w')
            fid.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'
                      % (headlineStr[0],headlineStr[1],headlineStr[2],headlineStr[3],headlineStr[4],
                         headlineStr[5],headlineStr[6],headlineStr[7],headlineStr[8],headlineStr[9]
                         ,headlineStr[10],headlineStr[11],headlineStr[12]))

            for repIndex in range(numberOfAllRepitionsParTable):
                fid.write('%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n'
                          % (cur_slc_i32entries[repIndex,0],cur_slc_i32entries[repIndex,1],
                             cur_slc_i32entries[repIndex,2],cur_slc_i32entries[repIndex,3],
                             parTable[repIndex,0],parTable[repIndex,1],parTable[repIndex,2],
                             parTable[repIndex,3],parTable[repIndex,4],parTable[repIndex,5],
                             driftTable[repIndex,0],driftTable[repIndex,1],driftTable[repIndex,2]))
            fid.close()
    return 0

if __name__ == "__main__":


    import argparse
    parser = argparse.ArgumentParser(description='Generate Regression Table out of par and I32')

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-i','--input', help='Path to input data',required=True)
    args = parser.parse_args()


    if args.input is not None and args.input is not None:
        input = args.input
    if not os.path.exists(input):
        sys.exit("Error: '%s' is not an existing directory or file %s is not in directory." % (input, args.file,))

    result = getRegrTable(input)
