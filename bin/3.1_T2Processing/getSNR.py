'''
Created on 10/08/2017

@author: Niklas Pallast
Neuroimaging & Neuroengineering
Department of Neurology
University Hospital Cologne

'''

import os
import changSNR as ch
import brummerSNR as bm
import sijbersSNR as sj
import numpy as np
import glob
import nibabel as nii


def snrCalclualtor(input_file):
    fileSNR = open(os.path.join(os.path.dirname(input_file), 'snr.txt'), 'w')

    data = nii.load(input_file)
    imgData = data.get_data()

    # nx = imgData.shape[0] # Images size in x - direction
    # ny = imgData.shape[1] # Images size in y - direction
    ns = imgData.shape[2]  # Number of slices

    noiseChSNR = np.zeros(ns)
    noiseBmSNR = np.zeros(ns)
    noiseSjSNR = np.zeros(ns)

    imgData = np.ndarray.astype(imgData, 'float64')
    for slc in range(ns):
        #   Print % of progress
        print('Slice: ' + str(slc + 1))

        # Temporal image containing all TE values for the selected slice
        slice = imgData[:, :, slc]

        curSnrCHMap, estStdChang, estStdChangNorm = ch.calcSNR(slice, 0, 1)

        curSnrBMMap, estStdBrummer, estStdBrummerNorm = bm.calcSNR(slice, 0, 1)

        curSnrSJMap, estStdSijbers, estStdSijbersNorm = sj.calcSNR(slice, 0, 1)

        noiseChSNR[slc] = estStdChang
        noiseBmSNR[slc] = estStdBrummer
        noiseSjSNR[slc] = estStdSijbers

    snrCh = 20 * np.log10(np.mean(imgData) / np.mean(noiseChSNR))
    snrBrum = 20 * np.log10(np.mean(imgData) / np.mean(noiseBmSNR))
    snrSij = 20 * np.log10(np.mean(imgData) / np.mean(noiseSjSNR))

    fileSNR.write("Mean of Chang: %0.3f dB \n" % snrCh)
    fileSNR.write("Mean of Brummer: %0.3f dB\n" % snrBrum)
    fileSNR.write("Mean of Sijbers: %0.3f dB\n" % snrSij)
    fileSNR.close()


def findRegisteredData(path, studyPrefix):
    regMR_list = []
    for filename in glob.iglob(path + '/' + studyPrefix + '*/T2w/*1.nii.gz', recursive=True):
        regMR_list.append(filename)

    return regMR_list


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Calculate SNR')
    requiredNamed = parser.add_argument_group('Required named arguments')
    requiredNamed.add_argument('-i', '--inputData', help='file name:folder of all files of one study')
    requiredNamed.add_argument('-s', '--studyPrefix', help='file name:first letter of files in input data')

    args = parser.parse_args()

    pathData = args.inputData
    studyPrefix = args.studyPrefix

    listMr = findRegisteredData(pathData, studyPrefix)

    for i in listMr:
        print(i)
        snrCalclualtor(i)
