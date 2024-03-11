"""
Created on 10/08/2017

@author: Niklas Pallast
Neuroimaging & Neuroengineering
Department of Neurology
University Hospital Cologne



"""

import os, sys
import changSNR as ch
import brummerSNR as bm
import sijbersSNR as sj
import numpy as np
import glob
import nibabel as nii


def snrCalclualtor(input_file, method):
    fileSNR = open(os.path.join(os.path.dirname(input_file),'snr.txt'), 'w')

    data = nii.load(input_file)
    imgData = data.get_data()

    #nx = imgData.shape[0] # Images size in x - direction
    #ny = imgData.shape[1] # Images size in y - direction
    ns = imgData.shape[2] # Number of slices

    noiseSNR = np.zeros(ns)

    imgData = np.ndarray.astype(imgData,'float64')
    for slc in range(ns):
        #   Print % of progress
        print('Slice: ' + str(slc + 1))


        # Temporal image containing all TE values for the selected slice
        slice = imgData[:, :, slc]

        if method == 1:
            curSnr, estStd, estStdNorm = ch.calcSNR(slice, 0, 1)

        elif method == 2:
            curSnr, estStd, estStdNorm = bm.calcSNR(slice, 0, 1)

        elif method == 3:
            curSnr, estStd, estStdNorm = sj.calcSNR(slice, 0, 1)
        else:
            sys.exit(
                "Error: '%i' is not an existing choice for a SNR method!" % (method,))

        noiseSNR[slc] = estStd

    snr = 20 * np.log10(np.mean(imgData) / np.mean(noiseSNR))

    fileSNR.write("SNR [dB]: %0.3f \n" % snr)
    fileSNR.close()


def findRegisteredData(path):
    regMR_list = []
    for filename in glob.iglob(path + '*/T2w/*1.nii.gz', recursive=True):
        regMR_list.append(filename)

    return regMR_list

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Calculates SNR and generates snr.txt in T2w files')
    parser.add_argument('-p', '--pathData', help='src path to all processed files')
    parser.add_argument('-f', '--filePrefix', help='file prefix in src path')
    parser.add_argument('-m', '--SNRmethod', help='1: Brummer(default) 2:Chang 3:Sijbers', nargs='?', type=int,
                        default=1)

    args = parser.parse_args()

    pathData = args.pathData + '/' + args.filePrefix

    listMr = findRegisteredData(pathData)
    method = args.SNRmethod

    for i in listMr:
        print(i)
        snrCalclualtor(i, method)
