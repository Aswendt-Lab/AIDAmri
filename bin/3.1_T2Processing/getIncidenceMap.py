'''
Created on 10/08/2017

@author: Niklas Pallast
Neuroimaging & Neuroengineering
Department of Neurology
University Hospital Cologne

'''

import os,sys
import nibabel as nii
import glob
import numpy as np
import scipy.io as sc
import scipy.ndimage as ndimage
import progressbar
import  getHeatmap as gHm


def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]


def thresholding(volumeMR,maskImg,thres):
    volumeMR=ndimage.gaussian_filter(volumeMR, sigma=(1.3, 1.3, 1))
    zvalues = volumeMR != 0


    if thres == 0:
        thres = np.mean(volumeMR[zvalues]) + 2*np.std(volumeMR[zvalues])

    bvalues = volumeMR < thres
    volumeMR[bvalues] = 0

    fvalues = volumeMR >= thres
    volumeMR[fvalues] = 1

    fvalues = volumeMR == 1
    return volumeMR,fvalues


def incidenceMap2(path_listInc, araTemplate):
    araDataTemplate = nii.load(araTemplate)
    realAraImg = araDataTemplate.get_data()
    overlazedInciedences = np.zeros([np.size(realAraImg, 0), np.size(realAraImg, 1), np.size(realAraImg, 2)])
    bar = progressbar.ProgressBar()
    for fileIndex in bar(range(len(path_listInc))):
        dataMRI = nii.load(path_listInc[fileIndex])
        volumeMRI = dataMRI.get_data()

        bvalues = volumeMRI <= 0
        volumeMRI[bvalues] = 0

        fvalues = volumeMRI > 0
        volumeMRI[fvalues] = 1

        overlazedInciedences = overlazedInciedences + volumeMRI

    gHm.heatMap(incidenceMap=overlazedInciedences, araVol=realAraImg)


def findIncData(path):
    regMR_list = []

    for filename in glob.iglob(path + '/T2w/*IncidenceData_mask.nii.gz', recursive=False):
        regMR_list.append(filename)

    return regMR_list



def findIncData(path):
    regMR_list = []

    for filename in glob.iglob(path + '/T2w/*IncidenceData_mask.nii.gz', recursive=False):
        regMR_list.append(filename)

    return regMR_list



if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Calculate an Incidence Map')
    parser.add_argument('-i','--inputFile', help='file name:Brain extracted input data')
    parser.add_argument('-s', '--studyname', help='prefix of the study in the input folder - for exmpale MA_30584')
    parser.add_argument('-t', '--threshold', help='threshold for stroke values ',  nargs='?', type=int,
                        default=0)
    parser.add_argument('-a', '--allenBrainTemplate', help='file name:Annotations of Allen Brain', nargs='?', type=str,
                        default=os.path.abspath(
                            os.path.join(os.getcwd(), os.pardir, os.pardir)) + '/lib/average_template_50.nii.gz')

    args = parser.parse_args()
    inputFile = None
    studyname = None
    allenBrainTemplate = None

    if args.inputFile is not None:
        inputFile = args.inputFile
    if not os.path.exists(inputFile):
        sys.exit("Error: '%s' is not an existing directory." % (inputFile,))

    if args.allenBrainTemplate is not None:
        allenBrainTemplate = args.allenBrainTemplate
    if not os.path.exists(allenBrainTemplate):
        sys.exit("Error: '%s' is not an existing directory." % (allenBrainTemplate,))

    studyname = args.studyname
    path = os.path.join(inputFile, studyname)
    regInc_list = findIncData(path)

    if len(regInc_list) < 1:
        sys.exit("Error: '%s' has no masked strokes." % (studyname,))

    print("'%i' folders are part of the incidence map." % (len(regInc_list),))
    incidenceMap2(regInc_list, allenBrainTemplate)
    sys.exit(0)
