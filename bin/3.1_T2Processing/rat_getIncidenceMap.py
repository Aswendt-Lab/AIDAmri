'''
Created on 10/08/2017
Updated on 18/12/2023

@author: Niklas Pallast, Markus Aswendt
Neuroimaging & Neuroengineering
Department of Neurology
University Hospital Cologne

'''

import os
import sys
import nibabel as nii
import glob
import numpy as np
import progressbar
import matplotlib.pyplot as plt


def heatMap(incidenceMap, araVol):
    maxV = int(np.max(incidenceMap))
    fig, axes = plt.subplots(nrows=3, ncols=4)
    t = 1
    for ax in axes.flat:
        im = ax.imshow(np.transpose(np.round(incidenceMap[:, :, t * 16])), cmap='gnuplot', vmin=0, vmax=maxV)
        ax.imshow(np.transpose(araVol[:, :, t * 16]), alpha=0.55, cmap='gray')
        ax.axis('off')
        t = t + 1

    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    bounds = np.linspace(0, maxV, maxV + 1)
    cbar = fig.colorbar(im, cax=cbar_ax, format='%1i', ticks=bounds)
    cbar.ax.tick_params(labelsize=14)
    plt.show()


def incidenceMap2(path_listInc, araTemplate, inputFile):
    araDataTemplate = nii.load(araTemplate)
    realAraImg = np.asanyarray(araDataTemplate.dataobj)
    overlazedInciedences = np.zeros([np.size(realAraImg, 0), np.size(realAraImg, 1), np.size(realAraImg, 2)])
    bar = progressbar.ProgressBar()
    for fileIndex in bar(range(len(path_listInc))):
        dataMRI = nii.load(path_listInc[fileIndex])
        volumeMRI = np.asanyarray(dataMRI.dataobj)

        bvalues = volumeMRI <= 0
        volumeMRI[bvalues] = 0

        fvalues = volumeMRI > 0
        volumeMRI[fvalues] = 1

        overlazedInciedences = overlazedInciedences + volumeMRI

    overlayNII = nii.Nifti1Image(overlazedInciedences, araDataTemplate.affine)
    output_file = os.path.join(inputFile, 'incMap.nii.gz')
    nii.save(overlayNII, output_file)
    heatMap(incidenceMap=overlazedInciedences, araVol=realAraImg)


def findIncData(path):
    regMR_list = []

    for filename in glob.iglob(path + '/anat/*IncidenceData_mask.nii.gz', recursive=False):
        regMR_list.append(filename)
    return regMR_list


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Calculate an Incidence Map')
    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument('-i', '--inputFile', help='File: Brain extracted input data')
    requiredNamed.add_argument('-s', '--studyname', help='Prefix of the study in the input folder - for example "Mouse"*')

    parser.add_argument('-a', '--sigmaBrainTemplate', help='File: Annotations of Sigma Brain', nargs='?', type=str,
                        default=os.path.abspath(os.path.join(os.getcwd(), os.pardir, os.pardir)) +
                        '/lib/SIGMA_InVivo_Brain_Template_Masked.nii.gz')

    args = parser.parse_args()
    inputFile = None
    studyname = None
    sigmaBrainTemplate = None

    if args.inputFile is not None:
        inputFile = args.inputFile
    if not os.path.exists(inputFile):
        sys.exit("Error: '%s' is not an existing directory." % (inputFile,))

    if args.sigmaBrainTemplate is not None:
        sigmaBrainTemplate = args.sigmaBrainTemplate
    if not os.path.exists(sigmaBrainTemplate):
        sys.exit("Error: '%s' is not an existing directory." % (sigmaBrainTemplate,))

    studyname = args.studyname
    path = os.path.join(inputFile, studyname)
    regInc_list = findIncData(path)

    if len(regInc_list) < 1:
        sys.exit("Error: '%s' has no masked strokes." % (studyname,))

    print("'%i' folders are part of the incidence map." % (len(regInc_list),))
    incidenceMap2(regInc_list, sigmaBrainTemplate, inputFile)
    sys.exit(0)
