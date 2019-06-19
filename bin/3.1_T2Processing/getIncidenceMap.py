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
    # uncomment the following lines to save incidence map
    # overlayNII = nii.Nifti1Image(overlazedInciedences, araDataTemplate.affine)
    # output_file = os.path.join(inputFile, 'incMap.nii.gz')
    # nii.save(overlayNII, output_file)
    heatMap(incidenceMap=overlazedInciedences, araVol=realAraImg)


def findIncData(path):
    regMR_list = []

    for filename in glob.iglob(path + '/T2w/*IncidenceData_mask.nii.gz', recursive=False):
        regMR_list.append(filename)
    return regMR_list

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Calculate an Incidence Map')
    requiredNamed = parser.add_argument_group('Required named arguments')
    requiredNamed.add_argument('-i', '--inputFile', help='file name:Brain extracted input data')
    requiredNamed.add_argument('-s', '--studyname', help='prefix of the study in the input folder - for exmpale S*')

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
    incidenceMap2(regInc_list, allenBrainTemplate, inputFile)
    sys.exit(0)
