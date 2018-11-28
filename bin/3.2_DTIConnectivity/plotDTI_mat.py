"""
Created on 10/08/2017

@author: Niklas Pallast
Neuroimaging & Neuroengineering
Department of Neurology
University Hospital Cologne



"""

import matplotlib.pyplot as plt
import os, sys
import numpy as np
import scipy.io as sio

np.seterr(divide='ignore', invalid='ignore')
import seaborn as sns


def intersect_mtlb(a, b):
    a1, ia = np.unique(a, return_index=True)
    b1, ib = np.unique(b, return_index=True)
    aux = np.concatenate((a1, b1))
    aux.sort()
    c = aux[:-1][aux[1:] == aux[:-1]]
    return ia[np.isin(a1, c)]


def getRefLabels(prefix):
    if "rsfMRISplit" in prefix:
        dataTemplate = np.loadtxt(
            os.path.abspath(os.path.join(os.getcwd(), os.pardir, os.pardir)) + '/lib/annoVolume+2000_rsfMRI.nii.txt',
            dtype=str)
        refLabels = dataTemplate[:, 1]

    elif "rsfMRI" in prefix:
        dataTemplate = np.loadtxt(
            os.path.abspath(os.path.join(os.getcwd(), os.pardir, os.pardir)) + '/lib/annoVolume.nii.txt', dtype=str)
        refLabels = dataTemplate[:, 1]

    else:
        dataTemplate = np.loadtxt(
            os.path.abspath(os.path.join(os.getcwd(), os.pardir, os.pardir)) + '/lib/ARA_changedAnnotatiosn2DTI.txt',
            dtype=str)
        refLabels = dataTemplate[:, 1]

    return refLabels


def matrixMaker(inputPath):
    # Read pass and end
    if "pass" in inputPath:
        matData = sio.loadmat(inputPath)
        connectivityPass = matData['connectivity']
        matData = sio.loadmat(inputPath.replace('.pass.', '.end.'))
        connectivityEnd = matData['connectivity']
        connectivity = connectivityEnd + connectivityPass

    elif "end" in inputPath:
        matData = sio.loadmat(inputPath)
        connectivityEnd = matData['connectivity']
        matData = sio.loadmat(inputPath.replace('.end.', '.pass.'))
        connectivityPass = matData['connectivity']
        connectivity = connectivityEnd + connectivityPass

    else:
        sys.exit("Error: %s path do not conatain path or end data." % (inputPath,))

    labels = matData['name']
    tempLabels = ""
    labels = tempLabels.join([chr(a) for a in labels[0]]).split('\n')

    # Get reference Labels
    refLabels = getRefLabels(os.path.basename(inputPath))

    # Intersection between ref and cur labels
    ia = intersect_mtlb(refLabels, labels)
    missingLabels = np.setdiff1d(np.arange(1, len(refLabels)), ia)

    # Adapt labels to pyplot
    labels = [s.replace('_', ' ') for s in labels]

    zeroVec = np.zeros([len(refLabels), len(refLabels)])
    zeroVec[np.ix_(np.sort(ia), np.sort(ia))] = connectivity

    connectivityFilled = zeroVec

    fig, ax = plt.subplots()

    sns.heatmap(connectivityFilled)
    ax.axis('tight')

    # Set labels
    ax.set(xticks=np.arange(len(labels)), xticklabels=labels,
           yticks=np.arange(len(labels)), yticklabels=labels)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")

    ax.set_title("DTI conncectivity between ARA regions")
    plt.show()

    return connectivity


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Visulize mat file of DTI ')
    parser.add_argument('-i', '--inputMat', help='file name:fMRI mat-File')

    args = parser.parse_args()

    inputPath = None
    if args.inputMat is not None and args.inputMat is not None:
        inputPath = args.inputMat
    if not os.path.exists(inputPath):
        sys.exit("Error: %s path is not an existing directory." % (args.inputPath,))

    inputPath = args.inputMat

    # generate Matrix
    matrixMaker(inputPath)
