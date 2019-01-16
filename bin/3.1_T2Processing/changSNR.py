""" Changs's method
{chang2005automatic,
title={An automatic method for estimating noise-induced signal variance in magnitude-reconstructed magnetic resonance images},
author={Chang, Lin-Ching and Rohde, Gustavo K and Pierpaoli, Carlo},
booktitle={Medical Imaging},
pages={1136--1142},
year={2005},
organization={International Society for Optics and Photonics}"""

from math import *
import numpy as np
import matplotlib.pyplot as plt


def calcSNR(img, show, fac):
    # Normalize input dataset and plot histogram
    # img = np.fliplr(img)

    img = img.astype(int)
    maxi = img.max()
    imgFlat = img.flatten(2)
    imgNorm = imgFlat / maxi
    bins = ceil(sqrt(imgNorm.size)) * fac
    binCount, binLoc = np.histogram(imgNorm, int(bins))
    n = len(imgNorm)

    estStd = np.argmax(binCount)
    estStd = (estStd) / binCount.shape

    x = np.linspace(0, 1, bins)
    fhat = np.zeros([1, len(x)])

    h = 1.06 * n ** (-1 / 5) * estStd

    # define function
    gauss = lambda x: gaussianFct(x)

    for i in range(n):
        # get each kernel function evaluated at x
        # centered at data
        f = gauss((x - imgNorm[i]) / h)
        # plot(x, f / (n * h))
        fhat = fhat + f

    fhat = fhat / (n * h)
    # SNR-Map
    maxPos = np.argmax(fhat)
    estStdNorm = binLoc[maxPos]
    estStd = (binLoc[maxPos] * maxi) / 10

    snrMap = np.sqrt(abs(np.square(img) - (np.square(estStd)))) / estStd

    if show > 0:
        if len(img.shape) == 2:
            figChang = plt.figure(3)
            plt.imshow(snrMap)
            plt.show()
        elif len(img.shape) == 3:
            figChang = plt.figure(3)
            plt.imshow(snrMap[:, :, int(np.ceil(len(img.shape) / 2))])
            plt.show()

    return snrMap, estStd, estStdNorm


def gaussianFct(x):
    y = 1 / sqrt(2 * pi) * np.exp((-(np.square(x))) / 2)
    return y
