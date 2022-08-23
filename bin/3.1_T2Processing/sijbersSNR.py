""""" Sijbers's method
sijbers2007automatic,
title={Automatic estimation of the noise variance from the histogram of a magnetic resonance image},
author={Sijbers, Jan and Poot, Dirk and den Dekker, Arnold J and Pintjens, Wouter},
journal={Physics in medicine and biology},
volume={52},
number={5},
pages={1335},
year={2007},
publisher={IOP Publishing} """

from math import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize


def calcSNR(img, show, fac):
    # Normalize input dataset and plot histogram
    # img = np.fliplr(img)

    img = img.astype(int)
    maxi = img.max()
    imgFlat = img.flatten(2)
    imgNorm = imgFlat / maxi
    bins = ceil(sqrt(imgNorm.size)) * fac
    binCount, binLoc = np.histogram(imgNorm, int(bins))

    estStd = np.argmax(binCount)
    fc = binLoc[2 * estStd]
    [n, l] = np.histogram(imgNorm[imgNorm <= fc], int(bins))

    Nk = np.sum(n)
    K = bins

    mlfunc = lambda x: maxLikelihood(x, Nk, K, l, n)
    sigma0 = binLoc[estStd]
    out = scipy.optimize.fmin(func=mlfunc, x0=sigma0, disp=False)

    estStdNorm = out
    estStd = (out * maxi) / 10

    snrMap = np.sqrt(abs(np.square(img) - (np.square(estStd)))) / estStd

    if show > 0:
        if len(img.shape) == 2:
            figSijbers = plt.figure(3)
            plt.imshow(snrMap)
            plt.show()
        elif len(img.shape) == 3:
            figSijbers = plt.figure(3)
            plt.imshow(snrMap[:, :, int(np.ceil(len(img.shape) / 2))])
            plt.show()

    return snrMap, estStd, estStdNorm


def maxLikelihood(x, Nk, K, l, n):
    y = Nk * np.log(np.exp(-l[0] ** 2 / (2 * x ** 2)) - np.exp(-l[K] ** 2 / (2 * x ** 2))) \
        - np.sum(n[1:K] * np.log(np.exp(-l[0:K - 1] ** 2 / (2 * x ** 2)) - np.exp(-l[1:K] ** 2. / (2 * x ** 2))))
    return y
