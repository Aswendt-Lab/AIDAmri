""" Brummer's Method
brummer1993automatic,
title={Automatic detection of brain contours in MRI data sets},
author={Brummer, Marijn E and Mersereau, Russell M and Eisner, Robert L and Lewine, Richard RJ},
journal={IEEE Transactions on medical imaging},
volume={12},
number={2},
pages={153--166},
year={1993},
publisher={IEEE}


Created on 10/08/2017

@author: Niklas Pallast
Neuroimaging & Neuroengineering
Department of Neurology
University Hospital Cologne



"""

from math import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize, scipy.signal
def calcSNR(img,show,fac):
    # Normalize input dataset and plot histogram
    #img = np.fliplr(img)
    img = img.astype(int)
    maxi =  img.max()
    imgFlat = img.flatten()
    imgNorm= imgFlat/maxi
    bins  = ceil(sqrt(imgNorm.size))*fac
    binCount, binLoc = np.histogram(imgNorm, int(bins))



    maxRayl = max(binCount)
    estStd = np.argmax(binCount)
    cutOff = 2 * estStd
    estStd = (estStd)/len(binCount)


    # define function
    raylfunc = lambda x: rayl_2p(x, binLoc[0:cutOff-1] , binCount[0:cutOff-1])
    yout = scipy.optimize.fmin(func=raylfunc, x0=[maxRayl,estStd],disp=False)

    estStdNorm = yout[1]
    estStd = (yout[1] * maxi)/10

    snrMap = np.sqrt(abs(np.square(img) - (np.square(estStd))))/estStd

    if show > 0:
        if len(img.shape) == 2:
            plt.figure(3)
            plt.imshow(snrMap)
            plt.show()
        elif len(img.shape) == 3:
            plt.figure(3)
            plt.imshow(snrMap[:,:,int(np.ceil(len(img.shape)/2))])
            plt.show()

    return snrMap, estStd, estStdNorm


def rayl_2p(fitPar,x,data):
    ray=x/(fitPar[1]**2)*np.exp(-np.square(x)/(2*fitPar[1]**2))
    err=sum((fitPar[0]*ray-data)**2)
    return err