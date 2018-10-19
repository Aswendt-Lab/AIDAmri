"""
Created on 10/08/2017

@author: Niklas Pallast
Neuroimaging & Neuroengineering
Department of Neurology
University Hospital Cologne

"""

import numpy as np
import scipy.ndimage

def applyFilter(im, num_iter, delta_t, kappa, option):


    # Convert input image to float.
    im.astype(float)

    # PDE(partial differential equation) initial condition.
    diff_im = im


    # Center pixel distances.
    dx = 1
    dy = 1
    dd = np.sqrt(2)

    # 2D convolution masks - finite differences.
    hN = np.array(([0, 1, 0],[0, -1, 0],[0, 0, 0]))

    hS = np.array(([0, 0, 0],[0, -1, 0],[0, 1, 0]))

    hE = np.array(([0, 0, 0],[0, -1, 1],[0, 0, 0]))

    hW = np.array(([0, 0, 0],[1, -1, 0],[0, 0, 0]))

    hNE = np.array(([0, 0, 1],[0, -1, 0],[0, 0, 0]))

    hSE = np.array(([0, 0, 0],[0, -1, 0],[0, 0, 1]))

    hSW = np.array(([0, 0, 0],[0, -1, 0],[1, 0, 0]))

    hNW = np.array(([1, 0, 0],[0, -1, 0],[0, 0, 0]))




    for t in range(num_iter):
        nablaN = scipy.ndimage.convolve(diff_im, hN, mode='nearest')
        nablaS = scipy.ndimage.convolve(diff_im, hS, mode='nearest')
        nablaE = scipy.ndimage.convolve(diff_im, hE, mode='nearest')
        nablaW = scipy.ndimage.convolve(diff_im, hW, mode='nearest')

        nablaNE = scipy.ndimage.convolve(diff_im, hNE, mode='nearest')
        nablaSE = scipy.ndimage.convolve(diff_im, hSE, mode='nearest')
        nablaSW = scipy.ndimage.convolve(diff_im, hSW, mode='nearest')
        nablaNW = scipy.ndimage.convolve(diff_im, hNW, mode='nearest')

        # Diffusion function.
        if option == 1:
            cN = np.exp(-(nablaN / kappa)**2)
            cS = np.exp(-(nablaS / kappa)**2)
            cW = np.exp(-(nablaW / kappa)**2)
            cE = np.exp(-(nablaE / kappa)**2)
            cNE = np.exp(-(nablaNE / kappa)**2)
            cSE = np.exp(-(nablaSE / kappa)**2)
            cSW = np.exp(-(nablaSW / kappa)**2)
            cNW = np.exp(-(nablaNW / kappa)**2)
        elif option == 2:
            cN = 1 / (1 + (nablaN / kappa)**2)
            cS = 1 / (1 + (nablaS / kappa)**2)
            cW = 1 / (1 + (nablaW / kappa)**2)
            cE = 1 / (1 + (nablaE / kappa)**2)
            cNE = 1 / (1 + (nablaNE / kappa)**2)
            cSE = 1 / (1 + (nablaSE / kappa)**2)
            cSW = 1 / (1 + (nablaSW / kappa)**2)
            cNW = 1 / (1 + (nablaNW / kappa)**2)

        #Discrete PDE solution
        diff_im = diff_im + delta_t * (
            (1 / (dy **2)) * cN * nablaN + (1 / (dy ** 2)) * cS * nablaS +
            (1 / (dx ** 2)) * cW * nablaW + (1 / (dx ** 2)) * cE * nablaE +
            (1 / (dd ** 2)) * cNE * nablaNE + (1 / (dd ** 2)) * cSE * nablaSE +
            (1 / (dd ** 2)) * cSW * nablaSW + (1 / (dd ** 2)) * cNW * nablaNW)

    return diff_im