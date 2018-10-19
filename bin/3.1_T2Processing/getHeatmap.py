"""
Created on 10/08/2017

@author: Niklas Pallast
Neuroimaging & Neuroengineering
Department of Neurology
University Hospital Cologne
"""

import matplotlib.pyplot as plt
import numpy as np

def heatMap(incidenceMap,araVol):
    fig = plt.figure(frameon=False)
    im = []
    for i in range(12):
        t = 1 + i
        fig.add_subplot(3, 4, t)
        plt.imshow(np.transpose(incidenceMap[:, :, t * 16]), cmap='gnuplot')

        if i == 8:
            im = plt.imshow(np.transpose(incidenceMap[:, :, t * 16]), cmap='gnuplot')
            plt.imshow(np.transpose(araVol[:, :, t * 16]), alpha=0.55, cmap='gray')
        else:
            plt.imshow(np.transpose(incidenceMap[:, :, t * 16]), cmap='gnuplot')
            plt.imshow(np.transpose(araVol[:, :, t * 16]), alpha=0.55, cmap='gray')
        plt.axis('off')
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(im, cax=cbar_ax)
    plt.show()