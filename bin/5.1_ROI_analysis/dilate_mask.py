'''
Created on 20.08.2020

Author:
Michael Diedenhofen
Max Planck Institute for Metabolism Research, Cologne
'''

from __future__ import print_function

try:
    zrange = xrange
except NameError:
    zrange = range

import os
import sys

import numpy as np

from scipy import ndimage

import proc_tools as pt

def circle_mask(n=8):
    nn = 2 * n + 1
    xx, yy = np.mgrid[:nn, :nn]
    circle = ((xx - n) ** 2 + (yy - n) ** 2) < (n * n + 2)
    #print(circle.astype(np.int16))

    return circle

def dilate_repeat(image, connectivity=1, n=8):
    dilated = np.copy(image).astype(np.bool)
    struct = ndimage.generate_binary_structure(2, connectivity)
    for _ in zrange(n):
        dilated = ndimage.binary_dilation(dilated, structure=struct)
    dilated = np.subtract(dilated.astype(image.dtype), image)

    return dilated

def dilate_struct(image, struct):
    dilated = np.copy(image).astype(np.bool)
    dilated = ndimage.binary_dilation(dilated, structure=struct)
    dilated = np.subtract(dilated.astype(image.dtype), image)

    return dilated

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Dilate input mask.')
    parser.add_argument('in_mask', help='input mask file name')
    parser.add_argument('-o', '--out_mask', help='output mask file name', required=True)
    args = parser.parse_args()

    # input mask file
    if not os.path.isfile(args.in_mask):
        sys.exit("Error: '%s' is not a regular file." % (args.in_mask,))

    # read input mask data
    data, voxel_dims = pt.read_data(args.in_mask)

    struct = circle_mask()
    for k in zrange(data.shape[2]):
        image = data[:, :, k]
        if np.any(image.astype(np.bool)):
            #data[:, :, k] = dilate_repeat(image)
            data[:, :, k] = dilate_struct(image, struct)

    # save mask data as NIfTI file
    pt.save_data(data, voxel_dims, args.out_mask, dtype=None)
