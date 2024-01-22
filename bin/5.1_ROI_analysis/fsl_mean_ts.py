'''
Created on 31.08.2020

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
import scipy.io as sio

import proc_tools as pt

def mean_ts(path_data, path_mask, path_out, label_names):
    # Read 4D data file (NIfTI)
    data, _ = pt.read_data(path_data)
    file_data = os.path.basename(path_data)
    if len(data.shape) != 4:
        sys.exit("Error: %s is not 4D shape %s." % (file_data, str(data.shape)))

    # Read 4D mask file (NIfTI)
    mask, _ = pt.read_data(path_mask)
    file_mask = os.path.basename(path_mask)
    if len(mask.shape) != 4:
        sys.exit("Error: %s is not 4D shape $s." % (file_mask, str(mask.shape)))

    if data.shape[:3] != mask.shape[:3]:
        sys.exit("Error: %s %s and %s %s are not the same shape." % (file_data, str(data.shape[:3]), file_mask, str(mask.shape[:3])))

    #path_out_mat = path_out + '.mat'
    path_out_mat = os.path.join(os.path.dirname(path_out), os.path.basename(path_out) + '.mat')

    m = np.zeros((mask.shape[3], data.shape[3]), dtype=data.dtype)
    for k in zrange(mask.shape[3]):
        msk = np.array(mask[:, :, :, k]) > 0
        if data[msk, :].size > 0:
            m[k] = np.mean(data[msk, :], 0)

    mT = np.transpose(m)

    #np.savetxt(path_out, mT, fmt='%.4f', delimiter=' ')

    #s = [['%.4f' % (x,) for x in line] for line in mT.tolist()]
    s = [map(lambda x: '%.4f' % (x,), line) for line in mT.tolist()]

    #pt.save_csv(path_out, s)
    pt.save_text(path_out, s)

    sio.savemat(path_out_mat, dict([('matrix', mT), ('label', label_names)]))
    print(path_out_mat)

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('in_data', help='input 4D data (x, y, slc, rep)')
    parser.add_argument('in_mask', help='input 4D mask (x, y, slc, msk)')
    parser.add_argument('-o', '--out_text', help='output text file', required=True)
    args = parser.parse_args()

    # input data (NIfTI)
    if not os.path.isfile(args.in_data):
        sys.exit("Error: '%s' is not a regular file." % (args.in_data,))

    # input mask (NIfTI)
    if not os.path.isfile(args.in_mask):
        sys.exit("Error: '%s' is not a regular file." % (args.in_mask,))

    # output text file
    path_out = args.out_text if args.out_text is not None else os.path.abspath(os.path.join(args.in_data, os.pardir, 'MasksTCs.' + os.path.split(args.in_data)[1]))
    if path_out.endswith('.nii.gz'):
        path_out = path_out[:-7]
    elif path_out.endswith('.nii'):
        path_out = path_out[:-4]
    path_out = path_out + '.txt'

    label_names_2000 = pt.read_text(pt.path_label_names_2000)

    #print(get_date())

    mean_ts(args.in_data, args.in_mask, path_out, label_names_2000)
