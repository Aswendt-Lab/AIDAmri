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

import proc_tools as pt

def get_mat_flip_x_z(data_dims, voxel_dims):
    mat = np.zeros((4, 4), dtype=np.float64)
    mat[0, 0] = -1
    mat[1, 1] = 1
    mat[2, 2] = -1
    mat[0, 3] = (data_dims[0] - 1) * voxel_dims[0]
    mat[2, 3] = (data_dims[2] - 1) * voxel_dims[2]
    mat[3, 3] = 1

    return mat

def get_mat_voxel_to_world(voxel_dims, origin=(0.0, 0.0, 0.0)):
    mat = np.zeros((4, 4), dtype=np.float64)
    mat[0, 0] = voxel_dims[0]
    mat[1, 1] = voxel_dims[1]
    mat[2, 2] = voxel_dims[2]
    mat[0, 3] = -origin[0] * voxel_dims[0] 
    mat[1, 3] = -origin[1] * voxel_dims[1]
    mat[2, 3] = -origin[2] * voxel_dims[2]
    mat[3, 3] = 1

    return mat

def get_mat_world_to_voxel(voxel_dims, origin=(0.0, 0.0, 0.0)):
    mat = np.zeros((4, 4), dtype=np.float64)
    mat[0, 0] = 1.0 / voxel_dims[0]
    mat[1, 1] = 1.0 / voxel_dims[1]
    mat[2, 2] = 1.0 / voxel_dims[2]
    mat[0, 3] = origin[0] * voxel_dims[0]
    mat[1, 3] = origin[1] * voxel_dims[1]
    mat[2, 3] = origin[2] * voxel_dims[2]
    mat[3, 3] = 1

    return mat

def make_matrix(lines):
    mat = np.zeros((4, 4), dtype=np.float64)
    for k in range(4):
        mat[k] = np.array(lines[k].split(), dtype=np.float64)

    return mat

def matrix_to_text(mat):
    return '\n'.join('  '.join(str(x) for x in mat[y]) + '  ' for y in range(mat.shape[0]))

def interp_nearest(data, v):
    vb = np.ones(v.shape[1], dtype=np.bool)
    v0 = np.int32(np.floor(v))
    for i in range(3):
        v1 = v0[i]
        n = data.shape[i]
        vb[np.logical_or(v1<0, v1>=n-1)] = 0
        v0[i, v1>=n-1] = n - 2
    v0[v0<0] = 0
    d0 = v - v0
    v0 = v0 + np.int32(d0 > 0.5)

    return data[v0[0], v0[1], v0[2]] * vb

def interp_trilinear(data, v):
    vb = np.ones(v.shape[1], dtype=np.bool)
    v0 = np.int32(np.floor(v))
    for i in range(3):
        v1 = v0[i]
        n = data.shape[i]
        vb[np.logical_or(v1<0, v1>=n-1)] = 0
        v0[i, v1>=n-1] = n - 2
    v0[v0<0] = 0
    v1 = v0 + 1
    # processing x
    d0 = v[0] - v0[0]
    d1 = 1 - d0
    c00 = data[v0[0], v0[1], v0[2]] * d1 + data[v1[0], v0[1], v0[2]] * d0
    c10 = data[v0[0], v1[1], v0[2]] * d1 + data[v1[0], v1[1], v0[2]] * d0
    c01 = data[v0[0], v0[1], v1[2]] * d1 + data[v1[0], v0[1], v1[2]] * d0
    c11 = data[v0[0], v1[1], v1[2]] * d1 + data[v1[0], v1[1], v1[2]] * d0
    # processing y
    d0 = v[1] - v0[1]
    c00 = (c10 - c00) * d0 + c00
    c01 = (c11 - c01) * d0 + c01
    # processing z
    d0 = v[2] - v0[2]
    c00 = (c01 - c00) * d0 + c00

    return c00 * vb

def xfm_serial(data, matrix, xfm_dims, voxel_dims_xfm, voxel_dims, interp=0, inverse=False, origin=(0.0, 0.0, 0.0)):
    interpolation = [interp_nearest, interp_trilinear][interp]
    if not inverse:
        matrix = np.linalg.inv(matrix)
        #print("inverse matrix:", matrix.shape); print(matrix)
    nx, ny = (xfm_dims[0], xfm_dims[1])
    v2w = get_mat_voxel_to_world(voxel_dims_xfm, origin)
    w2v = get_mat_world_to_voxel(voxel_dims, origin)
    #print("v2w:", v2w.shape); print(v2w)
    #print("w2v:", w2v.shape); print(w2v)
    v = np.flipud(np.indices((1, 1, ny, nx), dtype=np.float32)).reshape(4, nx * ny)
    v[3] = 1
    #print("v:", v.shape); print(v)
    #xfm = np.ndarray(shape=xfm_dims, dtype=np.float32)
    xfm = np.zeros(xfm_dims, dtype=data.dtype, order='F')
    for slc in zrange(xfm_dims[2]):
        #print("Slice:", slc)
        v[2] = slc
        w = np.dot(v2w, v) # voxel to world
        xfm_v = np.dot(w2v, np.dot(matrix, w)) # world to voxel
        #xfm[:, :, slc] = interpolation(data, xfm_v[:3]).reshape(ny, nx).T
        xfm[:, :, slc] = interpolation(data, xfm_v[:3]).reshape(nx, ny, order='F').astype(data.dtype)

    return xfm

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('in_data', help='input data file name')
    parser.add_argument('in_matrix', help='input matrix file name')
    parser.add_argument('-o', '--out_data', help='output data file name', required=True)
    parser.add_argument('-s', '--out_shape', help='output shape (e.g. "256,256,48")')
    parser.add_argument('-d', '--out_dims', help='output voxel dimensions (in mm)')
    parser.add_argument('-n', '--interp', type=int, choices=[0, 1], help='interpolation method (0: nearestneighbour, 1: trilinear)')
    parser.add_argument('-i', '--inverse', action='store_true', help='inverse transformation')
    args = parser.parse_args()

    # input data NIfTI file
    if not os.path.isfile(args.in_data):
        sys.exit("Error: '%s' is not a regular file." % (args.in_data,))

    # input matrix text file
    if not os.path.isfile(args.in_matrix):
        sys.exit("Error: '%s' is not a regular file." % (args.in_matrix,))

    # output shape
    out_shape = '' if args.out_shape is None else args.out_shape

    # output voxel dimensions
    out_dims = '' if args.out_dims is None else args.out_dims

    # interpolation method (0: nearestneighbour, 1: trilinear)
    interp = 1 if args.interp is None else args.interp

    # read input data
    data, voxel_dims = pt.read_data(args.in_data)

    # read matrix text file
    matrix = make_matrix(pt.read_text(args.in_matrix))

    data_dims_xfm = tuple(int(v) for v in out_shape.split(',')) if len(out_shape) > 0 else data.shape[:]
    voxel_dims_xfm = tuple(float(v) for v in out_dims.split(',')) if len(out_dims) > 0 else voxel_dims[:]

    # transform input data
    data_xfm = xfm_serial(data, matrix, data_dims_xfm, voxel_dims_xfm, voxel_dims, interp=interp, inverse=args.inverse)

    # save transformed data as NIfTI file
    pt.save_data(data_xfm, voxel_dims_xfm, args.out_data, dtype=None)
