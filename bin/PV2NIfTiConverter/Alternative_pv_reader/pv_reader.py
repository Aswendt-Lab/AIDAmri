'''
Created on 19.10.2020

Author:
Michael Diedenhofen
Max Planck Institute for Metabolism Research, Cologne

Read Bruker ParaVision data (2dseq) and save as NIfTI file.
Create a b-table text file with b-values and directions for diffusion data.
'''

from __future__ import print_function

try:
    zrange = xrange
except NameError:
    zrange = range

VERSION = 'pv_reader.py v 1.1.2 20201019'

import os
import sys

import numpy as np
import nibabel as nib
import nibabel.nifti1 as nii

import pv_parser as par

class ParaVision:
    """
    Read ParaVision data and save as NIfTI file
    """

    def __init__(self, procfolder, rawfolder, study, expno, procno):
        self.procfolder = procfolder
        self.rawfolder = rawfolder
        self.study = study
        self.expno = int(expno)
        self.procno = int(procno)
        self.name = '.'.join([study, str(expno), str(procno)])

    def __check_params(self, params_name, labels):
        misses = [label for label in labels if label not in getattr(self, params_name)]
        if len(misses) > 0:
            sys.exit("Missing labels in %s: %s" % (params_name, str(misses),))

    def __check_path(self, header_path):
        header_path = header_path.split('/')
        study, expno, procno = (header_path[-5], int(header_path[-4]), int(header_path[-2]))
        if self.study != study:
            print("Warning: Study '%s' differs from '%s' in the visu_pars header." % (self.study, study), file=sys.stderr)
        if self.expno != expno:
            print("Warning: Experiment number %s differs from %s in the visu_pars header." % (self.expno, expno), file=sys.stderr)
        if self.procno != procno:
            print("Warning: Processed images number %s differs from %s in the visu_pars header." % (self.procno, procno), file=sys.stderr)

    def __get_data_dims(self):
        labels_visu_pars = ['VisuCoreDim', 'VisuCoreSize', 'VisuCoreWordType', 'VisuCoreByteOrder']
        self.__check_params('visu_pars', labels_visu_pars)

        #VisuCoreFrameCount = self.visu_pars.get('VisuCoreFrameCount') # Number of frames
        VisuCoreDim = self.visu_pars.get('VisuCoreDim')
        VisuCoreSize = self.visu_pars.get('VisuCoreSize')
        VisuCoreDimDesc = self.visu_pars.get('VisuCoreDimDesc')
        VisuCoreWordType = self.visu_pars.get('VisuCoreWordType')
        #VisuCoreByteOrder = self.visu_pars.get('VisuCoreByteOrder')
        #VisuFGOrderDescDim = self.visu_pars.get('VisuFGOrderDescDim')
        VisuFGOrderDesc = self.visu_pars.get('VisuFGOrderDesc')

        dim_desc = None if VisuCoreDimDesc is None else VisuCoreDimDesc[0]

        # FrameGroup dimensions and names
        if (VisuFGOrderDesc is not None) and len(VisuFGOrderDesc) > 0:
            #fg_dims = list(map(lambda item: int(item[0]), VisuFGOrderDesc))
            #fg_names = list(map(lambda item: str(item[1]), VisuFGOrderDesc))
            fg_dims = [int(item[0]) for item in VisuFGOrderDesc]
            fg_names = [str(item[1]) for item in VisuFGOrderDesc]
            fg_names = [par.extract_jcamp_strings(item, get_all=False) for item in fg_names]
        else:
            fg_dims = []
            fg_names = []

        # Data dimensions
        data_dims = list(map(int, VisuCoreSize)) + fg_dims

        # FrameGroup FG_SLICE index
        fg_index, fg_slice = (None, None)
        if VisuCoreDim == 2:
            fg_slices = ('FG_SLICE', 'FG_IRMODE')
            fg_indices = [fg_names.index(x) for x in fg_slices if x in fg_names]
            if len(fg_indices) > 0:
                fg_index = fg_indices[0]
                fg_slice = fg_slices[fg_index]
                fg_index += VisuCoreDim

        # ParaVision to NumPy data-type conversion
        if VisuCoreWordType == '_8BIT_UNSGN_INT':
            data_type = 'uint8'
        elif VisuCoreWordType == '_16BIT_SGN_INT':
            data_type = 'int16'
        elif VisuCoreWordType == '_32BIT_SGN_INT':
            data_type = 'int32'
        elif VisuCoreWordType == '_32BIT_FLOAT':
            data_type = 'float32'
        else:
            sys.exit("The data format is not correct specified.")

        return (data_dims, data_type, dim_desc, fg_index, fg_slice)

    def __get_voxel_dims(self, data_dims, scale=1.0):
        labels_visu_pars = ['VisuCoreExtent']
        self.__check_params('visu_pars', labels_visu_pars)

        ACQ_slice_sepn = self.acqp.get('ACQ_slice_sepn')
        #PVM_SPackArrSliceGap = self.method.get('PVM_SPackArrSliceGap')
        PVM_SPackArrSliceDistance = self.method.get('PVM_SPackArrSliceDistance')
        VisuCoreExtent = self.visu_pars.get('VisuCoreExtent')
        VisuCoreFrameThickness = self.visu_pars.get('VisuCoreFrameThickness')
        VisuCoreUnits = self.visu_pars.get('VisuCoreUnits')
        VisuCoreSlicePacksSliceDist = self.visu_pars.get('VisuCoreSlicePacksSliceDist')
        VisuAcqRepetitionTime = self.visu_pars.get('VisuAcqRepetitionTime')

        nd = min(len(data_dims), 4)
        dims = [1] * 4
        dims[:nd] = data_dims
        nx, ny, nz, nt = dims

        # Voxel dimensions
        if len(VisuCoreExtent) > 1:
            dx = scale * float(VisuCoreExtent[0]) / nx
            dy = scale * float(VisuCoreExtent[1]) / ny
        else:
            dx = 1.0
            dy = 0.0
        if len(VisuCoreExtent) > 2:
            dz = scale * float(VisuCoreExtent[2]) / nz
        elif ACQ_slice_sepn is not None: # Slice thickness inclusive gap
            dz = scale * float(ACQ_slice_sepn[0])
        elif PVM_SPackArrSliceDistance is not None: # Slice thickness inclusive gap
            dz = scale * float(PVM_SPackArrSliceDistance[0])
        elif VisuCoreSlicePacksSliceDist is not None: # Slice thickness inclusive gap (PV6)
            dz = scale * float(VisuCoreSlicePacksSliceDist[0])
        elif VisuCoreFrameThickness is not None: # Slice thickness
            dz = scale * float(VisuCoreFrameThickness[0])
        else:
            dz = 0.0
        if (VisuAcqRepetitionTime is not None) and (nt > 1):
            dt = float(VisuAcqRepetitionTime[0]) / 1000.0
        else:
            dt = 0.0
        #voxel_dims = [dx, dy, dz, dt][:nd]
        voxel_dims = [dx, dy, dz, dt]

        # Units of the voxel dimensions
        voxel_unit = par.extract_unit_string(par.extract_jcamp_strings(VisuCoreUnits[0], get_all=False))

        return (voxel_dims, voxel_unit)

    def __map_data(self, data, map_pv6):
        VisuCoreExtent = self.visu_pars.get('VisuCoreExtent')
        VisuCoreDataOffs = self.visu_pars.get('VisuCoreDataOffs')
        VisuCoreDataSlope = self.visu_pars.get('VisuCoreDataSlope')

        n = min(len(VisuCoreExtent), 3)
        dims = data.shape[n:]

        if VisuCoreDataOffs.size > 1:
            VisuCoreDataOffs = VisuCoreDataOffs.reshape(dims, order='F').astype(np.float32)
        else:
            VisuCoreDataOffs = np.float32(VisuCoreDataOffs[0])

        if VisuCoreDataSlope.size > 1:
            VisuCoreDataSlope = VisuCoreDataSlope.reshape(dims, order='F').astype(np.float32)
        else:
            VisuCoreDataSlope = np.float32(VisuCoreDataSlope[0])

        if map_pv6:
            data = data.astype(np.float32) / VisuCoreDataSlope
            data = data + VisuCoreDataOffs
        else:
            data = data.astype(np.float32) * VisuCoreDataSlope
            data = data + VisuCoreDataOffs

        return (data, 'float32')

    def __make_subfolder(self, subfolder=''):
        procfolder = os.path.join(self.procfolder, self.study)
        if not os.path.isdir(procfolder):
            os.mkdir(procfolder)

        if len(subfolder) > 0:
            procfolder = os.path.join(self.procfolder, self.study, subfolder)
            if not os.path.isdir(procfolder):
                os.mkdir(procfolder)

        return procfolder

    def __get_matrix(self):
        VisuCoreOrientation = self.visu_pars.get('VisuCoreOrientation')
        VisuCoreOrientation = VisuCoreOrientation.flatten()[:9].astype(np.float32)
        VisuCoreOrientation = VisuCoreOrientation.reshape((3, 3), order='F')
        VisuCorePosition = self.visu_pars.get('VisuCorePosition')
        VisuCorePosition = VisuCorePosition.flatten()[:3].astype(np.float32)

        matrix = np.zeros((4, 4), dtype=np.float32)
        matrix[:3, :3] = VisuCoreOrientation
        matrix[:3, 3] = self.scale * VisuCorePosition
        matrix[3, 3] = 1

        return matrix

    def __save_matrix(self, matrix, procfolder, ext='mat'):
        lines = '\n'.join(('  '.join('%.12g' % (x,) for x in matrix[y]) + '  ') for y in range(matrix.shape[0]))

        fname = '.'.join([self.name, ext])
        fpath = os.path.join(procfolder, fname)
        print(fpath)

        # Open text file to write binary (Unix format)
        fid = open(fpath, 'wb')

        # Write text file
        for line in lines.splitlines():
            print(line, end=chr(10), file=fid)

        # Close text file
        fid.close()

    def read_2dseq(self, map_raw=False, map_pv6=False, roll_fg=False, squeeze=False, compact=False, swap_vd=False, scale=1.0):
        self.scale = float(scale)

        # Get acqp and method parameters
        datadir = os.path.join(self.rawfolder, self.study, str(self.expno))
        _header, self.acqp = par.read_param_file(os.path.join(datadir, 'acqp'))
        _header, self.method = par.read_param_file(os.path.join(datadir, 'method'))

        # Get visu_pars parameters
        datadir = os.path.join(self.rawfolder, self.study, str(self.expno), 'pdata', str(self.procno))
        #_header, self.d3proc = par.read_param_file(os.path.join(datadir, 'd3proc')) # Removed for PV6
        header, self.visu_pars = par.read_param_file(os.path.join(datadir, 'visu_pars'))

        self.__check_path(header['Path'])

        # Remove selected parameters from the visu_pars dictionary
        #if 'VisuCoreDataMin' in self.visu_pars: del self.visu_pars['VisuCoreDataMin']
        #if 'VisuCoreDataMax' in self.visu_pars: del self.visu_pars['VisuCoreDataMax']
        #if 'VisuCoreDataOffs' in self.visu_pars: del self.visu_pars['VisuCoreDataOffs']
        #if 'VisuCoreDataSlope' in self.visu_pars: del self.visu_pars['VisuCoreDataSlope']
        #if 'VisuAcqImagePhaseEncDir' in self.visu_pars: del self.visu_pars['VisuAcqImagePhaseEncDir']

        #VisuCoreFrameType = self.visu_pars.get('VisuCoreFrameType')
        #VisuCoreDiskSliceOrder = self.visu_pars.get('VisuCoreDiskSliceOrder')

        # Get data dimensions
        data_dims, data_type, dim_desc, fg_index, fg_slice = self.__get_data_dims()

        # Open 2dseq file
        path_2dseq = os.path.join(datadir, '2dseq')
        try:
            fid = open(path_2dseq, 'rb')
        except IOError as V:
            if V.errno == 2:
                sys.exit("Cannot open 2dseq file %s" % (path_2dseq,))
            else:
                raise

        # Read 2dseq file
        data = np.fromfile(fid, dtype=np.dtype(data_type)).reshape(data_dims, order='F')

        # Close 2dseq file
        fid.close()

        # Map to raw data range
        if map_raw:
            data, data_type = self.__map_data(data, map_pv6)

        # Move FrameGroup FG_SLICE axis to position 2
        self.roll_fg = False
        if roll_fg:
            if fg_index is None:
                print("Warning: Could not find FrameGroup.", file=sys.stderr)
            elif fg_index > 2:
                print("Warning: Move axis %d (FrameGroup %s) to position %d." % (fg_index, fg_slice, 2), file=sys.stderr)
                data = np.rollaxis(data, fg_index, 2)
                data_dims = list(data.shape)
                self.roll_fg = True
            else:
                print("Warning: Could not move FrameGroup %s." % (fg_slice,), file=sys.stderr)

        # Remove data dimensions of size 1
        if squeeze and (1 in data_dims):
            data = np.squeeze(data)
            data_dims = list(data.shape)

        # Reduce data dimensions to 4
        if compact and (len(data_dims) > 4):
            nt = int(np.prod(data_dims[3:]))
            data_dims[3:] = [nt]
            data = data.reshape(data_dims, order='F')

        # Get voxel dimensions
        voxel_dims, voxel_unit = self.__get_voxel_dims(data_dims, scale=self.scale)
        self.swap_vd = False
        if swap_vd:
            if (not self.roll_fg) and (fg_index != 2) and (len(data_dims) > 3):
                print("Warning: Swap third and fourth voxel dimension.", file=sys.stderr)
                voxel_dims[2:4] = voxel_dims[3:1:-1]
                self.swap_vd = True
            else:
                print("Warning: Could not swap third and fourth voxel dimension.", file=sys.stderr)

        pixdim = [0.0] * 8
        pixdim[0] = 1.0 # NIfTI qfac which is either -1 or 1
        pixdim[1:len(voxel_dims)+1] = voxel_dims

        # Info parameters
        self.data_dims = data_dims
        self.data_type = data_type
        self.voxel_dims = voxel_dims
        self.voxel_unit = voxel_unit

        # NIfTI image
        self.nifti_image = nii.Nifti1Image(data.reshape(data_dims, order='F'), None)

        # NIfTI header
        header = self.nifti_image.get_header()
        header.set_data_dtype(data.dtype)
        header.set_data_shape(data_dims)
        #header.set_zooms(voxel_dims)
        header['pixdim'] = pixdim
        if dim_desc != 'spectroscopic':
            header.set_xyzt_units(xyz=voxel_unit, t=None)
        #print("header:"); print(header)

    def save_nifti(self, ftype='NIFTI_GZ', subfolder=''):
        if ftype == 'NIFTI_GZ':
            ext = 'nii.gz'
        elif ftype == 'NIFTI':
            ext = 'nii'
        elif ftype == 'ANALYZE':
            ext = 'img'
        else:
            ext = 'nii.gz'

        fproc = self.__make_subfolder(subfolder=subfolder)

        fname = '.'.join([self.name, ext])
        fpath = os.path.join(fproc, fname)

        # Write NIfTI file
        nib.save(self.nifti_image, fpath)
        #self.nifti_image.to_filename(fpath)
        print(self.nifti_image.get_filename())

    def get_matrix(self):
        matrix = self.__get_matrix()

        return (matrix, np.linalg.inv(matrix))

    def save_matrix(self, subfolder=''):
        matrix = self.__get_matrix()

        fproc = self.__make_subfolder(subfolder=subfolder)

        self.__save_matrix(matrix, fproc, ext='omat')
        self.__save_matrix(np.linalg.inv(matrix), fproc, ext='imat')

    def save_table(self, eff_bval=False, subfolder=''):
        DwAoImages = int(self.method.get('PVM_DwAoImages'))
        DwNDiffDir = int(self.method.get('PVM_DwNDiffDir'))
        DwNDiffExpEach = int(self.method.get('PVM_DwNDiffExpEach'))
        #DwNDiffExp = int(self.method.get('PVM_DwNDiffExp'))
        #print("DwAoImages:", DwAoImages)
        #print("DwNDiffDir:", DwNDiffDir)
        #print("DwNDiffExpEach:", DwNDiffExpEach)
        #print("DwNDiffExp:", DwNDiffExp)

        nd = DwAoImages + DwNDiffDir * DwNDiffExpEach
        bvals = np.zeros(nd, dtype=np.float64)
        dwdir = np.zeros((nd, 3), dtype=np.float64)

        if eff_bval:
            DwEffBval = self.method.get('PVM_DwEffBval').astype(np.float64)
            #print("DwEffBval:"); print(DwEffBval)
            bvals[DwAoImages:] = DwEffBval[DwAoImages:]
        else:
            DwBvalEach = self.method.get('PVM_DwBvalEach').astype(np.float64)
            #print("DwBvalEach:", DwBvalEach)
            bvals[DwAoImages:] = np.tile(DwBvalEach, DwNDiffDir)

        DwDir = self.method.get('PVM_DwDir').astype(np.float64)
        #DwDir = DwDir.reshape((DwNDiffDir * DwNDiffExpEach, 3))
        #print("DwDir:"); print(DwDir)
        dwdir[DwAoImages:] = np.repeat(DwDir, DwNDiffExpEach, axis=0)

        fproc = self.__make_subfolder(subfolder=subfolder)

        fname = '.'.join([self.name, 'btable', 'txt'])
        fpath = os.path.join(fproc, fname)
        print(fpath)

        # Open btable file to write binary (Windows format)
        fid = open(fpath, 'wb')

        for i in zrange(nd):
            print("%.4f" % (bvals[i],) + " %.8f %.8f %.8f" % tuple(dwdir[i]), end="\r\n", file=fid)

        # Close file
        fid.close()

        fname = '.'.join([self.name, 'bvals', 'txt'])
        fpath = os.path.join(fproc, fname)
        print(fpath)

        # Open bvals file to write binary (Unix format)
        fid = open(fpath, 'wb')

        print(" ".join("%.4f" % (bvals[i],) for i in zrange(nd)), end=chr(10), file=fid)

        # Close bvals file
        fid.close()

        fname = '.'.join([self.name, 'bvecs', 'txt'])
        fpath = os.path.join(fproc, fname)
        print(fpath)

        # Open bvecs file to write binary (Unix format)
        fid = open(fpath, 'wb')

        for k in range(3):
            print(" ".join("%.8f" % (dwdir[i, k],) for i in zrange(nd)), end=chr(10), file=fid)

        # Close bvecs file
        fid.close()

def check_args(proc_folder, raw_folder, study, expno, procno):
    # processed data folder
    if not os.path.isdir(proc_folder):
        sys.exit("Error: '%s' is not an existing directory." % (proc_folder,))

    # raw data folder
    if not os.path.isdir(raw_folder):
        sys.exit("Error: '%s' is not an existing directory." % (raw_folder,))

    # study name
    path = os.path.join(raw_folder, study)
    if not os.path.isdir(path):
        sys.exit("Error: '%s' is not an existing directory." % (path,))

    # experiment number
    path = os.path.join(raw_folder, study, str(expno))
    if not os.path.isdir(path):
        sys.exit("Error: '%s' is not an existing directory." % (path,))

    # processed images number
    path = os.path.join(raw_folder, study, str(expno), 'pdata', str(procno))
    if not os.path.isdir(path):
        sys.exit("Error: '%s' is not an existing directory." % (path,))

def main():
    import argparse

    parser = argparse.ArgumentParser(description='Read ParaVision data and save as NIfTI file')
    parser.add_argument('proc_folder', help='processed data folder')
    parser.add_argument('raw_folder', help='raw data folder')
    parser.add_argument('study', help='study name')
    parser.add_argument('expno', help='experiment number')
    parser.add_argument('procno', help='processed (reconstructed) images number')
    parser.add_argument('-s', '--scale', default=1.0, help='voxel dimensions scale factor')
    parser.add_argument('-m', '--map_raw', action='store_true', help='map the data to get the real values')
    parser.add_argument('-p', '--map_pv6', action='store_true', help='map the data by dividing (ParaVision 6)')
    parser.add_argument('-r', '--roll_fg', action='store_true', help='move slice framegroup to third dimension')
    parser.add_argument('-q', '--squeeze', action='store_true', help='remove data dimensions of size 1')
    parser.add_argument('-c', '--compact', action='store_true', help='reduce data dimensions to 4')
    parser.add_argument('-v', '--swap_vd', action='store_true', help='swap third and fourth voxel dimension')
    parser.add_argument('-t', '--table', action='store_true', help='save b-values and diffusion directions')
    args = parser.parse_args()

    check_args(args.proc_folder, args.raw_folder, args.study, args.expno, args.procno)

    pv = ParaVision(args.proc_folder, args.raw_folder, args.study, args.expno, args.procno)
    pv.read_2dseq(map_raw=args.map_raw, map_pv6=args.map_pv6, roll_fg=args.roll_fg, squeeze=args.squeeze, compact=args.compact, swap_vd=args.swap_vd, scale=args.scale)
    pv.save_nifti(ftype='NIFTI_GZ')
    pv.save_matrix()
    if args.table:
        pv.save_table()

if __name__ == '__main__':
    main()
