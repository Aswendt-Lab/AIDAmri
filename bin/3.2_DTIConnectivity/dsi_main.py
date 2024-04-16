#!/opt/env/bin/python
"""
Created on 10/08/2017

@author: Niklas Pallast
Neuroimaging & Neuroengineering
Department of Neurology
University Hospital Cologne
"""
from __future__ import print_function

import argparse
import os
import glob
import dsi_tools
import shutil

if __name__ == '__main__':
    # default dsi studio directory
    f = open(os.path.join(os.getcwd(), "dsi_studioPath.txt"), "r")
    dsi_studio = f.read().split("\n")[0]
    f.close()

    # default b-table in input directory
    b_table = os.path.abspath(os.path.join(os.getcwd(), os.pardir,os.pardir)) + '/lib/DTI_Jones30.txt'

    # default connectivity directory relative to input directory
    dir_con = r'connectivity'

    # Defining CLI flags
    parser = argparse.ArgumentParser(description='Get connectivity of DTI dataset')
    requiredNamed = parser.add_argument_group('Required named arguments')
    requiredNamed.add_argument('-i',
                               '--file_in',
                               help = 'path to the raw NIfTI DTI file (ends with *dwi.nii.gz)',
                               required=True
                               )
    parser.add_argument('-b',
                        '--b_table',
                        default='auto',  # Default to 'auto' for automatic selection
                        help='Specify the b-table source: "auto" (will look for bvec and bval, create the btable. If val or vec can not be found, it uses the Jones30 file)'
                        )
    parser.add_argument('-o',
                        '--optional',
                        nargs = '*',
                        help = 'Optional arguments.\n\t"fa0": Renames the FA metric data to former DSI naming convention.\n\t"nii_gz": Converts ROI labeling relating files from .nii to .nii.gz format to match former data structures.'
                        )    
    args = parser.parse_args()
        
     # Determine the btable source based on the -b option
    if args.b_table.lower() == 'auto':
        # Use the merge_bval_bvec_to_btable function with folder_path as file_in
        b_table = dsi_tools.merge_bval_bvec_to_btable(os.path.dirname(args.file_in))
        if b_table is False:
        # Use the default "Jones30" btable
            b_table = os.path.abspath(os.path.join(os.getcwd(), os.pardir, os.pardir)) + '/lib/DTI_Jones30.txt'


    # Preparing directories
    file_cur = os.path.dirname(args.file_in)
    dsi_path = os.path.join(file_cur, 'DSI_studio')
    mcf_path = os.path.join(file_cur, 'mcf_Folder')
    dir_mask = glob.glob(os.path.join(dsi_path, '*BetMask_scaled.nii'))
    if not dir_mask:
        dir_mask = glob.glob(os.path.join(dsi_path, '*BetMask_scaled.nii.gz')) # check for ending (either .nii or .nii.gz)
    dir_mask = dir_mask[0]

    dir_out = args.file_in

    if os.path.exists(mcf_path):
        shutil.rmtree(mcf_path)
    os.mkdir(mcf_path)
    file_in = dsi_tools.fsl_SeparateSliceMoCo(args.file_in, mcf_path)
    dsi_tools.srcgen(dsi_studio, file_in, dir_mask, dir_out, b_table)
    file_in = os.path.join(file_cur,'fib_map')

    # Fiber tracking
    dir_out = os.path.dirname(args.file_in)
    dsi_tools.tracking(dsi_studio, file_in)

    # Calculating connectivity
    suffixes = ['*StrokeMask_scaled.nii', '*parental_Mask_scaled.nii', '*Anno_scaled.nii', '*AnnoSplit_parental_scaled.nii']
    for f in suffixes:
        dir_seeds = glob.glob(os.path.join(file_cur, 'DSI_studio', f))
        if not dir_seeds:
            dir_seeds = glob.glob(os.path.join(file_cur, 'DSI_studio', f + '.gz')) # check for ending (either .nii or .nii.gz)
        if not dir_seeds:
            continue
        dir_seeds = dir_seeds[0]
        dsi_tools.connectivity(dsi_studio, file_in, dir_seeds, dir_out, dir_con)

    # rename files to reduce path length
    confiles = os.path.join(file_cur,dir_con)
    data_list = os.listdir(confiles)
    for filename in data_list:
        splittedName = filename.split('.src.gz.dti.fib.gz.trk.gz.')
        if len(splittedName)>1:
            newName = splittedName[1]
            newName = os.path.join(confiles,newName)
            if os.path.isfile(newName):
                os.remove(newName)
            oldName = os.path.join(confiles,filename)
            os.rename(oldName,newName)

    # Including optional arguments regarding deprecated terminology
    if args.optional is not None:
        file_list = os.listdir(dsi_path)
        for f in file_list:

            # fa0 was a former term used in earlier DSI-studio versions; the '0' in fa0 referred to the first fiber track. However, DTI can only result in one track, therefore only one fractional anisotropy value per voxel is given, thus the collective values are referred to as fa. With the 'fa0' flag toggled on, the 'fa' data file is renamed to the former naming convention (fa0).
            if 'fa0' in [s.lower() for s in args.optional] and f.endswith('fa.nii.gz'):
                newName = f.split('fa.nii.gz')[0] + 'fa0.nii.gz'
                newName = os.path.join(dsi_path, newName)
                oldName = os.path.join(dsi_path, f)
                if os.path.isfile(newName):
                    os.remove(newName)
                os.rename(oldName, newName)
            
            # Due to changes in ROI annotations the corresponding files are saved as '.nii' files as opposed to '.nii.gz' files in earlier versions of DSI studio. With the 'nii_gz' flag toggled on, the '.nii' files are renamed to '.nii.gz'.
            if 'nii_gz' in args.optional and f.endswith('.nii'):
                newName = f + '.gz'
                newName = os.path.join(dsi_path, newName)
                oldName = os.path.join(dsi_path, f)
                if os.path.isfile(newName):
                    os.remove(newName)
                os.rename(oldName, newName)

