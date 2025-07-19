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
    parser.add_argument('-r',
                        '--recon_method',
                        default='dti',
                        help='Specify diffusion reconstruction method ("gqi" or default "dti").',
                        required=False
                       )
    parser.add_argument('-v',
                        '--vivo',
                        default='in_vivo',
                        help='Specify in vivo or ex vivo data to adjust sampling length ratio (param0). "in_vivo" param0=1.25 (default), "ex_vivo" param0=0.60.',
                        required=False
                       )
    parser.add_argument('-m',
                        '--make_isotropic',
                        default=0,
                        help='Specify an isotropic voxel size in mm for resampling. Default 0 = no resampling. "auto" uses nibabel to read the NIFTI header for the minimum voxel size',
                        required=False
                       )
    parser.add_argument('-t',
                        '--track_params',
                        default='default',
                        help='Specify tracking parameters from a pre-defined set ("aida_optimized", "rat", or "mouse") or as a list of values for fiber_count, interpolation, step_size, turning_angle, check_ending, fa_threshold, smoothing, min_length, and max_length.',
                        required=False
                       )
    parser.add_argument('-y',
                        '--flip_image_y',
                        default=False,
                        help='Specify whether to flip the image in the y-direction. Default is None (no flip). Set to "true" to flip the image.',
                        required=False
                       )
    parser.add_argument('-template',
                        '--template',
                        default=1,
                        help='Specify the template to use for the reconstruction. Default is 1 (mouse). Other options are "Rat" (5) or "Mouse" (1).',
                        required=False
                       )
    parser.add_argument('-thread_count',
                        '--thread_count',
                        default=1,
                        help='Specify the number of threads to use for fiber tracking. Default is 1.',
                        required=False
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

    make_isotropic=0
    if args.make_isotropic != 0:
        make_isotropic=args.make_isotropic
    
    flip_image_y = False
    if args.flip_image_y is None:
        flip_image_y = False
    elif str(args.flip_image_y).lower() == 'true':
        flip_image_y = True
    
    template = 1
    if args.template.lower() == 'rat':
        template = 5
    elif args.template.lower() == 'mouse':
        template = 1
    else:
        try:
            template = int(args.template)
        except ValueError:
            print(f"Invalid template value: {args.template}. Using default template 1 (mouse).")
            template = 1

    if os.path.exists(mcf_path):
        shutil.rmtree(mcf_path)
    os.mkdir(mcf_path)
    file_in = dsi_tools.fsl_SeparateSliceMoCo(args.file_in, mcf_path)
    voxel_size = dsi_tools.srcgen(dsi_studio, file_in, dir_mask, dir_out, b_table, args.recon_method, args.vivo, make_isotropic, flip_image_y, template)
    file_in = os.path.join(file_cur,'fib_map')

    track_param = args.track_params

    # Fiber tracking
    dir_out = os.path.dirname(args.file_in)
    dsi_tools.tracking(dsi_studio, file_in, track_param, voxel_size, args.thread_count)

    # Calculating connectivity
    suffixes = ['*StrokeMask_scaled.nii', '*parental_Mask_scaled.nii', '*Anno_scaled.nii', '*AnnoSplit_parental_scaled.nii']
    for f in suffixes:
        dir_seeds = glob.glob(os.path.join(file_cur, 'DSI_studio', f))
        if not dir_seeds:
            dir_seeds = glob.glob(os.path.join(file_cur, 'DSI_studio', f + '.gz')) # check for ending (either .nii or .nii.gz)
        if not dir_seeds:
            continue
        dir_seeds = dir_seeds[0]
        dsi_tools.connectivity(dsi_studio, file_in, dir_seeds, dir_out, dir_con, make_isotropic, flip_image_y)

    # rename files to reduce path length
    confiles = os.path.join(file_cur,dir_con)
    data_list = os.listdir(confiles)
    for filename in data_list:
        if args.recon_method == "dti":
            splittedName = filename.split('.src.gz.dti.fib.gz.trk.gz.')
        elif args.recon_method == "gqi":
            splittedName = filename.split('.src.gz.gqi.fib.gz.trk.gz.')
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

