"""
Created on 10/08/2017

@author: Niklas Pallast
Neuroimaging & Neuroengineering
Department of Neurology
University Hospital Cologne

Documentation preface, added 23/05/09 by Victor Vera Frazao:
This document is currently in revision for improvement and fixing.
Specifically changes are made to allow compatibility of the pipeline with Ubuntu 18.04 systems 
and Ubuntu 18.04 Docker base images, respectively, as well as adapting to appearent changes of 
DSI-Studio that were applied since the AIDAmri v.1.1 release. As to date the DSI-Studio version 
used is the 2022/08/03 Ubuntu 18.04 release.
All changes and additional documentations within this script carry a signature with the writer's 
initials (e.g. VVF for Victor Vera Frazao) and the date at application, denoted after '//' at 
the end of the comment line. If code segments need clearance the comment line will be prefaced 
by '#?'. Changes are prefaced by '#>' and other comments are prefaced ordinalrily 
by '#'.
"""
from __future__ import print_function

import argparse
import os
import glob
import dsi_tools_20170214
import shutil

if __name__ == '__main__':
    # default dsi studio directory
    f = open(os.path.join(os.getcwd(), "dsi_studioPath.txt"), "r")
    dsi_studio = f.read().split("\n")[0]
    f.close()

    #dsi_studio = os.path.abspath(os.path.join(os.getcwd(), os.pardir,os.pardir,os.pardir))+'/lib/dsi_studio'
    # default b-table in input directory
    b_table = os.path.abspath(os.path.join(os.getcwd(), os.pardir,os.pardir))+'/lib/DTI_Jones30.txt'
    # default connectivity directory relative to input directory
    dir_con = r'connectivity'

    parser = argparse.ArgumentParser(description='Get connectivity of DTI dataset')
    requiredNamed = parser.add_argument_group('Required named arguments')
    requiredNamed.add_argument('-i', '--file_in', help='path to the raw NIfTI DTI file (ends with *1.nii.gz)', required=True)

    parser.add_argument('-b', '--b_table', default=b_table, help='b-table in input directory: %s' % (b_table,))

    parser.add_argument('-o', '--optional', nargs = '*', help = 'Optional arguments (e.g. deprecated terminology)') #> New flag for depr. keywords // VVF 23/05/26
    args = parser.parse_args()


    file_cur = os.path.dirname(args.file_in)
    dsi_path = os.path.join(file_cur, 'DSI_studio')
    mcf_path = os.path.join(file_cur, 'mcf_Folder')
    dir_mask = glob.glob(os.path.join(dsi_path, '*BetMask_scaled.nii'))[0]
    dir_out = args.file_in

    if os.path.exists(mcf_path):
        shutil.rmtree(mcf_path)
    os.mkdir(mcf_path)
    file_in = dsi_tools_20170214.fsl_SeparateSliceMoCo(args.file_in, mcf_path)
    dsi_tools_20170214.srcgen(dsi_studio, file_in, dir_mask, dir_out, args.b_table)
    file_in = os.path.join(file_cur,'fib_map')

    dir_out = os.path.dirname(args.file_in)
    dsi_tools_20170214.tracking(dsi_studio, file_in)
    dir_seeds = glob.glob(os.path.join(file_cur, 'DSI_studio', '*StrokeMask_scaled.nii')) #> changed 'nii.gz' to '.nii' to correct ROI label transfer // VVF 23/05/10
    if len(dir_seeds)>0:
        dir_seeds = glob.glob(os.path.join(file_cur, 'DSI_studio', '*StrokeMask_scaled.nii'))[0] #> changed 'nii.gz' to '.nii' to correct ROI label transfer // VVF 23/05/10
        dsi_tools_20170214.connectivity(dsi_studio, file_in, dir_seeds, dir_out, dir_con)

        dir_seeds = glob.glob(os.path.join(file_cur, 'DSI_studio', '*rsfMRI_Mask_scaled.nii'))[0] #> changed 'nii.gz' to '.nii' to correct ROI label transfer // VVF 23/05/10
        dsi_tools_20170214.connectivity(dsi_studio, file_in, dir_seeds, dir_out, dir_con)



    dir_seeds = glob.glob(os.path.join(file_cur, 'DSI_studio', '*Anno_scaled.nii'))[0] #> changed 'nii.gz' to '.nii' to correct ROI label transfer // VVF 23/05/10
    dsi_tools_20170214.connectivity(dsi_studio, file_in, dir_seeds, dir_out, dir_con)

    dir_seeds = glob.glob(os.path.join(file_cur, 'DSI_studio', '*Anno_rsfMRISplit_scaled.nii'))[0] #> changed 'nii.gz' to '.nii' to correct ROI label transfer // VVF 23/05/10
    dsi_tools_20170214.connectivity(dsi_studio, file_in, dir_seeds, dir_out, dir_con)

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
    print('DTI Connectivity  \033[0;30;42m COMPLETED \33[0m')
