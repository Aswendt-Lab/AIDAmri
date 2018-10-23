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
import dsi_tools_20170214
import shutil

if __name__ == '__main__':
    # default dsi studio directory
    f = open(os.path.join(os.getcwd(), "dsi_studioPath.txt"), "r")
    dsi_studio = f.read()
    f.close()

    #dsi_studio = os.path.abspath(os.path.join(os.getcwd(), os.pardir,os.pardir,os.pardir))+'/lib/dsi_studio'
    # default b-table in input directory
    b_table = os.path.abspath(os.path.join(os.getcwd(), os.pardir,os.pardir))+'/lib/DTI_Jones30.txt'
    # default connectivity directory relative to input directory
    dir_con = r'connectivity'



    parser = argparse.ArgumentParser(description='dsi_srcgen')
    parser.add_argument('-i', '--dir_in', help='input directory')
    parser.add_argument('-b', '--b_table', default=b_table, help='b-table in input directory: %s' % (b_table,))
    args = parser.parse_args()


    dir_cur = os.path.dirname(args.dir_in)
    dsi_path = os.path.join(dir_cur, 'DSI_studio')
    mcf_path = os.path.join(dir_cur, 'mcf_Folder')
    dir_mask = glob.glob(os.path.join(dsi_path, '*BetMask_scaled.nii.gz'))[0]
    dir_out = args.dir_in

    if os.path.exists(mcf_path):
        shutil.rmtree(mcf_path)
    os.mkdir(mcf_path)
    dir_in = dsi_tools_20170214.fsl_SeparateSliceMoCo(args.dir_in, mcf_path)

    dsi_tools_20170214.srcgen(dsi_studio, dir_in, dir_mask, dir_out, args.b_table)
    dir_in = os.path.join(dir_cur,'fib_map')

    dir_out = os.path.dirname(args.dir_in)
    dsi_tools_20170214.tracking(dsi_studio, dir_in)
    dir_seeds = glob.glob(os.path.join(dir_cur, 'DSI_studio', '*StrokeMask_scaled.nii.gz'))
    if len(dir_seeds)>0:
        dir_seeds = glob.glob(os.path.join(dir_cur, 'DSI_studio', '*StrokeMask_scaled.nii.gz'))[0]
        dsi_tools_20170214.connectivity(dsi_studio, dir_in, dir_seeds, dir_out, dir_con)

        dir_seeds = glob.glob(os.path.join(dir_cur, 'DSI_studio', '*rsfMRI_Mask_scaled.nii.gz'))[0]
        dsi_tools_20170214.connectivity(dsi_studio, dir_in, dir_seeds, dir_out, dir_con)



    dir_seeds = glob.glob(os.path.join(dir_cur, 'DSI_studio', '*Anno_scaled.nii.gz'))[0]
    dsi_tools_20170214.connectivity(dsi_studio, dir_in, dir_seeds, dir_out, dir_con)

    dir_seeds = glob.glob(os.path.join(dir_cur, 'DSI_studio', '*Anno_rsfMRI_scaled.nii.gz'))[0]
    dsi_tools_20170214.connectivity(dsi_studio, dir_in, dir_seeds, dir_out, dir_con)

    # rename files to reduce path length
    confiles = os.path.join(dir_cur,dir_con)
    data_list = os.listdir(confiles)
    for filename in data_list:
        splittedName = filename.split('.src.gz.dti.fib.gz.')
        if len(splittedName)>1:
            newName = splittedName[1]
            newName = os.path.join(confiles,newName)
            if os.path.isfile(newName):
                os.remove(newName)
            oldName = os.path.join(confiles,filename)
            os.rename(oldName,newName)

    print('DTI Connectivity  \033[0;30;42m COMPLETED \33[0m')
