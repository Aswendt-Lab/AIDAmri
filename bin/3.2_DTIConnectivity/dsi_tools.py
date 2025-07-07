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

import os
import re
import sys
import time
import glob
import nibabel as nii
import numpy as np
import nipype.interfaces.fsl as fsl
import shutil
import subprocess
import pandas as pd

def scaleBy10(input_path, inv):
    data = nii.load(input_path)
    imgTemp = data.get_fdata()
    if inv is False:
        scale = np.eye(4) * 10
        scale[3][3] = 1
        scaledNiiData = nii.Nifti1Image(imgTemp, data.affine * scale)
        # overwrite old nifti
        fslPath = os.path.join(os.path.dirname(input_path), 'fslScaleTemp.nii.gz')
        nii.save(scaledNiiData, fslPath)
        return fslPath
    elif inv is True:
        scale = np.eye(4) / 10
        scale[3][3] = 1
        unscaledNiiData = nii.Nifti1Image(imgTemp, data.affine * scale)
        hdrOut = unscaledNiiData.header
        hdrOut.set_xyzt_units('mm')

        nii.save(unscaledNiiData, input_path)
        return input_path
    else:
        sys.exit("Error: inv - parameter should be a boolean.")


def findSlicesData(path, pre):
    regMR_list = []
    fileALL = glob.iglob(path + '/' + pre + '*.nii.gz', recursive=True)
    for filename in fileALL:
        regMR_list.append(filename)
    regMR_list.sort()
    return regMR_list


def fsl_SeparateSliceMoCo(input_file, par_folder):
    # scale Nifti data by factor 10
    dataName = os.path.basename(input_file).split('.')[0]
    fslPath = scaleBy10(input_file, inv=False)

    aidamri_dir = os.getcwd()
    temp_dir = os.path.join(os.path.dirname(input_file), "temp")
    if not os.path.exists(temp_dir):
        os.mkdir(temp_dir)

    os.chdir(temp_dir)
    mySplit = fsl.Split(in_file=fslPath, dimension='z', out_base_name=dataName)
    mySplit.run()
    os.remove(fslPath)

    # sparate ref and src volume in slices
    sliceFiles = findSlicesData(os.getcwd(), dataName)

    # start to correct motions slice by slice
    for i in range(len(sliceFiles)):
        slc = sliceFiles[i]
        output_file = os.path.join(par_folder, os.path.basename(slc))
        myMCFLIRT = fsl.preprocess.MCFLIRT(in_file=slc, out_file=output_file, save_plots=True, terminal_output='none')
        myMCFLIRT.run()
        os.remove(slc)

    # merge slices to a single volume
    mcf_sliceFiles = findSlicesData(par_folder, dataName)
    output_file = os.path.join(os.path.dirname(input_file),
                               os.path.basename(input_file).split('.')[0]) + '_mcf.nii.gz'
    myMerge = fsl.Merge(in_files=mcf_sliceFiles, dimension='z', merged_file=output_file)
    myMerge.run()

    for slc in mcf_sliceFiles: 
        os.remove(slc)

    # unscale result data by factor 10**(-1)
    output_file = scaleBy10(output_file, inv=True)
    
    os.chdir(aidamri_dir)

    return output_file


def make_dir(dir_out, dir_sub):
    """
    Creates new directory.
    """
    dir_out = os.path.normpath(os.path.join(dir_out, dir_sub))
    if not os.path.exists(dir_out):
        os.mkdir(dir_out)
        time.sleep(1.0)
        if not os.path.exists(dir_out):
            sys.exit("Could not create directory \"%s\"" % (dir_out,))
    return dir_out

def move_files(dir_in, dir_out, pattern):
    time.sleep(1.0)
    file_list = glob.glob(dir_in+pattern)
    file_list.sort()

    time.sleep(1.0)
    for file_mv in file_list: # move files from input to output directory
        file_in = os.path.join(dir_in, file_mv)
        shutil.copy(file_in, dir_out)

    for file_mv in file_list: # remove files in output directory
        file_out = os.path.join(dir_out, file_mv)
        if os.path.isfile(file_out):
            os.remove(file_out)

def connectivity(dsi_studio, dir_in, dir_seeds, dir_out, dir_con):
    """
    Calculates connectivity data (types: pass and end).
    """
    if not os.path.exists(dir_in):
        sys.exit("Input directory \"%s\" does not exist." % (dir_in,))

    dir_seeds = os.path.normpath(os.path.join(dir_in, dir_seeds))
    if not os.path.exists(dir_seeds):
        sys.exit("Seeds directory \"%s\" does not exist." % (dir_seeds,))

    if not os.path.exists(dir_out):
        sys.exit("Output directory \"%s\" does not exist." % (dir_out,))

    dir_con = make_dir(dir_out, dir_con)

    # change to input directory
    os.chdir(os.path.dirname(dir_in))
    cmd_ana = r'%s --action=%s --source=%s --tract=%s --connectivity=%s --connectivity_value=%s --connectivity_type=%s'

    filename = glob.glob(dir_in+'/*fib.gz')[0]
    file_trk = glob.glob(dir_in+'/*trk.gz')[0]
    file_seeds = dir_seeds

    # Performs analysis on every connectivity value within the list ('qa' may not be necessary; might be removed in the future.)
    connect_vals = ['qa', 'count']
    for i in connect_vals:
        parameters = (dsi_studio, 'ana', filename, file_trk, file_seeds, i, 'pass,end')
        os.system(cmd_ana % parameters)

    #move_files(dir_in, dir_con, re.escape(filename) + '\.' + re.escape(pre_seeds) + '.*(?:\.pass\.|\.end\.)')
    move_files(os.path.dirname(file_trk), dir_con, '/*.txt')
    move_files(os.path.dirname(file_trk), dir_con, '/*.mat')

def mapsgen(dsi_studio, dir_in, dir_msk, b_table, pattern_in, pattern_fib):
    """
    FUNCTION DEPRECATED. REMOVAL PENDING.
    """
    pre_msk = 'bet.bin.'

    ext_src = '.src.gz'
    ext_nii = '.nii.gz'

    if not os.path.exists(dir_in):
        sys.exit("Input directory \"%s\" does not exist." % (dir_in,))

    dir_msk = os.path.normpath(os.path.join(dir_in, dir_msk))
    if not os.path.exists(dir_msk):
        sys.exit("Masks directory \"%s\" does not exist." % (dir_msk,))

    b_table = os.path.join(dir_in, b_table)
    if not os.path.isfile(b_table):
        sys.exit("File \"%s\" does not exist." % (b_table,))

    # change to input directory
    os.chdir(dir_in)

    cmd_src = r'%s --action=%s --source=%s --output=%s --b_table=%s'
    # method: 0:DSI, 1:DTI, 4:GQI 7:QSDR, param0: 1.25 (in vivo) diffusion sampling lenth ratio for GQI and QSDR reconstruction, --thread_count: number of multi-threads used to conduct reconstruction 
    cmd_rec = r'%s --action=%s --source=%s --mask=%s --method=%d --param0=%s --thread_count=%d --check_btable=%d'

    file_list = [x for x in os.listdir(dir_in) if os.path.isfile(os.path.join(dir_in, x)) and re.match(pattern_in, x)]
    file_list.sort()

    for index, filename in enumerate(file_list):
        # create source files
        pos = filename.rfind('_')

        file_src = filename[:pos] + ext_src
        parameters = (dsi_studio, 'src', filename, file_src, b_table)
        subprocess.call(cmd_src % parameters)

        # create fib files
        file_msk = os.path.join(dir_msk, pre_msk + filename[:pos] + ext_nii)
        parameters = (dsi_studio, 'rec', file_src, file_msk, 3, '1.25', 2, 0)
        subprocess.call(cmd_rec % parameters)

    # extracts maps: 2 ways:
    cmd_exp = r'%s --action=%s --source=%s --export=%s'

    file_list = [x for x in os.listdir(dir_in) if os.path.isfile(os.path.join(dir_in, x)) and re.match(pattern_fib, x)]
    file_list.sort()

    for index, filename in enumerate(file_list):
        #file_fib = os.path.join(dir_in, filename)
        #parameters = (dsi_studio, 'exp', file_fib, 'fa')
        parameters = (dsi_studio, 'exp', filename, 'fa')
        print("%d of %d:" % (index + 1, len(file_list)), cmd_exp % parameters)
        subprocess.call(cmd_exp % parameters)

def srcgen(dsi_studio, dir_in, dir_msk, dir_out, b_table, recon_method='dti', vivo='in_vivo', make_isotropic=0):
    """
    Sources and creates fib files. Diffusivity and aniosotropy metrics are exported from data.
    """
    dir_src = r'src'
    dir_fib = r'fib_map'
    dir_qa  = r'DSI_studio'
    dir_con = r'connectivity'
    ext_src = '.src.gz'

    if not os.path.exists(dir_in):
        sys.exit("Input directory \"%s\" does not exist." % (dir_in,))

    dir_msk = os.path.normpath(os.path.join(dir_in, dir_msk))
    if not os.path.exists(dir_msk):
        sys.exit("Masks directory \"%s\" does not exist." % (dir_msk,))

    if not os.path.exists(dir_out):
        sys.exit("Output directory \"%s\" does not exist." % (dir_out,))

    b_table = os.path.join(dir_in, b_table)
    if not os.path.isfile(b_table):
        sys.exit("File \"%s\" does not exist." % (b_table,))

    dir_src = make_dir(os.path.dirname(dir_out), dir_src)
    dir_fib = make_dir(os.path.dirname(dir_out), dir_fib)
    dir_qa  = make_dir(os.path.dirname(dir_out), dir_qa)

    # change to input directory
    os.chdir(os.path.dirname(dir_in))

    cmd_src = r'%s --action=%s --source=%s --output=%s --b_table=%s'
    # method: 0:DSI, 1:DTI, 4:GQI 7:QSDR, param0: 1.25 (in vivo) diffusion sampling lenth ratio for GQI and QSDR reconstruction, 
    # check_btable: Set –check_btable=1 to test b-table orientation and apply automatic flippin, thread_count: number of multi-threads used to conduct reconstruction
    # flip image orientation in x, y or z direction !! needs to be adjusted according to your data, check fiber tracking result to be anatomically meaningful
    cmd_rec = r'%s --action=%s --source=%s --mask=%s --method=%d --param0=%s --other_output=all --check_btable=%d --half_sphere=%d --cmd=%s'

    # create source files
    filename = os.path.basename(dir_in)
    pos = filename.rfind('.')
    file_src = os.path.join(dir_src, filename[:pos] + ext_src)
    parameters = (dsi_studio, 'src', filename, file_src, b_table)
    os.system(cmd_src % parameters)

    # create fib files
    file_msk = dir_msk

    param_zero='1.25'
    # Select reconstruction parameters
    if vivo == "ex_vivo":
        param_zero='0.60'
        print(f'Using param0 value {param_zero} recommended for ex vivo data')
    elif vivo == "in_vivo":
        print(f'Using param0 value {param_zero} recommended for in vivo data')

    # get voxel size from nifti header to resample if make_isotropic == 'auto'
    img = nii.load(filename)
    header = img.header
    mm_voxel_size_arr = header.get_zooms()
    spatial_dims = mm_voxel_size_arr[:3]
    min_vox_size_mm = str(min(spatial_dims))
    if make_isotropic == "auto":
        make_isotropic = min_vox_size_mm
    
    additional_cmd='"[Step T2][B-table][flip by]+[Step T2][B-table][flip bz]"'
    if make_isotropic != 0:
        additional_cmd=f'"[Step T2][Resample]={make_isotropic}+[Step T2][B-table][flip by]+[Step T2][B-table][flip bz]"'
        print(f'Resampling to {make_isotropic} mm isotropic voxel size')

    # default method value for DTI 
    method_rec=1
    if recon_method == "gqi":
        method_rec=4

    parameters = (dsi_studio, 'rec', file_src, file_msk, method_rec, param_zero, 0, 1, additional_cmd)
    os.system(cmd_rec % parameters)

    # move fib to corresponding folders
    move_files(dir_src, dir_fib, '/*fib.gz')

    # extracts maps: 2 ways:
    cmd_exp = r'%s --action=%s --source=%s --export=%s'
    file_fib = glob.glob(dir_fib+'/*fib.gz')[0]
    parameters = (dsi_studio, 'exp', file_fib, 'fa')
    os.system(cmd_exp % parameters)

    # extracts maps: 2 ways:
    cmd_exp = r'%s --action=%s --source=%s --export=%s'
    file_fib = glob.glob(dir_fib + '/*fib.gz')[0]
    parameters = (dsi_studio, 'exp', file_fib, 'md')
    os.system(cmd_exp % parameters)

    # extracts maps: 2 ways:
    cmd_exp = r'%s --action=%s --source=%s --export=%s'
    file_fib = glob.glob(dir_fib + '/*fib.gz')[0]
    parameters = (dsi_studio, 'exp', file_fib, 'ad')
    os.system(cmd_exp % parameters)

    # extracts maps: 2 ways:
    cmd_exp = r'%s --action=%s --source=%s --export=%s'
    file_fib = glob.glob(dir_fib + '/*fib.gz')[0]
    parameters = (dsi_studio, 'exp', file_fib, 'rd')
    os.system(cmd_exp % parameters)

    move_files(dir_fib, dir_qa, '/*qa.nii.gz')
    move_files(dir_fib, dir_qa, '/*fa.nii.gz')
    move_files(dir_fib, dir_qa, '/*md.nii.gz')
    move_files(dir_fib, dir_qa, '/*ad.nii.gz')
    move_files(dir_fib, dir_qa, '/*rd.nii.gz')
    
    fa_file = nii.load(glob.glob(os.path.join(dir_qa,"*fa.nii*"))[0])
    fa_data = fa_file.get_fdata()
    fa_data_flipped = np.flip(fa_data,0)
    fa_data_flipped = np.flip(fa_data_flipped,1)
    fa_file_flipped = nii.Nifti1Image(fa_data_flipped, fa_file.affine)
    nii.save(fa_file_flipped,os.path.join(dir_qa,"fa_flipped.nii.gz"))
    
    md_file = nii.load(glob.glob(os.path.join(dir_qa,"*md.nii*"))[0])
    md_data = md_file.get_fdata()
    md_data_flipped = np.flip(md_data,0)
    md_data_flipped = np.flip(md_data_flipped,1)
    md_file_flipped = nii.Nifti1Image(md_data_flipped, md_file.affine)
    nii.save(md_file_flipped,os.path.join(dir_qa,"md_flipped.nii.gz"))
    
    ad_file = nii.load(glob.glob(os.path.join(dir_qa,"*ad.nii*"))[0])
    ad_data = ad_file.get_fdata()
    ad_data_flipped = np.flip(ad_data,0)
    ad_data_flipped = np.flip(ad_data_flipped,1)
    ad_file_flipped = nii.Nifti1Image(ad_data_flipped, ad_file.affine)
    nii.save(ad_file_flipped,os.path.join(dir_qa,"ad_flipped.nii.gz"))
    
    rd_file = nii.load(glob.glob(os.path.join(dir_qa,"*rd.nii*"))[0])
    rd_data = rd_file.get_fdata()
    rd_data_flipped = np.flip(rd_data,0)
    rd_data_flipped = np.flip(rd_data_flipped,1)
    rd_file_flipped = nii.Nifti1Image(rd_data_flipped, rd_file.affine)
    nii.save(rd_file_flipped,os.path.join(dir_qa,"rd_flipped.nii.gz"))

    return min_vox_size_mm
    
def tracking(dsi_studio, dir_in, track_param='default', min_voxel_size_mm=0.1):
    """
    Performs seed-based fiber-tracking. 
    Default parameters are used unless a custom parameter is specified.
    """
    if not os.path.exists(dir_in):
        sys.exit("Input directory \"%s\" does not exist." % (dir_in,))

    # Define parameter sets
    param_sets = {
        'default':        ['0AD7A33C9A99193FE8D5123F0AD7233CCDCCCC3D9A99993EbF04240420FdcaCDCC4C3Ec'],
        'aida_optimized': [1000000, 0, '.01', '55', 0, '.02', '.1', '.3', '120.0'],
        'rat':            [1000000, 0, '.01', '60', 0, '.02', '.1', '.3', '20.0'],
        'mouse':          [1000000, 0, '.01', '45', 0, '.02', '.1', '.3', '15.0'],
    }

    if isinstance(track_param, str):
        params = param_sets.get(track_param.lower())
        if params is None:
            sys.exit(f'Unknown track_param set: {track_param}')
    elif isinstance(track_param, list):
        params = track_param
    else:
        sys.exit('track_param must be "default", "aida_optimized", "rat", "mouse", or a list of parameter values for fiber_count, interpolation, step_size, turning_angle, check_ending, fa_threshold, smoothing, min_length, max_length.')

    # change to input directory
    os.chdir(os.path.dirname(dir_in))

    # Set tracking based on track_param:
    if track_param == 'default':
        print('Using DSI Studio default tracking parameters')
        # Use this tracking parameters in the form of parameter_id that you can get directly from the dsi_studio gui console. (this is here now the defualt mode)
        cmd_trk = r'%s --action=%s --source=%s --output=%s --parameter_id=%s'
        # Use this tracking parameters in the form of parameter_id that you can get directly from the dsi_studio gui console. (this is here now the defualt mode)
        parameters = (dsi_studio, 'trk', filename, os.path.join(dir_in, filename+'.trk.gz'), '0AD7A33C9A99193FE8D5123F0AD7233CCDCCCC3D9A99993EbF04240420FdcaCDCC4C3Ec')
    else:
        # Use this tracking parameters if you want to specify each tracking parameter separately.
        cmd_trk = r'%s --action=%s --source=%s --output=%s --fiber_count=%d --interpolation=%d --step_size=%s --turning_angle=%s --check_ending=%d --fa_threshold=%s --smoothing=%s --min_length=%s --max_length=%s'
        #parameters = (dsi_studio, 'trk', filename, os.path.join(dir_in, filename+'.trk.gz'), 1000000, 0, '.5', '55', 0, '.02', '.1', '.5', '12.0') #Our Old parameters
        #parameters = (dsi_studio, 'trk', filename, os.path.join(dir_in, filename+'.trk.gz'), 1000000, 0, '.01', '55', 0, '.02', '.1', '.3', '120.0') #Here are the optimized parameters (fatemeh)
        if track_param != 'aida_optimized':
            # step size = 1/2 (voxel size)
            params[2] = min_voxel_size_mm / 2
            # min streamline length = 2 * (voxel_size)
            params[7] = min_voxel_size_mm * 2
        parameters = (dsi_studio, 'trk', filename, os.path.join(dir_in, filename+'.trk.gz'), params)
    
    filename = glob.glob(dir_in+'/*fib.gz')[0]
    
    os.system(cmd_trk % parameters)

def merge_bval_bvec_to_btable(folder_path):
    # List files in the specified folder
    files = os.listdir(folder_path)

    # Find bval and bvec files in the folder
    bval_file = None
    bvec_file = None

    for file in files:
        if file.endswith(".bval"):
            bval_file = os.path.join(folder_path, file)
        elif file.endswith(".bvec"):
            bvec_file = os.path.join(folder_path, file)

    # Check if both bval and bvec files were found
    if bval_file is not None or bvec_file is not None:
        print("Both bval and bvec files must be present in the folder.")
        fileName = os.path.basename(bvec_file).replace(".bvec","")
       
    try:
        with open(bval_file, 'r') as bval_file:
            bval_contents = bval_file.read()
            # Split the content into a list of values (assuming it's space-separated)
            bval_values = bval_contents.strip().split()
            # Convert the list to a Pandas DataFrame and cast the 'bval' column to integers
            bval_table = pd.DataFrame({'bval': bval_values}).astype(float)

        with open(bvec_file, 'r') as bvec_file:
            # Read lines and split each line into values
            bvec_lines = bvec_file.readlines()
            bvec_values = [line.strip().split() for line in bvec_lines]

            # Create a Pandas DataFrame from the values
            bvec_table = pd.DataFrame(bvec_values, columns=[f'bvec_{i+1}' for i in range(len(bvec_values[0]))])
            # Transpose the bvec_table
            bvec_table = bvec_table.T

        # Merge bval_table and bvec_table
        merged_table = np.hstack((bval_table, bvec_table))
        # Convert the merged_table content to float
        merged_table = merged_table.astype(float)
        # Define the path for the final merged table
        final_path = os.path.join(folder_path, fileName + "_btable.txt")

        # Save the merged table to the final file
        np.savetxt(final_path, merged_table, fmt='%f', delimiter='\t')
        print(f"Merged table saved to {final_path}")
        return final_path
    except FileNotFoundError:
        print("One or both of the bval and bvec files were not found.")
        return False
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        return False

if __name__ == '__main__':
    pass
