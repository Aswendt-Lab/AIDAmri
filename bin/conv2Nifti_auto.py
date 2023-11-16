"""
Created on 18/10/2023

@author: Marc Schneider
AG Neuroimaging and Neuroengineering of Experimental Stroke
Department of Neurology, University Hospital Cologne

This script automates the conversion from the raw bruker data format to the NIfTI
format for the whole dataset using brkraw. The raw
data needs to be stored in one folder.
All the data which is contained in the input folder will be converted to nifti. During the processing a new folder called proc_data is being
created in the same directory where the raw data folder is located. If you wish to save the output elsewhere you can specify the output directory with the -o flag when starting the script.

Example:
python conv2Nifti_auto.py -i /Volumes/Desktop/MRI/raw_data -o /Volumes/Desktop/MRI/raw_data/proc_data
"""

import os
import csv
import json
import pandas as pd
import nibabel as nii
import glob as glob
from pathlib import Path
import numpy as np
import re
from PV2NIfTiConverter import P2_IDLt2_mapping


def create_slice_timings(method_file, scanid, out_file):
    # read in method file to search for parameters
    with open(method_file, "r") as infile:
        lines = infile.readlines()
        interleaved = False
        repetition_time = None
        slicepack_delay = None
        slice_order = []
        n_slices = 0
        reverse = False
        
        # iterate over line to find parameters
        for idx, line in enumerate(lines):
            if "RepetitionTime=" in line:
                repetition_time = int(float(line.split("=")[1]))
                repetition_time = int(repetition_time)
            if "PackDel=" in line:
                slicepack_delay = int(float(line.split("=")[1]))
            if "ObjOrderScheme=" in line:
                slice_order = line.split("=")[1]
            if slice_order == 'Sequential':
                interleaved = False
            else:
                interleaved = True
            if "ObjOrderList=" in line:    
                n_slices = re.findall(r'\d+', line)
                if len(n_slices) == 1:
                    n_slices = int(n_slices[0])
                if lines[idx+1]:
                    slice_order = [int(float(s)) for s in re.findall(r'\d+', lines[idx+1])]
                    if slice_order[0] > slice_order[-1]:
                        reverse = True

        # calculate actual slice timings
        slice_timings = calculate_slice_timings(n_slices, repetition_time, slicepack_delay, slice_order, reverse)

        # adjust slice order to start at 1
        slice_order = [x+1 for x in slice_order]
           
        #save metadata
        mri_meta_data = {}
        mri_meta_data["RepetitionTime"] = repetition_time
        mri_meta_data["ObjOrderList"] = slice_order
        mri_meta_data["n_slices"] = n_slices
        mri_meta_data["costum_timings"] = slice_timings
        mri_meta_data["ScanID"] = scanid
        
        if os.path.exists(out_file):
            with open(out_file, "r") as outfile:
                content = json.load(outfile)
                #update brkraw content with own slice timings
                content.update(mri_meta_data)
                with open(out_file, "w") as outfile:
                    json.dump(content, outfile)

        # if json has different naming than usual adjust path
        else:
            parent_path = Path(out_file).parent

            search_path = os.path.join(parent_path, "*.json")
            json_files = glob.glob(search_path)
            
            for json_file in json_files:
                if os.path.exists(json_file):
                    with open(json_file, "r") as outfile:
                        content = json.load(outfile)
                        #update brkraw content with own slice timings
                        content.update(mri_meta_data)
                        with open(json_file, "w") as outfile:
                            json.dump(content, outfile)
                 

def calculate_slice_timings(n_slices, repetition_time, slicepack_delay, slice_order, reverse=False):
    n_slices_2 = int(n_slices / 2)
    slice_spacing = float(repetition_time - slicepack_delay) / float(n_slices * repetition_time)
    if n_slices % 2 == 1: # odd
        slice_timings = list(range(n_slices_2, -n_slices_2 - 1, -1))
        slice_timings = list(map(float, slice_timings))
    else: # even
        slice_timings = list(range(n_slices_2, -n_slices_2, -1))
        slice_timings = list(map(lambda x: float(x) - 0.5, slice_timings))

    if reverse:
        slice_order.reverse()
    
    slice_timings = list(slice_timings[x] for x in slice_order)

    return list((slice_spacing * x) for x in slice_timings)

def get_visu_pars(path):
    echotimes = []
    if os.path.exists(path):
        with open(path, 'r') as infile:
            lines = infile.readlines()
            for idx, line in enumerate(lines):
                if "VisuAcqEchoTime=" in line:    
                    if lines[idx+1]:
                        echotimes = [float(s) for s in re.findall(r'\d+', lines[idx+1])]
                        echotimes = np.array(echotimes)
    return echotimes

def bids_convert(input_dir, out_path):
    ## rearrange proc data in BIDS-format   
    print('brkraw bids_helper ' + input_dir + ' ' + "dataset" + ' -j')
    os.system('brkraw bids_helper ' + input_dir + ' ' + "dataset" + ' -j')
    
    # # adjust dataset.json template
    dataset_json = glob.glob(os.path.join(os.getcwd(),"data*.json"))[0]
    dataset_csv = glob.glob(os.path.join(os.getcwd(),"data*.csv"))[0]
    if os.path.exists(dataset_json):
        with open(dataset_json, 'r') as infile:
            meta_data = json.load(infile)
            meta_data["common"]["RepititionTime"] = ""
            if meta_data["common"]["EchoTime"]:
                del meta_data["common"]["EchoTime"]
                
            with open(dataset_json, 'w') as outfile:
                json.dump(meta_data, outfile)
          
    ## convert to bids
    os.system('brkraw bids_convert ' + input_dir + ' ' + dataset_csv + ' -j ' + dataset_json + ' -o ' + out_path) 

def nifti_convert(input_dir, raw_data_list):
    # create list with full paths of raw data
    list_of_paths = []        
    for idx, raw_path in enumerate(raw_data_list):
        full_path = os.path.join(input_dir, raw_path)
        list_of_paths.append(full_path)
        
    # ## perform bruker to nifti conversion for all files    
    for idx, sub in enumerate(list_of_paths):
        os.system('cd ' + sub)
        os.system('brkraw tonii ' + sub + ' -o ' + sub)

def create_mems_and_map(mese_scan_ids, mese_scan_data):
    # iterate over every subject and ses to check if MEMS files are included
    for idx, sub in enumerate(mese_scan_ids):
        mese_scan_path = os.path.join(output_dir, "sub-" + sub)
        sessions = os.listdir(mese_scan_path)
        for ses in sessions:
            anat_data_path = os.path.join(mese_scan_path, ses, "anat", "*MESE.nii*")
            mese_data_paths = glob.glob(anat_data_path, recursive=True)
            
            #skip the subject if no MEMS files are found
            if not mese_data_paths:
                continue
            
            # collect data of all individual MEMS files of one subject and session
            img_array_data = {}
            for m_d_p in mese_data_paths:
                # find slice numer in path. e.g.: *echo-10_MESE.nii.gz, extract number 10
                slice_number = int(((Path(m_d_p).name).split('-')[-1]).split('_')[0])
            
                # load nifti image and save the array in a dict while key is the slice number
                data = nii.load(m_d_p)
                img_array = data.dataobj.get_unscaled()
                img_array_data[slice_number] = img_array

                # remove single mese file
                os.remove(m_d_p)
                os.remove(m_d_p.replace(".nii.gz", ".json"))
            
            # sort imgs into right order 
            sorted_imgs = []
            for key in sorted(img_array_data):
                sorted_imgs.append(img_array_data[key])
              
            # stack all map related niftis
            new_img = np.stack(sorted_imgs, axis=2)
            qform = data.header.get_qform()
            sform = data.header.get_sform()
            data.header.set_qform(None)
            data.header.set_sform(None)
            nii_img = nii.Nifti1Image(new_img, None, data.header)
            
            # save nifti file in anat folder
            img_name = "sub-" + sub + "_" + ses + "_T2w_MEMS.nii.gz"
            t2_mems_path = os.path.join(output_dir, "sub-" + sub, ses, "anat", img_name)
            nii.save(nii_img, t2_mems_path)

            # create t2 mapping
            # check visu_pars and echo time
            visu_pars_path = os.path.join(pathToRawData, mese_scan_data[sub]["RawData"], str(mese_scan_data[sub]["ScanID"]), "visu_pars")
            
            echotimes = get_visu_pars(visu_pars_path)
        
            if len(echotimes) > 3:
                img_name = "sub-" + sub + "_" + ses + "_T2w_MAP.nii.gz"
                t2map_path = os.path.join(output_dir,  "sub-" + sub, ses, "t2map", img_name)

                if not os.path.exists(os.path.join(output_dir,  "sub-" + sub, ses, "t2map")):
                    os.mkdri(os.path.join(output_dir,  "sub-" + sub, ses, "t2map"))

                P2_IDLt2_mapping.getT2mapping(t2_mems_path, 'T2_2p', 100, 1.5, 'Brummer', echotimes, t2map_path)

                correct_orientation(qform,sform,t2_mems_path,t2map_path)

            # generate transposed MEMS img for later registration
            org_mems_scan = nii.load(t2_mems_path)
            mems_data = org_mems_scan.dataobj.get_unscaled()
            
            mems_data_transposed = np.transpose(mems_data, axes=(0,1,3,2))
            mems_data_first_slice = mems_data_transposed[:,:,:,1]
            
            for i in range(mems_data_transposed.shape[3]):
                mems_data_transposed[:,:,:,i] = mems_data_first_slice
                
            transposed_copied_img = nii.Nifti1Image(mems_data_transposed, org_mems_scan.affine)
            
            img_name = "sub-" + sub + "_" + ses + "_T2w_transposed_MEMS.nii.gz"
            t2_mems_transposed_path = os.path.join(output_dir, "sub-" + sub, ses, "t2map", img_name)
            
            if not os.path.exists(os.path.join(output_dir, "sub-" + sub, ses, "t2map")):
                os.mkdir(os.path.join(output_dir, "sub-" + sub, ses, "t2map"))
            nii.save(transposed_copied_img, t2_mems_transposed_path)


def correct_orientation(qform,sform, t2_mems_img, t2_map_img):
    # overwrite img with correct orienation
    mems_img = nii.load(t2_mems_img)
    imgTemp = mems_img.dataobj.get_unscaled()

    mems_img.header.set_qform(qform)
    mems_img.header.set_sform(sform)

    new_img = nii.Nifti1Image(imgTemp, None, mems_img.header)
    nii.save(new_img, t2_mems_img)

    # overwrite img with correct orienation
    map_img = nii.load(t2_map_img)
    imgTemp = map_img.dataobj.get_unscaled()

    map_img.header.set_qform(qform)
    map_img.header.set_sform(sform)

    new_img = nii.Nifti1Image(imgTemp, None, map_img.header)
    nii.save(new_img, t2_map_img)



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='This script automates the conversion from the raw bruker data format to the NIfTI format using 1_PV2NIfTiConverter/pv_conv2Nifti.py. The raw data needs to be in the following structure: projectfolder/days/subjects/data/. For this script to work, the groupMapping.csv needs to be adjusted, where the group name of every subject''s folder in the raw data structure needs to be specified. This script computes the converison either for all data in the raw project folder or for certain days and/or groups specified through the optional arguments -d and -g. During the processing a new folder called proc_data is being created in the same directory where the raw data folder is located. Example: python conv2Nifti_auto.py -f /Volumes/Desktop/MRI/raw_data -d Baseline P1 P7 P14 P28')
    parser.add_argument('-i', '--input', required=True,
                        help='Path to the parent project folder of the dataset, e.g. raw_data', type=str)
    parser.add_argument('-o', '--output', type=str, required=False, help='Output directory where the results will be saved.')

    ## read out parameters
    args = parser.parse_args()
    pathToRawData = args.input
    if args.output == None:
        output_dir = os.path.join(pathToRawData, "proc_data")
    else:
        output_dir = args.output
    
    # get list of raw data in input folder
    list_of_raw = sorted([d for d in os.listdir(pathToRawData) if os.path.isdir(os.path.join(pathToRawData, d)) \
                               or (os.path.isfile(os.path.join(pathToRawData, d)) and (('zip' in d) or ('PvDataset' in d)))])

    # convert data into nifti format
    nifti_convert(pathToRawData, list_of_raw)
    
    # convert data into BIDS format
    bids_convert(pathToRawData, output_dir)
  
    # delete duplicated files in input folder
    all_files_input_folder = os.listdir(pathToRawData)
    del_file_ext = [".nii", ".bval", ".bvec"]
    
    for file in all_files_input_folder:
        if file not in list_of_raw and file != output_dir and any(ext in file for ext in del_file_ext):
            os.remove(os.path.join(pathToRawData,file))
    
    # find MEMS and fmri files 
    mese_scan_data = {}
    mese_scan_ids = []
    fmri_scan_ids = {}
    dataset_csv = glob.glob(os.path.join(os.getcwd(), "data*.csv"))[0]
    if os.path.exists(dataset_csv):
        with open(dataset_csv, 'r') as csvfile:
            df = pd.read_csv(csvfile, delimiter=',')
            for index, row in df.iterrows():
                # save every sub which has MEMS scans
                if row["modality"] == "MESE":
                    mese_scan_ids.append(row["SubjID"])
                    mese_scan_data[row["SubjID"]] = {}
                    mese_scan_data[row["SubjID"]]["ScanID"] = row["ScanID"]
                    mese_scan_data[row["SubjID"]]["RawData"] = row["RawData"]
                # save every sub and scanid wich is fmri scan
                if row["DataType"] == "func":
                    fmri_scan_ids[row["RawData"]] = {}
                    fmri_scan_ids[row["RawData"]]["ScanID"] = row["ScanID"] 
                    fmri_scan_ids[row["RawData"]]["SessID"] = row["SessID"]
                    fmri_scan_ids[row["RawData"]]["SubjID"] = row["SubjID"]
    
    create_mems_and_map(mese_scan_ids, mese_scan_data)            
    
    # iterate over all fmri scans to calculate and save costum slice timings
    for sub, data in fmri_scan_ids.items():
        scanid = str(data["ScanID"])
        sessid = str(data["SessID"])
        subjid = str(data["SubjID"])
        
        # determine method file path
        fmri_scan_method_file = os.path.join(pathToRawData, sub, scanid, "method")
        
        # determine output json file path
        out_file = os.path.join(output_dir, "sub-" + subjid, "ses-" + sessid, "func", "sub-" + subjid + "_ses-" + sessid + "_EPI.json")
        
        # calculate slice timings
        create_slice_timings(fmri_scan_method_file, scanid, out_file)
     
    print("\n")
    print("###")
    print("All duplicated files have been deleted")     

    print("\n")
    print("###")
    print("Finished converting raw data into nifti format!")
  


