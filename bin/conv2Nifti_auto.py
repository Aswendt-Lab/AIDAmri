"""
Created on 18/10/2023

@author: Marc Schneider
AG Neuroimaging and Neuroengineering of Experimental Stroke
Department of Neurology, University Hospital Cologne

This script automates the conversion from the raw bruker data format to the NIfTI
format for the whole dataset using brkraw. The raw
data needs to be stored in one folder.
All the data which is contained in the input folder will be converted to nifti. During the processing a new folder called proc_data is being
created next to the raw data folder. If you wish to save the output elsewhere you can specify the output directory with the -o flag when starting the script.

Example:
python conv2Nifti_auto.py -i /Volumes/Desktop/MRI/raw_data
"""

import os
import sys
import csv
import json
from calendar import month_name
from datetime import datetime
from zoneinfo import ZoneInfo
import pandas as pd
import nibabel as nii
import glob as glob
from pathlib import Path
import numpy as np
import re
import concurrent.futures
from PV2NIfTiConverter import P2_IDLt2_mapping
import functools
import subprocess
import shlex
import logging
import shutil
import openpyxl
import contextlib
import io

from helper_tools.plot_sourcedata_niftis import process_subject, write_html_report

REPORT_TIMEZONE = ZoneInfo("Europe/Berlin")


class BerlinTimeFormatter(logging.Formatter):
    def formatTime(self, record, datefmt=None):
        timestamp = datetime.fromtimestamp(record.created, REPORT_TIMEZONE)
        return (
            f"{timestamp.day:02d} {month_name[timestamp.month]} {timestamp.year} "
            f"{timestamp:%H:%M:%S} {timestamp.tzname()}"
        )

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


def extract_command_lines(output, tokens):
    return [
        line.strip()
        for line in output.splitlines()
        if any(token in line for token in tokens)
    ]


def log_command_output(label, output, returncode, warning_tokens=None, error_tokens=None):
    warning_tokens = warning_tokens or ["Warning", "FutureWarning", "UserWarning"]
    error_tokens = error_tokens or ["Traceback", "Error", "FAILED", "Failed"]
    warning_lines = extract_command_lines(output, warning_tokens)
    error_lines = extract_command_lines(output, error_tokens)
    message = f"{label}:\n{output}"

    if returncode != 0 or error_lines:
        logging.warning(message)
    elif warning_lines:
        logging.warning(message)
    else:
        logging.info(message)

    if returncode != 0:
        status = "failed"
    elif error_lines:
        status = "completed_with_issues"
    elif warning_lines:
        status = "completed_with_warnings"
    else:
        status = "completed"

    return {
        "status": status,
        "returncode": returncode,
        "warning_lines": warning_lines,
        "error_lines": error_lines,
    }


def extract_failed_scan_ids(issue_lines):
    scan_ids = []
    for line in issue_lines:
        match = re.search(r"Conversion failed:\s*ScanID:(\d+)", line)
        if match:
            scan_ids.append(int(match.group(1)))
    return scan_ids


def get_dataset_scan_ids(input_dir):
    dataset_csv_candidates = glob.glob(os.path.join(input_dir, "dataset*.csv"))
    if not dataset_csv_candidates:
        return set()
    df = pd.read_csv(dataset_csv_candidates[0])
    if "ScanID" not in df.columns:
        return set()
    return set(pd.to_numeric(df["ScanID"], errors="coerce").dropna().astype(int))


def classify_nifti_conversion_results(results, relevant_scan_ids):
    real_failures = []
    relevant_issue_count = 0
    ignored_issue_count = 0

    for result in results:
        if result["status"] == "failed":
            real_failures.append(result)
            continue

        failed_scan_ids = extract_failed_scan_ids(result.get("issue_lines", []))
        relevant_failed_ids = [scan_id for scan_id in failed_scan_ids if scan_id in relevant_scan_ids]
        ignored_failed_ids = [scan_id for scan_id in failed_scan_ids if scan_id not in relevant_scan_ids]
        relevant_issue_count += len(relevant_failed_ids)
        ignored_issue_count += len(ignored_failed_ids)

        if relevant_failed_ids:
            logging.warning(
                "NIfTI conversion failed for relevant dataset scan IDs %s in %s",
                relevant_failed_ids,
                result.get("dataset"),
            )
        if ignored_failed_ids:
            logging.info(
                "Ignoring NIfTI conversion failures for non-dataset/support scan IDs %s in %s",
                ignored_failed_ids,
                result.get("dataset"),
            )

    return real_failures, relevant_issue_count, ignored_issue_count

def bids_convert(input_dir, output_dir):
    ## rearrange proc data in BIDS-format    
    temp_dir = os.path.join(input_dir,"temp")   
    command = f"brkraw bids_helper {input_dir} dataset -j"
    command_args = shlex.split(command)
    
    os.chdir(input_dir)
    command_results = []
    
    try:
        result = subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        command_results.append(log_command_output("Output bids helper", result.stdout, result.returncode))
    except Exception as e:
        logging.error(f'Fehler bei der Ausführung des Befehls: {command_args}\nFehlermeldung: {str(e)}')
        raise
    
    # # adjust dataset.json template
    dataset_json = glob.glob(os.path.join(os.getcwd(),"data*.json"))[0]
    dataset_csv = glob.glob(os.path.join(os.getcwd(),"data*.csv"))[0]
    if os.path.exists(dataset_json):
        with open(dataset_json, 'r') as infile:
            meta_data = json.load(infile)
            if meta_data["common"]["EchoTime"]:
                del meta_data["common"]["EchoTime"]
                
            with open(dataset_json, 'w') as outfile:
                json.dump(meta_data, outfile)
          
    ## convert to bids
    command = f"brkraw bids_convert {input_dir} {dataset_csv} -j {dataset_json} -o {output_dir}"
    command_args = shlex.split(command)
    try:
        result = subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        command_results.append(log_command_output("Output bids convert", result.stdout, result.returncode))
    except Exception as e:
        logging.error(f'Fehler bei der Ausführung des Befehls: {command_args}\nFehlermeldung: {str(e)}')
        raise

    shutil.rmtree(temp_dir)
    return command_results


def nifti_convert(input_dir, raw_data_list, output_dir):
    # create list with full paths of raw data
    list_of_paths = []  
    aidamri_dir = os.getcwd()
    temp_dir = os.path.join(input_dir,"temp")
    if not os.path.exists(temp_dir):
        os.mkdir(temp_dir)
    os.chdir(temp_dir)
    conversion_results = []
    try:
        with concurrent.futures.ProcessPoolExecutor() as executor:
            futures = {executor.submit(brkraw_tonii, path): path for path in raw_data_list}
            for future in concurrent.futures.as_completed(futures):
                input_path = futures[future]
                try:
                    conversion_results.append(future.result())
                except Exception as e:
                    conversion_results.append({
                        "dataset": os.path.basename(input_path),
                        "status": "failed",
                        "returncode": None,
                        "issue_lines": [str(e)],
                    })
                    logging.error(f"NIfTI conversion failed for {input_path}: {e}")
    finally:
        os.chdir(aidamri_dir)
    return conversion_results
        
def brkraw_tonii(input_path):
    
    command = f"brkraw tonii {input_path}"
    command_args = shlex.split(command)
    try:
        result = subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        issue_lines = [
            line.strip()
            for line in result.stdout.splitlines()
            if "Conversion failed" in line or "Traceback" in line
        ]
        log_message = f"Output nifti conversion of dataset {os.path.basename(input_path)}:\n{result.stdout}"
        if result.returncode != 0 or issue_lines:
            logging.warning(log_message)
        else:
            logging.info(log_message)
        status = "failed" if result.returncode != 0 else "completed_with_issues" if issue_lines else "completed"
        return {
            "dataset": os.path.basename(input_path),
            "status": status,
            "returncode": result.returncode,
            "issue_lines": issue_lines,
        }
    except Exception as e:
        logging.error(f'Fehler bei der Ausführung des Befehls: {command_args}\nFehlermeldung: {str(e)}')
        raise

def create_mems_and_map(mese_scan_ses, mese_scan_data, output_dir):
    # iterate over every subject and ses to check if MEMS files are included
    
    sub = os.path.basename(os.path.dirname(mese_scan_ses))
    ses = os.path.basename(mese_scan_ses)
    result = {
        "session": mese_scan_ses,
        "subject": sub,
        "status": "skipped",
        "reason": "",
        "output_files": [],
    }

    anat_data_path = os.path.join(mese_scan_ses, "anat", "*MESE.nii*")
    mese_data_paths = glob.glob(anat_data_path, recursive=True)

    #skip the subject if no MEMS files are found
    if not mese_data_paths:
        result["reason"] = "no converted anat/*MESE.nii* files found"
        return result
    
    # collect data of all individual MEMS files of one subject and session
    img_array_data = {}
    for m_d_p in mese_data_paths:
        # find slice numer in path. e.g.: *echo-10_MESE.nii.gz, extract number 10
        slice_number = int(((Path(m_d_p).name).split('-')[-1]).split('_')[0])
    
        # load nifti image and save the array in a dict while key is the slice number
        data = nii.load(m_d_p)
        # Use scaled intensities; converted NIfTIs may store int16 data with
        # scl_slope/scl_inter, and get_unscaled() would discard that scaling.
        img_array = data.get_fdata(dtype=np.float32)
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
    img_name = sub + "_" + ses + "_T2w_MEMS.nii.gz"
    t2_mems_path = os.path.join(output_dir, sub, ses, "anat", img_name)
    nii.save(nii_img, t2_mems_path)

    # create t2 map
    sub_num = sub.split("-")[1]
    visu_pars_path = os.path.join(pathToRawData, mese_scan_data[sub_num]["RawData"], str(mese_scan_data[sub_num]["ScanID"]), "visu_pars")

    # get echotimes of scan
    echotimes = get_visu_pars(visu_pars_path)

    map_created = False
    if len(echotimes) > 3:
        img_name = sub + "_" + ses + "_T2w_MAP.nii.gz"
        t2map_path = os.path.join(output_dir, sub, ses, "t2map", img_name)

        if not os.path.exists(os.path.join(output_dir, sub, ses, "t2map")):
            os.mkdir(os.path.join(output_dir, sub, ses, "t2map"))
        try:
            P2_IDLt2_mapping.getT2mapping(t2_mems_path, 'T2_2p', 100, 1.5, 'Brummer', echotimes, t2map_path)
            logging.info(f"Map created for: {os.path.basename(t2_mems_path)}")
        except Exception as e:
            logging.error(f"Error while computing T2w Map:\n{e}")
            raise

        correct_orientation(qform,sform,t2_mems_path,t2map_path)
        result["output_files"].append(t2map_path)
        map_created = True

    # generate transposed MEMS img for later registration
    org_mems_scan = nii.load(t2_mems_path)
    mems_data = org_mems_scan.get_fdata(dtype=np.float32)
    
    mems_data_transposed = np.transpose(mems_data, axes=(0,1,3,2))
    mems_data_first_slice = mems_data_transposed[:,:,:,1]
    
    for i in range(mems_data_transposed.shape[3]):
        mems_data_transposed[:,:,:,i] = mems_data_first_slice
        
    transposed_copied_img = nii.Nifti1Image(
        mems_data_transposed.astype(np.float32),
        org_mems_scan.affine,
    )
    
    img_name = sub + "_" + ses + "_T2w_transposed_MEMS.nii.gz"
    t2_mems_transposed_path = os.path.join(output_dir, sub, ses, "t2map", img_name)
    
    if not os.path.exists(os.path.join(output_dir, sub, ses, "t2map")):
        os.mkdir(os.path.join(output_dir, sub, ses, "t2map"))
    nii.save(transposed_copied_img, t2_mems_transposed_path)
    result["output_files"].extend([t2_mems_path, t2_mems_transposed_path])
    if map_created:
        result["status"] = "created"
    else:
        result["status"] = "partial"
        result["reason"] = "MEMS image created, but no T2w map was created because fewer than 4 echo times were found"
    return result


def correct_orientation(qform,sform, t2_mems_img, t2_map_img):
    # overwrite img with correct orienation
    mems_img = nii.load(t2_mems_img)
    imgTemp = mems_img.get_fdata(dtype=np.float32)

    mems_header = mems_img.header.copy()
    mems_header.set_data_dtype(np.float32)
    mems_header.set_qform(qform)
    mems_header.set_sform(sform)

    new_img = nii.Nifti1Image(imgTemp, None, mems_header)
    nii.save(new_img, t2_mems_img)

    # overwrite img with correct orienation
    map_img = nii.load(t2_map_img)
    imgTemp = map_img.get_fdata(dtype=np.float32)

    map_header = map_img.header.copy()
    map_header.set_data_dtype(np.float32)
    map_header.set_qform(qform)
    map_header.set_sform(sform)

    new_img = nii.Nifti1Image(imgTemp, None, map_header)
    nii.save(new_img, t2_map_img)


#this is needed to be done for the bids converter to work correctly.
def fileCopy(list_of_data, input_path):
    for ll in list_of_data:
        if os.path.dirname(ll) != input_path:  # Use '!=' for inequality
            # Extract the filename from ll
            filename = os.path.basename(ll)
            # Create the destination path by combining input_path and the filename
            destination_path = os.path.join(input_path, filename)
            
            # Use shutil.copy to copy the file from ll to destination
            try:
                shutil.move(ll, destination_path)
                print(f"File '{filename}' moved to '{input_path}' successfully.")
            except Exception as e:
                print(f"Error moving '{filename}' to '{input_path}': {str(e)}")


def _count_files(pattern):
    return len(glob.glob(pattern, recursive=True))


def collect_output_counts(output_dir):
    counts = {
        "anat/T2w": _count_files(os.path.join(output_dir, "sub-*", "ses-*", "anat", "*_T2w.nii.gz")),
        "dwi": _count_files(os.path.join(output_dir, "sub-*", "ses-*", "dwi", "*_dwi.nii.gz")),
        "func/EPI": _count_files(os.path.join(output_dir, "sub-*", "ses-*", "func", "*_EPI.nii.gz")),
        "fmap/fieldmap": _count_files(os.path.join(output_dir, "sub-*", "ses-*", "fmap", "*_fieldmap.nii.gz")),
        "fmap/magnitude": _count_files(os.path.join(output_dir, "sub-*", "ses-*", "fmap", "*_magnitude.nii.gz")),
        "t2map/MEMS": _count_files(os.path.join(output_dir, "sub-*", "ses-*", "t2map", "*MEMS.nii.gz")),
        "t2map/T2w_MAP": _count_files(os.path.join(output_dir, "sub-*", "ses-*", "t2map", "*T2w_MAP.nii.gz")),
    }
    return counts


def print_output_counts(output_dir, title="Detected output files"):
    counts = collect_output_counts(output_dir)
    print(f"{title}:")
    for label, count in counts.items():
        print(f" - {label}: {count}")
    logging.info(f"{title}: {counts}")
    return counts


def _clean_values(values):
    cleaned = []
    for value in values:
        if pd.isna(value):
            continue
        value = str(value).strip()
        if value:
            cleaned.append(value)
    return sorted(set(cleaned))


def describe_dataset_csv(df):
    datatypes = _clean_values(df["DataType"]) if "DataType" in df.columns else []
    modalities = _clean_values(df["modality"]) if "modality" in df.columns else []
    datatype_text = ", ".join(datatypes) if datatypes else "none"
    modality_text = ", ".join(modalities) if modalities else "none"
    return datatype_text, modality_text


def print_conversion_summary(output_dir, postprocess_status):
    counts = collect_output_counts(output_dir)
    print("\nSummary:")
    print("Data type   Files   Status")
    print(f"anat        {counts['anat/T2w']}       converted")
    print(f"dwi         {counts['dwi']}       {postprocess_status.get('dwi', 'not checked')}")
    print(f"func        {counts['func/EPI']}       {postprocess_status.get('func', 'not checked')}")
    fmap_count = counts["fmap/fieldmap"] + counts["fmap/magnitude"]
    print(f"fmap        {fmap_count}       converted")
    t2map_count = counts["t2map/T2w_MAP"]
    print(f"t2map       {t2map_count}       {postprocess_status.get('t2map', 'not checked')}")


def run_with_logged_output(func, *args, log_prefix=None, **kwargs):
    stdout_buffer = io.StringIO()
    stderr_buffer = io.StringIO()
    with contextlib.redirect_stdout(stdout_buffer), contextlib.redirect_stderr(stderr_buffer):
        result = func(*args, **kwargs)

    stdout = stdout_buffer.getvalue().strip()
    stderr = stderr_buffer.getvalue().strip()
    prefix = f"{log_prefix}\n" if log_prefix else ""
    if stdout:
        logging.info("%s%s", prefix, stdout)
    if stderr:
        logging.warning("%s%s", prefix, stderr)
    return result




if __name__ == "__main__":
    import argparse
    from helper_tools.adjustbvecRep import adjust_bvec_rep

    parser = argparse.ArgumentParser(description='This script automates the conversion from the raw bruker data format to the NIfTI format using 1_PV2NIfTiConverter/pv_conv2Nifti.py. The raw data needs to be in the following structure: projectfolder/days/subjects/data/. For this script to work, the groupMapping.csv needs to be adjusted, where the group name of every subject''s folder in the raw data structure needs to be specified. This script computes the conversion either for all data in the raw project folder or for certain days and/or groups specified through the optional arguments -s. During the processing a new folder called proc_data is created next to the raw data folder unless an output directory is specified with -o. Example: python conv2Nifti_auto.py -i /Volumes/Desktop/MRI/raw_data -s Baseline P1 P7 P14 P28')
    parser.add_argument('-i', '--input', required=True,
                        help='Path to the parent project folder of the dataset, e.g. raw_data, WARNING:  all of the raw subjects have to be in one folder and not to have a subfolder structure. otherwise the conversion to bids wont work.', type=str)                 
    parser.add_argument('-s', '--sessions',
                        help='Select which sessions of your data should be processed, if no days are given all data will be used.', type=str, required=False)
    parser.add_argument('-o', '--output', type=str, required=False, help='Output directory where the results will be saved.')

    ## read out parameters
    args = parser.parse_args()
    pathToRawData = args.input
    if args.output == None:
        output_dir = str(Path(pathToRawData).expanduser().resolve().parent / "proc_data")
    else:
        output_dir = args.output

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # Create a new workbook and select the active sheet
    workbook = openpyxl.Workbook()
    sheet = workbook.active

    # Enter data in Row 1
    sheet['A1'] = "Subject"
    sheet['B1'] = "Group"

    #Create sourcedata folder
    sourcedata_dir = os.path.join(output_dir, "Sourcedata")
    os.makedirs(sourcedata_dir, exist_ok=True)
    # Save the workbook
    workbook.save(os.path.join(sourcedata_dir,"GroupMapping.xlsx"))
     
    # Configurate Logging-Modul
    log_file_path = os.path.join(sourcedata_dir, "conv2nifti_log.txt")
    log_handler = logging.FileHandler(log_file_path, mode='w')
    log_handler.setFormatter(BerlinTimeFormatter('%(asctime)s - %(levelname)s - %(message)s'))
    logging.basicConfig(level=logging.INFO, handlers=[log_handler], force=True)
    
    # get list of raw data in input folder
    #list_of_raw = sorted([d for d in os.listdir(pathToRawData) if os.path.isdir(os.path.join(pathToRawData, d)) \
    #                          or (os.path.isfile(os.path.join(pathToRawData, d)) and (('zip' in d) or ('PvDataset' in d)))])
    #list_of_raw = glob.glob(os.path.join(pathToRawData,"**","subject"),recursive=True)
    #list_of_data = []
    #for path in list_of_raw:
    #    list_of_data.append(os.path.dirname(path))

    #fileCopy(list_of_data,pathToRawData)

    list_of_raw = glob.glob(os.path.join(pathToRawData,"**","subject"),recursive=True)
    list_of_data = []
    for path in list_of_raw:
        list_of_data.append(os.path.dirname(path))

    logging.info(f"Converting following datasets: {list_of_data}")
    print(f"Converting following datasets: {list_of_data}")
    terminal_issues = []

    # convert data into nifti format
    print("Paravision to nifti conversion running \33[5m...\33[0m (wait!)")
    #nifti_convert(output_dir, list_of_data)
    nifti_results = nifti_convert(pathToRawData, list_of_data, output_dir)
    relevant_scan_ids = get_dataset_scan_ids(pathToRawData)
    nifti_failed, relevant_nifti_issue_count, ignored_nifti_issue_count = classify_nifti_conversion_results(
        nifti_results,
        relevant_scan_ids,
    )
    if nifti_failed:
        print(f"\rNifti conversion FAILED/PARTIAL: {len(nifti_results) - len(nifti_failed)} dataset(s) completed, {len(nifti_failed)} failed.                  ")
        terminal_issues.append(f"NIfTI conversion failed for {len(nifti_failed)} dataset(s)")
    elif relevant_nifti_issue_count > 0:
        print(f"\rNifti conversion COMPLETED WITH WARNINGS: {relevant_nifti_issue_count} relevant scan conversion issue(s) detected.                  ")
        terminal_issues.append(f"NIfTI conversion reported {relevant_nifti_issue_count} relevant scan conversion issue(s)")
    elif ignored_nifti_issue_count > 0:
        print(f"\rNifti conversion COMPLETED: relevant BIDS scans converted; ignored {ignored_nifti_issue_count} non-dataset/support scan conversion failure(s).                  ")
    else:
        print("\rNifti conversion \033[0;30;42m COMPLETED \33[0m                  ")
    if nifti_failed or relevant_nifti_issue_count > 0 or ignored_nifti_issue_count > 0:
        logging.warning("NIfTI conversion issue summary: %s", nifti_results)
    
    # convert data into BIDS format
    print("BIDS conversion running \33[5m...\33[0m (wait!)")
    bids_results = bids_convert(pathToRawData, output_dir)
    bids_failed = [r for r in bids_results if r["status"] == "failed"]
    bids_with_issues = [r for r in bids_results if r["status"] == "completed_with_issues"]
    bids_with_warnings = [r for r in bids_results if r["status"] == "completed_with_warnings"]
    if bids_failed:
        print(f"\rBIDS conversion FAILED/PARTIAL: {len(bids_failed)} command(s) failed.                   ")
        terminal_issues.append(f"BIDS conversion failed in {len(bids_failed)} command(s)")
    elif bids_with_issues:
        issue_count = sum(len(r["error_lines"]) for r in bids_with_issues)
        print(f"\rBIDS conversion COMPLETED WITH ISSUES: {issue_count} issue(s) detected.                   ")
        terminal_issues.append(f"BIDS conversion reported {issue_count} issue(s)")
    elif bids_with_warnings:
        warning_count = sum(len(r["warning_lines"]) for r in bids_with_warnings)
        print(f"\rBIDS conversion COMPLETED WITH WARNINGS: {warning_count} warning(s) logged.                   ")
    else:
        print("\rBIDS conversion \033[0;30;42m COMPLETED \33[0m                   ")
    print_output_counts(output_dir)
    postprocess_status = {}
    
    # adjust bvecs and bvals for diffusion data for each subject and session
    print("DWI post-processing: adjusting bvecs and bvals \33[5m...\33[0m (wait!)")
    dwi_processed = 0
    dwi_missing = 0
    dwi_failed = 0
    dwi_issue_details = []
    for subject_session_output_dir in glob.glob(os.path.join(output_dir, "sub-*", "ses-*")):
        # adjust bvecs and bvals for length of acquisition (number of repetitions)
        if os.path.exists(os.path.join(subject_session_output_dir, "dwi")):
            try:
                run_with_logged_output(
                    adjust_bvec_rep,
                    subject_session_output_dir,
                    log_prefix=f"DWI post-processing details for {subject_session_output_dir}:",
                )
                dwi_processed += 1
            except FileNotFoundError as e:
                logging.warning(f"Valid diffusion data files missing in {subject_session_output_dir}/dwi: {e}")
                dwi_issue_details.append(f"missing files in {subject_session_output_dir}/dwi: {e}")
                dwi_missing += 1
                continue
            except Exception as e:
                logging.error(f"Error processing diffusion data in {subject_session_output_dir}/dwi: {e}")
                dwi_issue_details.append(f"failed in {subject_session_output_dir}/dwi: {e}")
                dwi_failed += 1
                continue
        else:
            logging.info(f"No diffusion data found in {subject_session_output_dir}, skipping DWI post-processing.")
            continue
    if dwi_processed == 0 and dwi_missing == 0 and dwi_failed == 0:
        print("\rDWI post-processing SKIPPED: no dwi directories found.          ")
        postprocess_status["dwi"] = "skipped: no dwi data"
    elif dwi_failed > 0:
        print(f"\rDWI post-processing FAILED/PARTIAL: {dwi_processed} processed, {dwi_missing} missing files, {dwi_failed} failed.       ")
        postprocess_status["dwi"] = f"partial: {dwi_processed} processed, {dwi_failed} failed"
        terminal_issues.append(f"DWI post-processing had {dwi_failed} failed and {dwi_missing} skipped file set(s)")
    else:
        print(f"\rDWI post-processing COMPLETED: {dwi_processed} file set(s) adjusted, {dwi_missing} skipped.       ")
        postprocess_status["dwi"] = f"converted + bval/bvec adjusted ({dwi_processed})"
        if dwi_missing > 0:
            terminal_issues.append(f"DWI post-processing skipped {dwi_missing} file set(s) with missing files")
    if dwi_issue_details:
        logging.warning("DWI post-processing issue details:\n%s", "\n".join(dwi_issue_details))

    # plot QC images for nifti files
    print("Plotting Report images for nifti files \33[5m...\33[0m (wait!)")
    report_output_dir = os.path.join(output_dir, "Report", "Convert2Nifti")
    os.makedirs(report_output_dir, exist_ok=True)
    report_entries = []
    for subject_dir in glob.glob(os.path.join(output_dir, "sub-*")):
        subject_id = os.path.basename(subject_dir)
        print(f"Processing subject: {subject_id}")
        report_entries.extend(process_subject(subject_dir, report_output_dir, n_slices=10))
    if report_entries:
        write_html_report(report_entries, report_output_dir)
    else:
        print("No NIfTI files found for reporting.")
        logging.warning("No NIfTI files found for reporting.")

    # find MEMS and fmri files 
    mese_scan_data = {}
    mese_scan_ids = []
    fmri_scan_ids = {}
    dataset_csv_candidates = glob.glob(os.path.join(pathToRawData, "dataset*.csv"))
    if not dataset_csv_candidates:
        sys.exit(f"Error: no dataset*.csv found in {pathToRawData}")
    dataset_csv = dataset_csv_candidates[0]
    if os.path.exists(dataset_csv):
        with open(dataset_csv, 'r') as csvfile:
            df = pd.read_csv(csvfile, delimiter=',')
            for index, row in df.iterrows():
                # save every sub which has MEMS scans
                if row.get("modality") == "MESE":
                    mese_scan_ids.append(row["SubjID"])
                    mese_scan_data[row["SubjID"]] = {}
                    mese_scan_data[row["SubjID"]]["ScanID"] = row["ScanID"]
                    mese_scan_data[row["SubjID"]]["RawData"] = row["RawData"]
                # save every sub and scanid wich is fmri scan
                if row.get("DataType") == "func":
                    fmri_scan_ids[row["RawData"]] = {}
                    fmri_scan_ids[row["RawData"]]["ScanID"] = row["ScanID"] 
                    fmri_scan_ids[row["RawData"]]["SessID"] = row["SessID"]
                    fmri_scan_ids[row["RawData"]]["SubjID"] = row["SubjID"]           
    
    # iterate over all fmri scans to calculate and save custom slice timings
    print("Checking functional MRI metadata for slice timing updates ...")
    func_updated = 0
    func_failed = 0
    func_issue_details = []
    for sub, data in fmri_scan_ids.items():
        scanid = str(data["ScanID"])
        sessid = str(data["SessID"])
        subjid = str(data["SubjID"])
        
        # determine method file path
        fmri_scan_method_file = os.path.join(pathToRawData, sub, scanid, "method")
        
        # determine output json file path
        out_file = os.path.join(output_dir, "sub-" + subjid, "ses-" + sessid, "func", "sub-" + subjid + "_ses-" + sessid + "_EPI.json")
        
        # calculate slice timings
        try:
            create_slice_timings(fmri_scan_method_file, scanid, out_file)
            func_updated += 1
        except Exception as e:
            func_failed += 1
            logging.error(f"Error updating slice timing metadata for sub-{subjid} ses-{sessid}: {e}")
            func_issue_details.append(f"sub-{subjid} ses-{sessid}: {e}")
    if func_updated == 0 and func_failed == 0:
        print("FUNC metadata update SKIPPED: no func entries found in dataset.csv.")
        postprocess_status["func"] = "skipped: no func data"
    elif func_failed > 0:
        print(f"FUNC metadata update FAILED/PARTIAL: {func_updated} updated, {func_failed} failed.")
        postprocess_status["func"] = f"partial: {func_updated} updated, {func_failed} failed"
        terminal_issues.append(f"FUNC metadata update failed for {func_failed} file(s)")
    else:
        print(f"FUNC metadata update COMPLETED: {func_updated} EPI json file(s) updated.")
        postprocess_status["func"] = f"converted + slice timing updated ({func_updated})"
    if func_issue_details:
        logging.warning("FUNC metadata update issue details:\n%s", "\n".join(func_issue_details))
    
    ## use parallel computing for a faster generation of t2maps
    mese_scan_sessions = []
    for id in mese_scan_ids:
        mese_scan_path = os.path.join(output_dir, "sub-" + id)
        sessions = os.listdir(mese_scan_path)
        for ses in sessions:
            mese_scan_ses = os.path.join(mese_scan_path, ses)
            if mese_scan_ses not in mese_scan_sessions:
                mese_scan_sessions.append(os.path.join(mese_scan_path, ses))
   
    logging.info(f"Creating T2w maps for following datasets:\n{mese_scan_ids}")
    print("Checking for MESE/MEMS scans for T2 mapping ...")
    if not mese_scan_ids:
        datatypes, modalities = describe_dataset_csv(df)
        print(f"No MESE/MEMS scans found in {dataset_csv}.")
        print(f"Found DataTypes: {datatypes}. Found modalities: {modalities}.")
        print("T2 mapping SKIPPED: expected dataset.csv rows with modality == \"MESE\".")
        logging.info("T2 mapping skipped: no MESE/MEMS scans found in dataset.csv")
        postprocess_status["t2map"] = "skipped: no MESE/MEMS input"
    elif not mese_scan_sessions:
        print("T2 mapping SKIPPED: MESE entries were found, but no matching BIDS session directories exist.")
        logging.warning("T2 mapping skipped: MESE entries found but no matching session directories")
        postprocess_status["t2map"] = "skipped: no matching session directories"
    else:
        print(f"Found {len(mese_scan_ids)} MESE scan(s) for {len(mese_scan_sessions)} session(s).")
        print("T2 mapping running \33[5m...\33[0m (wait!)")
        t2_results = []
        t2_failures = []
        with concurrent.futures.ProcessPoolExecutor() as executor:
            futures = {
                executor.submit(create_mems_and_map, mese_scan_ses, mese_scan_data, output_dir): mese_scan_ses
                for mese_scan_ses in mese_scan_sessions
            }
            for future in concurrent.futures.as_completed(futures):
                session = futures[future]
                try:
                    t2_results.append(future.result())
                except Exception as e:
                    t2_failures.append((session, str(e)))
                    logging.error(f"T2 mapping failed for {session}: {e}")

        t2_created = [r for r in t2_results if r["status"] == "created"]
        t2_partial = [r for r in t2_results if r["status"] == "partial"]
        t2_skipped = [r for r in t2_results if r["status"] == "skipped"]
        if t2_failures:
            print(f"T2 mapping FAILED/PARTIAL: {len(t2_created)} map(s) created, {len(t2_partial)} partial, {len(t2_skipped)} skipped, {len(t2_failures)} failed.")
            postprocess_status["t2map"] = f"partial: {len(t2_created)} map(s) created, {len(t2_failures)} failed"
            terminal_issues.append(f"T2 mapping failed for {len(t2_failures)} session(s)")
            logging.warning(
                "T2 mapping failure details:\n%s",
                "\n".join([f"{session}: {reason}" for session, reason in t2_failures]),
            )
        elif t2_created:
            print(f"T2 mapping COMPLETED: {len(t2_created)} map(s) created, {len(t2_partial)} partial, {len(t2_skipped)} skipped.")
            postprocess_status["t2map"] = f"converted ({len(t2_created)} map(s) created)"
        else:
            print(f"T2 mapping SKIPPED: 0 maps created, {len(t2_partial)} partial, {len(t2_skipped)} skipped.")
            postprocess_status["t2map"] = "skipped: no T2w maps created"
            if t2_partial:
                terminal_issues.append(f"T2 mapping created partial outputs for {len(t2_partial)} session(s), but no T2w maps")
            logging.warning(
                "T2 mapping skipped/partial details:\n%s",
                "\n".join([f"{result['subject']}: {result['reason']}" for result in t2_partial + t2_skipped]),
            )
        logging.info(f"Finished T2 mapping with results: {t2_results}; failures: {t2_failures}")

    dataset_json = glob.glob(os.path.join(pathToRawData, "dataset*.json"))[0]

    #os.remove(dataset_csv)
    #os.remove(dataset_json)

    print_conversion_summary(output_dir, postprocess_status)

    if terminal_issues:
        print("\nIssues detected during conversion:")
        for issue in terminal_issues:
            print(f" - {issue}")
        print(f"See detailed log: {log_file_path}")
    else:
        print("\nNo conversion errors detected.")

    print("\n")
    print("###")
    print("Finished converting raw data into nifti format!")
    
    print("\n")
    print("###")
    print("For detailed information check logging file!")
    
    print("\n")
    print("###")
    print("Thank you for using AIDAmri!")
  
