import os
import glob
import subprocess
import shlex
import argparse
import nibabel as nib
import numpy as np
import shutil
from scipy.ndimage import zoom
from tqdm import tqdm

DATASET_PATH = "/mnt/data2/2024_Grandjean_Multiverse"  # Hardcoded dataset path
PROCESSED_LOG = os.path.join(DATASET_PATH, "processed_folders.log")

def ensure_processed_log_exists():
    """
    Ensure that the processed log file exists to prevent issues.
    """
    try:
        if not os.path.exists(PROCESSED_LOG):
            print(f"Creating processed log file at: {PROCESSED_LOG}")
            open(PROCESSED_LOG, "w").close()  # Create an empty log file
        else:
            print(f"Processed log file already exists at: {PROCESSED_LOG}")
    except Exception as e:
        print(f"Error creating processed log file: {e}")

def is_already_processed(folder_path):
    """
    Check if a folder has already been processed by looking at the processed log file.
    """
    ensure_processed_log_exists()
    try:
        with open(PROCESSED_LOG, "r") as log_file:
            processed_folders = log_file.read().splitlines()
        return folder_path in processed_folders
    except Exception as e:
        print(f"Error reading processed log file: {e}")
        return False

def mark_as_processed(folder_path):
    """
    Mark a folder as processed by appending its path to the processed log file.
    """
    ensure_processed_log_exists()
    try:
        with open(PROCESSED_LOG, "a") as log_file:
            log_file.write(folder_path + "\n")
        print(f"Marked folder as processed: {folder_path}")
    except Exception as e:
        print(f"Error marking folder as processed: {e}")

def execute_task(folder_path, script_path):
    """Execute the specified tasks on the given folder."""
    try:
        # Normalize folder path for consistent comparisons
        folder_path = os.path.abspath(folder_path)

        # Skip if the folder has already been processed
        if is_already_processed(folder_path):
            print(f"Skipping already processed folder: {folder_path}")
            return

        print(f"Processing folder: {folder_path}")

        # Step 1: Datalad get
        print(f"Running 'datalad get' for {folder_path}")
        result = subprocess.run(["datalad", "get", folder_path], text=True)
        if result.returncode != 0:
            print(f"Error during 'datalad get' for {folder_path}")
            return

        # Step 2: Datalad unlock
        print(f"Running 'datalad unlock' for {folder_path}")
        result = subprocess.run(["datalad", "unlock", folder_path], text=True)
        if result.returncode != 0:
            print(f"Error during 'datalad unlock' for {folder_path}")
            return

        # Step 3: Execute the Python script
        print(f"Executing script: {script_path} for folder: {folder_path}")
        result = subprocess.run(["python", script_path, "-i", folder_path], text=True)
        if result.returncode != 0:
            print(f"Error processing {folder_path}: {result.stderr}")
            return

        print(f"Task completed successfully for {folder_path}")

        # Step 4: Datalad save, push, and drop for the folder
        print(f"Saving and pushing changes for folder: {folder_path}")
        result = subprocess.run(["datalad", "save", "-m", f"add {os.path.basename(folder_path)} proc", folder_path], text=True)
        if result.returncode != 0:
            print(f"Error during 'datalad save' for {folder_path}: {result.stderr}")
            return

        result = subprocess.run(["datalad", "push", "--to", "origin"], text=True)
        if result.returncode != 0:
            print(f"Error during 'datalad push' for {folder_path}: {result.stderr}")
            return

        result = subprocess.run(["datalad", "drop", folder_path], text=True)
        if result.returncode != 0:
            print(f"Error during 'datalad drop' for {folder_path}: {result.stderr}")
            return

        # Step 5: Datalad save, push, and drop for results
        results_file = f"Multiverse_Results/task1/{os.path.basename(folder_path)}_task-rest_bold_mcf_st_f_registered_on_SIGMA_template_originated.nii"
        print(f"Saving and pushing changes for results file: {results_file}")
        result = subprocess.run(["datalad", "save", "-m", "update results", results_file], text=True)
        if result.returncode != 0:
            print(f"Error during 'datalad save' for results file {results_file}: {result.stderr}")
            return

        result = subprocess.run(["datalad", "push", "--to", "origin"], text=True)
        if result.returncode != 0:
            print(f"Error during 'datalad push' for results file {results_file}: {result.stderr}")
            return

        result = subprocess.run(["datalad", "drop", results_file], text=True)
        if result.returncode != 0:
            print(f"Error during 'datalad drop' for results file {results_file}: {result.stderr}")
            return

        # Mark the folder as processed
        mark_as_processed(folder_path)
    except Exception as e:
        print(f"Error during execution of task for {folder_path}: {e}")

def get_subfolders(input_path):
    """Generate a list of all folders starting with 'sub-' in the given input path."""
    return sorted([
        os.path.abspath(os.path.join(input_path, f))
        for f in os.listdir(input_path)
        if os.path.isdir(os.path.join(input_path, f)) and f.startswith("sub-")
    ])

def main(input_path, script_path):
    subfolders = get_subfolders(input_path)

    # Debugging output to confirm parsing of subfolders
    print(f"Subfolders to process: {subfolders}")

    # Process each folder exactly once
    for folder in subfolders:
        print(f"Processing folder: {folder}")
        execute_task(folder, script_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process folders starting with 'sub-'.")
    parser.add_argument("-i", "--input", required=True, help="Path to the input directory containing 'sub-' folders.")
    parser.add_argument("-s", "--script", required=True, help="Path to the Python script to execute.")

    args = parser.parse_args()

    main(args.input, args.script)
