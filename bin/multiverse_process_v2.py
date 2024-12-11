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
        # Skip if the folder has already been processed
        if is_already_processed(folder_path):
            print(f"Skipping already processed folder: {folder_path}")
            return

        print(f"Processing folder: {folder_path}")

        # Simulate processing tasks (replace with actual logic)
        print(f"Executing task logic for {folder_path}...")

        # Mark the folder as processed
        mark_as_processed(folder_path)
    except Exception as e:
        print(f"Error during execution of task for {folder_path}: {e}")

def get_subfolders(input_path):
    """Generate a list of all folders starting with 'sub-' in the given input path."""
    return sorted([
        os.path.join(input_path, f)
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
