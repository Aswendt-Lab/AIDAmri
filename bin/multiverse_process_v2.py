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
        print(f"Current log content: {open(PROCESSED_LOG, 'r').readlines()}")
    except Exception as e:
        print(f"Error marking folder as processed: {e}")

def task_results_exist(folder_path, task_folder):
    """
    Check if the results for a specific task already exist in the results folder.
    """
    outputfolder = os.path.join(os.path.dirname(folder_path), "Multiverse_Results", task_folder)
    return os.path.exists(outputfolder) and len(os.listdir(outputfolder)) > 0

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

def copy_files_to_results_folder(input_folder, new_epi_files, motion_parameters_list_of_folders, new_temporal_files):
    def move_with_unique_name(src, dst):
        """
        Move a file to a destination path. If a file with the same name already exists,
        append a unique counter to the filename to avoid overwriting.
        """
        base, ext = os.path.splitext(os.path.basename(src))
        counter = 1
        while os.path.exists(dst):
            dst = os.path.join(os.path.dirname(dst), f"{base}_{counter}{ext}")
            counter += 1
        shutil.move(src, dst)

    outputfolder = os.path.join(os.path.dirname(input_folder), "Multiverse_Results")
    os.makedirs(outputfolder, exist_ok=True)

    # Create task folders
    task_folders = ['task1', 'task2', 'task3']
    for task_folder in task_folders:
        task_path = os.path.join(outputfolder, task_folder)
        os.makedirs(task_path, exist_ok=True)

        if task_folder == "task1":
            # Move files to task1
            for epi_file in new_epi_files:
                destination = os.path.join(task_path, os.path.basename(epi_file))
                print(f"Moving {epi_file} to {destination}")
                move_with_unique_name(epi_file, destination)

        if task_folder == "task3":
            # Move temporal mean EPI files to task3
            for temporal_mean_epi in new_temporal_files:
                destination = os.path.join(task_path, os.path.basename(temporal_mean_epi))
                print(f"Moving {temporal_mean_epi} to {destination}")
                move_with_unique_name(temporal_mean_epi, destination)

        if task_folder == "task2":
            # Copy motion folders to task2
            for motion_folder in motion_parameters_list_of_folders:
                destination = os.path.join(task_path, os.path.basename(motion_folder))
                print(f"Copying motion folder {motion_folder} to {destination}")
                if os.path.exists(destination):
                    shutil.rmtree(destination)  # Clear any existing folder with the same name
                shutil.copytree(motion_folder, destination)

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