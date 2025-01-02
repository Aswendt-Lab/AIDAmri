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

# Set a custom temporary directory
os.environ['TMPDIR'] = '/mnt/DATA/tmp'

if not os.path.exists(os.environ['TMPDIR']):
    os.makedirs(os.environ['TMPDIR'])

def log_and_skip(message):
    """Log a message and skip the current iteration."""
    print(message)
    return None

def flip_file_in_z_dimension(input_file, header_from_smooth_bet):
    try:
        img = nib.load(input_file)
        data = img.get_fdata()
        data_flipped = np.flip(data, axis=2)

        affine = header_from_smooth_bet.get_best_affine()
        header = header_from_smooth_bet.copy()
        header.set_data_shape(data_flipped.shape)

        flipped_img = nib.Nifti1Image(data_flipped, affine, header)
        flipped_file = input_file.replace('.nii.gz', '_flippedinZ.nii.gz')
        nib.save(flipped_img, flipped_file)
        return flipped_file
    except Exception as e:
        log_and_skip(f"Error flipping file in Z dimension: {input_file} - {str(e)}")
        return None

def apply_affine_transformations(files_list, func_folder, anat_folder, sigma_template_address):
    transformed_files = []

    for file in tqdm(files_list, desc="Applying Transformations", unit="file"):
        try:
            output_file = file.replace('.nii.gz', '_registered_on_SIGMA_template.nii.gz')
            if os.path.exists(output_file):
                log_and_skip(f"Skipping file {file} as output already exists: {output_file}")
                transformed_files.append(output_file)
                continue

            bet_file = glob.glob(os.path.join(func_folder, "*SmoothBet.nii.gz"))
            if not bet_file:
                log_and_skip(f"Missing *SmoothBet.nii.gz in {func_folder}")
                continue
            bet_file = bet_file[0]

            func_trafo_files = glob.glob(os.path.join(func_folder, "*transMatrixAff.txt"))
            if not func_trafo_files:
                log_and_skip(f"Missing *transMatrixAff.txt in {func_folder}")
                continue
            func_trafo = func_trafo_files[0]  # Using file from func folder

            anat_trofo_inv_files = glob.glob(os.path.join(anat_folder, "*MatrixInv.txt"))
            if not anat_trofo_inv_files:
                log_and_skip(f"Missing *MatrixInv.txt in {anat_folder}")
                continue
            anat_trofo_inv = anat_trofo_inv_files[0]

            bet_img = nib.load(bet_file)
            bet_header = bet_img.header
            file_flipped = flip_file_in_z_dimension(file, bet_header)

            if file_flipped is None:
                continue

            func_trafo_inv = os.path.join(func_folder, "rs-fMRI_niiData", "func_trafo_inv.txt")
            merged_inverted = os.path.join(func_folder, "rs-fMRI_niiData", "merged_inverted.txt")

            subprocess.run(shlex.split(f"reg_transform -invAff {func_trafo} {func_trafo_inv}"), check=True)
            subprocess.run(shlex.split(f"reg_transform -comp {anat_trofo_inv} {func_trafo_inv} {merged_inverted}"), check=True)

            subprocess.run(shlex.split(
                f"reg_resample -ref {sigma_template_address} -flo {file_flipped} -trans {merged_inverted} -res {output_file} -inter 3"
            ), check=True)

            os.remove(file_flipped)
            transformed_files.append(output_file)
        except Exception as e:
            log_and_skip(f"Error applying transformations to {file} - {str(e)}")

    return transformed_files

def create_results_folder(input_folder):
    results_folder = os.path.join(input_folder, "Multiverse_Results")
    subfolders = ["task1", "task2", "task3"]
    if not os.path.exists(results_folder):
        os.makedirs(results_folder)
    for subfolder in subfolders:
        os.makedirs(os.path.join(results_folder, subfolder), exist_ok=True)
    return results_folder

def main(input_folder):
    parent_dir = os.path.dirname(os.path.abspath(__file__))
    sigma_template_address = os.path.join(parent_dir, "../lib/SIGMA_InVivo_Brain_Template_Masked.nii.gz")

    # Step 1: Iteratively search for *_bold_mcf_f.nii.gz files
    epi_files_list = []
    for root, dirs, files in os.walk(input_folder):
        if root.endswith("func/rs-fMRI_niiData"):
            for file in files:
                if file.endswith("_bold_mcf_f.nii.gz"):
                    epi_files_list.append(os.path.join(root, file))

    print(f"Found {len(epi_files_list)} EPI files to process.")
    if not os.path.exists(sigma_template_address):
        raise FileNotFoundError(f"Sigma template not found: {sigma_template_address}")

    # Create results folder
    results_folder = create_results_folder(input_folder)

    # Step 2: Process the _bold_mcf_f.nii.gz files
    for epi_file in tqdm(epi_files_list, desc="Processing EPI Files", unit="file"):
        func_folder = os.path.dirname(os.path.dirname(epi_file))
        anat_folder = os.path.join(os.path.dirname(func_folder), "anat")
        transformed_files = apply_affine_transformations([epi_file], func_folder, anat_folder, sigma_template_address)
        for i, transformed_file in enumerate(transformed_files):
            if i % 3 == 0:
                task_folder = os.path.join(results_folder, "task1")
            elif i % 3 == 1:
                task_folder = os.path.join(results_folder, "task2")
            else:
                task_folder = os.path.join(results_folder, "task3")
            shutil.move(transformed_file, task_folder)
        print(f"Transformed files moved to respective task folders: {transformed_files}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Apply affine transformations to imaging data.")
    parser.add_argument('-i', '--input_folder', required=True, help="Path to the input folder.")
    args = parser.parse_args()
    main(args.input_folder)
