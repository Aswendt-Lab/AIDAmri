import os
import glob
import subprocess
import shlex
import argparse
import nibabel as nib
import numpy as np
import shutil
from scipy.ndimage import zoom
from tqdm import tqdm  # Import tqdm for the progress bar

def flip_file_in_z_dimension(input_file, header_from_smooth_bet):
    """
    Flip the input file in the Z dimension and save it with '_flippedinZ' suffix.
    Uses header from SmoothBet file if provided.
    """
    # Load the image
    img = nib.load(input_file)
    data = img.get_fdata()
    
    # Flip the data in the Z dimension (axis=2)
    data_flipped = np.flip(data, axis=2)
    
    # Use the header from the SmoothBet file if provided
    affine = header_from_smooth_bet.get_best_affine()
    header = header_from_smooth_bet.copy()

    # Check if the data is 4D
    if len(data.shape) == 4:
        # Data is 4D; adjust the header to reflect the fourth dimension
        header.set_data_shape(data_flipped.shape)
    else:
        # Data is 3D; use the original header
        affine = header_from_smooth_bet.get_best_affine()
        header = header_from_smooth_bet.copy()

    # Create a new NIfTI image with the flipped data
    flipped_img = nib.Nifti1Image(data_flipped, affine, header)
    
    # Define output file path
    flipped_file = input_file.replace('.nii.gz', '_flippedinZ.nii.gz')
    
    # Save the flipped image
    nib.save(flipped_img, flipped_file)
    
    return flipped_file

def reorient_and_save(input_file):
    """
    Reorient the image, downsample it by a factor of 2, and save it with '_originated.nii' suffix.
    After saving, delete the input file.
    """
    # Load the image
    img = nib.load(input_file)
    data = img.get_fdata()

    # Check if the data is 4D or 3D
    if len(data.shape) == 4:
        # Data is 4D; permute and flip along the third axis for each 3D volume
        data_reoriented = np.transpose(data, (0, 2, 1, 3))  # permute axes for 4D
        data_reoriented = np.flip(data_reoriented, axis=2)    # flip along the third axis
    else:
        # Data is 3D; permute and flip the 3D volume
        data_reoriented = np.transpose(data, (0, 2, 1))      # permute axes for 3D
        data_reoriented = np.flip(data_reoriented, axis=2)    # flip along the third axis
    
    # Downsample the reoriented data by a factor of 2
    downsampling_factor = 0.5
    data_downsampled = zoom(data_reoriented, zoom=[downsampling_factor] * len(data_reoriented.shape), order=3)

    # Update the affine matrix to account for the new voxel sizes (adjust by a factor of 2)
    affine = img.affine.copy()
    affine[:3, :3] *= 2  # Adjust voxel scaling (double the voxel size because we halved the data size)

    # Update header with new data shape and zooms (voxel sizes)
    header = img.header.copy()
    header.set_data_shape(data_downsampled.shape)
    zooms = np.array(header.get_zooms())
    header.set_zooms(zooms * 2)  # Since the data is downsampled, voxel size doubles
    
    # Create a new NIfTI image with the reoriented and downsampled data
    reoriented_img = nib.Nifti1Image(data_downsampled, affine, header)
    
    # Define output file path
    reoriented_file = input_file.replace('.nii.gz', '_originated.nii')
    
    # Save the reoriented image
    nib.save(reoriented_img, reoriented_file)

    # Delete the original input file
    os.remove(input_file)
    
    return reoriented_file

def apply_affine_transformations(files_list, func_folder, anat_folder, sigma_template_address):
    transformed_files = []
    
    for file in tqdm(files_list, desc="Applying Transformations", unit="file"):  # Add progress bar for file processing
        # Find necessary transformation files
        bet_file = glob.glob(os.path.join(func_folder, "*SmoothBet.nii.gz"))[0]
        func_trafo = glob.glob(os.path.join(func_folder, "*transMatrixAff.txt"))[0]
        anat_trofo_inv = glob.glob(os.path.join(anat_folder, "*MatrixInv.txt"))[0]

        # Load header from SmoothBet file
        bet_img = nib.load(bet_file)
        bet_header = bet_img.header
        
        # Define file paths for new transformations
        func_trafo_inv = os.path.join(func_folder, "func_trafo_inv.txt")
        merged_inverted = os.path.join(func_folder, "merged_inverted.txt")
        file_dir = os.path.dirname(file)
        
        # Flip the file in the Z dimension
        file_flipped = flip_file_in_z_dimension(file, bet_header)
        
        # Output file paths in the same directory as the original file
        base_filename = os.path.basename(file)
        file_st_f_on_template = os.path.join(file_dir, base_filename.replace('.nii.gz', '_registered_on_SIGMA_template.nii.gz'))

        # Step 1: Invert the affine matrix
        command = f"reg_transform -invAff {func_trafo} {func_trafo_inv}"
        subprocess.run(shlex.split(command), check=True)

        # Step 2: Combine the two affine matrices
        command = f"reg_transform -comp {anat_trofo_inv} {func_trafo_inv} {merged_inverted}"
        subprocess.run(shlex.split(command), check=True)

        # Step 4: Resample the flipped file to the template
        command = f"reg_resample -ref {sigma_template_address} -flo {file_flipped} -trans {merged_inverted} -res {file_st_f_on_template} -inter 0"
        subprocess.run(shlex.split(command), check=True)

        # Reorient the transformed file
        reoriented_file = reorient_and_save(file_st_f_on_template)
        transformed_files.append(reoriented_file)
        os.remove(file_flipped)
    
    return transformed_files

def copy_files_to_results_folder(input_folder, new_epi_files, motion_parameters_list_of_folders, new_temporal_files):
    outputfolder = os.path.join(os.path.dirname(input_folder), "Multiverse_Results")
    os.makedirs(outputfolder, exist_ok=True)

    # Create task folders
    task_folders = ['task1', 'task2', 'task3']
    for task_folder in task_folders:
        task_path = os.path.join(outputfolder, task_folder)
        os.makedirs(task_path, exist_ok=True)
        
        if task_folder == "task1":
            # Copy files to each task folder
            for epi_file in new_epi_files:
                shutil.move(epi_file, task_path)
        if task_folder == "task3":        
            for temporal_mean_epi in new_temporal_files:
                shutil.move(temporal_mean_epi, task_path)
        if task_folder == "task2": 
            for motion_folder in motion_parameters_list_of_folders:
                # Copy the folder itself to task_path, keeping the folder name
                destination = os.path.join(task_path, os.path.basename(motion_folder))
                shutil.copytree(motion_folder, destination)
                
    # Print results
    print("EPI Files:", new_epi_files)
    print("Temporal Mean EPI Files:", new_temporal_files)

def main(input_folder):
    parent_dir = os.path.dirname(os.path.abspath(__file__))

    sigma_template_address = os.path.join(parent_dir, "..", "lib", "SIGMA_InVivo_Brain_Template_Masked.nii.gz")
    epi_processed_address = os.path.join(input_folder, "**", "rs-fMRI_niiData", "*_bold_mcf_st_f.nii.gz")
    motion_parameters_address = os.path.join(input_folder, "**", "rs-fMRI_mcf","*mcf.mat")
    temporal_mean_epi_single_frame_original_contrast = os.path.join(input_folder, "**", "rs-fMRI_niiData", "*boldmean.nii.gz")

    # Find files
    epi_files_list = glob.glob(epi_processed_address, recursive=True)
    motion_parameters_list_of_folders = glob.glob(motion_parameters_address, recursive=True)
    temporal_mean_epi_single_frame_original_contrast_list_files = glob.glob(temporal_mean_epi_single_frame_original_contrast, recursive=True)

    # Apply affine transformations to temporal mean EPI files
    new_temporal_files = []
    print("Transforming Temporal Mean EPI Files")
    for temporal_file in tqdm(temporal_mean_epi_single_frame_original_contrast_list_files, desc="Transforming Temporal Mean EPI", unit="file"):
        func_folder = os.path.dirname(os.path.dirname(temporal_file))
        anat_folder = os.path.join(os.path.dirname(func_folder), "anat")
        transformed_files = apply_affine_transformations([temporal_file], func_folder, anat_folder, sigma_template_address)
        new_temporal_files.extend(transformed_files)
        print(f"Transformed Temporal Mean EPI file: {transformed_files[0]}")

    # Apply affine transformations to EPI files
    new_epi_files = []
    print("Transforming EPI Files")
    for epi_file in tqdm(epi_files_list, desc="Transforming EPI", unit="file"):
        func_folder = os.path.dirname(os.path.dirname(epi_file))
        anat_folder = os.path.join(os.path.dirname(func_folder), "anat")
        transformed_files = apply_affine_transformations([epi_file], func_folder, anat_folder, sigma_template_address)
        new_epi_files.extend(transformed_files)
        print(f"Transformed EPI file: {transformed_files[0]}")

    # Copy the newly created files to the results folder
    copy_files_to_results_folder(input_folder, new_epi_files, motion_parameters_list_of_folders, new_temporal_files)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Apply affine transformations to imaging data.")
    parser.add_argument('-i', '--input_folder', required=True, help="Path to the input folder.")
    args = parser.parse_args()
    main(args.input_folder)
