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
import os
import tempfile

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
    Reorient only by swapping y and z axes.
    For 3D: in-memory.
    For 4D: chunked volume-wise to reduce RAM.
    Writes '*_originated.nii' and deletes input_file.
    """

    img = nib.load(input_file)
    hdr = img.header.copy()
    aff = img.affine.copy()

    out_file = input_file.replace('.nii.gz', '_originated.nii')
    shape = img.shape
    ndim = len(shape)

    out_dtype = np.float32
    hdr.set_data_dtype(out_dtype)

    if ndim == 3:
        # small enough: handle in RAM
        data = img.get_fdata(dtype=out_dtype)
        data = np.transpose(data, (0, 2, 1))  # swap y <-> z
        nib.save(nib.Nifti1Image(data, aff, hdr), out_file)

    elif ndim == 4:
        X, Y, Z, T = shape
        out_shape = (X, Z, Y, T)  # because we swap y <-> z

        # Create output memmap file
        with tempfile.TemporaryDirectory() as tmpd:
            mmap_path = os.path.join(tmpd, "reoriented.dat")
            mm = np.memmap(mmap_path, dtype=out_dtype, mode='w+', shape=out_shape)

            dataobj = img.dataobj  # lazy loader
            for t in range(T):
                vol = np.asanyarray(dataobj[..., t], dtype=out_dtype)  # only one volume loaded
                mm[..., t] = np.transpose(vol, (0, 2, 1))  # swap y <-> z
                if t % 50 == 0:
                    print(f"[reorient] processed volume {t + 1}/{T}")

            # Save without copying entire memmap into RAM
            img_out = nib.Nifti1Image(mm, aff, hdr)
            nib.save(img_out, out_file)

    else:
        raise ValueError(f"Unsupported ndim {ndim} for {input_file}")

    os.remove(input_file)
    return out_file


def resample_4d_in_batches(flo_file, ref_file, trans_file, out_file, batch_size=100):
    """
    Resample a 4D NIfTI file in smaller batches using NiftyReg (reg_resample).
    Each batch is temporarily saved to disk, processed, and then merged back into one 4D file.

    Parameters
    ----------
    flo_file : str
        Path to the 4D floating (input) file to be resampled.
    ref_file : str
        Path to the reference image (template).
    trans_file : str
        Path to the transformation matrix.
    out_file : str
        Path where the merged resampled 4D file will be saved.
    batch_size : int, optional
        Number of volumes per batch (default = 50).

    Returns
    -------
    str
        Path to the saved resampled 4D file.
    """
    import tempfile

    # Lade die 4D Datei
    img = nib.load(flo_file)
    data = img.get_fdata()
    affine = img.affine
    header = img.header.copy()

    if len(data.shape) != 4:
        raise ValueError(f"Expected 4D input for batch resampling, got shape {data.shape}")

    n_vols = data.shape[3]
    resampled_vols = []

    for start in range(0, n_vols, batch_size):
        end = min(start + batch_size, n_vols)
        print(f"[Batching] Resampling volumes {start}–{end-1}/{n_vols-1}")

        batch = data[..., start:end]

        # temporäres Verzeichnis für diesen Block
        with tempfile.TemporaryDirectory() as tmpdir:
            batch_file = os.path.join(tmpdir, "batch_in.nii.gz")
            nib.save(nib.Nifti1Image(batch, affine, header), batch_file)

            batch_out = os.path.join(tmpdir, "batch_out.nii.gz")

            # NiftyReg reg_resample aufrufen
            cmd = f"reg_resample -ref {ref_file} -flo {batch_file} -trans {trans_file} -res {batch_out} -inter 3"
            subprocess.run(shlex.split(cmd), check=True)

            # Ergebnis laden
            resampled_batch = nib.load(batch_out).get_fdata().astype(np.float32)
            print(f"[Batching] Finished batch {start}-{end-1}, shape {resampled_batch.shape}")

            resampled_vols.append(resampled_batch)

    # Überprüfen, ob wir überhaupt was gesammelt haben
    if not resampled_vols:
        raise RuntimeError(f"No batches were resampled for {flo_file}")

    # Debug: Shapes der Batches ausgeben
    for i, arr in enumerate(resampled_vols):
        print(f"[Debug] Batch {i} shape: {arr.shape}")

    # Alle Batches wieder zusammenfügen
    resampled_data = np.concatenate(resampled_vols, axis=3)
    print(f"[Merging] Final merged shape: {resampled_data.shape}")

    resampled_img = nib.Nifti1Image(resampled_data, affine, header)
    nib.save(resampled_img, out_file)
    print(f"[Done] Saved merged 4D result: {out_file}")

    return out_file

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
        img = nib.load(file_flipped)
        if len(img.shape) == 4 and img.shape[3] > 1:
            # 4D → Batch-Verarbeitung
            print('>> resampling 4D')
            file_st_f_on_template = resample_4d_in_batches(
                file_flipped, sigma_template_address, merged_inverted, file_st_f_on_template, batch_size=100
            )
        else:
            # 3D → normal resamplen
            print('>> resampling 3D')
            command = f"reg_resample -ref {sigma_template_address} -flo {file_flipped} -trans {merged_inverted} -res {file_st_f_on_template} -inter 3"
            subprocess.run(shlex.split(command), check=True)

        # Reorient the transformed file
        print(">> calling reorient_and_save:", file_st_f_on_template)
        reoriented_file = reorient_and_save(file_st_f_on_template)
        print(">> reoriented ->", reoriented_file)
        transformed_files.append(reoriented_file)
        print(">> removing flipped:", file_flipped)
        os.remove(file_flipped)
        print(">> flipped removed")
    
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
