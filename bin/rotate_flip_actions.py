# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 19:43:22 2024

@author: arefks
"""

import os
import glob
import nibabel as nib
import numpy as np

# Set the initial folder path
folder_path = r"E:/CRC_data/multiverse/proc_data/orient_3"  # Change this to your folder path
file_extension = '*.nii.gz'  # Adjust the file extension to match your files

# Use glob to find all files with the specified extension
file_list = glob.glob(os.path.join(folder_path, "**", file_extension), recursive=True)

# Iterate over all found files
for file_path in file_list:
    # Load the NIfTI data
    print("fuck")
    data = nib.load(file_path)
    imgTemp = data.get_fdata()

    # Determine if the file is 3D or 4D based on the number of dimensions
    if len(imgTemp.shape) == 3:
        print(f"Processing 3D file: {file_path}")

        # 3D file processing
        # Permute the image data dimensions (example permutation)
        imgTemp = np.transpose(imgTemp, (0, 2, 1))
        imgTemp = np.flip(imgTemp, 1)
        # Adjust the header fields based on the permutation
        header = data.header
        header['dim'][1:4] = [imgTemp.shape[0], imgTemp.shape[2], imgTemp.shape[1]]
        header['pixdim'][1:4] = [header['pixdim'][1], header['pixdim'][3], header['pixdim'][2]]

        # Apply additional flips or rotations if needed
        # imgTemp = np.flip(imgTemp, 2)
        
        # imgTemp = np.flip(imgTemp, 0)
        # imgTemp = np.rot90(imgTemp, 2)

    elif len(imgTemp.shape) == 4:
        print(f"Processing 4D file: {file_path}")

        # 4D file processing
        # Permute the image data dimensions (example permutation)
        imgTemp = np.transpose(imgTemp, (0, 2, 1, 3))
        imgTemp = np.flip(imgTemp, 1)
        # Adjust the header fields based on the permutation
        header = data.header
        header['dim'][1:5] = [imgTemp.shape[0], imgTemp.shape[2], imgTemp.shape[1], imgTemp.shape[3]]
        header['pixdim'][1:5] = [header['pixdim'][1], header['pixdim'][3], header['pixdim'][2], header['pixdim'][4]]

        # Apply additional flips or rotations if needed
        # imgTemp = np.flip(imgTemp, 2)
        #imgTemp = np.flip(imgTemp, 1)
        # imgTemp = np.flip(imgTemp, 0)
        # imgTemp = np.rot90(imgTemp, 2)

    else:
        print(f"File {file_path} is neither 3D nor 4D, skipping.")
        continue

    # Create a new NIfTI image with the updated data and header
    new_img = nib.Nifti1Image(imgTemp, data.affine, header)

    # Save the new NIfTI image to the same file path, overwriting the old file
    outPath = file_path.replace("orient_3","orient_3_new")
    nib.save(outPath, file_path)

    print(f"Processed and saved: {file_path}")
