# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 14:17:16 2024

@author: arefks
"""

import os
import numpy as np
import nibabel as nib

# Path to the folder containing registered MRI files
input_folder = r"E:\CRC_data\SP_allrestingstate_BiasBets_registered"

# Use the same folder for output
output_folder = input_folder

# Get list of registered MRI files
registered_files = [filename for filename in os.listdir(input_folder) if filename.startswith("registered_")]

# Load the first registered file to get image shape and affine matrix
first_file = os.path.join(input_folder, registered_files[0])
first_img = nib.load(first_file)
data_shape = first_img.shape
affine = first_img.affine

# Initialize an array to store voxel-wise sum
sum_data = np.zeros(data_shape)

# Iterate through each registered file and compute sum
for filename in registered_files:
    file_path = os.path.join(input_folder, filename)
    img = nib.load(file_path)
    data = img.get_fdata()
    sum_data += data

# Compute voxel-wise mean
average_data = sum_data / len(registered_files)

# Save the averaged image
output_file = os.path.join(output_folder, "AVERAGED_template.nii.gz")
averaged_img = nib.Nifti1Image(average_data, affine)
nib.save(averaged_img, output_file)
