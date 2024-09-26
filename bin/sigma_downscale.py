import nibabel as nib
import numpy as np
from scipy.ndimage import zoom

# Load the original NIfTI file
original_file_path = "E:\CRC_data\multiverse\multivereAIDAmri\lib\sigma\SIGMA_InVivo_Anatomical_Brain_Atlas.nii"
SigmaOriginal_img = nib.load(original_file_path)
SigmaOriginal_data = SigmaOriginal_img.get_fdata()
SigmaOriginal_affine = SigmaOriginal_img.affine
SigmaOriginal_header = SigmaOriginal_img.header

# Define the new voxel size and new image size
new_voxel_size = [0.3, 0.3, 0.3]  # New pixel dimensions in mm
new_image_size = [64, 109,64]    # New target image size

# Calculate the zoom factor for resizing
zoom_factors = np.array(new_image_size) / np.array(SigmaOriginal_data.shape)

# Resample the image to the new size
SigmaResized_data = zoom(SigmaOriginal_data, zoom_factors, order=0)  # Using cubic interpolation (order=3)

# Update the affine matrix: adjust the scaling based on new voxel size
scaling_factors = np.diag(new_voxel_size + [1])
new_affine = SigmaOriginal_affine @ scaling_factors  # Adjust affine matrix for voxel scaling

# Update the NIfTI header to reflect the new voxel size
SigmaOriginal_header['pixdim'][1:4] = new_voxel_size  # Update voxel size
SigmaOriginal_header.set_data_shape(new_image_size)   # Update the image size in the header

# Create a new NIfTI image with the resized data and updated affine and header
SigmaResized_img = nib.Nifti1Image(SigmaResized_data, new_affine, SigmaOriginal_header)

# Define the new file path for saving the modified NIfTI file
new_file_path = "E:\CRC_data\multiverse\multivereAIDAmri\lib\sigma/SIGMA_InVivo_Anatomical_Brain_Atlas_downsampled.nii"

# Save the resized NIfTI image with updated metadata
nib.save(SigmaResized_img, new_file_path)

# Display success message
print(f'Modified NIfTI file saved to: {new_file_path}')
