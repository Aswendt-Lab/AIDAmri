#!/usr/bin/env python
# Crop T2-weighted images to a usable FOV
#
# Usage: crop_T2.py -i <input_dir> -x_min <x_min> -x_max <x_max> -y_min <y_min> -y_max <y_max> -z_min <z_min> -z_max <z_max> -o <output_dir>

from nipype.interfaces.fsl import ExtractROI
import nibabel as nii
import argparse
import pathlib
import subprocess

def crop_T2(input_file, x_min, x_max, y_min, y_max, z_min, z_max, output_path):
    # Get dimensions (x_min, x_max, y_min, y_max, z_min, z_max) of input_file
    img = nii.load(input_file)
    img_data = img.get_fdata()
    img_dims = img_data.shape
    print(f"Input image dimensions: {img_dims}")
    crop_z = False
    if crop_z is False:
        z_min = 0
        z_size = img_dims[2]
    else:
        if z_max > img_dims[2]:
            z_size = img_dims[2]
    x_size = x_max - x_min
    y_size = y_max - y_min
    # Use FSL's fslroi to crop the image
    fslroi = ExtractROI()
    fslroi.inputs.in_file = input_file
    fslroi.inputs.roi_file = output_path
    fslroi.inputs.x_min = x_min
    fslroi.inputs.x_size = x_size
    fslroi.inputs.y_min = y_min
    fslroi.inputs.y_size = y_size
    fslroi.inputs.z_min = z_min
    fslroi.inputs.z_size = z_size
    fslroi.inputs.t_min = 0
    fslroi.inputs.t_size = 1
    print(f"Output image dimensions: ({fslroi.inputs.x_size}, {fslroi.inputs.y_size}, {fslroi.inputs.z_size})")
    fslroi.run()
    print(f"Output image saved to: {output_path}")
    print(f"Check {output_path} to ensure brain tissue was not removed in cropping")
    # Create QC images to check cropping using FSL's slicer via subprocess
    for nifti in [input_file, output_path]:
        base = pathlib.Path(nifti).with_suffix('').with_suffix('')  # Remove .nii.gz or .nii
        slicer_path = base.with_suffix('.png')
        cmd = [
            "slicer",
            nifti,
            "-L",  # show middle slices
            "-x", "0",  # no crosshairs
            "-w", "800",  # image width
            "-a",
            str(slicer_path)
        ]
        try:
            subprocess.run(cmd, check=True)
            print(f"Overlay image saved to: {slicer_path}")
        except subprocess.CalledProcessError as e:
            print(f"Failed to create overlay image for {nifti}: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Crop T2-weighted images to a usable FOV")
    parser.add_argument("-i", "--input_file", required=True, help="Input file (T2 image)")
    parser.add_argument("-x_min", type=int, required=True, help="Minimum x index")
    parser.add_argument("-x_max", type=int, required=True, help="Maximum x index")
    parser.add_argument("-y_min", type=int, required=True, help="Minimum y index")
    parser.add_argument("-y_max", type=int, required=True, help="Maximum y index")
    parser.add_argument("-z_min", type=int, required=True, help="Minimum z index")
    parser.add_argument("-z_max", type=int, required=True, help="Maximum z index")
    parser.add_argument("-o", "--output_file", required=True, help="Output file for cropped image")
    args = parser.parse_args()

    input_file = args.input_file
    output_path = args.output_file

    crop_T2(
        input_file,
        args.x_min,
        args.x_max,
        args.y_min,
        args.y_max,
        args.z_min,
        args.z_max,
        output_path
    )