import os
import argparse

def register_image(input_file, output_file, reference_image, affine_file):
    # Command to register MRI file using reg_aladin
    command = f"reg_aladin -ref {reference_image} -res {output_file} -flo {input_file} -aff {affine_file}"
    os.system(command)

def main(input_folder, output_folder):
    # Get list of MRI files
    mri_files = sorted([filename for filename in os.listdir(input_folder) if filename.endswith(".nii") or filename.endswith(".nii.gz")])

    # Reference image will be the first MRI file
    reference_image = os.path.join(input_folder, mri_files[0])

    # Create output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)

    # Iterate through each MRI file in the input folder (excluding the reference image)
    for i, filename in enumerate(mri_files[1:], start=1):
        # Path to the current MRI file
        input_file = os.path.join(input_folder, filename)
        
        # Output filename for the registered image
        output_file = os.path.join(output_folder, f"registered_{i}.nii.gz")
        
        # Output filename for the affine transformation
        affine_file = os.path.join(output_folder, f"affine_{i}.txt")
        
        # Register the image sequentially
        register_image(input_file, output_file, reference_image, affine_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Register MRI files using reg_aladin.")
    parser.add_argument("-i", "--input", required=True, help="Path to the folder containing MRI files.")
    parser.add_argument("-o", "--output", required=True, help="Path to the output folder where registered images will be saved.")
    args = parser.parse_args()
    
    main(args.input, args.output)
