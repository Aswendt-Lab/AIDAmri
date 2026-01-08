import os
import argparse
import nibabel as nib

def extracted_brain_output(directory, output_file):
    # Find files with "brain" in the name
    brain_files = [f for f in os.listdir(directory) if "brain" in f and "mask" not in f and os.path.isfile(os.path.join(directory, f))]
    if len(brain_files) < 2:
        raise ValueError("Less than two files with 'brain' in the name found in the directory.")

    # Take the first two different files
    file1, file2 = brain_files[:2]
    path1 = os.path.join(directory, file1)
    path2 = os.path.join(directory, file2)

    # Load images
    img1 = nib.load(path1)
    img2 = nib.load(path2)

    # Sum voxel values
    sum1 = img1.get_fdata().sum()
    sum2 = img2.get_fdata().sum()

    # Select file with lower sum
    if sum1 < sum2:
        selected_file = path1
    else:
        selected_file = path2

    # Rename selected file to output_file
    output_path = os.path.join(directory, output_file)
    os.rename(selected_file, output_path)
    print(f"Selected file: {selected_file} -> {output_path}")

    return output_path

def main():
    parser = argparse.ArgumentParser(description="Pick BET output with lower voxel sum.")
    parser.add_argument("--directory", help="Directory containing brain files")
    parser.add_argument("--output_file", help="Output file name for selected file")
    args = parser.parse_args()
    extracted_brain_output(args.directory, args.output_file)

if __name__ == "__main__":
    main()