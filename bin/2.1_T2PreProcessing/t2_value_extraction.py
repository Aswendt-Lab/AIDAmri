import atexit
import nibabel as nii
import numpy as np
import argparse
import os
import glob
import csv
import sys


def get_default_label_file(script_dir):
    repo_root = os.path.abspath(os.path.join(script_dir, os.pardir, os.pardir))
    candidates = [
        os.path.join(repo_root, "lib", "sigma", "SIGMA_InVivo_Anatomical_Brain_Atlas_Labels.txt"),
        os.path.join(repo_root, "lib", "SIGMA_InVivo_Anatomical_Brain_Atlas_Labels.txt"),
        os.path.join(repo_root, "lib", "SIGMA_InVivo_Anatomical_Brain_Atlas_Labels.txt"),
        os.path.join(script_dir, "SIGMA_InVivo_Anatomical_Brain_Atlas_Labels.txt"),
    ]
    for candidate in candidates:
        if os.path.exists(candidate):
            return candidate
    return candidates[0]


def load_label_lookup(label_file):
    labels = {}
    with open(label_file, "r") as label_handle:
        for line_number, line in enumerate(label_handle, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue

            parts = stripped.split()
            try:
                label_id = int(float(parts[0]))
            except (IndexError, ValueError):
                print(f"Warning: Skipping invalid label line {line_number} in '{label_file}'.")
                continue

            if '"' in stripped:
                label_name = stripped.split('"', 1)[1].rsplit('"', 1)[0]
            elif "\t" in stripped:
                tab_parts = stripped.split("\t")
                label_name = tab_parts[1].strip() if len(tab_parts) > 1 else ""
            else:
                label_name = " ".join(parts[1:])

            labels[label_id] = label_name

    if len(labels) == 0:
        sys.exit(f"Error: No labels could be read from '{label_file}'.")

    return labels

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))
from common.artifact_manifest import start_output_tracking
from common.script_logging import append_git_information

def getOutfile(atlas_type, img_file, suffix):
    imgName = os.path.basename(img_file)
    if imgName.endswith(".nii.gz"):
        t2map = imgName[:-7]
    elif imgName.endswith(".nii"):
        t2map = imgName[:-4]
    else:
        t2map = os.path.splitext(imgName)[0]
    acronym_name = os.path.basename(atlas_type).split('.')[0]
    outFile = os.path.join(os.path.dirname(img_file),"t2_values_extraction",f"{t2map}_T2values_{acronym_name}_{suffix}.csv")
    return outFile

def extractT2MapdataMean(img, rois, outfile, label_lookup):
    slices = np.unique(np.where(rois > 0)[2])
    regions = np.unique(rois)
    regions = regions[regions > 0]

    with open(outfile, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(["Slice", "Atlas IDs", "Names", "Mean T2 Values", "Region Sizes"])
        
        for s in slices:
            for r in regions:
                region_voxels = np.where((rois[:, :, s] == r) & (rois[:, :, s] > 0))
                if len(region_voxels[0]) == 0:
                    continue
                mean_value = np.mean(img[region_voxels])
                region_size = len(region_voxels[0])
                label_id = int(round(r))
                label_name = label_lookup.get(label_id, "")
                csv_writer.writerow([s, label_id, label_name, "%.2f" % mean_value, "%.2f" % region_size])

def extractT2MapdataPerRegion(img, rois, outfile, label_lookup):
    regions = np.unique(rois)
    regions = regions[regions > 0]

    with open(outfile, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(["Atlas IDs", "Names", "Mean T2 Values", "Region Sizes"])
        
        for r in regions:
            region_voxels = np.where((rois == r) & (rois > 0))
            if len(region_voxels[0]) == 0:
                continue
            mean_value = np.mean(img[region_voxels])
            region_size = len(region_voxels[0])
            label_id = int(round(r))
            label_name = label_lookup.get(label_id, "")
            csv_writer.writerow([label_id, label_name, "%.2f" % mean_value, "%.2f" % region_size])

if __name__ == '__main__':
    script_dir = os.path.dirname(os.path.abspath(__file__))
    parser = argparse.ArgumentParser(description='Extracts the T2w values for every atlas region')
    requiredNamed = parser.add_argument_group('Required named arguments')
    requiredNamed.add_argument('-i', '--input', help='Input T2w file, should be a nifti file', required=True)
    parser.add_argument(
        '-l',
        '--label-file',
        help='Atlas label lookup file. Defaults to the SIGMA rat atlas labels if available.',
        default=get_default_label_file(script_dir),
    )
    args = parser.parse_args()
    atexit.register(
        append_git_information,
        output_dir=os.path.dirname(os.path.abspath(args.input)) if args.input else None,
    )

    print(f"Extracting T2 values for: {args.input}")
    print(f"Label file: {args.label_file}")
    if not os.path.exists(args.label_file):
        sys.exit(f"Error: Label file does not exist: '{args.label_file}'.")
    label_lookup = load_label_lookup(args.label_file)

    image_file = args.input
    if not os.path.exists(image_file):
        sys.exit(f"Error: '{image_file}' is not an existing image nii-file.")
    start_output_tracking(os.path.dirname(image_file), "anat", "registration")

    img_data = nii.load(image_file)
    img = img_data.get_fdata()  # Using get_fdata() for compatibility
    
    image_dir = os.path.dirname(image_file)
    parental_atlases = sorted(glob.glob(os.path.join(image_dir, "*AnnoSplit_parental.nii*")))
    non_parental_atlases = sorted(glob.glob(os.path.join(image_dir, "*AnnoSplit.nii*")))

    if len(parental_atlases) == 0:
        sys.exit(f"Error: No parental atlas file '*AnnoSplit_parental.nii*' found in '{image_dir}'.")
    if len(non_parental_atlases) == 0:
        sys.exit(f"Error: No non-parental atlas file '*AnnoSplit.nii*' found in '{image_dir}'.")

    parental_atlas = parental_atlases[0]
    non_parental_atlas = non_parental_atlases[0]
    
    output_dir = os.path.join(os.path.dirname(image_file), "t2_values_extraction")
    os.makedirs(output_dir, exist_ok=True)

    atlas_jobs = [
        ("non-parental", non_parental_atlas),
        ("parental", parental_atlas),
    ]

    for atlas_type, atlas in atlas_jobs:
        try:
            roi_data = nii.load(atlas)
            rois = roi_data.get_fdata()  # Using get_fdata() for compatibility

            outFileMean = getOutfile(atlas_type, image_file, "SliceWisePerRegion")
            print(f"Outfile (Slice-wise mean per region): {outFileMean}")
            extractT2MapdataMean(img, rois, outFileMean, label_lookup)

            outFilePerRegion = getOutfile(atlas_type, image_file, "PerRegion")  # Fixed suffix to "PerRegion"
            print(f"Outfile (Per Region): {outFilePerRegion}")
            extractT2MapdataPerRegion(img, rois, outFilePerRegion, label_lookup)
        except Exception as e:
            print(f'Error while processing the T2 values: {str(e)}')  # Improved error message
            raise  # Raising the exception to halt execution

    print("Finished T2 map processing")
