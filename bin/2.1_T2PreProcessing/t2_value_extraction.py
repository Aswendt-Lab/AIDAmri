import nibabel as nii
import numpy as np
import argparse
import os
import glob
import csv
import sys  # Added import statement for sys module

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))
from common.artifact_manifest import start_output_tracking

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

def extractT2MapdataMean(img, rois, outfile, txt_file):
    slices = np.unique(np.where(rois > 0)[2])
    regions = np.unique(rois)
    regions = regions[regions > 0]
    
    indices = {}
    if txt_file is not None:
        ref_lines = open(txt_file).readlines()
        indices = {int(line.split('\t')[0]): line.split('\t')[1].strip() for line in ref_lines}
    
    with open(outfile, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(["Slice", "ARA IDs", "Names", "Mean T2 Values", "Region Sizes"])
        
        for s in slices:
            for r in regions:
                region_voxels = np.where((rois[:, :, s] == r) & (rois[:, :, s] > 0))
                if len(region_voxels[0]) == 0:
                    continue
                mean_value = np.mean(img[region_voxels])
                region_size = len(region_voxels[0])
                acro = indices.get(r, "")  # Using dict.get() to avoid KeyError
                csv_writer.writerow([s, r, acro, "%.2f" % mean_value, "%.2f" % region_size])

def extractT2MapdataPerRegion(img, rois, outfile, txt_file):
    regions = np.unique(rois)
    regions = regions[regions > 0]
    
    indices = {}
    if txt_file is not None:
        ref_lines = open(txt_file).readlines()
        indices = {int(line.split('\t')[0]): line.split('\t')[1].strip() for line in ref_lines}
    
    with open(outfile, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(["ARA IDs", "Names", "Mean T2 Values", "Region Sizes"])
        
        for r in regions:
            region_voxels = np.where((rois == r) & (rois > 0))
            if len(region_voxels[0]) == 0:
                continue
            mean_value = np.mean(img[region_voxels])
            region_size = len(region_voxels[0])
            acro = indices.get(r, "")  # Using dict.get() to avoid KeyError
            csv_writer.writerow([r, acro, "%.2f" % mean_value, "%.2f" % region_size])

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extracts the T2w values for every atlas region')
    requiredNamed = parser.add_argument_group('Required named arguments')
    requiredNamed.add_argument('-i', '--input', help='Input T2w file, should be a nifti file')
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    acronyms_files = sorted(glob.glob(os.path.join(script_dir, "*.txt")))
    print(f"Extracting T2 values for: {args.input}")
    print(f"Acronym files: {acronyms_files}")
    if len(acronyms_files) == 0:
        sys.exit(f"Error: No acronym text files '*.txt' found in '{script_dir}'.")

    # Checking if input file is provided
    if args.input is None:
        sys.exit("Error: No input file provided.")

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
    
    if not os.path.exists(os.path.join(os.path.dirname(image_file), "t2_values_extraction")):
        os.mkdir(os.path.join(os.path.dirname(image_file), "t2_values_extraction"))

    for acronyms in acronyms_files:  # Corrected variable name to acronyms
        try:
            if "parentARA_LR" in acronyms:
                atlas_type = "parental"
                atlas = parental_atlas
            else:
                atlas_type = "non-parental"
                atlas = non_parental_atlas

            roi_data = nii.load(atlas)
            rois = roi_data.get_fdata()  # Using get_fdata() for compatibility

            outFileMean = getOutfile(atlas_type, image_file, "SliceWisePerRegion")
            print(f"Outfile (Slice-wise mean per region): {outFileMean}")
            extractT2MapdataMean(img, rois, outFileMean, acronyms)

            outFilePerRegion = getOutfile(atlas_type, image_file, "PerRegion")  # Fixed suffix to "PerRegion"
            print(f"Outfile (Per Region): {outFilePerRegion}")
            extractT2MapdataPerRegion(img, rois, outFilePerRegion, acronyms)
        except Exception as e:
            print(f'Error while processing the T2 values: {str(e)}')  # Improved error message
            raise  # Raising the exception to halt execution

    print("Finished T2 map processing")
