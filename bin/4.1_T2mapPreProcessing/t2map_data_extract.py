import argparse
import csv
import glob
import os
import sys

import nibabel as nii
import numpy as np

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
DEFAULT_LABEL_FILE = os.path.join(
    REPO_ROOT,
    "lib",
    "sigma",
    "SIGMA_InVivo_Anatomical_Brain_Atlas_Labels.txt",
)


def find_single_file(input_folder, pattern, description):
    matches = sorted(glob.glob(os.path.join(input_folder, pattern)))
    if len(matches) == 0:
        sys.exit("Error: No %s found in '%s'." % (description, input_folder,))
    if len(matches) > 1:
        sys.exit("Error: Multiple %s found in '%s': %s" % (description, input_folder, ", ".join(matches),))
    return matches[0]


sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))
from common.artifact_manifest import start_output_tracking

def getOutfile(atlas_type, img_file, suffix):
    imgName = os.path.basename(img_file)
    if imgName.endswith(".nii.gz"):
        t2map = imgName[:-len(".nii.gz")]
    elif imgName.endswith(".nii"):
        t2map = imgName[:-len(".nii")]
    else:
        t2map = os.path.splitext(imgName)[0]
    outFile = os.path.join(os.path.dirname(img_file), f"{t2map}_T2values_{atlas_type}_{suffix}.csv")
    return outFile


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
                print("Warning: Skipping invalid label line %i in '%s'." % (line_number, label_file,))
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
        sys.exit("Error: No labels could be read from '%s'." % (label_file,))

    return labels


def extractT2MapdataMean(img, rois, outfile, label_lookup):
    slices = np.unique(np.where(rois > 0)[2])
    regions = np.unique(rois)
    regions = regions[regions > 0]

    with open(outfile, "w", newline="") as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(["Slice", "Atlas IDs", "Names", "T2 Values", "Region Sizes"])

        for s in slices:
            for r in regions:
                region_voxels = np.where((rois[:, :, s] == r) & (rois[:, :, s] > 0))
                if len(region_voxels[0]) == 0:
                    continue
                mean_value = np.mean(img[region_voxels])
                region_size = len(region_voxels[0])
                label_id = int(round(float(r)))
                label_name = label_lookup.get(label_id, "")
                csv_writer.writerow([s, label_id, label_name, "%.2f" % mean_value, "%.2f" % region_size])


def extractT2MapdataPerRegion(img, rois, outfile, label_lookup):
    regions = np.unique(rois)
    regions = regions[regions > 0]

    with open(outfile, "w", newline="") as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(["Atlas IDs", "Names", "T2 Values", "Region Sizes"])

        for r in regions:
            region_voxels = np.where((rois == r) & (rois > 0))
            if len(region_voxels[0]) == 0:
                continue
            mean_value = np.mean(img[region_voxels])
            region_size = len(region_voxels[0])
            label_id = int(round(float(r)))
            label_name = label_lookup.get(label_id, "")
            csv_writer.writerow([label_id, label_name, "%.2f" % mean_value, "%.2f" % region_size])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extracts the T2 values from the T2 map for every atlas region")
    requiredNamed = parser.add_argument_group("Required named arguments")
    requiredNamed.add_argument("-i", "--input", help="Input T2 map, should be a NIfTI file", required=True)
    parser.add_argument("-l", "--label-file", help="Atlas label lookup file", default=DEFAULT_LABEL_FILE)
    args = parser.parse_args()

    image_file = args.input
    if not os.path.exists(image_file):
        sys.exit("Error: '%s' is not an existing image NIfTI file." % (image_file,))

    if not os.path.exists(args.label_file):
        sys.exit("Error: '%s' is not an existing label file." % (args.label_file,))

    image_dir = os.path.dirname(image_file)
    parental_atlas = find_single_file(image_dir, "*AnnoSplit_parental.nii*", "parental atlas")
    non_parental_atlas = find_single_file(image_dir, "*AnnoSplit.nii*", "non-parental atlas")

    print("Extracting T2 values for: %s" % image_file)
    print("Using label file: %s" % args.label_file)
    start_output_tracking(os.path.dirname(image_file), "t2map", "processing")

    img_data = nii.load(image_file)
    img = img_data.get_fdata()
    label_lookup = load_label_lookup(args.label_file)

    for atlas_type, atlas in (("parental", parental_atlas), ("non-parental", non_parental_atlas)):
        roi_data = nii.load(atlas)
        rois = roi_data.get_fdata()

        if img.shape[:3] != rois.shape[:3]:
            sys.exit("Error: image shape %s and atlas shape %s do not match." % (img.shape[:3], rois.shape[:3],))

        outFileMean = getOutfile(atlas_type, image_file, "Mean")
        print("Outfile (Mean): %s" % outFileMean)
        extractT2MapdataMean(img, rois, outFileMean, label_lookup)

        outFilePerRegion = getOutfile(atlas_type, image_file, "PerRegion")
        print("Outfile (Per Region): %s" % outFilePerRegion)
        extractT2MapdataPerRegion(img, rois, outFilePerRegion, label_lookup)

    print("Finished T2 map processing")
