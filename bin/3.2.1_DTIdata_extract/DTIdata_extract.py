""""
Created on 06.04.2019

@authors: Niklas Pallast
latest Update: Julian Cramer


"""

import os
import sys
import argparse
import numpy as np
import nibabel as nii

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, os.pardir, os.pardir))
DEFAULT_SIGMA_LABEL_FILE = os.path.join(
    REPO_ROOT,
    "lib",
    "sigma",
    "SIGMA_InVivo_Anatomical_Brain_Atlas_Labels.txt",
)


def getOutfile(roi_file,img_file):
    imgName = os.path.basename(img_file)
    baseName = str.split(os.path.basename(roi_file),'.')[0]
    dtiParam = str.split(imgName,'.')[-3]

    print('\nStart processing DTI parameter: %s' % str.upper(dtiParam))
    outFile = os.path.join(os.path.dirname(img_file),baseName+'_'+str.split(imgName,'.')[-3])+'.txt'
    return outFile


def label_sidecar_candidates(roi_file):
    candidates = []
    if roi_file.endswith(".nii.gz"):
        candidates.append(roi_file[:-len(".nii.gz")] + ".txt")
        candidates.append(roi_file[:-len(".gz")] + ".txt")
    elif roi_file.endswith(".nii"):
        candidates.append(roi_file[:-len(".nii")] + ".txt")
        candidates.append(roi_file + ".txt")
    else:
        candidates.append(roi_file + ".txt")

    candidates.append(DEFAULT_SIGMA_LABEL_FILE)
    return candidates


def find_label_file(roi_file):
    for candidate in label_sidecar_candidates(roi_file):
        if os.path.exists(candidate):
            return candidate
    return None


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


def extractDTIData(img,rois,outfile,txt_file):
    if img.shape[:3] != rois.shape[:3]:
        sys.exit(f"Error: image shape {img.shape} and roi shape {rois.shape} do not match.")

    regions = np.unique(rois)
    regions = regions[regions > 0]

    mapping = None
    if txt_file is not None:
        mapping = load_label_lookup(txt_file)

    # robust parametername from outfile
    param_name = os.path.basename(outfile).replace(".txt", "").split("_")[-1].upper()

    with open(outfile, "w") as fileID:
        fileID.write(f"{param_name} values for {regions.size} given regions:\n\n")

        for r in regions:
            msk = (rois == r)
            if not np.any(msk):
                continue
            paramValue = float(np.mean(img[msk]))
            label_id = int(round(float(r)))

            if mapping is not None:
                acro = mapping.get(label_id, "NA")
                fileID.write(f"{label_id}\t{acro}\t{paramValue:.2f}\n")
            else:
                fileID.write(f"{label_id}\t{paramValue:.2f}\n")

    return outfile



if __name__ == '__main__':
    # default values


    parser = argparse.ArgumentParser(description='Extracts the major DTI parameters (apparent diffusion coefficients) '
                                                 'axial diffusivity (AD), fractional anisotropy (FA), mean diffusivity (MD), and radial diffusivity (RD)')
    requiredNamed = parser.add_argument_group('Required named arguments')
    requiredNamed.add_argument('image_file', help='Input file of AIDA pipeline with related folder')
    requiredNamed.add_argument('roi_file', help='Input file of related roi')
    parser.add_argument('-t', '--translatorTXT',
                        help='txt file to translate ROI Number to label names. If omitted, the ROI sidecar .txt or SIGMA labels are used.',
                        type=str)
    args = parser.parse_args()

    # read image data
    image_file=None
    if args.image_file is not None and args.image_file is not None:
        image_file = args.image_file
        if not os.path.exists(image_file):
            sys.exit("Error: '%s' is not an existing image nii-file." % (image_file))

    img_data=nii.load(image_file)
    img = img_data.get_fdata(dtype=np.float32)

    # read roi data
    roi_file = None
    if args.roi_file is not None and args.roi_file is not None:
        roi_file = args.roi_file
        if not os.path.exists(roi_file):
            sys.exit("Error: '%s' is not an existing roi file." % (roi_file))

    # read translation TXT file
    txt_file = None
    if args.translatorTXT is not None:
        txt_file = args.translatorTXT
        if  not os.path.exists(args.translatorTXT):
            sys.exit("Error: '%s' is not an existing translation txt file." % (txt_file))

    roi_data = nii.load(roi_file)
    rois = np.asanyarray(roi_data.dataobj).copy()

    if txt_file is None:
        txt_file = find_label_file(roi_file)
        if txt_file is not None:
            print("Using label file: %s" % txt_file)

    outFile = getOutfile(roi_file, image_file)
    file = extractDTIData(img,rois,outFile,txt_file)
    print("\033[0;30;42m Done \33[0m'  %s" % file)
    # save output image and txtFile
    #save_data(image_out, peaks)
