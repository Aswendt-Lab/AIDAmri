"""
Created on 11/09/2023, Updated on 02/09/2024

@author: Marc Schneider, Markus Aswendt
Neuroimaging & Neuroengineering
Department of Neurology
University Hospital Cologne
"""

import nibabel as nii
import numpy as np
import argparse
import os
import glob


def getOutfile(roi_file, img_file, acronmys):
    imgName = os.path.basename(img_file)

    t2map = str.split(imgName, '.')[-3]

    acronym_name = str.split(os.path.basename(acronmys), '.')[0]

    outFile = os.path.join(os.path.dirname(img_file), t2map + "_T2values_" + acronym_name + '.txt')

    return outFile


def extractT2Mapdata(img, rois, outfile, txt_file):
    regions = np.uint16(np.unique(rois))
    regions = np.delete(regions, 0)  # Exclude background label (0)

    indices = None
    region_sizes = np.zeros_like(regions, dtype=float)

    if txt_file is not None:
        ref_lines = open(txt_file).readlines()
        indices = np.zeros_like(ref_lines)
        for idx in range(np.size(ref_lines)):
            curNum = int(str.split(ref_lines[idx], '\t')[0])
            indices[idx] = curNum
        indices = np.uint16(indices)

    fileID = open(outfile, 'w')
    fileID.write("%s ARA IDs,  names,  T2 values, and regions sizes (separated by TAB) for %i given regions:\n\n" % (
    str.upper(outfile[-6:-4]), np.size(regions)))

    for idx, r in enumerate(regions):
        region_size = np.sum(rois == r)  # Calculate region size
        paramValue = np.mean(img[rois == r])  # Calculate T2 value

        if indices is not None:
            if len(np.argwhere(indices == r)) == 0:
                continue
            acro_idx = int(np.argwhere(indices == r)[0, 0])
            acro = str.split(ref_lines[acro_idx], '\t')[1][:-1]
            fileID.write("%i\t%s\t%.2f\t%.2f\n" % (r, acro, paramValue, region_size))
        else:
            fileID.write("%i\t%.2f\t%.2f\n" % (r, paramValue, region_size))

    fileID.close()
    return outfile


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extracts the T2 values from the T2 map for every atlas region')
    requiredNamed = parser.add_argument_group('Required named arguments')
    requiredNamed.add_argument('-i', '--input', help='Input T2 map, should be a nifti file')
    args = parser.parse_args()

    acronyms_files = glob.glob(os.path.join(os.getcwd(), "*.txt"))
    print(f"Extracting T2 values for: {args.input}")
    print(f"Acronym files: {acronyms_files}")

    # read image data
    if args.input is not None:
        image_file = args.input
        if not os.path.exists(image_file):
            sys.exit("Error: '%s' is not an existing image nii-file." % (image_file))

    img_data = nii.load(image_file)
    img = img_data.dataobj.get_unscaled()

    parental_atlas = glob.glob(os.path.join(os.path.dirname(image_file), "*AnnoSplit_par.nii*"))[0]
    non_parental_atlas = glob.glob(os.path.join(os.path.dirname(image_file), "*AnnoSplit.nii*"))[0]

    for acronmys in acronyms_files:
        try:
            if "parentARA_LR" in acronmys:
                atlas = parental_atlas
            else:
                atlas = non_parental_atlas

            roi_data = nii.load(atlas)
            rois = roi_data.dataobj.get_unscaled()

            outFile = getOutfile(atlas, image_file, acronmys)
            print(f"Outifle: {outFile}")
            file = extractT2Mapdata(img, rois, outFile, acronmys)
        except Exception as e:
            print(f'Error while processing the T2 values Errorcode: {str(e)}')
            raise

    print("Finished T2 map processing")
