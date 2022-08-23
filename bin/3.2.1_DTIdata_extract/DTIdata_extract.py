""""
Created on 06.04.2019

@authors: Niklas Pallast


"""

import os
import sys
import argparse
import numpy as np
import nibabel as nii

def getOutfile(roi_file,img_file):
    imgName = os.path.basename(img_file)
    baseName = str.split(os.path.basename(roi_file),'.')[0]
    dtiParam = str.split(imgName,'.')[-3]

    print('\nStart processing DTI parameter: %s' % str.upper(dtiParam))
    outFile = os.path.join(os.path.dirname(img_file),baseName+'_'+str.split(imgName,'.')[-3])+'.txt'
    return outFile

def extractDTIData(img,rois,outfile,txt_file):

    regions = np.uint16(np.unique(rois))
    regions=np.delete(regions,0)

    indices = None
    if txt_file is not None:

        ref_lines = open(txt_file).readlines()
        indices = np.zeros_like(ref_lines)
        for idx in range(np.size(ref_lines)):
            curNum = int(str.split(ref_lines[idx], '\t')[0])

            indices[idx] = curNum
        indices = np.uint16(indices)

    fileID = open(outfile, 'w')
    fileID.write("%s values for %i given regions:\n\n" % (str.upper(outfile[-6:-4]),np.size(regions)))

    for r in regions:
        paramValue = np.mean(img[rois==r])
        if indices is not None:
            str_idx = ref_lines[int(np.argwhere(indices == r)[0])]
            acro = str.split(str_idx,'\t')[1][:-1]
            fileID.write("%i\t%s\t%.2f\n" % (r,acro ,paramValue))
        else:
            fileID.write("%i\t%.2f\n" % (r, paramValue))

    fileID.close()
    return outfile



if __name__ == '__main__':
    # default values


    parser = argparse.ArgumentParser(description='Extracts the major DTI parameters (apparent diffusion coefficients) '
                                                 'axial diffusivity (AD), fractional anisotropy (FA), mean diffusivity (MD), and radial diffusivity (RD)')
    requiredNamed = parser.add_argument_group('Required named arguments')
    requiredNamed.add_argument('image_file', help='Input file of AIDA pipeline with related folder')
    requiredNamed.add_argument('roi_file', help='Input file of related roi')
    parser.add_argument('-t', '--translatorTXT',
                        help='txt file to translate ROI Number to acronyms',type=str)
    args = parser.parse_args()

    # read image data
    image_file=None
    if args.image_file is not None and args.image_file is not None:
        image_file = args.image_file
        if not os.path.exists(image_file):
            sys.exit("Error: '%s' is not an existing image nii-file." % (image_file))

    img_data=nii.load(image_file)
    img = img_data.get_data()

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
    rois = roi_data.get_data()

    outFile = getOutfile(roi_file, image_file)
    file = extractDTIData(img,rois,outFile,txt_file)
    print("\033[0;30;42m Done \33[0m'  %s" % file)
    # save output image and txtFile
    #save_data(image_out, peaks)
