import nibabel as nii
import numpy as np 
import argparse
import os
import glob
import logging



def getOutfile(roi_file,img_file, acronmys):
    imgName = os.path.basename(img_file)

    t2map = str.split(imgName,'.')[-3]

    acronym_name = str.split(os.path.basename(acronmys),'.')[0]

    outFile = os.path.join(os.path.dirname(img_file),t2map + "_T2values_" + acronym_name + '.txt')

    return outFile


def extractT2Mapdata(img,rois,outfile,txt_file):
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
        print(rois)
        print(r)
        test = test3234
        paramValue = np.mean(img[rois==r])
        if indices is not None:
            if len(np.argwhere(indices==r)) == 0:
                continue
            idx = int(np.argwhere(indices==r)[0])
            str_idx = ref_lines[idx]
            acro = str.split(str_idx,'\t')[1][:-1]
            fileID.write("%i\t%s\t%.2f\n" % (r,acro,paramValue))
        else:
            fileID.write("%i\t%.2f\n" % (r, paramValue))


    fileID.close()
    return outfile



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extracts the intensity of the t2map of every region')
    requiredNamed = parser.add_argument_group('Required named arguments')
    requiredNamed.add_argument('-i','--input', help='Input t2map, should be a nifti file')
    args = parser.parse_args()

    log_file_path = os.path.join(os.path.dirname(args.input), "process.txt")
    logging.basicConfig(filename=log_file_path, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    acronyms_files = glob.glob(os.path.join(os.getcwd(),"*.txt"))
    logging.info(f"Extracting T2values for: {args.input}")
    logging.info(f"Acronym files: {acronyms_files}")
  
    # read image data
    if args.input is not None:
        image_file = args.input
        if not os.path.exists(image_file):
            sys.exit("Error: '%s' is not an existing image nii-file." % (image_file))

    img_data=nii.load(image_file)
    img = img_data.dataobj.get_unscaled()

    parental_atlas = glob.glob(os.path.join(os.path.dirname(image_file), "*AnnoSplit_t2map.nii*"))[0]
    non_parental_atlas = glob.glob(os.path.join(os.path.dirname(image_file), "*AnnoSplit.nii*"))[0]

    print(parental_atlas)
    print(non_parental_atlas)
    

    for acronmys in acronyms_files:
        try:
            if "parentARA_LR" in acronmys:
                atlas = parental_atlas
            else:
                atlas = non_parental_atlas
            
            roi_data = nii.load(atlas)
            rois = roi_data.dataobj.get_unscaled()
               
            outFile = getOutfile(atlas, image_file, acronmys)
            logging.info(f"Outifle: {outFile}")
            file = extractT2Mapdata(img,rois,outFile,acronmys)
        except Exception as e:
            logging.error(f'Error while processing the T2values Errorcode: {str(e)}')
            raise

    logging.info("Finished t2map processing")
