"""
Created on 10/08/2017

@author: Niklas Pallast
Neuroimaging & Neuroengineering
Department of Neurology
University Hospital Cologne

"""


import nipype.interfaces.fsl as fsl
import os, sys
import nibabel as nii
import numpy as np
import applyMICO
import cv2
from pathlib import Path
import subprocess
import shutil


def reset_orientation(input_file):

    brkraw_dir = os.path.join(os.path.dirname(input_file), "brkraw")
    if os.path.exists(brkraw_dir):
        return 

    os.mkdir(brkraw_dir)
    dst_path = os.path.join(brkraw_dir, os.path.basename(input_file))

    shutil.copyfile(input_file, dst_path)

    data = nii.load(input_file)
    raw_img = data.dataobj.get_unscaled()

    raw_nii = nii.Nifti1Image(raw_img, data.affine)
    nii.save(raw_nii, input_file)

    delete_orient_command = f"fslorient -deleteorient {input_file}"
    subprocess.run(delete_orient_command, shell=True)

    # Befehl zum Festlegen der radiologischen Orientierung
    forceradiological_command = f"fslorient -forceradiological {input_file}"
    subprocess.run(forceradiological_command, shell=True)

def applyBET(input_file: str, frac: float, radius: int, output_path: str) -> str:
    """
    Performs brain extraction via the FSL Brain Extraction Tool (BET). Requires an appropriate input file (input_file), the fractional intensity threshold (frac), the head radius (radius) and the output path (output_path).
    """
    # scale Nifti data by factor 10
    data = nii.load(input_file)
    imgTemp = data.get_fdata()
    scale = np.eye(4)* 10
    scale[3][3] = 1
    imgTemp = np.flip(imgTemp, 2)

    scaledNiiData = nii.Nifti1Image(imgTemp, data.affine * scale)
    hdrIn = scaledNiiData.header
    hdrIn.set_xyzt_units('mm')
    scaledNiiData = nii.as_closest_canonical(scaledNiiData)

    fsl_path = os.path.join(os.path.dirname(input_file),'fslScaleTemp.nii.gz')
    nii.save(scaledNiiData, fsl_path)

    # extract brain
    output_file = os.path.join(output_path, os.path.basename(input_file).split('.')[0] + 'Bet.nii.gz')
    myBet = fsl.BET(in_file=fsl_path, out_file=output_file,frac=frac,radius=radius,robust=True, mask = True)
    myBet.run()
    os.remove(fsl_path)


    # unscale result data by factor 10ˆ(-1)
    dataOut = nii.load(output_file)
    imgOut = dataOut.get_fdata()
    scale = np.eye(4)/ 10
    scale[3][3] = 1

    unscaledNiiData = nii.Nifti1Image(imgOut, dataOut.affine * scale)
    hdrOut = unscaledNiiData.header
    hdrOut.set_xyzt_units('mm')
    nii.save(unscaledNiiData, output_file)
    return output_file

def smoothIMG(input_file, output_path):
    """
    Smoothes image via FSL. Only input and output has do be specified. Parameters are fixed to box shape and to the kernel size of 0.1 voxel.
    """
    data = nii.load(input_file)
    vol = data.get_fdata()
    ImgSmooth = np.min(vol, 3)

    unscaledNiiData = nii.Nifti1Image(ImgSmooth, data.affine)
    hdrOut = unscaledNiiData.header
    hdrOut.set_xyzt_units('mm')
    output_file = os.path.join(os.path.dirname(input_file),
                               os.path.basename(input_file).split('.')[0] + 'DN.nii.gz')
    nii.save(unscaledNiiData, output_file)
    input_file = output_file
    output_file = os.path.join(output_path, os.path.basename(input_file).split('.')[0] + 'Smooth.nii.gz')
    myGauss =  fsl.SpatialFilter(
        in_file = input_file,
        out_file = output_file, 
        operation = 'median',
        kernel_shape = 'box',
        kernel_size = 0.1
    )
    myGauss.run()
    return output_file

def thresh(input_file, output_path):
    #output_file = os.path.join(os.path.dirname(input_file),os.path.basename(input_file).split('.')[0]+ 'Thres.nii.gz')
    output_file = os.path.join(output_path, os.path.basename(input_file).split('.')[0] + 'Thres.nii.gz')
    myThres = fsl.Threshold(in_file=input_file,out_file=output_file,thresh=20)#,direction='above')
    myThres.run()
    return output_file

def cropToSmall(input_file,output_path):
    #output_file = os.path.join(os.path.dirname(input_file),os.path.basename(input_file).split('.')[0]  + 'Crop.nii.gz')
    output_file = os.path.join(output_path, os.path.basename(input_file).split('.')[0] + 'Crop.nii.gz')
    myCrop = fsl.ExtractROI(in_file=input_file,roi_file=output_file,x_min=40,x_size=130,y_min=50,y_size=110,z_min=0,z_size=12)
    myCrop.run()
    return  output_file


if __name__ == "__main__":
    import argparse


    parser = argparse.ArgumentParser(description='Preprocessing of DTI Data')

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-i', '--input', help='Path to the raw NIfTI DTI file', required=True)

    parser.add_argument('-f', '--frac', help='Fractional intensity threshold - default=0.26, smaller values give larger brain outline estimates', nargs='?', type=float,default=0.26)
    parser.add_argument('-r', '--radius', help='Head radius (mm not voxels) - default=55', nargs='?', type=int ,default=55)
    parser.add_argument('-g', '--vertical_gradient', help='Vertical gradient in fractional intensity threshold - default=0.07, positive values give larger brain outlines at bottom and smaller brain outlines at top', nargs='?',
                        type=float,default=0.07)
    args = parser.parse_args()

    # set parameters
    input_file = None
    if args.input is not None and args.input is not None:
        input_file = args.input

    if not os.path.exists(input_file):
        sys.exit("Error: '%s' is not an existing directory or file %s is not in directory." % (input_file, args.file,))
        
    frac = args.frac
    radius = args.radius
    vertical_gradient = args.vertical_gradient
    output_path = os.path.dirname(input_file)

    print(f"Frac: {frac} Radius: {radius} Gradient {vertical_gradient}")

    reset_orientation(input_file)
    print("Orientation resetted to RAS")

    try:
        output_smooth = smoothIMG(input_file = input_file, output_path = output_path)
        print("Smoothing completed")
    except Exception as e:
        print(f'Fehler in der Biasfieldcorrecttion\nFehlermeldung: {str(e)}')
        raise

    # intensity correction using non parametric bias field correction algorithm
    try:
        output_mico = applyMICO.run_MICO(output_smooth,output_path)
        print("Biasfieldcorrecttion was successful")
    except Exception as e:
        print(f'Fehler in der Biasfieldcorrecttion\nFehlermeldung: {str(e)}')
        raise

    # get rid of your skull         
    outputBET = applyBET(input_file = output_mico, frac = frac, radius = radius, output_path = output_path)
    print("Brainextraction was successful")



















