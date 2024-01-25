"""
Created on 10/08/2017

@author: Niklas Pallast
Neuroimaging & Neuroengineering
Department of Neurology
University Hospital Cologne

"""



import nipype.interfaces.fsl as fsl
import os,sys
import nibabel as nii
import numpy as np
import nipype.interfaces.ants as ants
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

def applyBET(input_file,frac,radius,outputPath):

    # scale Nifti data by factor 10
    data = nii.load(input_file)
    imgTemp = data.get_fdata()
    scale = np.eye(4)* 10
    scale[3][3] = 1
    #imgTemp = np.rot90(imgTemp,2)
    imgTemp = np.flip(imgTemp, 2)
    #imgTemp = np.flip(imgTemp, 0)
    scaledNiiData = nii.Nifti1Image(imgTemp, data.affine * scale)
    hdrIn = scaledNiiData.header
    hdrIn.set_xyzt_units('mm')
    scaledNiiData = nii.as_closest_canonical(scaledNiiData)
    print('Orientation:' + str(nii.aff2axcodes(scaledNiiData.affine)))

    fslPath = os.path.join(os.path.dirname(input_file),'fslScaleTemp.nii.gz')
    nii.save(scaledNiiData, fslPath)

    # extract brain
    output_file = os.path.join(outputPath, os.path.basename(input_file).split('.')[0] + 'Bet.nii.gz')
    myBet = fsl.BET(in_file=fslPath, out_file=output_file,frac=frac,radius=radius,robust=True, mask = True)
    myBet.run()
    os.remove(fslPath)


    # unscale result data by factor 10ˆ(-1)
    dataOut = nii.load(output_file)
    imgOut = dataOut.get_fdata()
    scale = np.eye(4)/ 10
    scale[3][3] = 1

    unscaledNiiData = nii.Nifti1Image(imgOut, dataOut.affine * scale)
    hdrOut = unscaledNiiData.header
    hdrOut.set_xyzt_units('mm')
    nii.save(unscaledNiiData, output_file)

    print("Brain extraction completed")
    return output_file

def biasfieldcorr(input_file,outputPath):
    output_file = os.path.join(outputPath, os.path.basename(input_file).split('.')[0] + 'Bias.nii.gz')
    myAnts = ants.N4BiasFieldCorrection(input_image=input_file,output_image=output_file,shrink_factor=4,dimension=3)
    myAnts.run()
    print("Biasfield correction completed")
    return output_file

def smoothIMG(input_file,outputPath):
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
    # hdrOut['sform_code'] = 1
    nii.save(unscaledNiiData, output_file)
    input_file = output_file
    #output_file =  os.path.join(os.path.dirname(input_file),os.path.basename(input_file).split('.')[0] + 'Smooth.nii.gz')
    output_file = os.path.join(outputPath, os.path.basename(inputFile).split('.')[0] + 'Smooth.nii.gz')
    myGauss =  fsl.SpatialFilter(in_file=input_file,out_file=output_file,operation='median',kernel_shape='box',kernel_size=0.1)
    myGauss.run()
    print("Smoothing completed")
    return output_file

def thresh(input_file,outputPath):
    #output_file = os.path.join(os.path.dirname(input_file),os.path.basename(input_file).split('.')[0]+ 'Thres.nii.gz')
    output_file = os.path.join(outputPath, os.path.basename(input_file).split('.')[0] + 'Thres.nii.gz')
    myThres = fsl.Threshold(in_file=input_file,out_file=output_file,thresh=20)#,direction='above')
    myThres.run()
    print("Thresholding completed")
    return output_file

def cropToSmall(input_file,outputPath):
    #output_file = os.path.join(os.path.dirname(input_file),os.path.basename(input_file).split('.')[0]  + 'Crop.nii.gz')
    output_file = os.path.join(outputPath, os.path.basename(input_file).split('.')[0] + 'Crop.nii.gz')
    myCrop = fsl.ExtractROI(in_file=input_file,roi_file=output_file,x_min=40,x_size=130,y_min=50,y_size=110,z_min=0,z_size=12)
    myCrop.run()
    print("Cropping done")
    return  output_file


if __name__ == "__main__":
    import argparse


    parser = argparse.ArgumentParser(description='Preprocessing of rsfMRI Data')

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-i', '--input', help='Path to the RAW data of rsfMRI NIfTI file', required=True)

    parser.add_argument('-f', '--frac',
                        help='Fractional intensity threshold - default=0.3, smaller values give larger brain outline estimates',
                        nargs='?', type=float, default=0.15)
    parser.add_argument('-r', '--radius', help='Head radius (mm not voxels) - default=45', nargs='?', type=int ,default=45)
    parser.add_argument('-g', '--vertical_gradient', help='Vertical gradient in fractional intensity threshold - default=0.0, positive values give larger brain outlines at bottom and smaller brain outlines at top', nargs='?',
                        type=float,default=0.0)
    args = parser.parse_args()

    # set parameters
    inputFile = None
    if args.input is not None and args.input is not None:
        inputFile = args.input

    if not os.path.exists(inputFile):
        sys.exit("Error: '%s' is not an existing directory or file %s is not in directory." % (inputFile, args.file,))

    frac = args.frac
    radius = args.radius
    vertical_gradient = args.vertical_gradient
    outputPath = os.path.dirname(inputFile)

    print(f"Frac: {frac} Radius: {radius} Gradient {vertical_gradient}")

    reset_orientation(inputFile)
    print("Orientation resetted to RAS")

    outputSmooth = smoothIMG(input_file=inputFile,outputPath=outputPath)

    # get rid of your skull
    outputBET = applyBET(input_file=outputSmooth,frac=frac,radius=radius,outputPath=outputPath)

    print("Preprocessing completed")










