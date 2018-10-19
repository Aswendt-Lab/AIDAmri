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
import applyMICO

def applyBET(input_file,frac,radius,vertical_gradient):
    """Apply BET"""
    # scale Nifti data by factor 10
    data = nii.load(input_file)
    imgTemp = data.get_data()
    scale = np.eye(4)* 10
    scale[3][3] = 1
    imgTemp = np.flip(imgTemp,2)
    #imgTemp = np.flip(imgTemp, 0)
    #imgTemp = np.rot90(imgTemp, 2)

    scaledNiiData = nii.Nifti1Image(imgTemp, data.affine * scale)
    hdrIn = scaledNiiData.header
    hdrIn.set_xyzt_units('mm')
    scaledNiiData = nii.as_closest_canonical(scaledNiiData)
    print('Orientation:' + str(nii.aff2axcodes(scaledNiiData.affine)))

    fslPath = os.path.join(os.getcwd(),'fslScaleTemp.nii.gz')
    nii.save(scaledNiiData, fslPath)

    # extract brain
    output_file = os.path.join(os.path.dirname(input_file),os.path.basename(input_file).split('.')[0] + 'Bet.nii.gz')

    myBet = fsl.BET(in_file=fslPath, out_file=output_file,frac=frac,radius=radius,
                    vertical_gradient=vertical_gradient,robust=True, mask = True)
    print(myBet.cmdline)
    myBet.run()
    os.remove(fslPath)

    # unscale result data by factor 10ˆ(-1)
    dataOut = nii.load(output_file)
    imgOut = dataOut.get_data()
    scale = np.eye(4)/ 10
    scale[3][3] = 1

    unscaledNiiData = nii.Nifti1Image(imgOut, dataOut.affine * scale)
    hdrOut = unscaledNiiData.header
    hdrOut.set_xyzt_units('mm')
    nii.save(unscaledNiiData, output_file)
    return output_file


if __name__ == "__main__":
    import argparse


    parser = argparse.ArgumentParser(description='Preprocessing of T2 Data')

    requiredNamed = parser.add_argument_group('Required named arguments')
    requiredNamed.add_argument('-i','--inputFile', help='file name of data',required=True)

    parser.add_argument('-f', '--frac', help='fractional intensity threshold - default=0.15  smaller values give larger brain outline estimates', nargs='?', type=float,default=0.15)
    parser.add_argument('-r', '--radius', help='head radius (mm not voxels) - default=45', nargs='?', type=int ,default=45)
    parser.add_argument('-g', '--vertical_gradient', help='vertical gradient in fractional intensity threshold - default=0.0   positive values give larger brain outline at bottom, smaller at top', nargs='?',
                        type=float,default=0.0)
    parser.add_argument('-b', '--bias_skip',
                        help='set value to 1 to skip bias field correction',
                        nargs='?',
                        type=float, default=0.0)


    args = parser.parse_args()

    # set Parameters
    inputFile = None
    if args.inputFile is not None and args.inputFile is not None:
        inputFile = args.inputFile
    if not os.path.exists(inputFile):
        sys.exit("Error: '%s' is not an existing directory of file %s is not in directory." % (inputFile, args.file,))

    frac = args.frac
    radius = args.radius
    vertical_gradient = args.vertical_gradient
    bias_skip = args.bias_skip

    # 1) Process MRI
    print("T2 Preprocessing  \33[5m...\33[0m (wait!)", end="\r")

    # generate log - file
    sys.stdout = open(os.path.join(os.path.dirname(inputFile), 'preprocess.log'), 'w')

    # print parameters
    print("Frac: %s" % frac)
    print("Radius: %s" % radius)
    print("Gradient: %s" %vertical_gradient)

    #intensity correction using non parametric bias field correction algorithm
    if bias_skip == 0:
        outputMICO = applyMICO.run_MICO(inputFile,os.path.dirname(inputFile))
    else:
        outputMICO = inputFile
    # get rid of your skull
    outputBET = applyBET(input_file=outputMICO,frac=frac,radius=radius,vertical_gradient=vertical_gradient)

    sys.stdout = sys.__stdout__
    print('T2 Preprocessing  \033[0;30;42m COMPLETED \33[0m')










