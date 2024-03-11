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
import subprocess
import shutil


def define_rodent_spezies():
    global rodent
    rodent = int(input("Select rodent: Mouse = 0 , Rat = 1 "))
    if rodent == 0 or rodent == 1:
        return rodent
    else:
        print("Invalid option. Enter 0 for mouse or 1 for rat.")
        return define_rodent_spezies()
        
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



def applyBET(input_file,frac,radius,vertical_gradient):
    """Apply BET"""
    # scale Nifti data by factor 10
    data = nii.load(input_file)
    imgTemp = data.get_fdata()
    scale = np.eye(4)* 10
    scale[3][3] = 1

    # this has to be adapted in the case the output image is not RAS orientated - Siding from feet to nose
    imgTemp = np.flip(imgTemp,2)
    imgTemp = np.flip(imgTemp,1)
    imgTemp = np.flip(imgTemp,0)
    #imgTemp = np.rot90(imgTemp, 2)

    scaledNiiData = nii.Nifti1Image(imgTemp, data.affine * scale)
    hdrIn = scaledNiiData.header
    hdrIn.set_xyzt_units('mm')
    scaledNiiData = nii.as_closest_canonical(scaledNiiData)

    fslPath = os.path.join(os.path.dirname(input_file),'fslScaleTemp.nii.gz')
    nii.save(scaledNiiData, fslPath)

    # extract brain
    output_file = os.path.join(os.path.dirname(input_file),os.path.basename(input_file).split('.')[0] + 'Bet.nii.gz')

    myBet = fsl.BET(in_file=fslPath, out_file=output_file,frac=frac,radius=radius,
                    vertical_gradient=vertical_gradient,robust=True, mask = True)
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
    return output_file

#%% Program

#specify default parameters by defining rodent spezies
define_rodent_spezies()

if rodent == 0:
    default_frac = 0.15
    default_rad  = 45
    default_vert = 0.0
elif rodent == 1:
    default_frac = 0.26
    default_rad  = 55
    default_vert = 0.07

if __name__ == "__main__":
    import argparse
    
    
    parser = argparse.ArgumentParser(description='Preprocessing of T2 Data')

    requiredNamed = parser.add_argument_group('Required named arguments')
    requiredNamed.add_argument('-i','--input_file', help='path to input file',required=True)

    parser.add_argument(
        '-f',
        '--frac',
        help='Fractional intensity threshold - default: Mouse=0.15 , Rat=0.26  smaller values give larger brain outline estimates',
        nargs='?',
        type=float,
        default=default_frac,
        )
    parser.add_argument(
        '-r', 
        '--radius',
        help='Head radius (mm not voxels) - default: Mouse=45 , Rat=55',
        nargs='?',
        type=int,
        default=default_rad,
        )
    parser.add_argument(
        '-g',
        '--vertical_gradient',
        help='Vertical gradient in fractional intensity threshold - default: Mouse=0.0 , Rat=0.07   positive values give larger brain outlines at bottom and smaller brain outlines at top',
        nargs='?',
        type=float,
        default=default_vert,
        )
    parser.add_argument(
        '-b',
        '--bias_skip',
        help='Set value to 1 to skip bias field correction',
        nargs='?',
        type=float, 
        default=0.0,
        )

    args = parser.parse_args()

    # set Parameters
    input_file = None
    if args.input_file is not None and args.input_file is not None:
        input_file = args.input_file
    if not os.path.exists(input_file):
        sys.exit("Error: '%s' is not an existing directory or file %s is not in directory." % (input_file, args.file,))

    frac = args.frac
    radius = args.radius
    vertical_gradient = args.vertical_gradient
    bias_skip = args.bias_skip

    print(f"Frac: {frac} Radius: {radius} Gradient {vertical_gradient}")

    reset_orientation(input_file)
    print("Orientation resetted to RAS")

    #intensity correction using non parametric bias field correction algorithm
    print("Starting Biasfieldcorrection:")
    if bias_skip == 0:
        try:
            outputMICO = applyMICO.run_MICO(input_file,os.path.dirname(input_file))
            print("Biasfieldcorrecttion was successful")
        except Exception as e:
            print(f'Fehler in der Biasfieldcorrecttion\nFehlermeldung: {str(e)}')
            raise
    else:
        outputMICO = input_file
        
    # brain extraction
    print("Starting brain extraction")
    try:
        outputBET = applyBET(input_file=outputMICO,frac=frac,radius=radius,vertical_gradient=vertical_gradient)
        print("Brain extraction was successful")
    except Exception as e:
        print(f'Error in brain extraction\nFehlermeldung: {str(e)}')
        raise
    
    print("Preprocessing completed")
 










