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
import nipype.interfaces.ants as ants

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


def n4biasfieldcorr(input_file):
    output_file = os.path.join(os.path.dirname(input_file), os.path.basename(input_file).split('.')[0] + 'Bias.nii.gz')
    # Note: shrink_factor is set to 4 to speed up the process, but can be adjusted
    myAnts = ants.N4BiasFieldCorrection(input_image=input_file, output_image=output_file,
                                        shrink_factor=2, bspline_fitting_distance=20,
                                        bspline_order=3, n_iterations=[50, 50, 50, 50, 0], dimension=3)
    myAnts.run()
    print("Biasfield correction completed")
    return output_file


def applyBET(input_file,frac,radius,vertical_gradient,use_bet4animal=False, species='mouse'):
    """Apply BET"""
    if use_bet4animal:
        # Use BET for animal brains
        print("Using BET for animal brains")
        species_id = 6 if species == 'mouse' else 5
        output_file = os.path.join(os.path.dirname(input_file), os.path.basename(input_file).split('.')[0] + 'Bet.nii.gz')
        command = f"/aida/bin/bet4animal {input_file} {output_file} -f {frac} -m -w {w_value} -z {species_id}"
        subprocess.run(command)
    else:
        # scale Nifti data by factor 10
        data = nii.load(input_file)
        imgTemp = data.get_fdata()
        scale = np.eye(4) * 10
        scale[3][3] = 1
        imgTemp = np.flip(imgTemp, 2)

        scaledNiiData = nii.Nifti1Image(imgTemp, data.affine * scale)
        hdrIn = scaledNiiData.header
        hdrIn.set_xyzt_units('mm')
        scaledNiiData = nii.as_closest_canonical(scaledNiiData)

        fslPath = os.path.join(os.path.dirname(input_file), 'fslScaleTemp.nii.gz')
        nii.save(scaledNiiData, fslPath)

        # extract brain
        output_file = os.path.join(os.path.dirname(input_file), os.path.basename(input_file).split('.')[0] + 'Bet.nii.gz')

        myBet = fsl.BET(in_file=fslPath, out_file=output_file, frac=frac, radius=radius,
                        vertical_gradient=vertical_gradient, robust=True, mask=True)
        myBet.run()
        os.remove(fslPath)

        # unscale result data by factor 10ˆ(-1)
        dataOut = nii.load(output_file)
        imgOut = dataOut.get_fdata()
        scale = np.eye(4) / 10
        scale[3][3] = 1

        unscaledNiiData = nii.Nifti1Image(imgOut, dataOut.affine * scale)
        hdrOut = unscaledNiiData.header
        hdrOut.set_xyzt_units('mm')
        nii.save(unscaledNiiData, output_file)
    print(f"Brain extraction completed, output saved to {output_file}")
    return output_file

#%% Program

if __name__ == "__main__":
    import argparse


    parser = argparse.ArgumentParser(description='Preprocessing of T2 Data')

    requiredNamed = parser.add_argument_group('Required named arguments')
    requiredNamed.add_argument('-i','--input_file', help='path to input file',required=True)

    parser.add_argument(
        '-f',
        '--frac',
        help='Fractional intensity threshold - default=0.15  smaller values give larger brain outline estimates',
        nargs='?',
        type=float,
        default=0.15,
        )
    parser.add_argument(
        '-r', 
        '--radius',
        help='Head radius (mm not voxels) - default=45',
        nargs='?',
        type=int,
        default=45,
        )
    parser.add_argument(
        '-g',
        '--vertical_gradient',
        help='Vertical gradient in fractional intensity threshold - default=0.0   positive values give larger brain outlines at bottom and smaller brain outlines at top',
        nargs='?',
        type=float,
        default=0.0,
        )
    parser.add_argument(
        '-b',
        '--bias_skip',
        help='Set value to 1 to skip bias field correction',
        nargs='?',
        type=float, 
        default=0.0,
        )
    parser.add_argument(
        '-bias_method',
        '--bias_method',
        help='Biasfield correction method - default="mico", other options are "mico" or "ants"',
        nargs='?',
        type=str,
        default="mico",
        )
    parser.add_argument(
        '-bet4animal',
        '--use_bet4animal',
        help='Set value to 1 to use BET for animal brains',
        nargs='?',
        type=bool,
        default=False,
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
    bias_method = args.bias_method

    print(f"Frac: {frac} Radius: {radius} Gradient {vertical_gradient}")

    reset_orientation(input_file)
    print("Orientation resetted to RAS")

    #intensity correction using non parametric bias field correction algorithm
    print("Starting Biasfieldcorrection:")
    if bias_skip == 0:
        if bias_method == "mico":
            try:
                outputBiasCorr = applyMICO.run_MICO(input_file, os.path.dirname(input_file))
                print("Biasfield correction was successful")
            except Exception as e:
                print(f'Error in bias field correction\nError message: {str(e)}')
                raise
        elif bias_method == "ants":
            try:
                outputBiasCorr = n4biasfieldcorr(input_file=input_file)
                print("Biasfield correction was successful")
            except Exception as e:
                print(f'Error in bias field correction\nError message: {str(e)}')
                raise
    else:
        outputBiasCorr = input_file
    
    print(os.path.exists(outputBiasCorr))

    use_bet4animal = args.use_bet4animal

    # brain extraction
    print("Starting brain extraction")
    try:
        outputBET = applyBET(input_file=outputBiasCorr,frac=frac,radius=radius,vertical_gradient=vertical_gradient,use_bet4animal=use_bet4animal)
        print("Brain extraction was successful")
    except Exception as e:
        print(f'Error in brain extraction\nFehlermeldung: {str(e)}')
        raise
    
    print("Preprocessing completed")
 










