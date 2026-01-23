"""
Created on 10/08/2017

@author: Niklas Pallast
Neuroimaging & Neuroengineering
Department of Neurology
University Hospital Cologne

"""


import nipype.interfaces.fsl as fsl
import os,sys
import nibabel as nib
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

    data = nib.load(input_file)
    raw_img = data.dataobj.get_unscaled()

    raw_nii = nib.Nifti1Image(raw_img, data.affine)
    nib.save(raw_nii, input_file)

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

def copy_xform(ref_file, dst_file):
    ref = nib.load(ref_file)
    dst = nib.load(dst_file)

    data = dst.get_fdata(dtype=np.float32)

    new = nib.Nifti1Image(data, ref.affine, header=dst.header)

    # Copy qform from reference
    qaff, qcode = ref.get_qform(coded=True)
    if qaff is not None:
        new.set_qform(qaff, int(qcode))

    # Set sform explicitly → Code = 2 (aligned anatomical)
    saff, _ = ref.get_sform(coded=True)
    if saff is None:
        saff = ref.affine

    new.set_sform(saff, code=2)

    nib.save(new, dst_file)

def estimate_center_intensity_based(nifti, percentile=60):
    """
    Estimate BET center (-c) using intensity-weighted center-of-gravity (fslstats -C),
    but excluding low-intensity voxels using a data-adaptive threshold (-l = P{percentile}).
    """
    # 1) Get intensity percentile
    p = subprocess.check_output(
        ["fslstats", nifti, "-P", str(percentile)]
    ).decode().strip()

    # 2) Compute center-of-gravity using only voxels > P{percentile}
    center = subprocess.check_output(
        ["fslstats", nifti, "-l", p, "-C"]
    ).decode().strip().split()

    cx, cy, cz = [int(round(float(v))) for v in center]
    return [cx, cy, cz], float(p)


def applyBET(input_file,frac,radius,vertical_gradient,use_bet4animal=False, species='mouse', center= None):
    """Apply BET"""
    if use_bet4animal == True:
        # Use BET for animal brains
        print("Using BET for animal brains")
        print("Note: bet4animal requires that the AC-PC line of brain is parallel to Y-axis")
        w_value = 2 #smooth the surface (lissencephalic weighting)
        species_id = 6 if species == 'mouse' else 5
        output_file = os.path.join(os.path.dirname(input_file), os.path.basename(input_file).split('.')[0] + 'Bet.nii.gz')
        if center is None:
            center, p = estimate_center_intensity_based(input_file)
        cx, cy, cz = center

        command = f"/aida/bin/bet4animal {input_file} {output_file} -f {frac} -m -w {w_value} -z {species_id} -c {cx} {cy} {cz}"  # m = binary mask output
        subprocess.run(command, shell=True, check=True)
    else:
        data = nib.load(input_file)
        imgTemp = data.get_fdata()
        # create 4x4 scaling matrix
        scale = np.eye(4) * 10
        #Set last element to 1 (important for affine matrix)
        scale[3][3] = 1
        imgTemp = np.flip(imgTemp, 2)

        #Create new Nifti image with flipped data and scaled affine
        scaledNiiData = nib.Nifti1Image(imgTemp, data.affine * scale)
        hdrIn = scaledNiiData.header
        hdrIn.set_xyzt_units('mm')
        scaledNiiData = nib.as_closest_canonical(scaledNiiData)

        fslPath = os.path.join(os.path.dirname(input_file), 'fslScaleTemp.nii.gz')
        nib.save(scaledNiiData, fslPath)

        # set output file path
        output_file = os.path.join(os.path.dirname(input_file), os.path.basename(input_file).split('.')[0] + 'Bet.nii.gz')
        # run BET with robust option and mask output
        myBet = fsl.BET(in_file=fslPath, out_file=output_file, frac=frac, radius=radius,
                        vertical_gradient=vertical_gradient, robust=True, mask=True)
        myBet.run()
        os.remove(fslPath) # remove temporary scaled file

        # unscale result data by factor 10ˆ(-1)
        dataOut = nib.load(output_file)
        imgOut = dataOut.get_fdata()
        scale = np.eye(4) / 10
        scale[3][3] = 1
        #create unscaled Nifti image with unscaled affine and flip
        unscaledNiiData = nib.Nifti1Image(imgOut, dataOut.affine * scale)
        hdrOut = unscaledNiiData.header
        hdrOut.set_xyzt_units('mm')
        nib.save(unscaledNiiData, output_file)
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
        '-c', '--center',
        help='Brain center in voxel coordinates: x y z',
        nargs=3,
        type=int,
        default=None
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
        '--use_bet4animal',
        help='Use BET for animal brains',
        action='store_true'
    )

    args = parser.parse_args()

    # set Parameters
    input_file = None
    if args.input_file is not None and args.input_file is not None:
        input_file = args.input_file
    if not os.path.exists(input_file):
        sys.exit(f"Error: input file does not exist: {input_file}")


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
                copy_xform(input_file, outputBiasCorr)
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
        outputBET = applyBET(input_file=outputBiasCorr,frac=frac,radius=radius,vertical_gradient=vertical_gradient,use_bet4animal=use_bet4animal, center=args.center)
        print("Brain extraction was successful")
    except Exception as e:
        print(f'Error in brain extraction\nFehlermeldung: {str(e)}')
        raise
    
    print("Preprocessing completed")
 










