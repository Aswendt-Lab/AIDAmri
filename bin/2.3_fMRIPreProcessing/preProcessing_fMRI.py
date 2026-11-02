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
import nipype.interfaces.ants as ants
import shutil
#makes sure to import bet.py
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))
from common.bet import applyBET, skip_bet_function
from common.artifact_manifest import start_output_tracking
from common.script_logging import setup_script_logging

FATAL_LIP_HEADER_EXIT_CODE = 86


def create_brkraw_backup(input_file):

    brkraw_dir = os.path.join(os.path.dirname(input_file), "brkraw")
    if os.path.exists(brkraw_dir):
        return 

    os.mkdir(brkraw_dir)
    dst_path = os.path.join(brkraw_dir, os.path.basename(input_file))

    # Keep the original NIfTI in brkraw and process a copy at the input path.
    shutil.move(input_file, dst_path)
    shutil.copyfile(dst_path, input_file)

    data = nib.load(input_file)
    # Preserve nibabel scaling (scl_slope/scl_inter). Using get_unscaled()
    # would write raw int16 scanner values into the working file.
    raw_img = data.get_fdata(dtype=np.float32)

    hdr = data.header.copy()
    hdr.set_data_dtype(np.float32)
    set_default_xyzt_units_if_unknown(hdr)

    raw_nii = nib.Nifti1Image(raw_img, data.affine, header=hdr)
    raw_nii.set_qform(data.affine, code=1)
    raw_nii.set_sform(data.affine, code=1)
    nib.save(raw_nii, input_file)


def header_check(input_file):
    img = nib.load(input_file)
    axcodes = nib.aff2axcodes(img.affine)

    if axcodes != ("L", "I", "P"):
        print(
            f"Fatal header check failure: expected LIP orientation, found {axcodes} in {input_file}",
            file=sys.stderr,
        )
        sys.exit(FATAL_LIP_HEADER_EXIT_CODE)

    data = img.get_fdata(dtype=np.float32)

    out = nib.Nifti1Image(data, img.affine, header=img.header.copy())
    out.set_qform(img.affine, code=1)
    out.set_sform(img.affine, code=1)

    hdr = out.header
    hdr.set_data_dtype(np.float32)
    set_default_xyzt_units_if_unknown(hdr)

    nib.save(out, input_file)
    return input_file


def set_default_xyzt_units_if_unknown(target):
    if isinstance(target, str):
        img = nib.load(target)
        hdr = img.header
        set_default_xyzt_units_if_unknown(hdr)
        nib.save(img, target)
        return target

    if hasattr(target, "get_xyzt_units") and hasattr(target, "set_xyzt_units"):
        space_unit, time_unit = target.get_xyzt_units()
        if not space_unit or space_unit.lower() == "unknown":
            space_unit = "mm"
        if not time_unit or time_unit.lower() == "unknown":
            time_unit = "sec"

        target.set_xyzt_units(space_unit, time_unit)
        return target

    raise TypeError("Expected a NIfTI file path or nibabel header object")


def set_xform_codes_to_one(input_file):
    img = nib.load(input_file)
    img.set_qform(img.affine, code=1)
    img.set_sform(img.affine, code=1)
    nib.save(img, input_file)
    return input_file


def n4biasfieldcorr(input_file):
    output_file = os.path.join(os.path.dirname(input_file), os.path.basename(input_file).split('.')[0] + 'AntsBias.nii.gz')
    myAnts = ants.N4BiasFieldCorrection(
        input_image=input_file,
        output_image=output_file,
        shrink_factor=4,
        bspline_fitting_distance=20,
        bspline_order=3,
        n_iterations=[1000, 0],
        dimension=3,
    )
    myAnts.run()

    img = nib.load(output_file)
    hdr = img.header
    hdr["pixdim"][4:8] = 1
    nib.save(img, output_file)

    print("Biasfield correction completed")
    return output_file


def smoothIMG(input_file, outputPath, skip_smoothing=False):
    """
    Prepare a 3D reference image and optionally apply FSL's median spatial filter.
    For 4D inputs, a voxel-wise median projection across the 4th dimension is
    written as *MP.nii.gz. For 3D inputs, the MP image is just a
    float32/header-normalized copy.
    """
    source_base = os.path.basename(input_file).split('.')[0]
    data = nib.load(input_file)
    vol = data.get_fdata(dtype=np.float32)
    if vol.ndim == 4:
        ImgSmooth = np.median(vol, axis=3).astype(np.float32)
        source_base = source_base + 'MP'
    elif vol.ndim == 3:
        ImgSmooth = vol.astype(np.float32)
    else:
        raise ValueError(f"Unsupported image dimensionality: {vol.ndim}")

    unscaledNiiData = nib.Nifti1Image(ImgSmooth, data.affine)
    unscaledNiiData.set_qform(data.affine, code=1)
    unscaledNiiData.set_sform(data.affine, code=1)
    hdrOut = unscaledNiiData.header
    hdrOut.set_data_dtype(np.float32)
    set_default_xyzt_units_if_unknown(hdrOut)
    output_file = os.path.join(os.path.dirname(input_file),
                               os.path.basename(input_file).split('.')[0] + 'MP.nii.gz')
    # hdrOut['sform_code'] = 1
    nib.save(unscaledNiiData, output_file)
    input_file = output_file

    if skip_smoothing:
        print("Spatial smoothing skipped")
        return input_file

    #output_file =  os.path.join(os.path.dirname(input_file),os.path.basename(input_file).split('.')[0] + 'Smooth.nii.gz')
    output_file = os.path.join(outputPath, source_base + 'Smooth.nii.gz')
    myGauss =  fsl.SpatialFilter(in_file=input_file,out_file=output_file,operation='median',kernel_shape='box',kernel_size=0.1)
    myGauss.run()
    set_xform_codes_to_one(output_file)
    set_default_xyzt_units_if_unknown(output_file)
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
    parser.add_argument(
        '-g',
        '--horizontal-gradient',
        help='Horizontal gradient in fractional intensity threshold - default=0.0. Not for bet4animals! Higher positive values make the BET stricter posterior and less stricter anterior (snout)',
        nargs='?',
        type=float,
        default=0.0,
    )
    parser.add_argument(
        '-c', '--center',
        help='Brain center in voxel coordinates: x y z',
        nargs=3,
        type=float,
        default=None
    )
    parser.add_argument(
        '--bet',
        choices=["skip", "bet", "bet4animal"],
        type=str.lower,
        default="bet",
        help='Brain extraction method for fMRI preprocessing: skip, bet or bet4animal. Default: bet'
    )
    parser.add_argument(
        '-b',
        '--bias-method',
        help='Biasfield correction method - default=None, other options are "ants" or "skip"',
        choices=["skip", "ants"],
        type=str.lower,
        default=None,
    )
    parser.add_argument(
        '--skip-smoothing',
        action='store_true',
        help='Skip the FSL spatial median smoothing step; still creates the 3D median reference image'
    )
    args = parser.parse_args()

    # set parameters
    inputFile = None
    if args.input is not None and args.input is not None:
        inputFile = args.input

    if not os.path.exists(inputFile):
        sys.exit(f"Error: input file does not exist: {inputFile}")

    frac = args.frac
    radius = args.radius
    horizontal_gradient = args.horizontal_gradient
    bias_method = args.bias_method
    outputPath = os.path.dirname(inputFile)
    start_output_tracking(outputPath, "func", "preprocessing")
    setup_script_logging(outputPath, "preprocess.log")

    if args.bet == "bet":
        print(f"Frac: {frac} Radius: {radius} Gradient {horizontal_gradient}")

    create_brkraw_backup(inputFile)
    header_check(inputFile)

    outputSmooth = smoothIMG(input_file=inputFile, outputPath=outputPath, skip_smoothing=args.skip_smoothing)

    if bias_method is None or bias_method == "skip":
        print("No bias field correction applied")
        outputBiasCorr = outputSmooth
    elif bias_method == "ants":
        print("Starting Biasfieldcorrection with ANTS:")
        try:
            outputBiasCorr = n4biasfieldcorr(input_file=outputSmooth)
            set_xform_codes_to_one(outputBiasCorr)
            set_default_xyzt_units_if_unknown(outputBiasCorr)
            print("Biasfield correction was successful")
        except Exception as e:
            print(f'Error in bias field correction\nError message: {str(e)}')
            raise

    if args.bet == "skip":
        print("Skipping brain extraction.")
        outputBET = skip_bet_function(outputBiasCorr)
    else:
        # get rid of your skull
        outputBET = applyBET(
            input_file=outputBiasCorr,
            frac=frac,
            radius=radius,
            horizontal_gradient=horizontal_gradient,
            use_bet4animal=args.bet == "bet4animal",
            center=args.center,
        )

    print("Preprocessing completed")
