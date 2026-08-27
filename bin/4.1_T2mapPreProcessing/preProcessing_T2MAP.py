"""
Created on 11/09/2023

@author: Marc Schneider
Neuroimaging & Neuroengineering
Department of Neurology
University Hospital Cologne

"""


import nipype.interfaces.fsl as fsl
import os, sys
import nibabel as nib
import numpy as np
import applyMICO
import shutil
#makes sure to import bet.py
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))
from common.bet import applyBET, skip_bet_function
from common.artifact_manifest import start_output_tracking

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

def smoothIMG(input_file, output_path):
    """
    Prepare a 3D image for smoothing and apply FSL's median spatial filter.
    For 4D inputs, a voxel-wise minimum projection across the 4th dimension is
    written as *MP.nii.gz before smoothing. For 3D inputs, the MP image is just
    a float32/header-normalized copy.
    """
    data = nib.load(input_file)
    vol = data.get_fdata()
    if vol.ndim == 4:
        img_smooth = np.min(vol, axis=3).astype(np.float32)
    elif vol.ndim == 3:
        img_smooth = vol.astype(np.float32)
    else:
        raise ValueError(f"Unsupported image dimensionality: {vol.ndim}")
    unscaledNiiData = nib.Nifti1Image(img_smooth, data.affine)
    unscaledNiiData.set_qform(data.affine, code=1)
    unscaledNiiData.set_sform(data.affine, code=1)
    hdrOut = unscaledNiiData.header
    hdrOut.set_data_dtype(np.float32)
    space_unit, time_unit = hdrOut.get_xyzt_units()

    if not space_unit or space_unit.lower() == "unknown":
        space_unit = "mm"
    if not time_unit or time_unit.lower() == "unknown":
        time_unit = "sec"
    hdrOut.set_xyzt_units(space_unit, time_unit)
    output_file = os.path.join(os.path.dirname(input_file),
                               os.path.basename(input_file).split('.')[0] + '_MP.nii.gz')
    nib.save(unscaledNiiData, output_file)
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

    img = nib.load(output_file)
    hdr = img.header
    hdr["pixdim"][4:8] = 1
    nib.save(img, output_file)

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


    parser = argparse.ArgumentParser(description='Preprocessing of T2map Data')

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-i', '--input', help='Path to the raw NIfTI T2map file', required=True)

    parser.add_argument('-f', '--frac', help='Fractional intensity threshold - default=0.3, smaller values give larger brain outline estimates', nargs='?', type=float,default=0.3)
    parser.add_argument('-r', '--radius', help='Head radius (mm not voxels) - default=45', nargs='?', type=int ,default=45)
    parser.add_argument('-g', '--horizontal-gradient', help='Horizontal gradient in fractional intensity threshold - default=0.0, positive values give larger brain outlines at bottom and smaller brain outlines at top', nargs='?',
                        type=float,default=0.0)
    parser.add_argument(
        '--bet',
        choices=["skip", "bet", "bet4animal"],
        type=str.lower,
        default="bet",
        help='Brain extraction method for T2map: skip, bet or bet4animal. Default: bet'
    )
    parser.add_argument(
        '-b',
        '--bias-method',
        choices=["skip", "mico"],
        type=str.lower,
        default="mico",
        help='Biasfield correction method for T2map: skip or mico. Default: mico'
    )
    parser.add_argument('-c', '--center', nargs=3, type=float, default=None, help='BET center as x y z')
    args = parser.parse_args()

    # set parameters
    input_file = None
    if args.input is not None and args.input is not None:
        input_file = args.input

    if not os.path.exists(input_file):
        sys.exit(f"Error: input file does not exist: {input_file}")

    frac = args.frac
    radius = args.radius
    horizontal_gradient = args.horizontal_gradient
    bias_method = args.bias_method
    output_path = os.path.dirname(input_file)
    start_output_tracking(output_path, "t2map", "preprocessing")

    if args.bet == "bet":
        print(f"Frac: {frac} Radius: {radius} Gradient {horizontal_gradient}")

    create_brkraw_backup(input_file)
    header_check(input_file)
    
    try:
        output_smooth = smoothIMG(input_file = input_file, output_path = output_path)
        print("Smoothing completed")
    except Exception as e:
        print(f'Error in smoothing\nError message: {str(e)}')
        raise

    # intensity correction using non parametric bias field correction algorithm
    if bias_method == "skip":
        print("No bias field correction applied")
        outputBiasCorr = output_smooth
    elif bias_method == "mico":
        try:
            outputBiasCorr = applyMICO.run_MICO(output_smooth,output_path)
            print("Biasfieldcorrecttion was successful")
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
            output_path=output_path,
            use_bet4animal=args.bet == "bet4animal",
            center=args.center,
        )
    print("Brainextraction was successful")
