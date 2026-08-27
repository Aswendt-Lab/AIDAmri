"""
Created on 10/08/2017

@author: Niklas Pallast
Neuroimaging & Neuroengineering
Department of Neurology
University Hospital Cologne

Edited by Paul Camacho 2025

"""


import nipype.interfaces.fsl as fsl
import os, sys
import nibabel as nib
import numpy as np
import applyMICO
import nipype.interfaces.ants as ants
import subprocess
import shutil
import averageb0
import dipy.denoise.patch2self as patch2self
import itertools
import sys
import threading
import time

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))
from common.bet import applyBET, skip_bet_function
from common.artifact_manifest import start_output_tracking
from common.script_logging import setup_script_logging

FATAL_LIP_HEADER_EXIT_CODE = 86

def creat_brkraw_backup(input_file):

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
    space_unit, time_unit = hdr.get_xyzt_units()

    if not space_unit or space_unit.lower() == "unknown":
        space_unit = "mm"
    if not time_unit or time_unit.lower() == "unknown":
        time_unit = "sec"

    hdr.set_xyzt_units(space_unit, time_unit)

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

    nib.save(out, input_file)
    return input_file

def set_default_xyzt_units_if_unknown(target):
    """
       Set NIfTI space/time units to mm/sec only if they are unknown.

       Accepts either:
       - a file path to a NIfTI file, or
       - a nibabel header object
       """
    if isinstance(target, str):
        img = nib.load(target)
        hdr = img.header

        space_unit, time_unit = hdr.get_xyzt_units()
        if not space_unit or space_unit.lower() == "unknown":
            space_unit = "mm"
        if not time_unit or time_unit.lower() == "unknown":
            time_unit = "sec"

        hdr.set_xyzt_units(space_unit, time_unit)
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

def n4biasfieldcorr(input_file):
    output_file = os.path.join(os.path.dirname(input_file), os.path.basename(input_file).split('.')[0] + 'AntsBias.nii.gz')
    # Note: shrink_factor is set to 4 to speed up the process, but can be adjusted
    myAnts = ants.N4BiasFieldCorrection(input_image=input_file,output_image=output_file,
                                        shrink_factor=4,bspline_fitting_distance=20,
                                        bspline_order=3,n_iterations=[1000,0],dimension=3)
    myAnts.run()

    img = nib.load(output_file)
    hdr = img.header
    hdr["pixdim"][4:8] = 1
    nib.save(img, output_file)

    print("Biasfield correction completed")
    return output_file

def spinner(stop_event, message="Working"):
    """
    Displays a simple terminal spinner while a long-running processing step is active.
    Does not report the actual progress of external tools such as FSL BET or ANTs.
    """
    if os.environ.get("AIDAMRI_DISABLE_SPINNER") == "1" or not sys.stdout.isatty():
        return

    for ch in itertools.cycle("|/-\\"):
        if stop_event.is_set():
            break
        sys.stdout.write(f"\r{message}... {ch}")
        sys.stdout.flush()
        time.sleep(0.1)
    sys.stdout.write(f"\r{message}... done\n")
    sys.stdout.flush()


def set_xform_codes_to_one(input_file):
    img = nib.load(input_file)
    img.set_qform(img.affine, code=1)
    img.set_sform(img.affine, code=1)
    nib.save(img, input_file)
    return input_file

def denoise_patch2self(input_file, output_path, b0_thresh=100):
    """
    Denoises the input DTI image using Patch2Self from DIPY.
    Requires an appropriate input file (input_file) and the output path (output_path).
    """
    bvalsname = input_file.replace(".nii.gz", ".bval")
    if not os.path.exists(bvalsname):
        try:
            bvalsname = input_file.replace(".nii.gz", ".btable")
            btable = np.loadtxt(bvalsname, dtype=float)
            bvalsname = os.path.splitext(bvalsname)[0] + ".bval"
            np.savetxt(bvalsname, btable[0, :], fmt='%.6f')
        except:
            sys.exit(f"Error: bvals file {bvalsname} not found.")
    bvals = np.loadtxt(bvalsname, dtype=float)
    data = nib.load(input_file)
    img = data.get_fdata()
    affine = data.affine
    debug = False
    if debug:
        print("Debugging information:")
        print("Image header:", data.header)
        print("Affine matrix:", affine)
        print("Image sform:", data.header.get_sform())
    if img.ndim != 4:
        raise ValueError("Input image must be a 4D NIfTI file.")
    
    # Apply Patch2Self denoising
    denoised_img = patch2self.patch2self(img, bvals, b0_threshold=b0_thresh, model='ols', out_dtype=np.float32)
    
    # Save the denoised image
    output_file = os.path.join(output_path, os.path.basename(input_file).split('.')[0] + 'Patch2SelfDenoised.nii.gz')
    denoised_nii = nib.Nifti1Image(denoised_img, affine)
    if debug:
        print("Denoised image header:", denoised_nii.header)
        print("Denoised affine matrix:", denoised_nii.affine)
        print("Denoised image sform:", denoised_nii.header.get_sform())
    nib.save(denoised_nii, output_file)

    # # Copy header from original image to denoised image using fslcpgeom
    # myFslCpGeom = fsl.utils.CopyGeom(dest_file=output_file, in_file=input_file)
    # myFslCpGeom.run()
    # print(f"Denoising completed, output saved to {output_file}")
    # if debug is True:
    #     output = nii.load(output_file)
    #     print("Final denoised image header after copying geometry:", output.header)
    #     print("Final denoised affine matrix after copying geometry:", output.affine)
    #     print("Final denoised image sform after copying geometry:", output.header.get_sform())
    return output_file

def smoothIMG(input_file, output_path, skip_smoothing=False):
    """
    Prepare a 3D reference image and optionally apply FSL's median spatial filter.
    For 4D inputs, a voxel-wise median projection across the 4th dimension is
    written as *MP.nii.gz. For 3D inputs, the MP image is just a
    float32/header-normalized copy.
    """
    source_base = os.path.basename(input_file).split('.')[0]
    data = nib.load(input_file)
    vol = data.get_fdata()
    if vol.ndim == 4:
        img_smooth = np.median(vol, axis=3).astype(np.float32)
        source_base = source_base + 'MP'
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
                               os.path.basename(input_file).split('.')[0] + 'MP.nii.gz')
    nib.save(unscaledNiiData, output_file)
    input_file = output_file

    if skip_smoothing:
        print("Spatial smoothing skipped")
        return input_file

    output_file = os.path.join(output_path, source_base + 'Smooth.nii.gz')
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


    parser = argparse.ArgumentParser(description='Preprocessing of DTI Data')

    requiredNamed = parser.add_argument_group('Required named arguments')
    requiredNamed.add_argument(
        '-i',
        '--input-file',
        help='Path to the raw NIfTI DTI file',
        required=True,
    )

    parser.add_argument(
        '-f',
        '--frac',
        help='Fractional intensity threshold - default=0.4, smaller values give larger brain outline estimates',
        type=float,
        default=0.4,
    )
    parser.add_argument(
        '-r',
        '--radius',
        help='Head radius (mm not voxels) - default=45',
        type=int,
        default=45,
    )
    parser.add_argument(
        '-g',
        '--horizontal-gradient',
        help='Horizontal gradient in fractional intensity threshold - default=0.0. Not for bet4animals! Higher positive values make the BET stricter posterior and less stricter anterior (snout)',
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
        '-b',
        '--bias-method',
        help='Biasfield correction method - default=skip, other options are "mico", "ants"',
        choices = ["skip", "mico", "ants"],
        type=str.lower,
        default="skip",
    )

    parser.add_argument(
        '--bet',
        choices=["skip", "bet", "bet4animal"],
        type=str.lower,
        default="bet",
        help='Brain extraction method for DTI: skip, bet or bet4animal. Default: bet'
    )

    parser.add_argument(
        '-d',
        '--denoiser',
        help='Denoising method - default=None, other option is "patch2self"',
        choices = ["patch2self"],
        type=str.lower,
        default=None
    )

    parser.add_argument(
        '--average-b0',
        help='Average the b0 volumes',
        action='store_true'
    )
    parser.add_argument(
        '--skip-smoothing',
        action='store_true',
        help='Skip the FSL spatial median smoothing step; still creates the 3D median reference image'
    )
    args = parser.parse_args()

    # set Parameters
    input_file = args.input_file
    if not os.path.exists(input_file):
        sys.exit(f"Error: input file does not exist: {input_file}")

    frac = args.frac
    radius = args.radius
    horizontal_gradient = args.horizontal_gradient
    bias_method = args.bias_method
    output_path = os.path.dirname(input_file)
    start_output_tracking(output_path, "dwi", "preprocessing")
    setup_script_logging(output_path, "preprocess.log")
    b0_thresh=100

    if args.bet == "bet":
        print(f"Frac: {frac} Radius: {radius} Gradient {horizontal_gradient}")

    creat_brkraw_backup(input_file)
    header_check(input_file)
    
    if args.denoiser == "patch2self":
        # Denoising using Patch2Self
        print("Starting denoising using patch2self")
        try:
            # start spinner
            stop_event = threading.Event()
            thread = threading.Thread(
                target=spinner,
                args=(stop_event, "Running Denoising with Patch2Self")
            )
            thread.start()

            try:
                denoised_image = denoise_patch2self(input_file, output_path, b0_thresh)
                set_xform_codes_to_one(denoised_image)
                set_default_xyzt_units_if_unknown(denoised_image)
            finally:
                stop_event.set()
                thread.join()

            print("Denoising completed, output saved to", denoised_image)
        except Exception as e:
            print(f'Error in Patch2Self denoising\nError message: {str(e)}')
            raise


        input_file = denoised_image

    if args.average_b0:
        # Average b0 volumes
        print("Starting averaging b0 volumes")
        try:
            # start spinner
            stop_event = threading.Event()
            thread = threading.Thread(
                target=spinner,
                args=(stop_event,)
            )
            thread.start()

            try:
                b0image = averageb0.averageb0(input_file,b0_thresh)
                set_xform_codes_to_one(b0image)
                set_default_xyzt_units_if_unknown(b0image)
            finally:
                stop_event.set()
                thread.join()
        # # Copy header with fslcopygeom
        # myFslCpGeom = fsl.utils.CopyGeom(dest_file=b0image, in_file=input_file)
        # myFslCpGeom.run()
            input_file = b0image
            print("Averaging b0 volumes completed, output saved to", input_file)
        except Exception as e:
            print(f'Error in averaging b0 volumes\nError message: {str(e)}')
            raise


    try:
        # start spinner
        stop_event = threading.Event()
        thread = threading.Thread(
            target=spinner,
            args=(stop_event, "Running smoothing")
        )
        thread.start()

        try:
            output_smooth = smoothIMG(input_file = input_file, output_path = output_path, skip_smoothing=args.skip_smoothing)
        finally:
            stop_event.set()
            thread.join()
        print(f"Smoothing completed, output saved to {output_path}")
    except Exception as e:
        print(f'Error in smoothing\nError message: {str(e)}')
        raise

    # intensity correction using non parametric bias field correction algorithm
    if bias_method == "skip":
        print("No bias field correction applied")
        outputBiasCorr = output_smooth
    elif bias_method == "mico":
        print("Starting Biasfieldcorrection with MICO:")
        try:
            outputBiasCorr = applyMICO.run_MICO(output_smooth, output_path)
            set_xform_codes_to_one(outputBiasCorr)
            set_default_xyzt_units_if_unknown(outputBiasCorr)
            print("Biasfield correction was successful")
        except Exception as e:
            print(f'Error in bias field correction\nError message: {str(e)}')
            raise
    elif bias_method == "ants":
        # intensity correction using ANTs N4BiasFieldCorrection
        print("Starting Biasfieldcorrection with ANTS:")
        try:
            #start spinner
            stop_event = threading.Event()
            thread = threading.Thread(
                target=spinner,
                args=(stop_event, "Running N4 ANTS bias correction")
            )
            thread.start()

            try:
                outputBiasCorr = n4biasfieldcorr(input_file=output_smooth)
            finally:
                stop_event.set()
                thread.join()
            print("Biasfield correction was successful")
        except Exception as e:
            print(f'Error in bias field correction\nError message: {str(e)}')
            raise
    #print(os.path.exists(outputBiasCorr))

    if args.bet == "skip":
        print("Skipping brain extraction.")
        outputBET = skip_bet_function(outputBiasCorr)
    else:
        # brain extraction
        print("Starting brain extraction")
        try:
            stop_event = threading.Event()
            thread = threading.Thread(
                target=spinner,
                args=(stop_event, "Running Brain extraction")
            )
            thread.start()
            try:
                outputBET = applyBET(
                    input_file=outputBiasCorr,
                    frac=frac,
                    radius=radius,
                    horizontal_gradient=horizontal_gradient,
                    use_bet4animal=args.bet == "bet4animal",
                    center=args.center)
            finally:
                stop_event.set()
                thread.join()
            print("Brain extraction was successful")
        except Exception as e:
            print(f'Error in brain extraction\nError messsage: {str(e)}')
            raise

    print("Preprocessing completed")
