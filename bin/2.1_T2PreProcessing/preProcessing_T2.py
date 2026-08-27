"""
Created on 10/08/2017

@author: Niklas Pallast
Neuroimaging & Neuroengineering
Department of Neurology
University Hospital Cologne

"""


import nipype.interfaces.fsl as fsl
import os, sys
import nibabel as nib
import numpy as np
import applyMICO
import nipype.interfaces.ants as ants
import subprocess
import shutil
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
    myAnts = ants.N4BiasFieldCorrection(input_image=input_file, output_image=output_file,
                                        shrink_factor=2, bspline_fitting_distance=20,
                                        bspline_order=3, n_iterations=[50, 50, 50, 50, 0], dimension=3)
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

#%% Program

if __name__ == "__main__":
    import argparse


    parser = argparse.ArgumentParser(description='Preprocessing of T2 Data')

    requiredNamed = parser.add_argument_group('Required named arguments')
    requiredNamed.add_argument('-i','--input-file', help='path to input file',required=True)

    parser.add_argument(
        '-f',
        '--frac',
        help='Fractional intensity threshold - default=0.15  smaller values give larger brain outline estimates',
        type=float,
        default=0.15,
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
        help='Biasfield correction method - default="mico", other options are "ants" or "skip"',
        choices = ["skip", "mico", "ants"],
        type=str.lower,
        default="mico",
    )

    parser.add_argument(
        '--bet',
        choices=["skip", "bet", "bet4animal"],
        type=str.lower,
        default="bet",
        help='Brain extraction method for T2: skip, bet or bet4animal. Default: bet'
    )

    args = parser.parse_args()

    # set Parameters
    input_file = args.input_file
    if not os.path.exists(input_file):
        sys.exit(f"Error: input file does not exist: {input_file}")
    start_output_tracking(os.path.dirname(input_file), "anat", "preprocessing")
    setup_script_logging(os.path.dirname(input_file), "preprocess.log")

    frac = args.frac
    radius = args.radius
    horizontal_gradient = args.horizontal_gradient
    bias_method = args.bias_method

    if args.bet == "bet":
        print(f"Frac: {frac} Radius: {radius} Gradient {horizontal_gradient}")

    creat_brkraw_backup(input_file)
    header_check(input_file)

    #intensity correction using non parametric bias field correction algorithm
    if bias_method == "skip":
        print("No bias field correction applied")
        outputBiasCorr = input_file
    elif bias_method == "mico":
        print("Starting Biasfieldcorrection with MICO:")
        try:
            outputBiasCorr = applyMICO.run_MICO(input_file, os.path.dirname(input_file))
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
                outputBiasCorr = n4biasfieldcorr(input_file=input_file)
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
