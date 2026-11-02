"""
Created on 10/08/2017

@author: Niklas Pallast
Neuroimaging & Neuroengineering
Department of Neurology
University Hospital Cologne

"""

import sys,os
import nibabel as nib
import numpy as np
import shutil
import glob
import subprocess
import shlex

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))
from common.artifact_manifest import start_output_tracking
from common.script_logging import setup_script_logging


def regABA2DTI(inputVolume,stroke_mask,refStroke_mask,T2data, splitAnno,splitAnno_rsfMRI,anno_rsfMRI,bsplineMatrix,outfile,allenBrain_anno):
    if not os.path.isfile(allenBrain_anno):
        raise FileNotFoundError(f"Original Allen Brain annotation not found: {allenBrain_anno}")
    outputT2w = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_T2w.nii.gz')
    outputAff = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + 'transMatrixAff.txt')

    # NiftyReg directly registers the floating T2 image to the DTI (BET)
    # reference grid and writes a NiftyReg affine
    command = f"reg_aladin -ref {inputVolume} -flo {T2data} -res {outputT2w} -rigOnly -aff {outputAff}"
    command_args = shlex.split(command)
    try:
        result = subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,text=True)
        print(f"Output of {command}:\n{result.stdout}")
    except Exception as e:
        print(f'Error while executing the command: {command_args}Errorcode: {str(e)}')
        raise
    # Check for errors in reg_aladin
    if result.returncode != 0:
        print(f"\nCommand failed: {command}\n")
        print("STDOUT:\n", result.stdout)
        print("STDERR:\n", result.stderr)
        raise RuntimeError(f"Command failed: {command}")

    # Compose transformation: DTI -> T2 -> atlas, i.e. bsplineMatrix(outputAff(x)).
    outputComposite = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_AtlasToDTI_deformation.nii.gz')
    command_args = [
        "reg_transform", "-ref", inputVolume, "-ref2", T2data,
        "-comp", outputAff, bsplineMatrix, outputComposite,
    ]
    command = shlex.join(command_args)
    try:
        result = subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        print(f"Output of {command}:\n{result.stdout}")
    except Exception as e:
        print(f'Error while executing the command: {command_args}Errorcode: {str(e)}')
        raise
    if result.returncode != 0:
        print(f"\nCommand failed: {command}\n")
        print("STDOUT:\n", result.stdout)
        print("STDERR:\n", result.stderr)
        raise RuntimeError(f"Command failed: {command}")

    # Resample the original unsplit Allen atlas once onto the DTI BET grid.
    outputAnno = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_Anno.nii.gz')
    command_args = [
        "reg_resample", "-ref", inputVolume, "-flo", allenBrain_anno,
        "-trans", outputComposite, "-inter", "0", "-res", outputAnno,
    ]
    command = shlex.join(command_args)
    try:
        result = subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        print(f"Output of {command}:\n{result.stdout}")
    except Exception as e:
        print(f'Error while executing the command: {command_args}Errorcode: {str(e)}')
        raise
    if result.returncode != 0:
        print(f"\nCommand failed: {command}\n")
        print("STDOUT:\n", result.stdout)
        print("STDERR:\n", result.stderr)
        raise RuntimeError(f"Command failed: {command}")

    # resample split  Annotation
    outputAnnoSplit = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_AnnoSplit.nii.gz')

    # Resample the original atlas once onto the DTI BET grid, preserving region labels.
    command = f"reg_resample -ref {inputVolume} -flo {splitAnno} -trans {outputComposite} -inter 0 -res {outputAnnoSplit}"
    command_args = shlex.split(command)
    try:
        result = subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,text=True)
        print(f"Output of {command}:\n{result.stdout}")
    except Exception as e:
        print(f'Error while executing the command: {command_args}Errorcode: {str(e)}')
        raise
    # Check for errors in reg_resample
    if result.returncode != 0:
        print(f"\nCommand failed: {command}\n")
        print("STDOUT:\n", result.stdout)
        print("STDERR:\n", result.stderr)
        raise RuntimeError(f"Command failed: {command}")

    # resample split par Annotation
    outputAnnoSplit_par = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_AnnoSplit_parental.nii.gz')

    # Apply the composed atlas-to-DTI transform with a single label resampling.
    command = f"reg_resample -ref {inputVolume} -flo {splitAnno_rsfMRI} -trans {outputComposite} -inter 0 -res {outputAnnoSplit_par}"
    command_args = shlex.split(command)
    try:
        result = subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,text=True)
        print(f"Output of {command}:\n{result.stdout}")
    except Exception as e:
        print(f'Error while executing the command: {command_args}Errorcode: {str(e)}')
        raise
    # Check for errors in reg_resample
    if result.returncode != 0:
        print(f"\nCommand failed: {command}\n")
        print("STDOUT:\n", result.stdout)
        print("STDERR:\n", result.stderr)
        raise RuntimeError(f"Command failed: {command}")

    # resample par Annotation
    outputAnno_par = os.path.join(outfile,os.path.basename(inputVolume).split('.')[0] + '_Anno_parental.nii.gz')

    # Apply the same composite directly to the original parental atlas.
    command = f"reg_resample -ref {inputVolume} -flo {anno_rsfMRI} -trans {outputComposite} -inter 0 -res {outputAnno_par}"
    command_args = shlex.split(command)
    try:
        result = subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,text=True)
        print(f"Output of {command}:\n{result.stdout}")
    except Exception as e:
        print(f'Error while executing the command: {command_args}Errorcode: {str(e)}')
        raise

    # Check for errors in reg_resample
    if result.returncode != 0:
        print(f"\nCommand failed: {command}\n")
        print("STDOUT:\n", result.stdout)
        print("STDERR:\n", result.stderr)
        raise RuntimeError(f"Command failed: {command}")

    # Some scaled data for DSI Studio
    outfileDSI = os.path.join(os.path.dirname(inputVolume), 'DSI_studio')
    if os.path.exists(outfileDSI):
        shutil.rmtree(outfileDSI) #? script-based removal of directories not recommended. Maybe change? // VVF 23/10/05
    os.makedirs(outfileDSI)
    outputRefStrokeMaskAff = None
    #only done if a reference stroke mask is provided
    if refStroke_mask is not None and len(refStroke_mask) > 0 and os.path.exists(refStroke_mask):
        refMatrix = find_RefAff(inputVolume)[0]
        refMTemplate = find_RefTemplate(inputVolume)[0]
        outputRefStrokeMaskAff = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_refStrokeMaskAff.nii.gz')

        #transform the reference stroke mask from atlas space to T2 space
        #refMTemplate: anat/*TemplateAff.nii.gz
        #refStroke_mask: Stroke mask in atlas space
        #refMatrix: affine matrix atlas to T2 space (anat/*MatrixAff.txt)
        command = f"reg_resample -ref {refMTemplate} -flo {refStroke_mask} -trans {refMatrix} -inter 0 -res {outputRefStrokeMaskAff}"
        command_args = shlex.split(command)
        try:
            result = subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,text=True)
            print(f"Output of {command}:\n{result.stdout}")
        except Exception as e:
            print(f'Error while executing the command: {command_args}Errorcode: {str(e)}')
            raise

        # Check for errors in reg_resample
        if result.returncode != 0:
            print(f"\nCommand failed: {command}\n")
            print("STDOUT:\n", result.stdout)
            print("STDERR:\n", result.stderr)
            raise RuntimeError(f"Command failed: {command}")

        stroke_mask = outputRefStrokeMaskAff



    if stroke_mask is not None and len(stroke_mask) > 0 and os.path.exists(stroke_mask):
        outputStrokeMask = os.path.join(outfile,
                                        os.path.basename(inputVolume).split('.')[0] + 'Stroke_mask.nii.gz')

        #resample the stroke mask to the DTI space using the affine T2-DTI transformation
        command = f"reg_resample -ref {inputVolume} -flo {stroke_mask} -inter 0 -trans {outputAff} -res {outputStrokeMask}"
        command_args = shlex.split(command)
        try:
            result = subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,text=True)
            print(f"Output of {command}:\n{result.stdout}")
        except Exception as e:
            print(f'Error while executing the command: {command_args}Errorcode: {str(e)}')
            raise

        # Check for errors in reg_resample
        if result.returncode != 0:
            print(f"\nCommand failed: {command}\n")
            print("STDOUT:\n", result.stdout)
            print("STDERR:\n", result.stderr)
            raise RuntimeError(f"Command failed: {command}")

        # Binary mask of the split annotation.
        dataAnno = nib.load(outputAnnoSplit)
        imgAnno_mask = dataAnno.get_fdata()
        imgAnno_mask[imgAnno_mask > 0] = 1
        imgAnno_mask[imgAnno_mask == 0] = 0
        imgAnno_mask = imgAnno_mask.astype(np.uint8)

        unscaledNiiData = nib.Nifti1Image(imgAnno_mask, dataAnno.affine)
        hdrOut = unscaledNiiData.header
        hdrOut.set_xyzt_units('mm')
        nib.save(unscaledNiiData,
                 os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + 'Anno_mask.nii.gz'))

        # Labelled stroke ROI for DSI Studio connectivity.
        dataAnno = nib.load(outputAnnoSplit_par)
        dataStroke = nib.load(outputStrokeMask)
        imgAnno = dataAnno.get_fdata()
        imgStroke = dataStroke.get_fdata()
        imgStroke[imgStroke > 0] = 1
        imgStroke[imgStroke == 0] = 0

        superPosAnnoStroke = imgStroke * imgAnno
        outputStrokeMaskAnno = os.path.join(
            outfile,
            os.path.basename(inputVolume).split('.')[0] + 'Stroke_mask_anno.nii.gz'
        )

        unscaledNiiDataMask = nib.Nifti1Image(superPosAnnoStroke, dataStroke.affine)

        hdrOut = unscaledNiiDataMask.header
        hdrOut.set_xyzt_units('mm')
        nib.save(unscaledNiiDataMask, outputStrokeMaskAnno)
    # --- Safety checks for DSI Studio inputs ---
    base = os.path.basename(inputVolume).split('.')[0]
    #os.makedirs(outfileDSI, exist_ok=True)

    bet_mask_path = os.path.join(outfile, f"{base}_mask.nii.gz")

    #Textfiles for DSI Studio lookup
    script_dir = os.path.dirname(os.path.abspath(__file__))
    lib_dir = os.path.abspath(os.path.join(script_dir, os.pardir, os.pardir, "lib"))
    anno_lut_src = os.path.join(lib_dir, "ARA_annotationR+2000.nii.txt")
    annop_lut_src = os.path.join(lib_dir, "annoVolume+2000_rsfMRI.nii.txt")

    if not os.path.exists(bet_mask_path):
        raise RuntimeError(
            f"Required BET brain mask is missing:\n  {bet_mask_path}\n"
            "BET mask is mandatory for DSI Studio reconstruction."
        )

    # --- DSI Studio LUTs for original DWI-space annotation files ---
    missing_core = [p for p in [outputAnnoSplit, outputAnnoSplit_par] if not os.path.exists(p)]
    if missing_core:
        print("Notice: Missing annotations for DSI Studio connectivity:")
        for p in missing_core:
            print("  -", p)
    else:
        # Copy LUTs next to the original DWI-space NIfTIs so DSI Studio can
        # find labels without requiring redundant NIfTI copies in DSI_studio.
        if os.path.exists(anno_lut_src):
            shutil.copyfile(anno_lut_src, os.path.join(outfile, f"{base}_AnnoSplit.txt"))
        else:
            print(f"Notice: LUT missing: {anno_lut_src} (DSI will still load NIfTI, but labels may be missing)")

        if os.path.exists(annop_lut_src):
            shutil.copyfile(annop_lut_src, os.path.join(outfile, f"{base}_AnnoSplit_parental.txt"))
        else:
            print(f"Notice: LUT missing: {annop_lut_src}")

    if outputRefStrokeMaskAff is not None:
        os.remove(outputRefStrokeMaskAff)

    return outputAnnoSplit

def find_RefStroke(refStrokePath,inputVolume):
    search_patterns = [
        os.path.join(refStrokePath, os.path.basename(inputVolume)[0:9], '*', 'anat', 'IncidenceData', '*IncidenceData_Lesion_mask.nii.gz'),
        os.path.join(refStrokePath, os.path.basename(inputVolume)[0:9], '*', 'anat', '*', '*IncidenceData_mask.nii.gz'),
        os.path.join(refStrokePath, os.path.basename(inputVolume)[0:9], '*', 'anat', '*IncidenceData_mask.nii.gz'),
    ]
    path = []
    for pattern in search_patterns:
        path.extend(glob.glob(pattern, recursive=False))
    return path

def find_RefAff(inputVolume):
    parent_dir = os.path.dirname(os.path.dirname(inputVolume))
    path = glob.glob(os.path.join(parent_dir, 'anat', '*MatrixAff.txt'))
    return path

def find_RefTemplate(inputVolume):
    parent_dir = os.path.dirname(os.path.dirname(inputVolume))
    path = glob.glob(os.path.join(parent_dir, 'anat', '*TemplateAff.nii.gz'))
    return path


def find_relatedData(pathBase):
    pathT2 = glob.glob(pathBase+'*/anat/*Bet.nii.gz', recursive=False)
    pathStroke_mask = glob.glob(pathBase + '*/anat/*Stroke_mask.nii.gz', recursive=False)
    bsplineMatrix = glob.glob(pathBase + '*/anat/*MatrixBspline.nii', recursive=False)
    return pathT2, pathStroke_mask, bsplineMatrix


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Registration Allen Brain to DTI')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-i', '--inputVolume', help='Path to the BET file of DTI data after preprocessing',
                               required=True)

    parser.add_argument('-r', '--referenceDay', help='Reference Stroke mask (for example: P5)', nargs='?', type=str,
                        default=None)
    parser.add_argument('-s', '--splitAnno', help='Split annotations atlas', nargs='?', type=str,
                        default=os.path.abspath(os.path.join(os.getcwd(), os.pardir,os.pardir))+'/lib/ARA_annotationR+2000.nii.gz')
    parser.add_argument('-f', '--splitAnno_rsfMRI', help='Split annotations atlas for rsfMRI/DTI', nargs='?', type=str,
                        default=os.path.abspath(os.path.join(os.getcwd(), os.pardir,os.pardir))+'/lib/annoVolume+2000_rsfMRI.nii.gz')
    parser.add_argument('-a', '--anno_rsfMRI', help='Parental Annotations atlas for rsfMRI/DTI', nargs='?', type=str,
                        default=os.path.abspath(os.path.join(os.getcwd(), os.pardir,os.pardir))+'/lib/annoVolume.nii.gz')
    parser.add_argument('--allenBrain_anno', help='Original unsplit Allen Brain annotation in atlas space', type=str,
                        default=os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir, "lib", "annotation_50CHANGEDanno.nii.gz"))
    )

    args = parser.parse_args()

    stroke_mask = None
    inputVolume = None
    splitAnno = None
    splitAnno_rsfMRI = None
    anno_rsfMRI = None
        
    if args.inputVolume is not None:
        inputVolume = args.inputVolume
    if not os.path.exists(inputVolume):
        sys.exit("Error: '%s' is not an existing directory." % (inputVolume,))

    outfile = os.path.join(os.path.dirname(inputVolume)) #this will be something like E:\CRC_data\proc_data\sub-GVsT3c3m2\ses-Baseline
    if not os.path.exists(outfile):
        os.makedirs(outfile)
    start_output_tracking(outfile, "dwi", "registration")
    setup_script_logging(outfile, "registration.log")

    # find related  data
    pathT2, pathStroke_mask, bsplineMatrix = find_relatedData(os.path.dirname(outfile)) #this will be something like E:\CRC_data\proc_data\sub-GVsT3c3m2
    if len(pathT2) == 0:
        T2data = []
        sys.exit("Error: %s' has no reference T2 template." % (os.path.basename(inputVolume),))
    else:
        T2data = pathT2[0]

    if len(pathStroke_mask) == 0:
        pathStroke_mask = []
        print("Notice: '%s' has no defined reference (stroke) mask - will proceed without." % (os.path.basename(inputVolume),))
    else:
        stroke_mask = pathStroke_mask[0]

    if len(bsplineMatrix) == 0:
        bsplineMatrix = []
        sys.exit("Error: %s' has no bspline Matrix." % (os.path.basename(inputVolume),))
    else:
        bsplineMatrix = bsplineMatrix[0]


    # finde reference stroke mask
    refStroke_mask = None
    if args.referenceDay is not None:
        referenceDay = args.referenceDay
        refStrokePath = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(outfile))), referenceDay)

        if not os.path.exists(refStrokePath):
            sys.exit("Error: '%s' is not an existing directory." % (refStrokePath,))
        refStroke_mask = find_RefStroke(refStrokePath, inputVolume)
        if len(refStroke_mask) == 0:
            refStroke_mask = []
            print("Notice: '%s' has no defined reference (stroke) mask - will proceed without." % (os.path.basename(inputVolume),))
        else:
            refStroke_mask = refStroke_mask[0]

    if args.splitAnno is not None:
        splitAnno = args.splitAnno
    if not os.path.exists(splitAnno):
        sys.exit("Error: '%s' is not an existing directory." % (splitAnno,))

    if args.splitAnno_rsfMRI is not None:
        splitAnno_rsfMRI = args.splitAnno_rsfMRI
    if not os.path.exists(splitAnno_rsfMRI):
        sys.exit("Error: '%s' is not an existing directory." % (splitAnno_rsfMRI,))

    if args.anno_rsfMRI is not None:
        anno_rsfMRI = args.anno_rsfMRI
    if not os.path.exists(anno_rsfMRI):
        sys.exit("Error: '%s' is not an existing directory." % (anno_rsfMRI,))

    regABA2DTI(inputVolume, stroke_mask, refStroke_mask, T2data, splitAnno,splitAnno_rsfMRI,anno_rsfMRI,bsplineMatrix,outfile,allenBrain_anno=args.allenBrain_anno)

    print("Registration completed")
