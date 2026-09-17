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
import shutil as sh
import logging
from calendar import month_name
from datetime import datetime
from zoneinfo import ZoneInfo

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))
from common.artifact_manifest import start_output_tracking
from common.script_logging import setup_script_logging

LOGGER = logging.getLogger(__name__)
DISABLE_LOG_ENV = "AIDAMRI_DISABLE_SCRIPT_LOG"
REPORT_TIMEZONE = ZoneInfo("Europe/Berlin")
POSTERIOR_CROP_FRACTION = 0.50


class BerlinTimeFormatter(logging.Formatter):
    def formatTime(self, record, datefmt=None):
        timestamp = datetime.fromtimestamp(record.created, REPORT_TIMEZONE)
        return (
            f"{timestamp.day:02d} {month_name[timestamp.month]} {timestamp.year} "
            f"{timestamp:%H:%M:%S} {timestamp.tzname()}"
        )

def setup_logging(outfile):
    handlers = [logging.StreamHandler()]
    if os.environ.get(DISABLE_LOG_ENV) != "1":
        handlers.append(logging.FileHandler(os.path.join(outfile, "registration.log"), mode="w"))
    formatter = BerlinTimeFormatter("%(asctime)s %(levelname)s: %(message)s")
    for handler in handlers:
        handler.setFormatter(formatter)
    logging.basicConfig(level=logging.INFO, handlers=handlers, force=True)

def regSIG2DTI(inputVolume,stroke_mask,refStroke_mask,T2data, bsplineMatrix, outfile, sigBrain_anno):
    if not os.path.isfile(sigBrain_anno):
        raise FileNotFoundError(f"Original Sigma Brain annotation not found: {sigBrain_anno}")
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
        "reg_resample", "-ref", inputVolume, "-flo", sigBrain_anno,
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

    ''' Atlas files are all the same the following steps are redundant but are kept for clarity and future flexibility
        
    command = f"reg_resample -ref {inputVolume} -flo {outputAnnoSplit_par} -trans {outputAff} -inter 0 -res {outputAnnoSplit_par}"
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
    outputAnno_par = os.path.join(outfile,
                                          os.path.basename(inputVolume).split('.')[0] + '_Anno_parental.nii.gz')
        
    command = f"reg_resample -ref {brain_anno} -flo {anno_rsfMRI} -trans {bsplineMatrix} -inter 0 -res {outputAnno_par}"
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
       
    command = f"reg_resample -ref {inputVolume} -flo {outputAnno_par} -trans {outputAff} -inter 0 -res {outputAnno_par}"
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

    # resample Template
    outputTemplate = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_Template.nii.gz')
        
    command = f"reg_resample -ref {inputVolume} -flo {brain_template} -trans {outputAff} -res {outputTemplate}"
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
    '''
    #Create pipline downstream compatible output (copies of fMRI anno)
    prefix = os.path.basename(inputVolume).split('.')[0]
    outputAnnoSplit = os.path.join(outfile, prefix + '_AnnoSplit.nii.gz')
    outputAnnoSplit_par = os.path.join(outfile, prefix + '_AnnoSplit_parental.nii.gz')
    outputAnno_par = os.path.join(outfile, prefix + '_Anno_parental.nii.gz')

    # The split/parental atlas inputs are currently identical to outputAnno, so
    # create downstream-compatible filenames as direct copies.
    for target in [outputAnnoSplit, outputAnnoSplit_par, outputAnno_par]:
        sh.copyfile(outputAnno, target)
        LOGGER.info("Copied %s to %s", outputAnno, target)

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
    anno_lut_src = os.path.join(lib_dir, "sigma", "SIGMA_InVivo_Anatomical_Brain_Atlas_Labels.txt")
    annop_lut_src = os.path.join(lib_dir, "sigma", "SIGMA_InVivo_Anatomical_Brain_Atlas_Labels.txt")

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

    parser = argparse.ArgumentParser(description='Registration SIGMA Brain to DTI')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-i', '--inputVolume', help='Path to the BET file of DTI data after preprocessing',
                               required=True)

    parser.add_argument('-r', '--referenceDay', help='Reference Stroke mask (for example: P5)', nargs='?', type=str,
                        default=None)
    parser.add_argument('-s', '--splitAnno', help='Split annotations atlas', nargs='?', type=str,
                        default=os.path.abspath(os.path.join(os.getcwd(), os.pardir,os.pardir))+'/lib/sigma/SIGMA_InVivo_Anatomical_Brain_Atlas.nii.gz')
    parser.add_argument('-f', '--splitAnno_rsfMRI', help='Split annotations atlas for rsfMRI/DTI', nargs='?', type=str,
                        default=os.path.abspath(os.path.join(os.getcwd(), os.pardir,os.pardir))+'/lib/sigma/SIGMA_InVivo_Anatomical_Brain_Atlas.nii.gz')
    parser.add_argument('-a', '--anno_rsfMRI', help='Parental Annotations atlas for rsfMRI/DTI', nargs='?', type=str,
                        default=os.path.abspath(os.path.join(os.getcwd(), os.pardir,os.pardir))+'/lib/sigma/SIGMA_InVivo_Anatomical_Brain_Atlas.nii.gz')
    parser.add_argument('--sigBrain_anno', help='Original unsplit SIGMA Brain annotation in atlas space', type=str,
                        default=os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir, "lib", "sigma", "SIGMA_InVivo_Anatomical_Brain_Atlas.nii.gz"))
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

    regSIG2DTI(inputVolume, stroke_mask, refStroke_mask, T2data, bsplineMatrix,outfile, sigBrain_anno=args.sigBrain_anno)

    print("Registration completed")
