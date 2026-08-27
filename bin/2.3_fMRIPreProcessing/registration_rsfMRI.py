"""
Created on 10/08/2017

@author: Niklas Pallast
Neuroimaging & Neuroengineering
Department of Neurology
University Hospital Cologne

"""

import sys,os
import glob
import shutil as sh
import subprocess
import shlex
import logging
from calendar import month_name
from datetime import datetime
from zoneinfo import ZoneInfo

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))
from common.artifact_manifest import start_output_tracking


LOGGER = logging.getLogger(__name__)
DISABLE_LOG_ENV = "AIDAMRI_DISABLE_SCRIPT_LOG"
REPORT_TIMEZONE = ZoneInfo("Europe/Berlin")


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


def require_single_match(matches, description):
    if len(matches) != 1:
        raise RuntimeError(
            "Expected exactly one %s, found %d: %s"
            % (description, len(matches), matches)
        )
    return matches[0]


def run_command(command):
    if isinstance(command, str):
        command_args = shlex.split(command)
        command_display = command
    else:
        command_args = command
        command_display = subprocess.list2cmdline(command_args)

    try:
        result = subprocess.run(
            command_args,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        LOGGER.info("Output of %s:\n%s", command_display, result.stdout)
        if result.stderr:
            LOGGER.info("Stderr of %s:\n%s", command_display, result.stderr)
        return result
    except Exception as e:
        LOGGER.error("Error while executing the command: %s\nErrorcode: %s", command_display, e)
        raise

# Parameters passed to regSIG2rsfMRI:
# inputVolume: Preprocessed rsfMRI reference image and final target space.
# T2data: Individual brain-extracted T2 image (*/anat/*Bet.nii.gz).
# brain_template: Allen intensity template already warped into T2 space.
# brain_anno: Allen label image already warped into T2 space.
# splitAnno: Detailed label atlas in the original Allen space.
# splitAnno_rsfMRI: Parental split-label atlas in the original Allen space.
# anno_rsfMRI: Parental label atlas in the original Allen space.
# bsplineMatrix: Existing non-linear Allen-to-T2 transformation.
# dref: If True, use existing results from the related DWI directory.
# outfile: Directory for registered images, matrix and log file.
# use_atlas_mask: If True, mask T2data with a dilated brain_anno mask.

def regABA2rsfMRI(inputVolume, T2data, brain_template, brain_anno, splitAnno, splitAnno_rsfMRI, anno_rsfMRI,
                  bsplineMatrix, dref, outfile):
    outputT2w = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_T2w.nii.gz')
    outputAff = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + 'transMatrixAff.txt')

    if dref:
        pathT2 = glob.glob(os.path.dirname(outfile) + '*/dwi/*T2w.nii.gz', recursive=False)
        sh.copy(require_single_match(pathT2, "DTI T2 image"), outputT2w)
    else:
        if use_atlas_mask:
            prefix = os.path.basename(inputVolume).split('.')[0]
            atlasMask = os.path.join(outfile, prefix + '_T2AtlasMask.nii.gz')
            maskedT2 = os.path.join(outfile, prefix + '_T2AtlasMasked.nii.gz')

            # Create a slightly dilated mask from the well-registered atlas
            # annotation and apply it to the imperfect T2 BET.
            for command_args in [
                ["fslmaths", brain_anno, "-bin", "-dilM", atlasMask], #binarize atlas and dilate it by one voxel layer for tolerance
                ["fslmaths", T2data, "-mas", atlasMask, maskedT2], #apply mask to the T2 BET
            ]:
                run_command(command_args)
            T2data = maskedT2

        # Rigidly registers the individual T2 BET image T2data (moving image,
        # -flo) to the rsfMRI space defined by inputVolume (reference, -ref).
        # outputT2w (-res) is the resampled T2 image in the rsfMRI grid.
        # At the same time, -aff creates outputAff, the new T2-to-rsfMRI
        # transformation matrix. Despite its name, this matrix contains only
        # rotations and translations because -rigOnly is specified.
        command = f"reg_aladin -ref {inputVolume} -flo {T2data} -res {outputT2w} -rigOnly -aff {outputAff}"
        run_command(command)
        outputAnno = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_Anno.nii.gz')
        # Transforms brain_anno from the individual T2 space (moving image,
        # -flo) into the rsfMRI grid of inputVolume (-ref). It applies
        # outputAff, the rigid T2-to-rsfMRI matrix created by reg_aladin above.
        # outputAnno (-res) is the individual annotation in rsfMRI space (*_Anno.nii.gz);
        # -inter 0 uses nearest-neighbour interpolation to preserve label IDs.
        command = f"reg_resample -ref {inputVolume} -flo {brain_anno} -cpp {outputAff} -inter 0 -res {outputAnno}"
        run_command(command)

    # resample split annotation
    outputAnnoSplit = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_AnnoSplit.nii.gz')
    if dref:
        pathT2 = glob.glob(os.path.dirname(outfile) + '*/dwi/*AnnoSplit.nii.gz', recursive=False)
        sh.copy(require_single_match(pathT2, "DTI split annotation"), outputAnnoSplit)
    else:
        outputAnnoSplit_T2 = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_AnnoSplit_T2w.nii.gz')
        # Transforms splitAnno from Allen atlas space (moving image, -flo) into
        # the individual T2 space, using the grid of brain_anno as the target
        # grid (-ref). bsplineMatrix is the existing non-linear Allen-to-T2
        # B-spline transformation created by the anatomical T2 workflow.
        # outputAnnoSplit_T2 (-res) is the split annotation in individual T2
        # space; -inter 0 preserves the discrete atlas label IDs.
        command = f"reg_resample -ref {brain_anno} -flo {splitAnno} -trans {bsplineMatrix} -inter 0 -res {outputAnnoSplit_T2}"
        run_command(command)
        # Transforms outputAnnoSplit_T2 from individual T2 space (moving image,
        # -flo) into the rsfMRI grid of inputVolume (-ref). It applies
        # outputAff, the rigid T2-to-rsfMRI matrix created by reg_aladin.
        # outputAnnoSplit (-res) is the split annotation in rsfMRI space;
        # -inter 0 preserves the discrete atlas label IDs.
        command = f"reg_resample -ref {inputVolume} -flo {outputAnnoSplit_T2} -trans {outputAff} -inter 0 -res {outputAnnoSplit}"
        run_command(command)
        os.remove(outputAnnoSplit_T2)

    # resample split parental annotation
    outputAnnoSplit_rsfMRI = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_AnnoSplit_parental.nii.gz')
    if dref:
        pathT2 = glob.glob(os.path.dirname(outfile) + '*/dwi/*AnnoSplit_parental.nii.gz', recursive=False)
        sh.copy(require_single_match(pathT2, "DTI parental split annotation"), outputAnnoSplit_rsfMRI)
    else:
        outputAnnoSplit_rsfMRI_T2 = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_AnnoSplit_parental_T2w.nii.gz')
        # Transforms the parental split annotation splitAnno_rsfMRI from Allen
        # atlas space (moving image, -flo) into the individual T2 grid defined
        # by brain_anno (-ref). It applies bsplineMatrix, the existing
        # non-linear Allen-to-T2 B-spline transformation from the anatomical
        # workflow. outputAnnoSplit_rsfMRI_T2 (-res) is the temporary parental
        # split annotation in T2 space; -inter 0 preserves its label IDs.
        command = f"reg_resample -ref {brain_anno} -flo {splitAnno_rsfMRI} -trans {bsplineMatrix} -inter 0 -res {outputAnnoSplit_rsfMRI_T2}"
        run_command(command)
        # Transforms outputAnnoSplit_rsfMRI_T2 from individual T2 space
        # (moving image, -flo) into the rsfMRI grid of inputVolume (-ref).
        # It applies outputAff, the rigid T2-to-rsfMRI matrix created above.
        # outputAnnoSplit_rsfMRI (-res) is the parental split annotation in
        # rsfMRI space; -inter 0 preserves the discrete atlas label IDs.
        command = f"reg_resample -ref {inputVolume} -flo {outputAnnoSplit_rsfMRI_T2} -trans {outputAff} -inter 0 -res {outputAnnoSplit_rsfMRI}"
        run_command(command)
        os.remove(outputAnnoSplit_rsfMRI_T2)

    # resample parental annotation
    outputAnno_rsfMRI = os.path.join(outfile,
                                          os.path.basename(inputVolume).split('.')[0] + '_Anno_parental.nii.gz')
    if dref:
        pathT2 = glob.glob(os.path.dirname(outfile) + '*/dwi/*Anno_parental.nii.gz', recursive=False)
        sh.copy(require_single_match(pathT2, "DTI parental annotation"), outputAnno_rsfMRI)
    else: 
        outputAnno_rsfMRI_T2 = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_Anno_parental_T2w.nii.gz')
        # Transforms the parental annotation anno_rsfMRI from Allen atlas space
        # (moving image, -flo) into the individual T2 grid defined by
        # brain_anno (-ref). It applies bsplineMatrix, the existing non-linear
        # Allen-to-T2 B-spline transformation from the anatomical workflow.
        # outputAnno_rsfMRI_T2 (-res) is the temporary parental annotation in
        # T2 space; -inter 0 preserves the discrete atlas label IDs.
        command = f"reg_resample -ref {brain_anno} -flo {anno_rsfMRI} -trans {bsplineMatrix} -inter 0 -res {outputAnno_rsfMRI_T2}"
        run_command(command)
        # Transforms outputAnno_rsfMRI_T2 from individual T2 space (moving
        # image, -flo) into the rsfMRI grid of inputVolume (-ref). It applies
        # outputAff, the rigid T2-to-rsfMRI matrix created by reg_aladin.
        # outputAnno_rsfMRI (-res) is the parental annotation in rsfMRI space;
        # -inter 0 preserves the discrete atlas label IDs.
        command = f"reg_resample -ref {inputVolume} -flo {outputAnno_rsfMRI_T2} -trans {outputAff} -inter 0 -res {outputAnno_rsfMRI}"
        run_command(command)
        os.remove(outputAnno_rsfMRI_T2)
        # resample in-house developed tempalate
        outputTemplate = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_Template.nii.gz')
        
        command = f"reg_resample -ref {inputVolume} -flo {brain_template} -trans {outputAff} -res {outputTemplate}"
        run_command(command)
    
    return outputAnnoSplit

def find_RefStroke(refStrokePath,inputVolume):
    search_patterns = [
        os.path.join(refStrokePath, os.path.basename(inputVolume)[0:9] + '*', 'anat', 'IncidenceData', 'IncidenceData_Lesion_mask.nii.gz'),
        os.path.join(refStrokePath, os.path.basename(inputVolume)[0:9] + '*', 'anat', '*IncidenceData_mask.nii.gz'),
    ]
    path = []
    for pattern in search_patterns:
        path.extend(glob.glob(pattern, recursive=False))
    return path

def find_RefAff(inputVolume):
    path =  glob.glob(os.path.dirname(os.path.dirname(inputVolume))+'/anat/*MatrixAff.txt', recursive=False)
    return path

def find_RefTemplate(inputVolume):
    path =  glob.glob(os.path.dirname(os.path.dirname(inputVolume))+'/anat/*TemplateAff.nii.gz', recursive=False)
    return path


def find_relatedData(pathBase):
    pathT2 =  glob.glob(pathBase+'*/anat/*Bet.nii.gz', recursive=False)
    pathStroke_mask = glob.glob(pathBase + '*/anat/*Stroke_mask.nii.gz', recursive=False)
    pathAnno = glob.glob(pathBase + '*/anat/*Anno.nii.gz', recursive=False)
    pathAllen = glob.glob(pathBase + '*/anat/*Allen.nii.gz', recursive=False)
    bsplineMatrix =  glob.glob(pathBase + '*/anat/*MatrixBspline.nii', recursive=False)
    return pathT2,pathStroke_mask,pathAnno,pathAllen,bsplineMatrix


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Registration of Allen Brain Atlas to rsfMRI')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-i', '--inputVolume', help='Path to rsfMRI data after preprocessing', required=True)
    parser.add_argument('-d', '--dtiasRef', action='store_true', help='use DTI as reference if data quality is low')
    parser.add_argument('-r', '--referenceDay', help='Reference Stroke mask', nargs='?', type=str,
                        default=None)
    parser.add_argument('-s', '--splitAnno', help='Split annotations atlas', nargs='?', type=str,
                        default=os.path.abspath(os.path.join(os.getcwd(), os.pardir,os.pardir))+'/lib/ARA_annotationR+2000.nii.gz')
    parser.add_argument('-f', '--splitAnno_rsfMRI', help='Split annotations atlas for rsfMRI', nargs='?', type=str,
                        default=os.path.abspath(os.path.join(os.getcwd(), os.pardir,os.pardir))+'/lib/annoVolume+2000_rsfMRI.nii.gz')
    parser.add_argument('-a', '--anno_rsfMRI', help='Annotations atlas for rsfMRI', nargs='?', type=str,
                        default=os.path.abspath(os.path.join(os.getcwd(), os.pardir,os.pardir))+'/lib/annoVolume.nii.gz')



    args = parser.parse_args()



    stroke_mask = None
    inputVolume = None
    refStrokePath = None
    splitAnno = None
    splitAnno_rsfMRI = None
    anno_rsfMRI = None

    if args.inputVolume is not None:
        inputVolume = args.inputVolume
    if not os.path.exists(inputVolume):
        sys.exit("Error: '%s' is not an existing directory." % (inputVolume,))

    outfile = os.path.join(os.path.dirname(inputVolume))
    if not os.path.exists(outfile):
        os.makedirs(outfile)
    start_output_tracking(outfile, "func", "registration")
    setup_logging(outfile)


    # find related  data
    pathT2, pathStroke_mask, pathAnno, pathTemplate, bsplineMatrix = find_relatedData(os.path.dirname(outfile))
    if len(pathT2) == 0:
        T2data = []
        sys.exit("Error: %s' has no reference T2 template." % (os.path.basename(inputVolume),))
    else:
        T2data = require_single_match(pathT2, "reference T2 template")

    if len(pathStroke_mask) == 0:
        pathStroke_mask = []
        LOGGER.info("Notice: '%s' has no defined reference (stroke) mask - will proceed without.", os.path.basename(inputVolume))
    else:
        stroke_mask = require_single_match(pathStroke_mask, "stroke mask")

    if len(pathAnno) == 0:
        pathAnno = []
        sys.exit("Error: %s' has no reference annotations." % (os.path.basename(inputVolume),))
    else:
        brain_anno = require_single_match(pathAnno, "reference annotation")

    if len(pathTemplate) == 0:
        pathTemplate = []
        sys.exit("Error: %s' has no reference template." % (os.path.basename(inputVolume),))
    else:
        brain_template = require_single_match(pathTemplate, "reference template")

    if len(bsplineMatrix) == 0:
        bsplineMatrix = []
        sys.exit("Error: %s' has no bspline Matrix." % (os.path.basename(inputVolume),))
    else:
        bsplineMatrix = require_single_match(bsplineMatrix, "B-spline matrix")


    # find reference stroke mask
    refStroke_mask = None
    if args.referenceDay is not None:
        referenceDay = args.referenceDay
        refStrokePath = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(outfile))), referenceDay)

        if not os.path.exists(refStrokePath):
            sys.exit("Error: '%s' is not an existing directory." % (refStrokePath,))
        refStroke_mask = find_RefStroke(refStrokePath, inputVolume)
        if len(refStroke_mask) == 0:
            refStroke_mask = []
            LOGGER.info("Notice: '%s' has no defined reference (stroke) mask - will proceed without.", os.path.basename(inputVolume))
        else:
            refStroke_mask = require_single_match(refStroke_mask, "reference stroke mask")

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

    output = regABA2rsfMRI(inputVolume, T2data, brain_template, brain_anno, splitAnno, splitAnno_rsfMRI,
                           anno_rsfMRI, bsplineMatrix, args.dtiasRef, outfile)
    sys.stdout = sys.__stdout__

    current_dir = os.path.dirname(inputVolume)
    search_string = os.path.join(current_dir, "*EPI.nii.gz")
    currentFile = glob.glob(search_string)

    search_string = os.path.join(current_dir, "*.nii*")
    created_imgs = glob.glob(search_string, recursive=True)

    os.chdir(os.path.dirname(os.getcwd()))
    for idx, img in enumerate(created_imgs):
        if img == None:
            continue
        #os.system('python adjust_orientation.py -i '+ str(img) + ' -t ' + currentFile[0])

    LOGGER.info("Registration done")
