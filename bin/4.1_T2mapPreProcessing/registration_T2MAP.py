"""
Created on 11/09/2023

@author: Marc Schneider
Neuroimaging & Neuroengineering
Department of Neurology
University Hospital Cologne


Documentation preface, added 23/05/09 by Victor Vera Frazao:
This document is currently in revision for improvement and fixing.
Specifically changes are made to allow compatibility of the pipeline with Ubuntu 18.04 systems 
and Ubuntu 18.04 Docker base images, respectively, as well as adapting to appearent changes of 
DSI-Studio that were applied since the AIDAmri v.1.1 release. As to date the DSI-Studio version 
used is the 2022/08/03 Ubuntu 18.04 release.
All changes and additional documentations within this script carry a signature with the writer's 
initials (e.g. VVF for Victor Vera Frazao) and the date at application, denoted after '//' at 
the end of the comment line. If code segments need clearance the comment line will be prefaced 
by '#?'. Changes are prefaced by '#>' and other comments are prefaced ordinalrily 
by '#'.
"""

import sys,os
import nibabel as nii
import numpy as np
import shutil
import glob
import subprocess
import shlex

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
DEFAULT_SIGMA_ATLAS = os.path.join(REPO_ROOT, 'lib', 'sigma', 'SIGMA_InVivo_Anatomical_Brain_Atlas.nii.gz')
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))
from common.artifact_manifest import start_output_tracking
from common.script_logging import setup_script_logging

def regRef2T2map(inputVolume,stroke_mask,refStroke_mask,T2data, brain_template,brain_anno, splitAnno,splitAnno_parental,anno_parental,bsplineMatrix,outfile):

    outputT2w = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_T2w.nii.gz')
    outputAff = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + 'transMatrixAff.txt')
    
    
    command = f"reg_aladin -ref {inputVolume} -flo {T2data} -res {outputT2w} -rigOnly -aff {outputAff}"
    command_args = shlex.split(command)
    try:
        result = subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,text=True)
        print(f"Output of {command}:\n{result.stdout}")
    except Exception as e:
        print(f'Error while executing the command: {command_args} Errorcode: {str(e)}')
        raise

    # resample split  Annotation
    outputAnnoSplit = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_AnnoSplit.nii.gz')
    
    command = f"reg_resample -ref {brain_anno} -flo {splitAnno} -trans {bsplineMatrix} -inter 0 -res {outputAnnoSplit}"
    command_args = shlex.split(command)
    try:
        result = subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,text=True)
        print(f"Output of {command}:\n{result.stdout}")
    except Exception as e:
        print(f'Error while executing the command: {command_args}Errorcode: {str(e)}')
        raise
        
    command = f"reg_resample -ref {inputVolume} -flo {outputAnnoSplit} -trans {outputAff} -inter 0 -res {outputAnnoSplit}"
    command_args = shlex.split(command)
    try:
        result = subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,text=True)
        print(f"Output of {command}:\n{result.stdout}")
    except Exception as e:
        print(f'Error while executing the command: {command_args} Errorcode: {str(e)}')
        raise    
        

    # resample split parental annotation
    outputAnnoSplit_parental = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_AnnoSplit_parental.nii.gz')
    
    command = f"reg_resample -ref {brain_anno} -flo {splitAnno_parental} -trans {bsplineMatrix} -inter 0 -res {outputAnnoSplit_parental}"
    command_args = shlex.split(command)
    try:
        result = subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,text=True)
        print(f"Output of {command}:\n{result.stdout}")
    except Exception as e:
        print(f'Error while executing the command: {command_args} Errorcode: {str(e)}')
        raise 
        
    command = f"reg_resample -ref {inputVolume} -flo {outputAnnoSplit_parental} -trans {outputAff} -inter 0 -res {outputAnnoSplit_parental}"
    command_args = shlex.split(command)
    try:
        result = subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,text=True)
        print(f"Output of {command}:\n{result.stdout}")
    except Exception as e:
        print(f'Error while executing the command: {command_args} Errorcode: {str(e)}')
        raise 

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
    pathTemplate = glob.glob(pathBase + '*/anat/*TemplateAff.nii.gz', recursive=False)
    if len(pathTemplate) == 0:
        pathTemplate = glob.glob(pathBase + '*/anat/*Template.nii.gz', recursive=False)
    bsplineMatrix =  glob.glob(pathBase + '*/anat/*MatrixBspline.nii', recursive=False)
    return pathT2,pathStroke_mask,pathAnno,pathTemplate,bsplineMatrix



if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Registration reference brain to T2map')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-i', '--inputVolume', help='Path to the BET file of T2map data after preprocessing',
                               required=True)

    parser.add_argument('-r', '--referenceDay', help='Reference Stroke mask (for example: P5)', nargs='?', type=str,
                        default=None)
    parser.add_argument('-s', '--splitAnno', help='Split annotations atlas', nargs='?', type=str,
                        default=DEFAULT_SIGMA_ATLAS)
    parser.add_argument('-f', '--splitAnno_parental', '--splitAnno_rsfMRI', dest='splitAnno_parental', help='Split annotations atlas for parental/T2map output', nargs='?', type=str,
                        default=DEFAULT_SIGMA_ATLAS)
    parser.add_argument('-a', '--anno_parental', '--anno_rsfMRI', dest='anno_parental', help='Parental annotations atlas for T2map output', nargs='?', type=str,
                        default=DEFAULT_SIGMA_ATLAS)

    args = parser.parse_args()

    stroke_mask = None
    inputVolume = None
    refStrokePath = None
    splitAnno = None
    splitAnno_parental = None
    anno_parental = None

    if args.inputVolume is not None:
        inputVolume = args.inputVolume
    if not os.path.exists(inputVolume):
        sys.exit("Error: '%s' is not an existing directory." % (inputVolume,))

    outfile = os.path.join(os.path.dirname(inputVolume))
    if not os.path.exists(outfile):
        os.makedirs(outfile)
    start_output_tracking(outfile, "t2map", "registration")
    setup_script_logging(outfile, "registration.log")

    # find related  data
    pathT2, pathStroke_mask, pathAnno, pathTemplate, bsplineMatrix = find_relatedData(os.path.dirname(outfile))
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

    if len(pathAnno) == 0:
        pathAnno = []
        sys.exit("Error: %s' has no reference annotations." % (os.path.basename(inputVolume),))
    else:
        brain_anno = pathAnno[0]

    if len(pathTemplate) == 0:
        pathTemplate = []
        sys.exit("Error: %s' has no reference template." % (os.path.basename(inputVolume),))
    else:
        brain_template = pathTemplate[0]

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

    if args.splitAnno_parental is not None:
        splitAnno_parental = args.splitAnno_parental
    if not os.path.exists(splitAnno_parental):
        sys.exit("Error: '%s' is not an existing directory." % (splitAnno_parental,))

    if args.anno_parental is not None:
        anno_parental = args.anno_parental
    if not os.path.exists(anno_parental):
        sys.exit("Error: '%s' is not an existing directory." % (anno_parental,))

    output = regRef2T2map(inputVolume, stroke_mask, refStroke_mask, T2data, brain_template, brain_anno, splitAnno,splitAnno_parental,anno_parental,bsplineMatrix,outfile)

    current_dir = os.path.dirname(inputVolume)
    search_string = os.path.join(current_dir, "*t2map.nii.gz")
    currentFile = glob.glob(search_string)

    search_string = os.path.join(current_dir, "*.nii*")
    created_imgs = glob.glob(search_string, recursive=True)

    os.chdir(os.path.dirname(os.getcwd()))
    for idx, img in enumerate(created_imgs):
        if img == None:
            continue
        #os.system('python adjust_orientation.py -i '+ str(img) + ' -t ' + currentFile[0])
        continue

    print("Registration completed")
