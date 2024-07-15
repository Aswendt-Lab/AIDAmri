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

        
def regABA2rsfMRI(inputVolume, T2data, brain_template, brain_anno, splitAnno, splitAnno_rsfMRI, anno_rsfMRI,
                  bsplineMatrix, dref, outfile):
    outputT2w = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_T2w.nii.gz')
    outputAff = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + 'transMatrixAff.txt')

    if dref:
        pathT2 = glob.glob(os.path.dirname(outfile) + '*/dwi/*T2w.nii.gz', recursive=False)
        sh.copy(pathT2[0], outputT2w)
    else:
        command = f"reg_aladin -ref {inputVolume} -flo {T2data} -res {outputT2w} -rigOnly -aff {outputAff}"
        command_args = shlex.split(command)
        try:
            result = subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,text=True)
            print(f"Output of {command}:\n{result.stdout}")
        except Exception as e:
            print(f'Error while executing the command: {command_args}\Errorcode: {str(e)}')
            raise
        #  resample Annotation
        outputAnno = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_Anno.nii.gz')

        command = f"reg_resample -ref {inputVolume} -flo {brain_anno} -cpp {outputAff} -inter 0 -res {outputAnno}"
        command_args = shlex.split(command)
        try:
            result = subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,text=True)
            print(f"Output of {command}:\n{result.stdout}")
        except Exception as e:
            print(f'Error while executing the command: {command_args}\Errorcode: {str(e)}')
            raise

    # resample split annotation
    outputAnnoSplit = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_AnnoSplit.nii.gz')
    if dref:
        pathT2 = glob.glob(os.path.dirname(outfile) + '*/dwi/*AnnoSplit.nii.gz', recursive=False)
        sh.copy(pathT2[0], outputAnnoSplit)
    else:
        command = f"reg_resample -ref {brain_anno} -flo {splitAnno} -trans {bsplineMatrix} -inter 0 -res {outputAnnoSplit}"
        command_args = shlex.split(command)
        try:
            result = subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,text=True)
            print(f"Output of {command}:\n{result.stdout}")
        except Exception as e:
            print(f'Error while executing the command: {command_args}\Errorcode: {str(e)}')
            raise
      
        command = f"reg_resample -ref {inputVolume} -flo {outputAnnoSplit} -trans {outputAff} -inter 0 -res {outputAnnoSplit}"
        command_args = shlex.split(command)
        try:
            result = subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,text=True)
            print(f"Output of {command}:\n{result.stdout}")
        except Exception as e:
            print(f'Error while executing the command: {command_args}\Errorcode: {str(e)}')
            raise

    # resample split parental annotation
    outputAnnoSplit_rsfMRI = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_AnnoSplit_parental.nii.gz')
    if dref:
        pathT2 = glob.glob(os.path.dirname(outfile) + '*/dwi/*AnnoSplit_parental.nii.gz', recursive=False)
        sh.copy(pathT2[0], outputAnnoSplit_rsfMRI)
    else:
        command = f"reg_resample -ref {brain_anno} -flo {splitAnno_rsfMRI} -trans {bsplineMatrix} -inter 0 -res {outputAnnoSplit_rsfMRI}"
        command_args = shlex.split(command)
        try:
            result = subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,text=True)
            print(f"Output of {command}:\n{result.stdout}")
        except Exception as e:
            print(f'Error while executing the command: {command_args}\Errorcode: {str(e)}')
            raise
        
        command = f"reg_resample -ref {inputVolume} -flo {outputAnnoSplit_rsfMRI} -trans {outputAff} -inter 0 -res {outputAnnoSplit_rsfMRI}"
        command_args = shlex.split(command)
        try:
            result = subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,text=True)
            print(f"Output of {command}:\n{result.stdout}")
        except Exception as e:
            print(f'Error while executing the command: {command_args}\Errorcode: {str(e)}')
            raise

    # resample parental annotation
    outputAnno_rsfMRI = os.path.join(outfile,
                                          os.path.basename(inputVolume).split('.')[0] + '_Anno_parental.nii.gz')
    if dref:
        pathT2 = glob.glob(os.path.dirname(outfile) + '*/dwi/*Anno_parental.nii.gz', recursive=False)
        sh.copy(pathT2[0], outputAnno_rsfMRI)
    else: 
        command = f"reg_resample -ref {brain_anno} -flo {anno_rsfMRI} -trans {bsplineMatrix} -inter 0 -res {outputAnno_rsfMRI}"
        command_args = shlex.split(command)
        try:
            result = subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,text=True)
            print(f"Output of {command}:\n{result.stdout}")
        except Exception as e:
            print(f'Error while executing the command: {command_args}\Errorcode: {str(e)}')
            raise
        
        command = f"reg_resample -ref {inputVolume} -flo {outputAnno_rsfMRI} -trans {outputAff} -inter 0 -res {outputAnno_rsfMRI}"
        command_args = shlex.split(command)
        try:
            result = subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,text=True)
            print(f"Output of {command}:\n{result.stdout}")
        except Exception as e:
            print(f'Error while executing the command: {command_args}\Errorcode: {str(e)}')
            raise
        
        # resample in-house developed tempalate
        outputTemplate = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_Template.nii.gz')
        
        command = f"reg_resample -ref {inputVolume} -flo {brain_template} -trans {outputAff} -res {outputTemplate}"
        command_args = shlex.split(command)
        try:
            result = subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,text=True)
            print(f"Output of {command}:\n{result.stdout}")
        except Exception as e:
            print(f'Error while executing the command: {command_args}\Errorcode: {str(e)}')
            raise


    return outputAnnoSplit

def find_RefStroke(refStrokePath,inputVolume):
    path =  glob.glob(refStrokePath+'/' + os.path.basename(inputVolume)[0:9]+'*/anat/*IncidenceData_mask.nii.gz', recursive=False)
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
    pathReference = glob.glob(pathBase + '*/anat/*Reference.nii.gz', recursive=False)
    bsplineMatrix =  glob.glob(pathBase + '*/anat/*MatrixBspline.nii', recursive=False)
    return pathT2,pathStroke_mask,pathAnno,pathReference,bsplineMatrix


#%% Program


#mice
#default_splitAnno = os.path.abspath(os.path.join(os.getcwd(), os.pardir,os.pardir))+'/lib/ARA_annotationR+2000.nii.gz'
#default_splitAnno_rsfMRI  = os.path.abspath(os.path.join(os.getcwd(), os.pardir,os.pardir))+'/lib/annoVolume+2000_rsfMRI.nii.gz'
#default_anno_rsfMRI = os.path.abspath(os.path.join(os.getcwd(), os.pardir,os.pardir))+'/lib/annoVolume.nii.gz'

default_splitAnno = os.path.abspath(os.path.join(os.getcwd(), os.pardir,os.pardir))+'/lib/SIGMA_InVivo_Anatomical_Brain_Atlas.nii.gz'
default_splitAnno_rsfMRI  = os.path.abspath(os.path.join(os.getcwd(), os.pardir,os.pardir))+'/lib/SIGMA_InVivo_Anatomical_Brain_Atlas.nii.gz'
default_anno_rsfMRI = os.path.abspath(os.path.join(os.getcwd(), os.pardir,os.pardir))+'/lib/SIGMA_InVivo_Anatomical_Brain_Atlas.nii.gz'
    
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Registration of Reference Brain Atlas to rsfMRI')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-i', '--inputVolume', help='Path to rsfMRI data after preprocessing', required=True)
    parser.add_argument('-d', '--dtiasRef', action='store_true', help='use DTI as reference if data quality is low')
    parser.add_argument('-r', '--referenceDay', help='Reference Stroke mask', nargs='?', type=str,
                        default=None)
    parser.add_argument('-s', '--splitAnno', help='Split annotations atlas', nargs='?', type=str,
                        default=default_splitAnno)
    parser.add_argument('-f', '--splitAnno_rsfMRI', help='Split annotations atlas for rsfMRI', nargs='?', type=str,
                        default=default_splitAnno_rsfMRI)
    parser.add_argument('-a', '--anno_rsfMRI', help='Annotations atlas for rsfMRI', nargs='?', type=str,
                        default=default_anno_rsfMRI)


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


    # find related  data
    pathT2, pathStroke_mask, pathAnno, pathTemplate, bsplineMatrix = find_relatedData(os.path.dirname(outfile))
    if len(pathT2) is 0:
        T2data = []
        sys.exit("Error: %s' has no reference T2 template." % (os.path.basename(inputVolume),))
    else:
        T2data = pathT2[0]

    if len(pathStroke_mask) is 0:
        pathStroke_mask = []
        print("Notice: '%s' has no defined reference (stroke) mask - will proceed without." % (os.path.basename(inputVolume),))
    else:
        stroke_mask = pathStroke_mask[0]

    if len(pathAnno) is 0:
        pathAnno = []
        sys.exit("Error: %s' has no reference annotations." % (os.path.basename(inputVolume),))
    else:
        brain_anno = pathAnno[0]

    if len(pathTemplate) is 0:
        pathTemplate = []
        sys.exit("Error: %s' has no reference template." % (os.path.basename(inputVolume),))
    else:
        brain_template = pathTemplate[0]

    if len(bsplineMatrix) is 0:
        bsplineMatrix = []
        sys.exit("Error: %s' has no bspline Matrix." % (os.path.basename(inputVolume),))
    else:
        bsplineMatrix = bsplineMatrix[0]


    # find reference stroke mask
    refStroke_mask = None
    if args.referenceDay is not None:
        referenceDay = args.referenceDay
        refStrokePath = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(outfile))), referenceDay)

        if not os.path.exists(refStrokePath):
            sys.exit("Error: '%s' is not an existing directory." % (refStrokePath,))
        refStroke_mask = find_RefStroke(refStrokePath, inputVolume)
        if len(refStroke_mask) is 0:
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

    print("Registration done")



