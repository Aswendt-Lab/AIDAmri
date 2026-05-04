"""
Created on 10/08/2017

@author: Niklas Pallast
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
import nibabel as nib
import numpy as np
import shutil
import glob
import subprocess
import shlex

def regABA2DTI(inputVolume,stroke_mask,refStroke_mask,T2data, brain_template,brain_anno, splitAnno,splitAnno_rsfMRI,anno_rsfMRI,bsplineMatrix,outfile):
    outputT2w = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_T2w.nii.gz')
    outputAff = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + 'transMatrixAff.txt')
    
    
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

    # resample Annotation
    #outputAnno = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_Anno.nii.gz')
    #os.system(
    #    'reg_resample -ref ' + inputVolume + ' -flo ' + brain_anno +
    #    ' -cpp ' + outputAff + ' -inter 0 -res ' + outputAnno)

    # resample split  Annotation
    outputAnnoSplit = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_AnnoSplit.nii.gz')
    
    command = f"reg_resample -ref {inputVolume} -flo {splitAnno} -trans {outputAff} -inter 0 -res {outputAnnoSplit}"
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
    
    command = f"reg_resample -ref {brain_anno} -flo {splitAnno_rsfMRI} -trans {bsplineMatrix} -inter 0 -res {outputAnnoSplit_par}"
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

    # Some scaled data for DSI Studio
    outfileDSI = os.path.join(os.path.dirname(inputVolume), 'DSI_studio')
    if os.path.exists(outfileDSI):
        shutil.rmtree(outfileDSI) #? script-based removal of directories not recommended. Maybe change? // VVF 23/10/05
    os.makedirs(outfileDSI)
    outputRefStrokeMaskAff = None
    if refStroke_mask is not None and len(refStroke_mask) > 0 and os.path.exists(refStroke_mask):
        refMatrix = find_RefAff(inputVolume)[0]
        refMTemplate = find_RefTemplate(inputVolume)[0]
        outputRefStrokeMaskAff = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_refStrokeMaskAff.nii.gz')
            
        command = f"reg_resample -ref {refMTemplate} -flo {refStroke_mask} -trans {refMatrix} -res {outputRefStrokeMaskAff}"
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

        # Superposition of annotations and mask
        dataAnno = nib.load(outputAnnoSplit_par)
        dataStroke = nib.load(outputStrokeMask)
        imgAnno = dataAnno.get_fdata()
        imgStroke = dataStroke.get_fdata()
        imgStroke[imgStroke > 0] = 1
        imgStroke[imgStroke == 0] = 0

        superPosAnnoStroke = imgStroke * imgAnno
        unscaledNiiData = nib.Nifti1Image(superPosAnnoStroke, dataAnno.affine)
        hdrOut = unscaledNiiData.header
        hdrOut.set_xyzt_units('mm')
        nib.save(unscaledNiiData,
                 os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + 'Anno_mask.nii.gz'))

        # Stroke Mask
        outputStrokeMaskScaled = os.path.join(outfileDSI,
                                        os.path.basename(inputVolume).split('.')[0] + 'StrokeMask_scaled.nii') #> removed '.gz' ending to correct atlas implementation // VVF 23/05/10
        # superPosAnnoStroke = np.flip(superPosAnnoStroke, 2)
        # uperPosAnnoStroke = np.rot90(superPosAnnoStroke, 2)
        # superPosAnnoStroke = np.flip(superPosAnnoStroke, 0)
        # scale = np.eye(4) * 10
        # scale[3][3] = 1
        # unscaledNiiDataMask = nii.Nifti1Image(superPosAnnoStroke, dataStroke.affine * scale)
        unscaledNiiDataMask = nib.Nifti1Image(superPosAnnoStroke, dataStroke.affine)

        hdrOut = unscaledNiiDataMask.header
        hdrOut.set_xyzt_units('mm')
        nib.save(unscaledNiiDataMask, outputStrokeMaskScaled)
        #src_file = os.path.join(os.path.abspath(os.path.join(os.getcwd(), os.pardir,os.pardir))+'/lib/', 'ARA_annotationR+2000.nii.txt')
        #dst_file = os.path.join(outfileDSI, os.path.basename(inputVolume).split('.')[0] + 'StrokeMask_scaled.txt')#> removed '.nii.' ending to correct atlas implementation // VVF 23/05/10
        # superPosAnnoStroke = np.flip(superPosAnnoStroke, 2)
        #shutil.copyfile(src_file, dst_file)

        # Superposition of rsfMRI annotations and mask
        dataAnno = nib.load(outputAnnoSplit_par)
        dataStroke = nib.load(outputStrokeMask)
        imgAnno = dataAnno.get_fdata()
        imgStroke = dataStroke.get_fdata()
        imgStroke[imgStroke > 0] = 1
        imgStroke[imgStroke == 0] = 0

        superPosAnnoStroke = imgStroke * imgAnno
        unscaledNiiData = nib.Nifti1Image(superPosAnnoStroke, dataAnno.affine)
        hdrOut = unscaledNiiData.header
        hdrOut.set_xyzt_units('mm')
        nib.save(unscaledNiiData,
                 os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + 'Anno_parental_mask.nii.gz'))
        #superPosAnnoStroke = np.flip(superPosAnnoStroke, 2)

        # Stroke Mask
        outputParentalStrokeLabelsScaled = os.path.join(outfileDSI,
                                        os.path.basename(inputVolume).split('.')[0] + 'parental_Mask_scaled.nii') #> removed '.gz' ending to correct atlas implementation // VVF 23/05/10
        #superPosAnnoStroke = np.flip(superPosAnnoStroke, 2)
        # superPosAnnoStroke = np.rot90(superPosAnnoStroke, 2)
        #superPosAnnoStroke = np.flip(superPosAnnoStroke, 0)

        # scale = np.eye(4) * 10
        # scale[3][3] = 1
        # unscaledNiiDataMask = nii.Nifti1Image(superPosAnnoStroke, dataStroke.affine * scale)
        unscaledNiiDataMask = nib.Nifti1Image(superPosAnnoStroke, dataStroke.affine)
        hdrOut = unscaledNiiDataMask.header
        hdrOut.set_xyzt_units('mm')
        nib.save(unscaledNiiDataMask, outputParentalStrokeLabelsScaled)
        #src_file = os.path.join(os.path.abspath(os.path.join(os.getcwd(),os.pardir,os.pardir))+'/lib/annoVolume+2000_rsfMRI.nii.txt')
        #dst_file = os.path.join(outfileDSI, os.path.basename(inputVolume).split('.')[0] + 'parental_Mask_scaled.txt') #> removed '.nii.' ending to correct atlas implementation // VVF 23/05/10
        # superPosAnnoStroke = np.flip(superPosAnnoStroke, 2)
        #shutil.copyfile(src_file, dst_file)

    # --- Safety checks for DSI export inputs ---
    base = os.path.basename(inputVolume).split('.')[0]
    #os.makedirs(outfileDSI, exist_ok=True)

    bet_mask_path = os.path.join(outfile, f"{base}_mask.nii.gz")
    anno_path = os.path.join(outfile, f"{base}_AnnoSplit.nii.gz")
    annop_path = os.path.join(outfile, f"{base}_AnnoSplit_parental.nii.gz")
    templ_path = os.path.join(outfile, f"{base}_Template.nii.gz")

    #Textfiles for DSI Studio lookup
    script_dir = os.path.dirname(os.path.abspath(__file__))
    lib_dir = os.path.abspath(os.path.join(script_dir, os.pardir, os.pardir, "lib"))
    anno_lut_src = os.path.join(lib_dir, "ARA_annotationR+2000.nii.txt")
    annop_lut_src = os.path.join(lib_dir, "annoVolume+2000_rsfMRI.nii.txt")

    needed = [anno_path, annop_path, templ_path]
    missing = [p for p in needed if not os.path.exists(p)]
    if missing:
        print("Notice: Skipping DSI_studio export, missing:")
        for p in missing:
            print("  -", p)

    # --- Mask export ---
    if not os.path.exists(bet_mask_path):
        raise RuntimeError(
            f"Required BET brain mask is missing:\n  {bet_mask_path}\n"
            "BET mask is mandatory for DSI Studio export."
        )

    outputBrainMaskScaled = os.path.join(outfileDSI, f"{base}Mask_scaled.nii")
    dataMask = nib.load(bet_mask_path)
    imgMask = dataMask.get_fdata()

    # imgMask = np.flip(imgMask, 2)
    # scale = np.eye(4) * 10
    # scale[3][3] = 1

    unscaledNiiDataMask = nib.Nifti1Image(imgMask, dataMask.affine)
    unscaledNiiDataMask.header.set_xyzt_units('mm')
    nib.save(unscaledNiiDataMask, outputBrainMaskScaled)

    # --- Anno/Template export (optional, but independent of mask) ---
    missing_core = [p for p in [anno_path, annop_path, templ_path] if not os.path.exists(p)]
    if missing_core:
        print("Notice: Missing Anno/Template for DSI export, skipping these:")
        for p in missing_core:
            print("  -", p)
    else:
        outputAnnoScaled = os.path.join(outfileDSI, f"{base}Anno_scaled.nii")
        outputAnnorparScaled = os.path.join(outfileDSI, f"{base}AnnoSplit_parental_scaled.nii")
        outputAllenBScaled = os.path.join(outfileDSI, f"{base}Allen_scaled.nii")

        # copy LUTs (only if present)
        if os.path.exists(anno_lut_src):
            shutil.copyfile(anno_lut_src, os.path.join(outfileDSI, f"{base}Anno_scaled.txt"))
        else:
            print(f"Notice: LUT missing: {anno_lut_src} (DSI will still load NIfTI, but labels may be missing)")

        if os.path.exists(annop_lut_src):
            shutil.copyfile(annop_lut_src, os.path.join(outfileDSI, f"{base}AnnoSplit_parental_scaled.txt"))
        else:
            print(f"Notice: LUT missing: {annop_lut_src}")

        dataAnno = nib.load(anno_path)
        dataAnnorspar = nib.load(annop_path)
        dataAllen = nib.load(templ_path)

        imgTempAnno = dataAnno.get_fdata()
        imgTempAnnorspar = dataAnnorspar.get_fdata()
        imgTempAllen = dataAllen.get_fdata()

        # pipeline expects these flips for DSI export (keep consistent with your existing behavior)
        # imgTempAllen = np.flip(imgTempAllen, 2)
        # imgTempAnno = np.flip(imgTempAnno, 2)
        # imgTempAnnorspar = np.flip(imgTempAnnorspar, 2)

        # scale = np.eye(4) * 10
        # scale[3][3] = 1

        unscaledNiiDataAnno = nib.Nifti1Image(imgTempAnno, dataAnno.affine)
        unscaledNiiDataAnnorspar = nib.Nifti1Image(imgTempAnnorspar, dataAnnorspar.affine)
        unscaledNiiDataAllen = nib.Nifti1Image(imgTempAllen, dataAllen.affine)

        unscaledNiiDataAnno.header.set_xyzt_units('mm')
        unscaledNiiDataAnnorspar.header.set_xyzt_units('mm')
        unscaledNiiDataAllen.header.set_xyzt_units('mm')

        nib.save(unscaledNiiDataAnno, outputAnnoScaled)
        nib.save(unscaledNiiDataAnnorspar, outputAnnorparScaled)
        nib.save(unscaledNiiDataAllen, outputAllenBScaled)

    if outputRefStrokeMaskAff is not None:
        os.remove(outputRefStrokeMaskAff)

    return outputAnnoSplit

def find_RefStroke(refStrokePath,inputVolume):
    path =  glob.glob(os.path.join(refStrokePath, os.path.basename(inputVolume)[0:9],'*',"anat","*","IncidenceData_mask.nii.gz"), recursive=False)
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
    pathAnno = glob.glob(pathBase + '*/anat/*Anno.nii.gz', recursive=False)

    pathTemplate = glob.glob(pathBase + '*/anat/*TemplateAff.nii.gz', recursive=False)
    if len(pathTemplate) == 0:
        pathTemplate = glob.glob(pathBase + '*/anat/*Template.nii.gz', recursive=False)

    bsplineMatrix = glob.glob(pathBase + '*/anat/*MatrixBspline.nii', recursive=False)
    return pathT2, pathStroke_mask, pathAnno, pathTemplate, bsplineMatrix


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

    outfile = os.path.join(os.path.dirname(inputVolume)) #this will be something like E:\CRC_data\proc_data\sub-GVsT3c3m2\ses-Baseline
    if not os.path.exists(outfile):
        os.makedirs(outfile)

    # find related  data
    pathT2, pathStroke_mask, pathAnno, pathTemplate, bsplineMatrix = find_relatedData(os.path.dirname(outfile)) #this will be something like E:\CRC_data\proc_data\sub-GVsT3c3m2
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

    if args.splitAnno_rsfMRI is not None:
        splitAnno_rsfMRI = args.splitAnno_rsfMRI
    if not os.path.exists(splitAnno_rsfMRI):
        sys.exit("Error: '%s' is not an existing directory." % (splitAnno_rsfMRI,))

    if args.anno_rsfMRI is not None:
        anno_rsfMRI = args.anno_rsfMRI
    if not os.path.exists(anno_rsfMRI):
        sys.exit("Error: '%s' is not an existing directory." % (anno_rsfMRI,))

    output = regABA2DTI(inputVolume, stroke_mask, refStroke_mask, T2data, brain_template, brain_anno, splitAnno,splitAnno_rsfMRI,anno_rsfMRI,bsplineMatrix,outfile)

    current_dir = os.path.dirname(inputVolume)
    # look for denoised data and register if found
    currentFile = glob.glob(os.path.join(current_dir, "*Patch2SelfDenoised.nii.gz"))
    if len(currentFile) == 0:
        currentFile = glob.glob(os.path.join(current_dir, "*dwi.nii.gz"))

    search_string = os.path.join(current_dir, "*.nii*")
    created_imgs = glob.glob(search_string, recursive=True)

    os.chdir(os.path.dirname(os.getcwd()))
    for idx, img in enumerate(created_imgs):
        if img == None:
            continue
        #os.system('python adjust_orientation.py -i '+ str(img) + ' -t ' + currentFile[0])
        continue

    print("Registration completed")

