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

def regABA2T2map(inputVolume,stroke_mask,refStroke_mask,T2data, brain_template,brain_anno, splitAnno,splitAnno_rsfMRI,anno_rsfMRI,bsplineMatrix,outfile):

    outputT2w = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_T2w.nii.gz')
    outputAff = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + 'transMatrixAff.txt')
    
    
    command = f"reg_aladin -ref {inputVolume} -flo {T2data} -res {outputT2w} -rigOnly -aff {outputAff}"
    command_args = shlex.split(command)
    try:
        result = subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,text=True)
        print(f"Output of {command}:\n{result.stdout}")
    except Exception as e:
        print(f'Error while executing the command: {command_args}\Errorcode: {str(e)}')
        raise

    # resample Annotation
    #outputAnno = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_Anno.nii.gz')
    #os.system(
    #    'reg_resample -ref ' + inputVolume + ' -flo ' + brain_anno +
    #    ' -cpp ' + outputAff + ' -inter 0 -res ' + outputAnno)

    # resample split  Annotation
    outputAnnoSplit = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_AnnoSplit.nii.gz')
    
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
        

    # resample split rsfMRI Annotation
    outputAnnoSplit_rsfMRI = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_AnnoSplit_parental.nii.gz')
    
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


    # resample rsfMRI Annotation
    outputAnno_rsfMRI = os.path.join(outfile,
                                          os.path.basename(inputVolume).split('.')[0] + '_Anno_parental.nii.gz')
        
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



    # resample Template
    outputTemplate = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_Template.nii.gz')
        
    command = f"reg_resample -ref {inputVolume} -flo {brain_template} -cpp {outputAff} -res {outputTemplate}"
    command_args = shlex.split(command)
    try:
        result = subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,text=True)
        print(f"Output of {command}:\n{result.stdout}")
    except Exception as e:
        print(f'Error while executing the command: {command_args}\Errorcode: {str(e)}')
        raise  

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
            
        command = f"reg_resample -ref {refMTemplate} -flo {refStroke_mask} -cpp {refMatrix} -res {outputRefStrokeMaskAff}"
        command_args = shlex.split(command)
        try:
            result = subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,text=True)
            print(f"Output of {command}:\n{result.stdout}")
        except Exception as e:
            print(f'Error while executing the command: {command_args}\Errorcode: {str(e)}')
            raise 

        stroke_mask = outputRefStrokeMaskAff



    if stroke_mask is not None and len(stroke_mask) > 0 and os.path.exists(stroke_mask):
        outputStrokeMask = os.path.join(outfile,
                                        os.path.basename(inputVolume).split('.')[0] + 'Stroke_mask.nii.gz')
         
        command = f"reg_resample -ref {inputVolume} -flo {stroke_mask} -inter 0 -cpp {outputAff} -res {outputStrokeMask}"
        command_args = shlex.split(command)
        try:
            result = subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,text=True)
            print(f"Output of {command}:\n{result.stdout}")
        except Exception as e:
            print(f'Error while executing the command: {command_args}\Errorcode: {str(e)}')
            raise 

        # Superposition of annotations and mask
        dataAnno = nii.load(outputAnnoSplit)
        dataStroke = nii.load(outputStrokeMask)
        imgAnno = dataAnno.get_fdata()
        imgStroke = dataStroke.get_fdata()
        imgStroke[imgStroke > 0] = 1
        imgStroke[imgStroke == 0] = 0

        superPosAnnoStroke = imgStroke * imgAnno
        unscaledNiiData = nii.Nifti1Image(superPosAnnoStroke, dataAnno.affine)
        hdrOut = unscaledNiiData.header
        hdrOut.set_xyzt_units('mm')
        nii.save(unscaledNiiData,
                 os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + 'Anno_mask.nii.gz'))

        # Stroke Mask
        outputMaskScaled = os.path.join(outfileDSI,
                                        os.path.basename(inputVolume).split('.')[0] + 'StrokeMask_scaled.nii') #> removed '.gz' ending to correct atlas implementation // VVF 23/05/10
        superPosAnnoStroke = np.flip(superPosAnnoStroke, 2)
        # uperPosAnnoStroke = np.rot90(superPosAnnoStroke, 2)
        # superPosAnnoStroke = np.flip(superPosAnnoStroke, 0)
        scale = np.eye(4) * 10
        scale[3][3] = 1
        unscaledNiiDataMask = nii.Nifti1Image(superPosAnnoStroke, dataStroke.affine * scale)

        hdrOut = unscaledNiiDataMask.header
        hdrOut.set_xyzt_units('mm')
        nii.save(unscaledNiiDataMask, outputMaskScaled)
        src_file = os.path.join(os.path.abspath(os.path.join(os.getcwd(), os.pardir,os.pardir))+'/lib/', 'ARA_annotationR+2000.nii.txt')
        dst_file = os.path.join(outfileDSI, os.path.basename(inputVolume).split('.')[0] + 'StrokeMask_scaled.txt')#> removed '.nii.' ending to correct atlas implementation // VVF 23/05/10
        superPosAnnoStroke = np.flip(superPosAnnoStroke, 2)
        shutil.copyfile(src_file, dst_file)

        # Superposition of rsfMRI annotations and mask
        dataAnno = nii.load(outputAnnoSplit_rsfMRI)
        dataStroke = nii.load(outputStrokeMask)
        imgAnno = dataAnno.get_fdata()
        imgStroke = dataStroke.get_fdata()
        imgStroke[imgStroke > 0] = 1
        imgStroke[imgStroke == 0] = 0

        superPosAnnoStroke = imgStroke * imgAnno
        unscaledNiiData = nii.Nifti1Image(superPosAnnoStroke, dataAnno.affine)
        hdrOut = unscaledNiiData.header
        hdrOut.set_xyzt_units('mm')
        nii.save(unscaledNiiData,
                 os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + 'Anno_parental_mask.nii.gz'))
        superPosAnnoStroke = np.flip(superPosAnnoStroke, 2)

        # Stroke Mask
        outputMaskScaled = os.path.join(outfileDSI,
                                        os.path.basename(inputVolume).split('.')[0] + 'parental_Mask_scaled.nii') #> removed '.gz' ending to correct atlas implementation // VVF 23/05/10
        superPosAnnoStroke = np.flip(superPosAnnoStroke, 2)
        # superPosAnnoStroke = np.rot90(superPosAnnoStroke, 2)
        #superPosAnnoStroke = np.flip(superPosAnnoStroke, 0)

        scale = np.eye(4) * 10
        scale[3][3] = 1
        unscaledNiiDataMask = nii.Nifti1Image(superPosAnnoStroke, dataStroke.affine * scale)
        hdrOut = unscaledNiiDataMask.header
        hdrOut.set_xyzt_units('mm')
        nii.save(unscaledNiiDataMask, outputMaskScaled)
        src_file = os.path.join(os.path.abspath(os.path.join(os.getcwd(),os.pardir,os.pardir))+'/lib/annoVolume+2000_rsfMRI.nii.txt') 
        dst_file = os.path.join(outfileDSI, os.path.basename(inputVolume).split('.')[0] + 'parental_Mask_scaled.txt') #> removed '.nii.' ending to correct atlas implementation // VVF 23/05/10
        superPosAnnoStroke = np.flip(superPosAnnoStroke, 2)
        shutil.copyfile(src_file, dst_file)

    # Mask
    outputMaskScaled = os.path.join(outfileDSI, os.path.basename(inputVolume).split('.')[0] + 'Mask_scaled.nii') #> removed '.gz' ending to correct atlas implementation // VVF 23/05/10
    dataMask = nii.load(os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_mask.nii.gz'))
    imgMask = dataMask.get_fdata()

    imgMask = np.flip(imgMask, 2)
    # imgMask = np.rot90(imgMask, 2)
    # imgMask = np.flip(imgMask, 0)
    scale = np.eye(4) * 10
    scale[3][3] = 1

    unscaledNiiDataMask = nii.Nifti1Image(imgMask, dataMask.affine * scale)
    hdrOut = unscaledNiiDataMask.header
    hdrOut.set_xyzt_units('mm')
    nii.save(unscaledNiiDataMask, outputMaskScaled)

    # Allen Brain
    outputAnnoScaled = os.path.join(outfileDSI, os.path.basename(inputVolume).split('.')[0] + 'Anno_scaled.nii') #> removed '.gz' ending to correct atlas implementation // VVF 23/05/10
    outputAnnorsfMRIScaled = os.path.join(outfileDSI, os.path.basename(inputVolume).split('.')[
        0] + 'AnnoSplit_parental_scaled.nii')  #> removed '.gz' ending to correct atlas implementation // VVF 23/05/10
    outputAllenBScaled = os.path.join(outfileDSI, os.path.basename(inputVolume).split('.')[0] + 'Allen_scaled.nii') #> removed '.gz' ending to correct atlas implementation // VVF 23/05/10

    src_file = os.path.join(os.path.abspath(os.path.join(os.getcwd(), os.pardir,os.pardir))+'/lib/', 'ARA_annotationR+2000.nii.txt')
    dst_file = os.path.join(outfileDSI, os.path.basename(inputVolume).split('.')[0] + 'Anno_scaled.txt') #> removed '.nii.' ending to correct atlas implementation // VVF 23/05/10
    shutil.copyfile(src_file, dst_file)

    src_file = os.path.join(os.path.abspath(os.path.join(os.getcwd(), os.pardir,os.pardir))+'/lib/', 'annoVolume+2000_rsfMRI.nii.txt')
    dst_file = os.path.join(outfileDSI, os.path.basename(inputVolume).split('.')[0] + 'AnnoSplit_parental_scaled.txt') #> removed '.nii.' ending to correct atlas implementation // VVF 23/05/10
    shutil.copyfile(src_file, dst_file)


    dataAnno = nii.load(os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_AnnoSplit.nii.gz'))
    dataAnnorsfMRI = nii.load(os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_AnnoSplit_parental.nii.gz'))
    dataAllen = nii.load(os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_Template.nii.gz'))

    imgTempAnno = dataAnno.get_fdata()
    imgTempAnnorsfMRI = dataAnnorsfMRI.get_fdata()
    imgTempAllen = dataAllen.get_fdata()

    imgTempAllen = np.flip(imgTempAllen, 2)
    imgTempAnno = np.flip(imgTempAnno, 2)
    imgTempAnnorsfMRI = np.flip(imgTempAnnorsfMRI, 2)
    # imgTempAllen = np.rot90(imgTempAllen, 2)
    # imgTempAnno = np.rot90(imgTempAnno, 2)
    # imgTempAnnorsfMRI = np.rot90(imgTempAnnorsfMRI, 2)
    # imgTempAllen = np.flip(imgTempAllen, 0)
    # imgTempAnno = np.flip(imgTempAnno, 0)
    #imgTempAnnorsfMRI = np.flip(imgTempAnnorsfMRI, 0)
    scale = np.eye(4) * 10
    scale[3][3] = 1

    unscaledNiiDataAnno = nii.Nifti1Image(imgTempAnno, dataAnno.affine * scale)
    unscaledNiiDataAnnorsfMRI = nii.Nifti1Image(imgTempAnnorsfMRI, dataAnnorsfMRI.affine * scale)
    unscaledNiiDataAllen = nii.Nifti1Image(imgTempAllen, dataAllen.affine * scale)
    hdrOut = unscaledNiiDataAnno.header
    hdrOut.set_xyzt_units('mm')
    hdrOut = unscaledNiiDataAnnorsfMRI.header
    hdrOut.set_xyzt_units('mm')
    hdrOut = unscaledNiiDataAllen.header
    hdrOut.set_xyzt_units('mm')
    nii.save(unscaledNiiDataAnno, outputAnnoScaled)
    nii.save(unscaledNiiDataAnnorsfMRI, outputAnnorsfMRIScaled)
    nii.save(unscaledNiiDataAllen, outputAllenBScaled)

    if outputRefStrokeMaskAff is not None:
        os.remove(outputRefStrokeMaskAff)

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
    pathAllen = glob.glob(pathBase + '*/anat/*Allen.nii.gz', recursive=False)
    bsplineMatrix =  glob.glob(pathBase + '*/anat/*MatrixBspline.nii', recursive=False)
    return pathT2,pathStroke_mask,pathAnno,pathAllen,bsplineMatrix



if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Registration Allen Brain to T2map')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-i', '--inputVolume', help='Path to the BET file of T2map data after preprocessing',
                               required=True)

    parser.add_argument('-r', '--referenceDay', help='Reference Stroke mask (for example: P5)', nargs='?', type=str,
                        default=None)
    parser.add_argument('-s', '--splitAnno', help='Split annotations atlas', nargs='?', type=str,
                        default=os.path.abspath(os.path.join(os.getcwd(), os.pardir,os.pardir))+'/lib/ARA_annotationR+2000.nii.gz')
    parser.add_argument('-f', '--splitAnno_rsfMRI', help='Split annotations atlas for rsfMRI/T2map', nargs='?', type=str,
                        default=os.path.abspath(os.path.join(os.getcwd(), os.pardir,os.pardir))+'/lib/annoVolume+2000_rsfMRI.nii.gz')
    parser.add_argument('-a', '--anno_rsfMRI', help='Parental Annotations atlas for rsfMRI/T2map', nargs='?', type=str,
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


    # finde reference stroke mask
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

    output = regABA2T2map(inputVolume, stroke_mask, refStroke_mask, T2data, brain_template, brain_anno, splitAnno,splitAnno_rsfMRI,anno_rsfMRI,bsplineMatrix,outfile)

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



