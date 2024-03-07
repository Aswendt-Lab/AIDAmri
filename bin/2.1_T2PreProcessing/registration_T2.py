"""
Created on 10/08/2017

@author: Niklas Pallast
Neuroimaging & Neuroengineering
Department of Neurology
University Hospital Cologne

"""


import sys,os
import numpy as np
import nibabel as nii
import glob
import subprocess
import shlex

def define_rodent_spezies():
    global rodent
    rodent = int(input("Select rodent: Mouse = 0 , Rat = 1 "))
    if rodent == 0 or rodent == 1:
        return rodent
    else:
        print("Invalid option. Enter 0 for mouse or 1 for rat.")
        return define_rodent_spezies()
        
def BET_2_MPIreg(inputVolume, stroke_mask,brain_template, ReferenceBrain_template,ReferenceBrain_anno,split_anno,anno_rsfMRI,split_ReferenceBrain_annorsfMRI,outfile,opt):
    output = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_TemplateAff.nii.gz')
    outputCPPAff = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + 'MatrixAff.txt')

    command = f"reg_aladin -ref {inputVolume} -flo {brain_template} -res {output} -aff {outputCPPAff}" 
    command_args = shlex.split(command)
    try:
        result = subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,text=True)
        print(f"Output of {command}:\n{result.stdout}")
    except Exception as e:
        print(f'Error while executing the command: {command_args}\Errorcode: {str(e)}')
        raise

    # Inverse registration
    outputInc = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_IncidenceData.nii.gz')
    outputIncAff = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + 'MatrixInv.txt')

    command = f"reg_aladin -ref {ReferenceBrain_template} -flo {inputVolume} -res {outputInc} -aff {outputIncAff}"
    command_args = shlex.split(command)
    try:
        result = subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,text=True)
        print(f"Output of {command}:\n{result.stdout}")
    except Exception as e:
        print(f'Error while executing the command: {command_args}\Errorcode: {str(e)}')
        raise

    # if region such as stroke_mask is defined
    if len(stroke_mask) > 0:
        outputIncStrokeMask = os.path.join(outfile, os.path.basename(outputInc).split('.')[0] + '_mask.nii.gz')

        command = f"reg_resample -ref {ReferenceBrain_template} -flo {stroke_mask} -trans {outputIncAff} -res {outputIncStrokeMask}"
        command_args = shlex.split(command)
        try:
            result = subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,text=True)
            print(f"Output of {command}:\n{result.stdout}")
        except Exception as e:
            print(f'Error while executing the command: {command_args}\Errorcode: {str(e)}')
            raise


    jac = 0.3
    # minimum defomraiton field in mm
    if opt == 1:
        s = [1, 1, 2]
    elif opt == 2:  s = [2,2,2]
    elif opt == 3:  s = [3,3,3]
    else:           s = [5,5,5]


    outputCPP = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + 'MatrixBspline.nii')

    # resample in-house developed template
    output = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_Template.nii.gz')

    command = f"reg_f3d -ref {inputVolume} -flo {brain_template} -sx {s[0]} -sy {s[1]} -sz {s[2]} -jl {jac} -res {output} -cpp {outputCPP} -aff {outputCPPAff}"
    command_args = shlex.split(command)
    try:
        result = subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,text=True)
        print(f"Output of {command}:\n{result.stdout}")
    except Exception as e:
        print(f'Error while executing the command: {command_args}\Errorcode: {str(e)}')
        raise

    # resmaple Reference Brain Template
    outputAnno = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_TemplateReference.nii.gz')

    command = f"reg_resample -ref {inputVolume} -flo {ReferenceBrain_template} -cpp {outputCPP} -res {outputAnno}"
    command_args = shlex.split(command)
    try:
        result = subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        print(f"Output of {command}:\n{result.stdout}")
    except Exception as e:
        print(f'Error while executing the command: {command_args}\Errorcode: {str(e)}')
        raise
        
     # resample parental annotations
    outputAnnorsfMRI = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_Anno_parental.nii.gz')

    command = f"reg_resample -ref {inputVolume} -flo {anno_rsfMRI} -inter 0 -cpp {outputCPP} -res {outputAnnorsfMRI}"
    command_args = shlex.split(command)
    try:
        result = subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,text=True)
        print(f"Output of {command}:\n{result.stdout}")
    except Exception as e:
        print(f'Error while executing the command: {command_args}\Errorcode: {str(e)}')
        raise    

    # resample parental split annotations
    outputAnnorsfMRI_split = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_AnnoSplit_parental.nii.gz')

    command = f"reg_resample -ref {inputVolume} -flo {split_ReferenceBrain_annorsfMRI} -inter 0 -cpp {outputCPP} -res {outputAnnorsfMRI_split}"
    command_args = shlex.split(command)
    try:
        result = subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,text=True)
        print(f"Output of {command}:\n{result.stdout}")
    except Exception as e:
        print(f'Error while executing the command: {command_args}\Errorcode: {str(e)}')
        raise

    # resample annotations
    outputAnno = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_Anno.nii.gz')

    command = f"reg_resample -ref {inputVolume} -flo {ReferenceBrain_anno} -inter 0 -cpp {outputCPP} -res {outputAnno}"
    command_args = shlex.split(command)
    try:
        result = subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,text=True)
        print(f"Output of {command}:\n{result.stdout}")
    except Exception as e:
        print(f'Error while executing the command: {command_args}\Errorcode: {str(e)}')
        raise
        
    # resample parental split annotations
    outputAnnoSplit = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_AnnoSplit.nii.gz')

    command = f"reg_resample -ref {inputVolume} -flo {split_anno} -inter 0 -cpp {outputCPP} -res {outputAnnoSplit}"
    command_args = shlex.split(command)
    try:
        result = subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,text=True)
        print(f"Output of {command}:\n{result.stdout}")
    except Exception as e:
        print(f'Error while executing the command: {command_args}\Errorcode: {str(e)}')
        raise    
    

    return outputAnno

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

def clearAnno(araAnno,realBrain_anno,outfile):
    araData = nii.load(araAnno)
    araVol = araData.get_data()
    nullValues = araVol < 0.0
    araVol[nullValues] = 0.0
    araVol = np.memmap.round(araVol)

    realData = nii.load(realBrain_anno)
    realVal = realData.get_data()
    realVal = realVal.tolist()
    uniqueList = np.unique(realVal)

    for i in np.nditer(araVol,op_flags=['readwrite']):
        i[...] = find_nearest(uniqueList, i)

    scaledNiiData = nii.Nifti1Image(araVol, araData.affine)
    hdrIn = scaledNiiData.header
    hdrIn.set_xyzt_units('mm')
    output_file = os.path.join(outfile, 'reconstructedAnno.nii.gz')
    nii.save(scaledNiiData, output_file)

    return outfile

def find_mask(inputVolume):
    return glob.glob(os.path.dirname(inputVolume)+'/*Stroke_mask.nii.gz', recursive=False)


#%% Program

#specify default Arguments by defining rodent spezies
define_rodent_spezies()

if rodent == 0:
    default_template = os.path.abspath(os.path.join(os.getcwd(), os.pardir,os.pardir))+'/lib/NP_template_sc0.nii.gz'
    default_ReferenceBrain_template  = os.path.abspath(
                            os.path.join(os.getcwd(), os.pardir, os.pardir)) + '/lib/average_template_50.nii.gz'
    default_ReferenceBrain_anno = os.path.abspath(
                            os.path.join(os.getcwd(), os.pardir, os.pardir)) + '/lib/annotation_50CHANGEDanno.nii.gz'
    default_splitAnno = os.path.abspath(os.path.join(os.getcwd(), os.pardir,os.pardir))+'/lib/ARA_annotationR+2000.nii.gz'
    default_anno_rsfMRI = os.path.abspath(os.path.join(os.getcwd(), os.pardir,os.pardir))+'/lib/annoVolume.nii.gz'
    default_split_annorsfMRI = os.path.abspath(
                            os.path.join(os.getcwd(), os.pardir, os.pardir)) + '/lib/annoVolume+2000_rsfMRI.nii.gz'
elif rodent == 1:
    default_template = os.path.abspath(os.path.join(os.getcwd(), os.pardir,os.pardir))+'/lib/SIGMA_InVivo_Brain_Template_Masked.nii.gz'
    default_ReferenceBrain_template  = os.path.abspath(
                            os.path.join(os.getcwd(), os.pardir, os.pardir)) + '/lib/SIGMA_InVivo_Brain_Template_Masked.nii.gz'
    default_ReferenceBrain_anno = os.path.abspath(
                            os.path.join(os.getcwd(), os.pardir, os.pardir)) + '/lib/SIGMA_InVivo_Anatomical_Brain_Atlas.nii.gz'
    default_splitAnno = os.path.abspath(os.path.join(os.getcwd(), os.pardir,os.pardir))+'/lib/SIGMA_InVivo_Anatomical_Brain_Atlas.nii.gz'
    default_anno_rsfMRI = os.path.abspath(os.path.join(os.getcwd(), os.pardir,os.pardir))+'/lib/SIGMA_InVivo_Anatomical_Brain_Atlas.nii.gz'
    default_split_annorsfMRI = os.path.abspath(
                            os.path.join(os.getcwd(), os.pardir, os.pardir)) + '/lib/SIGMA_InVivo_Anatomical_Brain_Atlas.nii.gz'
    
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Registration from ABA to T2 Data')

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-i', '--inputVolume', help='Path to input file', required=True)

    parser.add_argument('-s', '--deformationStrength', help='integer: 1 - very strong deformation, 2 - strong deformation, 3 - medium deformation, 4 - weak deformation ', nargs='?', type=int,
                        default=3)
    parser.add_argument('-g', '--template', help='File: Templates for Reference Brain', nargs='?', type=str,
                        default=default_template)
    parser.add_argument('-t','--ReferenceBrain_template', help='File: Templates of Reference Brain', nargs='?', type=str,
                        default=default_ReferenceBrain_template)
    parser.add_argument('-a','--ReferenceBrain_anno', help='File: Annotations of Reference Brain', nargs='?', type=str,
                        default=default_ReferenceBrain_anno)
    parser.add_argument('-sa', '--splitAnno', help='Split annotations atlas', nargs='?', type=str,
                        default=default_splitAnno)
    parser.add_argument('-f', '--anno_rsfMRI', help='Parental Annotations atlas', nargs='?', type=str,
                        default=default_anno_rsfMRI)
    parser.add_argument('-sf', '--split_annorsfMRI', help='File: Annotations of split Reference Brain', nargs='?',
                        type=str,
                        default=default_split_annorsfMRI)

    args = parser.parse_args()

    inputVolume = None
    ReferenceBrain_template = None
    ReferenceBrain_anno = None
    split_anno = None
    brain_template = None
    split_ReferenceBrain_annorsfMRI = None
    anno_rsfMRI = None
    deformationStrength = args.deformationStrength

    if args.inputVolume is not None:
        inputVolume = args.inputVolume
    if not os.path.exists(inputVolume):
        sys.exit("Error: '%s' is not an existing directory." % (inputVolume,))

    if args.ReferenceBrain_template is not None:
        ReferenceBrain_template = args.ReferenceBrain_template
    if not os.path.exists(ReferenceBrain_template):
        sys.exit("Error: '%s' is not an existing directory." % (ReferenceBrain_template,))
        
    if args.ReferenceBrain_anno is not None:
        ReferenceBrain_anno = args.ReferenceBrain_anno
    if not os.path.exists(ReferenceBrain_anno):
        sys.exit("Error: '%s' is not an existing directory." % (ReferenceBrain_anno,))
        
    if args.splitAnno is not None:
        split_anno = args.splitAnno
    if not os.path.exists(split_anno):
        sys.exit("Error: '%s' is not an existing directory." % (split_anno,))

    if args.split_annorsfMRI is not None:
        split_ReferenceBrain_annorsfMRI = args.split_annorsfMRI
    if not os.path.exists(split_ReferenceBrain_annorsfMRI):
        sys.exit("Error: '%s' is not an existing directory." % (split_ReferenceBrain_annorsfMRI,))
    
    if args.anno_rsfMRI is not None:
        anno_rsfMRI = args.anno_rsfMRI
    if not os.path.exists(anno_rsfMRI):
        sys.exit("Error: '%s' is not an existing directory." % (anno_rsfMRI,))

    if args.template is not None:
        brain_template = args.template
    if not os.path.exists(brain_template):
        sys.exit("Error: '%s' is not an existing directory." % (brain_template,))
        

    outfile = os.path.join(os.path.dirname(inputVolume))
    if not os.path.exists(outfile):
        os.makedirs(outfile)

    stroke_mask = find_mask(inputVolume)
    if len(stroke_mask) is 0:
        stroke_mask = []
        print("Notice: '%s' has no defined reference (stroke) mask - will proceed without." % (inputVolume,))
    else:
        stroke_mask = stroke_mask[0]

    transInput = BET_2_MPIreg(inputVolume, stroke_mask,brain_template,ReferenceBrain_template,ReferenceBrain_anno,split_anno,anno_rsfMRI,split_ReferenceBrain_annorsfMRI,outfile,deformationStrength)

    current_dir = os.path.dirname(inputVolume)
    search_string = os.path.join(current_dir, "*T2w.nii.gz")
    currentFile = glob.glob(search_string)

    search_string = os.path.join(current_dir, "*.nii*")
    created_imgs = glob.glob(search_string, recursive=True)

    os.chdir(os.path.dirname(os.getcwd()))
    for idx, img in enumerate(created_imgs):
        if img == None:
            continue
        #os.system('python adjust_orientation.py -i '+ str(img) + ' -t ' + currentFile[0])
        
    print("Registration completed")



