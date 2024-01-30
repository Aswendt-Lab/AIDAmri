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

def BET_2_MPIreg(inputVolume, stroke_mask,brain_template, sigmaBrain_template,sigmaBrain_anno,sigmaBrain_annorsfMRI,outfile,opt):
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

    command = f"reg_aladin -ref {sigmaBrain_template} -flo {inputVolume} -res {outputInc} -aff {outputIncAff}"
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

        command = f"reg_resample -ref {sigmaBrain_template} -flo {stroke_mask} -trans {outputIncAff} -res {outputIncStrokeMask}"
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

    # resmaple Sigma Brain Reference Template
    outputAnno = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_TemplateSigma.nii.gz')

    command = f"reg_resample -ref {inputVolume} -flo {sigmaBrain_template} -cpp {outputCPP} -res {outputAnno}"
    command_args = shlex.split(command)
    try:
        result = subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        print(f"Output of {command}:\n{result.stdout}")
    except Exception as e:
        print(f'Error while executing the command: {command_args}\Errorcode: {str(e)}')
        raise

    # resample parental annotations
    outputAnnorsfMRI = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_AnnorsfMRI.nii.gz')

    command = f"reg_resample -ref {inputVolume} -flo {sigmaBrain_annorsfMRI} -inter 0 -cpp {outputCPP} -res {outputAnnorsfMRI}"
    command_args = shlex.split(command)
    try:
        result = subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,text=True)
        print(f"Output of {command}:\n{result.stdout}")
    except Exception as e:
        print(f'Error while executing the command: {command_args}\Errorcode: {str(e)}')
        raise

    # resample annotations
    outputAnno = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_Anno.nii.gz')

    command = f"reg_resample -ref {inputVolume} -flo {sigmaBrain_anno} -inter 0 -cpp {outputCPP} -res {outputAnno}"
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



if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Registration from ABA to T2 Data')

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-i', '--inputVolume', help='Path to input file', required=True)

    parser.add_argument('-s', '--deformationStrength', help='integer: 1 - very strong deformation, 2 - strong deformation, 3 - medium deformation, 4 - weak deformation ', nargs='?', type=int,
                        default=3)
    parser.add_argument('-g', '--template', help='File: Templates for Sigma Brain', nargs='?', type=str,
                        default=os.path.abspath(os.path.join(os.getcwd(), os.pardir,os.pardir))+'/lib/SIGMA_InVivo_Brain_Template_Masked.nii.gz')
    parser.add_argument('-t','--sigmaBrain_template', help='File: Templates of Sigma Brain', nargs='?', type=str,
                        default=os.path.abspath(
                            os.path.join(os.getcwd(), os.pardir, os.pardir)) + '/lib/SIGMA_InVivo_Brain_Template_Masked.nii.gz')
    parser.add_argument('-a','--sigmaBrain_anno', help='File: Annotations of Sigma Brain', nargs='?', type=str,
                        default=os.path.abspath(
                            os.path.join(os.getcwd(), os.pardir, os.pardir)) + '/lib/SIGMA_InVivo_Anatomical_Brain_Atlas.nii.gz')
    parser.add_argument('-f', '--sigmaBrain_annorsfMRI', help='File: Annotations of Sigma Brain', nargs='?',
                        type=str,
                        default=os.path.abspath(
                            os.path.join(os.getcwd(), os.pardir, os.pardir)) + '/lib/SIGMA_InVivo_Anatomical_Brain_Atlas.nii.gz')

    args = parser.parse_args()

    inputVolume = None
    sigmaBrain_template = None
    sigmaBrain_anno = None
    brain_template = None
    sigmaBrain_annorsfMRI = None
    deformationStrength = args.deformationStrength

    if args.inputVolume is not None:
        inputVolume = args.inputVolume
    if not os.path.exists(inputVolume):
        sys.exit("Error: '%s' is not an existing directory." % (inputVolume,))

    if args.sigmaBrain_template is not None:
        sigmaBrain_template = args.sigmaBrain_template
    if not os.path.exists(sigmaBrain_template):
        sys.exit("Error: '%s' is not an existing directory." % (sigmaBrain_template,))
        
    if args.sigmaBrain_anno is not None:
        sigmaBrain_anno = args.sigmaBrain_anno
    if not os.path.exists(sigmaBrain_anno):
        sys.exit("Error: '%s' is not an existing directory." % (sigmaBrain_anno,))

    if args.sigmaBrain_annorsfMRI is not None:
        sigmaBrain_annorsfMRI = args.sigmaBrain_annorsfMRI
    if not os.path.exists(sigmaBrain_annorsfMRI):
        sys.exit("Error: '%s' is not an existing directory." % (sigmaBrain_annorsfMRI,))

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

    transInput = BET_2_MPIreg(inputVolume, stroke_mask,brain_template,sigmaBrain_template,sigmaBrain_anno,sigmaBrain_annorsfMRI,outfile,deformationStrength)

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



