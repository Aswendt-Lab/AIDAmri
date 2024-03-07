"""
Created on 10/08/2017

@author: Niklas Pallast
Neuroimaging & Neuroengineering
Department of Neurology
University Hospital Cologne

"""


import sys, os
import nipype.interfaces.fsl as fsl
import nibabel as nii
import numpy as np
import glob
import shutil
import regress
import getSingleRegTable
import scipy.misc as mc
import create_seed_rois
import fsl_mean_ts
from pathlib import Path 
import json

def define_rodent_spezies():
    global rodent
    rodent = int(input("Select rodent: Mouse = 0 , Rat = 1 "))
    if rodent == 0 or rodent == 1:
        return rodent
    else:
        print("Invalid option. Enter 0 for mouse or 1 for rat.")
        return define_rodent_spezies()
        
def copyAtlasOfData(path,post,labels):
    fileALL = glob.glob(path + '/*' + post + '.nii.gz')
    if fileALL.__len__()>1:
        sys.exit("Error: '%s' has no related Atlas File." % (path,))
    else:
        fileALL = fileALL[0]

    print("Copy Atlas Data and generate seed ROIs")
    #pathfMRI = os.path.join(os.path.dirname(path),'fMRI')
    outputRois = create_seed_rois.startSeedPoint(in_atlas=os.path.join(path, os.path.basename(fileALL)),in_labels=labels)
    return outputRois

def imgScaleResize(img):
    newImg = np.zeros([128,128,20,355])
    for i in range(img.shape[3]):
        for j in range(img.shape[2]):
            newImg[:,:,j,i]=mc.imresize(img[:,:,j,i],1.34,interp='nearest')

    return newImg

def scaleBy10(input_path,inv):
    data = nii.load(input_path)
    imgTemp = data.get_fdata()
    if inv is False:
        scale = np.eye(4) * 10
        scale[3][3] = 1
        scaledNiiData = nii.Nifti1Image(imgTemp, data.affine * scale)
        fslPath = os.path.join(os.path.dirname(input_path), 'fslScaleTemp.nii.gz')
        nii.save(scaledNiiData, fslPath)
        return fslPath
    elif inv is True:
        scale = np.eye(4) / 10
        scale[3][3] = 1
        unscaledNiiData = nii.Nifti1Image(imgTemp, data.affine * scale)
        hdrOut = unscaledNiiData.header
        hdrOut.set_xyzt_units('mm')

        # hdrOut['sform_code'] = 1
        nii.save(unscaledNiiData, input_path)
        return input_path
    else:
        sys.exit("Error: inv - parameter should be a boolean.")

def findSlicesData(path,pre):
    regMR_list = []
    fileALL = glob.iglob(path+'/'+pre+'*.nii.gz',recursive=True)
    for filename in fileALL:
        regMR_list.append(filename)
    regMR_list.sort()
    return regMR_list

def getRASorientation(file_name,proc_Path):
    data = nii.load(file_name)
    imgData = data.get_fdata()

    imgData = np.flip(imgData, 2)
    imgData = np.flip(imgData, 0)

    epiData = nii.Nifti1Image(imgData, data.affine)
    hdrIn = epiData.header
    hdrIn.set_xyzt_units('mm')
    epiData_RAS = nii.as_closest_canonical(epiData)
    print('Orientation:' + str(nii.aff2axcodes(epiData_RAS.affine)))
    output_file = os.path.join(proc_Path, os.path.basename(file_name))
    nii.save(epiData, output_file)
    return output_file

def getEPIMean(file_name,proc_Path):
    output_file = os.path.join(proc_Path, os.path.basename(file_name).split('.')[0]) + 'mean.nii.gz'
    myMean = fsl.MeanImage(in_file=file_name, out_file=output_file)
    print(myMean.cmdline)
    myMean.run()
    return output_file

def applyBET(input_file,frac,radius,vertical_gradient):

    # scale Nifti data by factor 10
    fslPath = scaleBy10(input_file,inv=False)
    # extract brain
    output_file = os.path.join(os.path.dirname(input_file),os.path.basename(input_file).split('.')[0]) + 'Bet.nii.gz'
    maskFile = os.path.join(os.path.dirname(input_file), os.path.basename(input_file).split('.')[0]) + 'Bet_mask.nii.gz'
    myBet = fsl.BET(in_file=fslPath, out_file=output_file,frac=frac,radius=radius,
                    vertical_gradient=vertical_gradient,robust=True, mask = True)
    print(myBet.cmdline)
    myBet.run()
    os.remove(fslPath)
    # unscale result data by factor 10ˆ(-1)
    output_file = scaleBy10(output_file,inv=True)
    return output_file,maskFile

def applyMask(input_file,mask_file):
    fslPath = scaleBy10(input_file, inv=False)
    # maks apply
    output_file = os.path.join(os.path.dirname(input_file), os.path.basename(input_file).split('.')[0]) + 'BET.nii.gz'
    myMaskapply = fsl.ApplyMask(in_file=fslPath, out_file=output_file, mask_file=mask_file)
    print(myMaskapply.cmdline)
    myMaskapply.run()
    os.remove(fslPath)
    # unscale result data by factor 10ˆ(-1)
    output_file = scaleBy10(output_file, inv=True)
    return output_file

def fsl_SeparateSliceMoCo(input_file,par_folder):
    # scale Nifti data by factor 10
    dataName = os.path.basename(input_file).split('.')[0]

    aidamri_dir = os.getcwd()
    temp_dir = os.path.join(os.path.dirname(input_file), "temp")
    if not os.path.exists(temp_dir):
        os.mkdir(temp_dir)

    fslPath = scaleBy10(input_file, inv=False)
    os.chdir(temp_dir)
    mySplit= fsl.Split(in_file=fslPath,dimension='z',out_base_name = dataName)
    print(mySplit.cmdline)
    mySplit.run()
    os.remove(fslPath)


    # sparate ref and src volume in slices
    sliceFiles = findSlicesData(os.getcwd(),dataName)

    #start to correct motions slice by slice
    for i in  range(len(sliceFiles)):
        slc = sliceFiles[i]
        # ref = refFiles[i]
        # take epi as ref
        output_file = os.path.join(par_folder,os.path.basename(slc))
        myMCFLIRT = fsl.preprocess.MCFLIRT(in_file=slc,out_file=output_file,save_plots=True,terminal_output='none')
        myMCFLIRT.run()
        os.remove(slc)
        # os.remove(ref)

    # merge slices to a single volume

    mcf_sliceFiles = findSlicesData(par_folder,dataName)
    output_file = os.path.join(os.path.dirname(input_file),
                               os.path.basename(input_file).split('.')[0]) + '_mcf.nii.gz'
    myMerge = fsl.Merge(in_files=mcf_sliceFiles,dimension='z',merged_file =output_file)
    print(myMerge.cmdline)
    myMerge.run()

    for slc in mcf_sliceFiles: os.remove(slc)

    # unscale result data by factor 10ˆ(-1)
    output_file = scaleBy10(output_file, inv=True)
    
    #os.remove(temp_dir)
    os.chdir(aidamri_dir)

    return output_file

def copyRawPhysioData(file_name,i32_Path):
    img_name = Path(file_name).name
    json_name = img_name.replace(".nii.gz", ".json")
    
    json_file = os.path.join(os.path.dirname(file_name), json_name)
    sub_name = (Path(file_name).name.split("_")[0]).split("-")[1]
    studyName = (Path(os.path.dirname(os.path.dirname(file_name))).name).split("-")[1]
    
    relatedPhysioData = []
    if os.path.exists(json_file):
        with open(json_file, 'r') as infile:
            content = json.load(infile)
            scanid = str(content["ScanID"]) + ".I32"

        physioPath=os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(file_name)))),'Physio')
        
        conditions = [sub_name , studyName]
        for file in glob.iglob(os.path.join(physioPath, "**", "*" + scanid), recursive=True):
            filename = os.path.basename(file)
            if all(condition in filename for condition in conditions):
                relatedPhysioData.append(file)

    if len(relatedPhysioData)>1:
        sys.exit("Warning: '%s' has no unique physio data for scan %s." % (physioPath, scanid,))
    if len(relatedPhysioData) is 0:
        print("Error: '%s' has no related physio data for scan %s." % (physioPath, scanid,))
        return []

    physioFile_name = relatedPhysioData[0]
    print('Copy related physio data %s to rawMoData' % (physioFile_name,))
    shutil.copyfile(physioFile_name, os.path.join(i32_Path,os.path.basename(physioFile_name)))

    return physioFile_name

def create_txt_file(file, data):
    import time
    with open(file, "w") as outfile:
        for data_point in data:
            outfile.write(''.join([str(data_point), '\n']))

    time.sleep(1.0)

def delete_txt_file(file):
    os.remove(file)

def startProcess(Rawfile_name):
    # generate folder for images
    origin_Path = os.path.dirname(Rawfile_name)
    proc_Path = os.path.join(origin_Path, 'rs-fMRI_niiData')
    if os.path.exists(proc_Path):
        shutil.rmtree(proc_Path)
    os.mkdir(proc_Path)

    # generate folder for motion correction files
    par_Path = os.path.join(origin_Path, 'rs-fMRI_mcf')
    if os.path.exists(par_Path):
        shutil.rmtree(par_Path)
    os.mkdir(par_Path)

    # generate folder for motion correction files
    subFile= os.path.basename(Rawfile_name).split('.')[0]
    subFile = '%s_mcf.mat' % subFile
    par_Path = os.path.join(par_Path,subFile)
    if os.path.exists(par_Path):
        shutil.rmtree(par_Path)
    os.mkdir(par_Path)

    # generate folder for physio data
    i32_Path = os.path.join(origin_Path, 'rawMonData')
    if  os.path.exists(i32_Path):
        shutil.rmtree(i32_Path)
    os.mkdir(i32_Path)

    # bring dataset to RAS orientation
    file_name = getRASorientation(Rawfile_name,proc_Path)

    # calculate EPIMean
    file_nameEPI = getEPIMean(file_name,proc_Path)

    # apply BET on EPImean
    if rodent == 0:
        file_nameEPI_BET,mask_file = applyBET(file_nameEPI,frac=0.35,radius=45,vertical_gradient=0.1)
    elif rodent == 1:
        file_nameEPI_BET,mask_file = applyBET(file_nameEPI,frac=0.26,radius=55,vertical_gradient=0.07)
    #apply Mask on original dataset
    maskedFile_data = applyMask(file_name,mask_file)

    # apply motion correction on original dataset with EPImean as reference
    mcfFile_name=fsl_SeparateSliceMoCo(file_name,par_Path)

    # apply mean on motion corrected data
    meanMcfFile_name = getEPIMean(mcfFile_name, proc_Path)

    # copy physio data to rawMonData-Folder
    relatedPhysioFolder = copyRawPhysioData(Rawfile_name,i32_Path)

    # get Regression Values
    if len(relatedPhysioFolder) is not 0:
        getSingleRegTable.getRegrTable(os.path.dirname(Rawfile_name),relatedPhysioFolder,par_Path)
    else:
        print("Error: Processing not possible, because either there is no folder called Physio or the related physio data for the scan is missing there.")

    return mcfFile_name


#%% Program

#specify default Arguments by defining rodent spezies
define_rodent_spezies()

if __name__ == "__main__":

    TR = 1.42
    cutOff_sec = 100.0
    FWHM = 3.0

    import argparse
    parser = argparse.ArgumentParser(description='Process fMRI data')
    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument('-i', '--input', help='Path to the RAW data of rsfMRI NIfTI file', required=True)

    parser.add_argument('-t', '--TR', default=TR, help='Current TR value')
    parser.add_argument('-c', '--cutOff_sec', default=cutOff_sec, help='High-pass filter cutoff sec')
    parser.add_argument('-f', '--FWHM', default=FWHM, help='Full width at half maximum')
    parser.add_argument('-stc', '--slicetimecorrection', default="False", type=str, help='choose to perform slice time correction or not')

    args = parser.parse_args()

    if args.slicetimecorrection == "True":
        stc = True
    else:
        stc = False
    if rodent == 0:
        labels = os.path.abspath(
            os.path.join(os.getcwd(), os.pardir, os.pardir)) + '/lib/annotation_50CHANGEDanno_label_IDs.txt'
        labelNames = os.path.abspath(os.path.join(os.getcwd(), os.pardir, os.pardir)) + '/lib/annoVolume.nii.txt'
        labels2000 = os.path.abspath(
            os.path.join(os.getcwd(), os.pardir, os.pardir)) + '/lib/annotation_50CHANGEDanno_label_IDs+2000.txt'
        labelNames2000 = os.path.abspath(
            os.path.join(os.getcwd(), os.pardir, os.pardir)) + '/lib/annoVolume+2000_rsfMRI.nii.txt'
    elif rodent == 1:
        labels = os.path.abspath(
            os.path.join(os.getcwd(), os.pardir, os.pardir)) + '/lib/SIGMA_InVivo_Anatomical_Brain_Atlas_Labels_fMRI.txt'
        labelNames = os.path.abspath(os.path.join(os.getcwd(), os.pardir, os.pardir)) + '/lib/SIGMA_InVivo_Anatomical_Brain_Atlas_Labels.txt'
        labels2000 = os.path.abspath(
            os.path.join(os.getcwd(), os.pardir, os.pardir)) + '/lib/SIGMA_InVivo_Anatomical_Brain_Atlas_Labels_fMRI.txt'
        labelNames2000 = os.path.abspath(
            os.path.join(os.getcwd(), os.pardir, os.pardir)) + '/lib/SIGMA_InVivo_Anatomical_Brain_Atlas_Labels.txt'    
    input_file = None
    if args.input is not None and args.input is not None:
        input_file = args.input
    if not os.path.exists(input_file):
        sys.exit("Error: '%s' is not an existing directory or file %s is not in directory." % (input_file, args.file,))

    mcfFile_name = startProcess(input_file)

    
    # if stc is activated find parameters
    if stc: 
        print("Starting Regression with slice time correction:")
        # find meta data json file
        meta_data_file_name = Path(args.input).name.replace(".nii.gz", ".json")
        meta_data_file = os.path.join(Path(args.input).parent, meta_data_file_name)

        with open(meta_data_file, "r") as infile:
            meta_data = json.load(infile)
            
        TR = meta_data["RepetitionTime"] / 1000
        slice_order = meta_data["ObjOrderList"]
        n_slices = meta_data["n_slices"]
        costum_timings = meta_data["costum_timings"]

        # create costum timings txt file
        costum_timings_path = os.path.join(Path(meta_data_file).parent, "tcostum.txt")
        create_txt_file(costum_timings_path, costum_timings)
        
        # create slice order txt file
        slice_order_path = os.path.join(Path(meta_data_file).parent, "slice_order.txt")
        create_txt_file(slice_order_path, slice_order)
        
        rgr_file, srgr_file, sfrgr_file = regress.startRegression(mcfFile_name, FWHM, cutOff_sec, TR, stc, slice_order_path, costum_timings_path)

        # delete temp txt files
        delete_txt_file(costum_timings_path)
        delete_txt_file(slice_order_path)

    else:
        print("Starting Regression without slice time correction:")
        rgr_file, srgr_file, sfrgr_file = regress.startRegression(mcfFile_name, FWHM, cutOff_sec, TR, stc)
        print(f"sfrgr_file {sfrgr_file}")

    

    atlasPath = os.path.dirname(input_file)
    roisPath = copyAtlasOfData(atlasPath,'Anno_parental',labels)

    fslMeantsFile = fsl_mean_ts.start_fsl_mean_ts(sfrgr_file, roisPath, labelNames, 'MasksTCs.')

    roisPath = copyAtlasOfData(atlasPath, 'AnnoSplit_parental', labels2000)

    fslMeantsFile = fsl_mean_ts.start_fsl_mean_ts(sfrgr_file, roisPath, labelNames2000, 'MasksTCsSplit.')
