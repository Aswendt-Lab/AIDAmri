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

def copyAtlasOfData(path,post,labels):
    fileALL = glob.glob(path + '/*' + post + '.nii.gz')
    if fileALL.__len__()>1:
        sys.exit("Error: '%s' has no related Atlas File." % (path,))
    else:
        fileALL = fileALL[0]

    print("Copy Atlas Data an generate seed ROIs")
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
    imgTemp = data.get_data()
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
    file_data = nii.load(file_name)
    imgData = file_data.get_data()
    imgData = np.flip(imgData, 2)

    imgData = np.flip(imgData, 0)
    y = imgData
    #y = imgScaleResize(imgData)


    epiData = nii.Nifti1Image(y, file_data.affine)
    hdrIn = epiData.header
    hdrIn.set_xyzt_units('mm')
    epiData_RAS = nii.as_closest_canonical(epiData)
    print('Orientation:' + str(nii.aff2axcodes(epiData_RAS.affine)))
    output_file = os.path.join(proc_Path, os.path.basename(file_name))
    nii.save(epiData, output_file)
    return output_file

def getEPIMean(file_name,proc_Path):
    output_file = os.path.join(proc_Path, os.path.basename(file_name).split('.')[0]) + 'EPI.nii.gz'
    myMean = fsl.MeanImage(in_file=file_name, out_file=output_file)
    print(myMean.cmdline)
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
    fslPath = scaleBy10(input_file, inv=False)
    mySplit= fsl.Split(in_file=fslPath,dimension='z',out_base_name = dataName)
    print(mySplit.cmdline)
    mySplit.run()
    os.remove(fslPath)


    # sparate ref and src volume in slices
    sliceFiles = findSlicesData(os.getcwd(),dataName)
    # refFiles = findSlicesData(os.getcwd(),'ref')
    print('For all slices ... ')


    #start to correct motions slice by slice
    for i in  range(len(sliceFiles)):
        slc = sliceFiles[i]
        # ref = refFiles[i]
        # take epi as ref
        output_file = os.path.join(par_folder,os.path.basename(slc))
        myMCFLIRT = fsl.preprocess.MCFLIRT(in_file=slc,out_file=output_file,save_plots=True,terminal_output='none')
        print(myMCFLIRT.cmdline)
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

    return output_file

def copyRawPhysioData(file_name,i32_Path):
    studyName = file_name.split('/')[-3]
    scanNo = file_name.split('.')[-4]+'.I32'
    physioPath=os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(file_name))),'Physio',studyName)

    relatedPhysioData = []
    fileALL = glob.glob(physioPath+'/'+studyName+'*'+scanNo)
    for filename in fileALL:
        relatedPhysioData.append(filename)

    if len(relatedPhysioData)>1:
        sys.exit("Error: '%s' has no unique physio data for scan %s." % (physioPath, scanNo,))
    if len(relatedPhysioData) is 0:
        print("Warning: '%s' has no related physio data for scan %s." % (physioPath, scanNo,))
        return []

    physioFile_name = relatedPhysioData[0]
    print('Copy related physio data %s to rawMoData' % (physioFile_name,))
    shutil.copyfile(physioFile_name, os.path.join(i32_Path,os.path.basename(physioFile_name)))

    return physioFile_name

def startProcess(Rawfile_name):
    # generate folder for images
    origin_Path = os.path.dirname(Rawfile_name)
    proc_Path = os.path.join(origin_Path, 'rs-fMRI_niiData')
    if os.path.exists(proc_Path):
        shutil.rmtree(proc_Path)
    os.mkdir(proc_Path)

    # generate folder for motion corrction files
    par_Path = os.path.join(origin_Path, 'rs-fMRI_mcf')
    if os.path.exists(par_Path):
        shutil.rmtree(par_Path)
    os.mkdir(par_Path)

    # generate folder for motion corrction files
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

    print("fMRI Processing \33[5m...\33[0m (wait!)", end="\r")

    # generate log - file
    sys.stdout = open(os.path.join(os.path.dirname(Rawfile_name), 'process.log'), 'w')

    # bring dataset to RAS orientation
    file_name = getRASorientation(Rawfile_name,proc_Path)

    # calculate EPIMean
    file_nameEPI = getEPIMean(file_name,proc_Path)

    # apply BET on EPImean
    file_nameEPI_BET,mask_file = applyBET(file_nameEPI,frac=0.35,radius=45,vertical_gradient=0.1)

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

    sys.stdout = sys.__stdout__
    print('fMRI Processing  \033[0;30;42m COMPLETED \33[0m')

    return mcfFile_name

if __name__ == "__main__":


    import argparse
    parser = argparse.ArgumentParser(description='Process fMRI data')

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-i','--input', help='file name of data',required=True)
    parser.add_argument('-l','--labels',help='labels for fMRI',
                        default=os.path.abspath(os.path.join(os.getcwd(), os.pardir,os.pardir))+'/lib/annotation_50CHANGEDanno_label_IDs.txt')
    parser.add_argument('-n', '--labelNames', help='labels names for fMRI',
                        default=os.path.abspath(os.path.join(os.getcwd(), os.pardir,os.pardir))+'/lib/annoVolume.nii.txt')
    parser.add_argument('-s', '--labels2000', help='labels for fMRI',
                        default=os.path.abspath(os.path.join(os.getcwd(), os.pardir,os.pardir)) + '/lib/annotation_50CHANGEDanno_label_IDs+2000.txt')
    parser.add_argument('-t', '--labelNames2000', help='labels names for fMRI',
                        default=os.path.abspath(os.path.join(os.getcwd(), os.pardir, os.pardir)) + '/lib/annoVolume+2000_rsfMRI.nii.txt')
    args = parser.parse_args()

    TR = 2.84
    cutOff_sec = 100
    FWHM = 3.0

    labels = args.labels
    labelNames = args.labelNames
    labels2000 = args.labels2000
    labelNames2000 = args.labelNames2000
    input = None
    if args.input is not None and args.input is not None:
        input = args.input
    if not os.path.exists(input):
        sys.exit("Error: '%s' is not an existing directory of file %s is not in directory." % (input, args.file,))

    mcfFile_name = startProcess(input)
    rgr_file, srgr_file, sfrgr_file = regress.startRegression(mcfFile_name, FWHM, cutOff_sec, TR)


    atlasPath = os.path.dirname(input)
    roisPath = copyAtlasOfData(atlasPath,'Anno_rsfMRI',labels)

    fslMeantsFile = fsl_mean_ts.start_fsl_mean_ts(sfrgr_file, roisPath, labelNames, 'MasksTCs.')

    roisPath = copyAtlasOfData(atlasPath, 'AnnoSplit_rsfMRI', labels2000)

    fslMeantsFile = fsl_mean_ts.start_fsl_mean_ts(sfrgr_file, roisPath, labelNames2000, 'MasksTCsSplit.')
