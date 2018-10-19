"""
Created on 10/08/2017

@author: Niklas Pallast
Neuroimaging & Neuroengineering
Department of Neurology
University Hospital Cologne

"""

import sys,os
import nibabel as nii
import numpy as np
import nipype.interfaces.fsl as fsl
import glob
import shutil



def scaleBy10(input_path,inv):
    data = nii.load(input_path)
    imgTemp = data.get_data()
    if inv == False:
        scale = np.eye(4) * 10
        scale[3][3] = 1
        scaledNiiData = nii.Nifti1Image(imgTemp, data.affine * scale)
        fslPath = os.path.join(os.path.dirname(input_path), os.path.basename(input_path).split('.')[0]+'_fslScaleTemp.nii.gz')
        nii.save(scaledNiiData, fslPath)
        return fslPath
    elif inv == True:
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

def findRegData(path):
    regMR_list = []
    fileALL = glob.iglob(path+'/*.txt',recursive=True)
    for filename in fileALL:
        regMR_list.append(filename)
    regMR_list.sort()
    return regMR_list

def findSlicesData(path,pre):
    regMR_list = []
    fileALL = glob.iglob(path+'/'+pre+'*.nii.gz',recursive=True)
    for filename in fileALL:
        regMR_list.append(filename)
    regMR_list.sort()
    return regMR_list

def delete5Slides(input_file,regr_Path):
    # scale Nifti data by factor 10
    fslPath = scaleBy10(input_file, inv=False)
    # delete 5 slides
    output_file = os.path.join(os.path.dirname(input_file), os.path.basename(input_file).split('.')[0]) + '_f.nii.gz'
    myROI = fsl.ExtractROI(in_file=fslPath, roi_file=output_file, t_min=5, t_size=-1)
    print(myROI.cmdline)
    myROI.run()
    os.remove(fslPath)
    # unscale result data by factor 10ˆ(-1)
    output_file = scaleBy10(output_file, inv=True)
    return output_file

def fsl_RegrSliceWise(input_file,txtregr_Path,regr_Path):
    # scale Nifti data by factor 10
    dataName = os.path.basename(input_file).split('.')[0]

    # proof  data existence
    regrTextFiles = findRegData(txtregr_Path)
    if len(regrTextFiles) == 0:
        print('No regression with physio data!')
        output_file = os.path.join(regr_Path,
                                   os.path.basename(input_file).split('.')[0]) + '_RGR.nii.gz'
        shutil.copyfile(input_file, output_file)
        return output_file


    fslPath = scaleBy10(input_file, inv=False)
    # split input_file in slices
    mySplit = fsl.Split(in_file=fslPath, dimension='z', out_base_name=dataName)
    print(mySplit.cmdline)
    mySplit.run()
    os.remove(fslPath)

    # sparate ref and src volume in slices
    sliceFiles = findSlicesData(os.getcwd(), dataName)




    if not len(regrTextFiles) == len(sliceFiles):
        sys.exit('Error: Not enough txt.Files in %s' % txtregr_Path)

    print('Start separate slice Regression ... ')

    # start to regression slice by slice
    print('For all Sices ...')
    for i in range(len(sliceFiles)):
        slc = sliceFiles[i]
        regr = regrTextFiles[i]
        # only take the columns [1,2,7,9,11,12,13] of the reg-.txt Files
        output_file = os.path.join(regr_Path, os.path.basename(slc))
        myRegr = fsl.FilterRegressor(in_file=slc,design_file=regr,out_file=output_file,filter_columns=[1,2,7,9,11,12,13])
        print(myRegr.cmdline)
        myRegr.run()
        os.remove(slc)


    # merge slices to a single volume
    mcf_sliceFiles = findSlicesData(regr_Path, dataName)
    output_file = os.path.join(regr_Path,
                               os.path.basename(input_file).split('.')[0]) + '_RGR.nii.gz'
    myMerge = fsl.Merge(in_files=mcf_sliceFiles, dimension='z', merged_file=output_file)
    print(myMerge.cmdline)
    myMerge.run()

    for slc in mcf_sliceFiles: os.remove(slc)

    # unscale result data by factor 10ˆ(-1)
    output_file = scaleBy10(output_file, inv=True)

    return output_file

def getMask(input_file,threshold):
    threshold = threshold/10
    output_file = os.path.join(os.path.dirname(input_file), 'mask.nii.gz')
    thres = fsl.Threshold(in_file=input_file,thresh=threshold,out_file=output_file,output_datatype='char',args='-Tmin -bin')
    print(thres.cmdline)
    thres.run()
    return output_file

def dilF(input_file):
    output_file = os.path.join(os.path.dirname(input_file), 'mask.nii.gz')
    mydilf = fsl.DilateImage(in_file=input_file, operation='max', out_file=output_file)
    print(mydilf.cmdline)
    mydilf.run()
    return output_file

def applyMask(input_file,mask_file,appendix):
    output_file = os.path.join(os.path.dirname(input_file), os.path.basename(input_file).split('.')[0]) + appendix+ '.nii.gz'
    myMaskapply = fsl.ApplyMask(in_file=input_file, out_file=output_file, mask_file=mask_file)
    print(myMaskapply.cmdline)
    myMaskapply.run()
    return output_file

def getMean(input_file,appendix):
    output_file = os.path.join(os.path.dirname(input_file), appendix+'.nii.gz')
    myMean = fsl.MeanImage(in_file=input_file, out_file=output_file)
    print(myMean.cmdline)
    myMean.run()
    return output_file

def applySusan(input_file,meanintensity,FWHM,mean_func):
    # scale Nifti data by factor 10
    fslPath = scaleBy10(input_file, inv=False)
    meanPath = scaleBy10(mean_func, inv=False)
    output_file = os.path.join(os.path.dirname(input_file),
                               os.path.basename(input_file).split('RGR')[0]) + 'SRGR.nii.gz'

    mySusan = fsl.SUSAN(in_file=fslPath, brightness_threshold=meanintensity, fwhm=FWHM, dimension=2,
                        use_median=1, usans=[(meanPath, meanintensity), ], out_file=output_file)
    print(mySusan.cmdline)
    mySusan.run()
    os.remove(fslPath)
    os.remove(meanPath)
    output_file = scaleBy10(output_file, inv=True)
    return  output_file

def mathOperation(input_file,scale_factor):
    output_file = os.path.join(os.path.dirname(input_file),  os.path.basename(input_file).split('.')[0])+'_intnorm.nii.gz'
    myMath = fsl.BinaryMaths(in_file=input_file,operand_value =scale_factor,operation='mul',out_file=output_file)
    print(myMath.cmdline)
    myMath.run()
    return output_file

def filterFSL(input_file,highpass,tempMean):
    outputSFRGR = os.path.join(os.path.dirname(input_file), os.path.basename(input_file).split('SRGR')[0])+'SFRGR.nii.gz'
    myHP = fsl.TemporalFilter(in_file = input_file,highpass_sigma=highpass, args='-add '+tempMean,out_file=outputSFRGR)
    print(myHP.cmdline)
    myHP.run()
    input_file = outputSFRGR
    output_file = os.path.join(os.path.dirname(input_file),  os.path.basename(input_file).split('.')[0])+'_thres_mask.nii.gz'
    #input_file = getMean(input_file,'HPmean')
    thres = fsl.Threshold(in_file=input_file, thresh=17, out_file=output_file, output_datatype='float',use_robust_range=True,args='-Tmean -bin')
    print(thres.cmdline)
    thres.run()

    return outputSFRGR


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
    maskFile = scaleBy10(maskFile,inv=True)
    return output_file,maskFile


def startRegression(input_File):
    # generate folder regr images
    print("Regression \33[5m...\33[0m (wait!)", end="\r")
    origin_Path = os.path.dirname(os.path.dirname(input_File))
    regr_Path = os.path.join(origin_Path, 'regr')

    if os.path.exists(regr_Path):
        shutil.rmtree(regr_Path)
    os.mkdir(regr_Path)

    # generatre log-File
    sys.stdout = open(os.path.join(regr_Path,'regress.log'),'w')

    # delete the first slides
    input_File5Sub = delete5Slides(input_File, regr_Path)

    # proof regression files
    txtregr_Path = os.path.join(origin_Path, 'txtRegrPython')

    # slive wise regression with physio data
    regr_FileReal = fsl_RegrSliceWise(input_File5Sub, txtregr_Path, regr_Path)

    # get mean
    meanRegr_File = getMean(regr_FileReal,'mean2')
    file_nameEPI_BET, mask_file = applyBET(meanRegr_File, frac=0.35, radius=45, vertical_gradient=0.1)
    os.remove(meanRegr_File)
    regr_File = applyMask(regr_FileReal,mask_file,'')


    #  "robust intensity range" which calculates values similar to the 98% percentiles
    myStat = fsl.ImageStats(in_file=regr_File,op_string='-p 98',terminal_output='allatonce')
    print(myStat.cmdline)
    stat_result = myStat.run()
    upperp = stat_result.outputs.out_stat

    # get binary mask
    mask = getMask(regr_File, upperp)

    # "robust intensity range" which calculates values similar to the 50% percentiles with mask
    myStat = fsl.ImageStats(in_file=regr_File, op_string=' -k ' +mask+ ' -p 50 ',mask_file=mask ,terminal_output='allatonce')
    print(myStat.cmdline)
    stat_result = myStat.run()
    meanintensity = stat_result.outputs.out_stat
    meanintensity = meanintensity*0.75

    # maxmium filter of mask
    mask = dilF(mask)

    # apply mask on regrFile
    thresRegr_file = applyMask(regr_File,mask,'thres')

    # get mean of masked regr-Dataset
    mean_func = getMean(thresRegr_file,'mean_func')

    FWHM = 3.0
    #sigma = FWHM/(2 * np.sqrt(2 * np.log(2)))
    srgr_file = applySusan(thresRegr_file,meanintensity,FWHM,mean_func)

    # apply mask on srgr_file
    smmothSRegr_file = applyMask(srgr_file,mask,'_smooth')
    inscalefactor = 10000.0/meanintensity


    # multiply image with inscalefactor
    intnormSrgr_file = mathOperation(smmothSRegr_file,inscalefactor)

    # mean of scaled Dataset
    tempMean  =  getMean(intnormSrgr_file,'tempMean')

    # filter image
    highpass = 17.6056338028
    #lowpass = 2.20070422535
    filtered_image = filterFSL(intnormSrgr_file,highpass,tempMean)

    sys.stdout = sys.__stdout__
    print('Regression  \033[0;30;42m COMPLETED \33[0m')
    return regr_FileReal, srgr_file ,filtered_image

if __name__ == "__main__":


    import argparse
    parser = argparse.ArgumentParser(description='Regression of fMRI data')

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-i','--input', help='file name of data',required=True)
    args = parser.parse_args()


    if args.input is not None and args.input is not None:
        input = args.input
    if not os.path.exists(input):
        sys.exit("Error: '%s' is not an existing directory of file %s is not in directory." % (input, args.file,))

    result = startRegression(input)
