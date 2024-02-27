
"""
Created on 10/08/2017

@author: Niklas Pallast
Neuroimaging & Neuroengineering
Department of Neurology
University Hospital Cologne


"""


import os,sys
import nibabel as nii
import glob
import numpy as np
import scipy.io as sc
import scipy.ndimage as ndimage


def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

def thresholdingSlc(volumeMR,maskImg,thres):
    volumeMR=ndimage.gaussian_filter(volumeMR, sigma=(1.2, 1.2, 1))

    volumeMR = volumeMR * maskImg[:, :, :, 0]

    scaledNiiData = nii.Nifti1Image(volumeMR, np.eye(4))
    hdrIn = scaledNiiData.header
    hdrIn.set_xyzt_units('mm')
    output_file = os.path.join(outfile,'maskedVolume.nii.gz')
    nii.save(scaledNiiData, output_file)

    for i in range(len(volumeMR[1,1,:])-1):
        voSlc = volumeMR[:,:,i]

        uvalues = voSlc >= thres
        voSlc[uvalues] = 0
        zvalues = voSlc != 0
        thresF = np.mean(voSlc[zvalues]) + 1.5*np.std(voSlc[zvalues])

        bvalues = voSlc < thresF
        voSlc[bvalues] = 0
        voSlc[bvalues] = 0
        fvalues = voSlc >= thresF
        voSlc[fvalues] = 1
        volumeMR[:, :, i] = voSlc

    fvalues = volumeMR == 1
    return volumeMR,fvalues



def thresholding(volumeMR,maskImg,thres,k):
    volumeMR=ndimage.gaussian_filter(volumeMR, sigma=(1.3, 1.3, 1))
    zvalues = volumeMR != 0



    if k==1:
        volumeMR = volumeMR * maskImg[:, :, :]#, 0]

    if thres == 0:
        thres = np.mean(volumeMR[zvalues]) + 2*np.std(volumeMR[zvalues])

    bvalues = volumeMR < thres
    volumeMR[bvalues] = 0

    fvalues = volumeMR >= thres
    volumeMR[fvalues] = 1

    return volumeMR


def incidenceMap(path_listInc,path_listMR ,path_listAnno, araDataTemplate,incidenceMask ,thres, outfile,labels):

    araDataTemplate  = nii.load(araDataTemplate)
    realAraImg = araDataTemplate.get_data()
    coloredAraLabels = np.zeros([np.size(realAraImg, 0), np.size(realAraImg, 1), np.size(realAraImg, 2)])

    matFile = sc.loadmat(labels)
    labMat = matFile['ABLAbelsIDsParental']

    maskData = nii.load(incidenceMask)
    maskImg = maskData.get_data()
    oneValues = maskImg > 0.0
    maskImg[oneValues] = 1.0
    fileIndex = 0

    # get warped annos of the current mr
    dataAnno = nii.load(path_listAnno[fileIndex])
    volumeAnno = np.round(dataAnno.get_data())
    dataMR = nii.load(path_listInc[fileIndex])
    volumeMR = dataMR.get_data()

    strokeVolume = thresholding(volumeMR, maskImg, thres,1)

    fValues_Anno = volumeAnno*strokeVolume

    scaledNiiData = nii.Nifti1Image(fValues_Anno, dataAnno.affine)
    hdrIn = scaledNiiData.header
    hdrIn.set_xyzt_units('mm')
    output_file =  os.path.join(outfile,os.path.basename(path_listMR[fileIndex]).split('.')[0]+ 'Anno_parmask.nii.gz')
    nii.save(scaledNiiData, output_file)

    ref_Image = fValues_Anno
    fValues_Anno = np.unique(fValues_Anno)
    nullValues = np.argwhere(fValues_Anno<=0.0)
    fValues_Anno = np.delete(fValues_Anno, nullValues)

    regionAffectPercent = np.zeros(np.size(fValues_Anno))
    for i in range(np.size(fValues_Anno)):
        regionAffectPercent[i] = (np.sum(ref_Image == fValues_Anno[i]) / np.sum(volumeAnno == fValues_Anno[i])) * 200
        regionAffectPercent[regionAffectPercent > 100] = 100
    labCounterList = np.isin(labMat[:, 0], fValues_Anno)
    labMat = labMat[labCounterList,0]
    labCounterColor = np.isin(realAraImg, fValues_Anno)
    coloredAraLabels[labCounterColor] = realAraImg[labCounterColor]
    xdim = np.size(coloredAraLabels, 0)
    coloredAraLabels[int(xdim / 2):xdim, :, :] = coloredAraLabels[int(xdim / 2):xdim, :, :] + 2000
    coloredAraLabels[coloredAraLabels == 2000] = 0
    scaledNiiData = nii.Nifti1Image(coloredAraLabels, araDataTemplate.affine)
    hdrIn = scaledNiiData.header
    hdrIn.set_xyzt_units('mm')
    output_file = os.path.join(outfile, 'affectedRegions_Parental.nii.gz')
    nii.save(scaledNiiData, output_file)

    # Stroke volume calculation
    betMask = nii.load(os.path.join(outfile,os.path.basename(path_listInc[fileIndex]).split('.')[0]+'_mask.nii.gz'))
    betMaskImg = betMask.get_data()
    oneValues = betMaskImg > 0.0
    betMaskImg[oneValues] = 1.0
    strokeVolumeInCubicMM = np.sum(maskImg * (dataMR.affine[0, 0] * dataMR.affine[1, 1] * dataMR.affine[2, 2]))
    brainVolumeInCubicMM = np.sum(betMaskImg * (dataMR.affine[0, 0] * dataMR.affine[1, 1] * dataMR.affine[2, 2]))

    lines =open(os.path.abspath(os.path.join(os.getcwd(), os.pardir,os.pardir))+ '/lib/annoVolume.nii.txt').readlines()
    o=open(os.path.join(outfile, 'affectedRegions_Parental.txt'), 'w')
    o.write("Stroke: %0.2f %% - Stroke Volume: %0.2f mm^3\n"  % (((strokeVolumeInCubicMM/brainVolumeInCubicMM)*100),strokeVolumeInCubicMM,))
    matIndex = 0
    labelNamesAffected = ["" for x in range(np.size(fValues_Anno))]
    labelNames = ["" for x in range(np.size(lines))]
    for i in range(len(lines)):
        labelNames[i] = lines[i].split('\t')[1]
        if  np.isin(int(lines[i].split('\t')[0]),labMat):
            o.write(lines[i][:-1] + "\t %0.2f %%\n" % regionAffectPercent[matIndex])
            labelNamesAffected[matIndex] = lines[i].split('\t')[1]
            matIndex = matIndex + 1

            #o.write(str(int(lines[i].split('	')[0]) + 2000) + '	R_' + lines[i].split('	')[1])
    o.close()
    rows = np.shape(labMat)[0]
    labMat = np.stack((labMat, regionAffectPercent[0:rows]))
    matFile['ABLAbelsIDsParental'] = labMat
    matFile['ABANamesPar'] = labelNamesAffected
    matFile['ABAlabels'] = labelNames
    matFile['volumePer'] = (strokeVolumeInCubicMM / brainVolumeInCubicMM) * 100
    matFile['volumeMM'] = strokeVolumeInCubicMM
    sc.savemat(os.path.join(outfile, 'labelCount_par.mat'), matFile)






def findIncData(path):
    regMR_list = []

    for filename in glob.iglob(path+'*/*IncidenceData.nii.gz', recursive=False):
        regMR_list.append(filename)

    return regMR_list

def findBETData(path):
    regMR_list = []

    for filename in glob.iglob(path+'*/*Bet.nii.gz', recursive=False):
        regMR_list.append(filename)

    return regMR_list



def findRegisteredData(path):
    regMR_list = []

    for filename in glob.iglob(path+'*/*_Template.nii.gz', recursive=True):
        regMR_list.append(filename)

    return regMR_list

def findRegisteredAnno(path):
    regANNO_list = []

    for filename in glob.iglob(path + '*/*_AnnoSplit_parental.nii.gz', recursive=True):
        regANNO_list.append(filename)

    return regANNO_list

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Calculate incidence sizes of parental regions. You do not need to enter single files, but the path to the .../T2w folder')
    requiredNamed = parser.add_argument_group('Required named arguments')
    requiredNamed.add_argument('-i', '--inputFolder', help='.../T2w')

    parser.add_argument('-t', '--threshold', help='Threshold for stroke values ',  nargs='?', type=int,
                        default=0)
    parser.add_argument('-a', '--allenBrain_anno', help='File: Annotations of Allen Brain', nargs='?', type=str,
                        default=os.path.abspath(
                            os.path.join(os.getcwd(), os.pardir, os.pardir)) + '/lib/average_template_50.nii.gz')

    inputFolder = None
    allenBrain_template = None
    allenBrain_anno = None
    outfile = None

    args = parser.parse_args()

    if args.inputFolder is not None:
        inputFolder = args.inputFolder
        outfile = args.inputFolder
    if not os.path.exists(inputFolder):
        sys.exit("Error: '%s' is not an existing directory." % (inputFolder,))


    if args.allenBrain_anno is not None:
        allenBrain_anno = args.allenBrain_anno
    if not os.path.exists(allenBrain_anno):
        sys.exit("Error: '%s' is not an existing directory." % (allenBrain_anno,))

    thres = args.threshold
    labels = os.path.abspath(os.path.join(os.getcwd(), os.pardir,os.pardir))+ '/lib/rsfMRILablelID.mat'
    araDataTemplate = os.path.abspath(os.path.join(os.getcwd(), os.pardir,os.pardir))+ '/lib/annoVolume.nii.gz'

    if len(glob.glob(inputFolder+'/*Stroke_mask.nii.gz')) > 0:
        incidenceMask = glob.glob(inputFolder+'/*Stroke_mask.nii.gz')[0]
    else:
        sys.exit("Error: '%s' has no affected or masked regions." % (inputFolder,))

    path = os.path.join(inputFolder)
    regMR_list = findBETData(path)
    regInc_list = findIncData(path)
    regANNO_list = findRegisteredAnno(path)


    print("'%i' folder will be proccessed..." % (len(regMR_list),))

    if not len(regANNO_list) == len(regMR_list):
        sys.exit("Error: For one or more annotations is no corresponding MR file defined in '%s'." % (inputFolder,))

    incidenceMap(regMR_list,regInc_list,regANNO_list,araDataTemplate,incidenceMask,thres,outfile,labels)
