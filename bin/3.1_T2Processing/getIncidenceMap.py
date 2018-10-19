'''
Created on 10/08/2017

@author: Niklas Pallast
Neuroimaging & Neuroengineering
Department of Neurology
University Hospital Cologne

'''

import os,sys
import nibabel as nii
import glob
import numpy as np
import scipy.io as sc
import scipy.ndimage as ndimage
import progressbar
import  getHeatmap as gHm




def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]


def thresholdingSlc(volumeMR,maskImg,thres):
    volumeMR=ndimage.gaussian_filter(volumeMR, sigma=(1.2, 1.2, 1))
    zvalues = volumeMR != 0

    volumeMR = volumeMR * maskImg[:, :, :, 0]

    scaledNiiData = nii.Nifti1Image(volumeMR, np.eye(4))
    hdrIn = scaledNiiData.header
    hdrIn.set_xyzt_units('mm')
    output_file = os.path.join(outfile,'maskedVolume.nii.gz')
    nii.save(scaledNiiData, output_file)

    for i in range(len(volumeMR[1,1,:])-1):
        voSlc = volumeMR[:,:,i]
        zvalues = voSlc != 0
        #if thres == 0:
        thresU = np.mean(voSlc[zvalues]) +2.5*np.std(voSlc[zvalues])
        #print("U: '%f'." % (thresU,))
        uvalues = voSlc >= thres
        voSlc[uvalues] = 0
        zvalues = voSlc != 0
        thresF = np.mean(voSlc[zvalues]) + 1.5*np.std(voSlc[zvalues])
        #print("F: '%f'." % (thresF,))

        bvalues = voSlc < thresF
        voSlc[bvalues] = 0
        voSlc[bvalues] = 0
        fvalues = voSlc >= thresF
        voSlc[fvalues] = 1
        volumeMR[:, :, i] = voSlc

    fvalues = volumeMR == 1
    return volumeMR,fvalues

def thresholding(volumeMR,maskImg,thres):
    volumeMR=ndimage.gaussian_filter(volumeMR, sigma=(1.3, 1.3, 1))
    zvalues = volumeMR != 0

    #volumeMR = volumeMR * maskImg[:, :, :, 0]

    # scaledNiiData = nii.Nifti1Image(volumeMR, np.eye(4))
    # hdrIn = scaledNiiData.header
    # hdrIn.set_xyzt_units('mm')
    # output_file = os.path.join(outfile,'maskedVolume.nii.gz')
    # nii.save(scaledNiiData, output_file)

    if thres == 0:
        thres = np.mean(volumeMR[zvalues]) + 2*np.std(volumeMR[zvalues])

    bvalues = volumeMR < thres
    volumeMR[bvalues] = 0

    fvalues = volumeMR >= thres
    volumeMR[fvalues] = 1

    fvalues = volumeMR == 1
    return volumeMR,fvalues


def incidenceMap(path_listInc,path_listMR ,path_listAnno, araDataTemplate,incidenceMask ,thres, outfile, araAnno,labels):

    araData = nii.load(araAnno)
    araVol = araData.get_data()
    incidenceMap = np.zeros([np.size(araVol, 0), np.size(araVol, 1), np.size(araVol, 2)])

    araDataTemplate  = nii.load(araDataTemplate)
    realAraImg = araDataTemplate.get_data()
    coloredAraLabels = np.zeros([np.size(realAraImg, 0), np.size(realAraImg, 1), np.size(realAraImg, 2)])



    matFile = sc.loadmat(labels)
    labMat = matFile['ABALabelIDs']

    maskData = nii.load(incidenceMask)
    maskImg = maskData.get_data()
    oneValues = maskImg > 0.0
    maskImg[oneValues] = 1.0

    bar = progressbar.ProgressBar()
    for fileIndex in bar(range(len(path_listMR))):
        #print(path_listMR[fileIndex])
        dataMR = nii.load(path_listMR[fileIndex])
        volumeMR = dataMR.get_data()

        # Thresholding
        [strokeVolume,fvalues] = thresholding(volumeMR,maskImg,thres)

        incidenceMap = incidenceMap + strokeVolume

        scaledNiiData = nii.Nifti1Image(strokeVolume, dataMR.affine)
        hdrIn = scaledNiiData.header
        hdrIn.set_xyzt_units('mm')
        output_file = os.path.join(outfile,os.path.basename(path_listMR[fileIndex]).split('.')[0]+ 'StrokeSeg.nii.gz')
        nii.save(scaledNiiData, output_file)


        # get warped annos of the current mr
        dataAnno = nii.load(path_listAnno[fileIndex])
        volumeAnno = np.round(dataAnno.get_data())
        dataMR = nii.load(path_listInc[fileIndex])
        volumeMR = dataMR.get_data()

        [strokeVolume, fvalues] = thresholding(volumeMR, maskImg, thres)

        fValues_Anno = volumeAnno*strokeVolume

        scaledNiiData = nii.Nifti1Image(fValues_Anno, dataAnno.affine)
        hdrIn = scaledNiiData.header
        hdrIn.set_xyzt_units('mm')
        output_file =  os.path.join(outfile,os.path.basename(path_listMR[fileIndex]).split('.')[0]+ 'AnnoValues.nii.gz')
        nii.save(scaledNiiData, output_file)


        fValues_Anno = np.unique(fValues_Anno)
        nullValues = np.argwhere(fValues_Anno<=0.0)
        fValues_Anno = np.delete(fValues_Anno, nullValues)

        labCounterList = np.isin(labMat[:, 0], fValues_Anno)
        labMat[labCounterList,1] = labMat[labCounterList,1]+1

        labCounterColor = np.isin(realAraImg, fValues_Anno)
        coloredAraLabels[labCounterColor] = coloredAraLabels[labCounterColor]+1


    matFile['ABALabelIDs'] = labMat
    sc.savemat(os.path.join(outfile,'labelCount.mat'),matFile)


    scaledNiiData = nii.Nifti1Image(incidenceMap, araData.affine)
    hdrIn = scaledNiiData.header
    hdrIn.set_xyzt_units('mm')
    output_file = os.path.join(outfile,'incidenceMap.nii.gz')
    nii.save(scaledNiiData, output_file)

    gHm.heatMap(incidenceMap=incidenceMap,araVol=araVol)

    scaledNiiData = nii.Nifti1Image(coloredAraLabels, araData.affine)
    hdrIn = scaledNiiData.header
    hdrIn.set_xyzt_units('mm')
    output_file = os.path.join(outfile, 'coloredAnnotationMap.nii.gz')
    nii.save(scaledNiiData, output_file)


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

    for filename in glob.iglob(path+'*/*linTransInput.nii.gz', recursive=True):
        regMR_list.append(filename)

    return regMR_list

def findRegisteredAnno(path):
    regANNO_list = []

    for filename in glob.iglob(path + '*/*transAnno.nii.gz', recursive=True):
        regANNO_list.append(filename)

    return regANNO_list

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Calculate an Incidence Map')
    parser.add_argument('-i','--inputFile', help='file name:Brain extracted input data')
    parser.add_argument('-s', '--studyname', help='prefix of the study in the input folder - for exmpale MA_30584')
    parser.add_argument('-t', '--threshold', help='threshold for stroke values ',  nargs='?', type=int,
                        default=0)
    parser.add_argument('-o','--outfile', help='file name')
    parser.add_argument('-a', '--allenBrain_anno', help='file name:Annotations of Allen Brain', nargs='?', type=str,
                        default=os.path.abspath(os.path.join(os.getcwd(), os.pardir, os.pardir,
                                                             os.pardir)) + '/lib/average_template_50.nii.gz')



    args = parser.parse_args()


    if args.inputFile is not None:
        inputFile = args.inputFile
    if not os.path.exists(inputFile):
        sys.exit("Error: '%s' is not an existing directory." % (inputFile,))

    if args.outfile is not None:
        outfile = args.outfile
    if not os.path.isdir(outfile):
        sys.exit("Error: '%s' is not an existing directory." % (outfile,))

    if args.allenBrain_anno is not None:
        allenBrain_anno = args.allenBrain_anno
    if not os.path.exists(allenBrain_anno):
        sys.exit("Error: '%s' is not an existing directory." % (allenBrain_anno,))

    thres = args.threshold
    study = args.studyname
    labels = os.path.abspath(os.path.join(os.getcwd(), os.pardir,os.pardir,os.pardir))+'/lib/ABALabelIDs.mat'
    araDataTemplate = os.path.abspath(
        os.path.join(os.getcwd(), os.pardir, os.pardir, os.pardir)) + '/lib/annotation_50CHANGEDanno.nii.gz'
    incidenceMask = os.path.abspath(
        os.path.join(os.getcwd(), os.pardir, os.pardir, os.pardir)) + '/lib/MPI_maskBIG_for_incidence2.nii.gz'

    path = os.path.join(inputFile,study)
    regMR_list = findBETData(path)
    regInc_list = findIncData(path)
    regANNO_list = findRegisteredAnno(path)


    print("'%i' folders are part of the incidence map." % (len(regMR_list),))

    if not len(regANNO_list) == len(regMR_list):
        sys.exit("Error: For one or more annotations is no corresponding MR file defined in '%s'." % (inputFile,))

    incidenceMap(regMR_list,regInc_list,regANNO_list,araDataTemplate,incidenceMask,thres,outfile,allenBrain_anno,labels)