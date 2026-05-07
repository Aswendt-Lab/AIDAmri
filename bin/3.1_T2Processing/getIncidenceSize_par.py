
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
    # Return the array element closest to value. This helper is currently unused.
    idx = (np.abs(array-value)).argmin()
    return array[idx]

def thresholdingSlc(volumeMR,maskImg,thres):
    # Smooth the incidence volume before slice-wise thresholding.
    volumeMR=ndimage.gaussian_filter(volumeMR, sigma=(1.2, 1.2, 1))

    # Restrict the incidence volume to the provided 4D mask.
    volumeMR = volumeMR * maskImg[:, :, :, 0]

    # Save the masked volume for inspection.
    scaledNiiData = nii.Nifti1Image(volumeMR, np.eye(4))
    hdrIn = scaledNiiData.header
    hdrIn.set_xyzt_units('mm')
    output_file = os.path.join(outfile,'maskedVolume.nii.gz')
    nii.save(scaledNiiData, output_file)

    # Threshold each slice independently using mean + 1.5 * std of non-zero voxels.
    for i in range(len(volumeMR[1,1,:])-1):
        voSlc = volumeMR[:,:,i]

        # Remove values above the user-provided threshold before deriving a slice threshold.
        uvalues = voSlc >= thres
        voSlc[uvalues] = 0
        zvalues = voSlc != 0
        thresF = np.mean(voSlc[zvalues]) + 1.5*np.std(voSlc[zvalues])

        # Binarize the slice: below slice threshold becomes 0, above becomes 1.
        bvalues = voSlc < thresF
        voSlc[bvalues] = 0
        voSlc[bvalues] = 0
        fvalues = voSlc >= thresF
        voSlc[fvalues] = 1
        volumeMR[:, :, i] = voSlc

    fvalues = volumeMR == 1
    return volumeMR,fvalues



def thresholding(volumeMR,maskImg,thres,k):
    # Smooth the incidence image before global thresholding.
    volumeMR=ndimage.gaussian_filter(volumeMR, sigma=(1.3, 1.3, 1))

    # Use only non-zero voxels to estimate an automatic threshold.
    zvalues = volumeMR != 0



    if k==1:
        # Restrict thresholding to the binary incidence/stroke mask.
        volumeMR = volumeMR * maskImg[:, :, :]#, 0]

    if thres == 0:
        # If no threshold is provided, derive one from the masked non-zero signal.
        thres = np.mean(volumeMR[zvalues]) + 2*np.std(volumeMR[zvalues])

    # Remove values below threshold.
    bvalues = volumeMR < thres
    volumeMR[bvalues] = 0

    # Binarize remaining values as stroke/incidence voxels.
    fvalues = volumeMR >= thres
    volumeMR[fvalues] = 1

    return volumeMR


def incidenceMap(path_listInc,path_listMR ,path_listAnno, araDataTemplate,incidenceMask ,thres, outfile,labels):

    # Load the reference Allen/ARA template used for the affected-regions overlay.
    araDataTemplate  = nii.load(araDataTemplate)
    realAraImg = araDataTemplate.get_fdata()
    coloredAraLabels = np.zeros([np.size(realAraImg, 0), np.size(realAraImg, 1), np.size(realAraImg, 2)])

    # Load the parental label IDs from the MATLAB label file.
    matFile = sc.loadmat(labels)
    labMat = matFile['ABLAbelsIDsParental']

    # Load and binarize the externally provided stroke/incidence mask.
    maskData = nii.load(incidenceMask)
    maskImg = maskData.get_fdata()
    oneValues = maskImg > 0.0
    maskImg[oneValues] = 1.0

    # The script processes one subject folder at a time and uses the first matched file.
    fileIndex = 0

    # Load the parental annotation image and incidence data for the current subject.
    dataAnno = nii.load(path_listAnno[fileIndex])
    volumeAnno = np.round(dataAnno.get_fdata())
    dataMR = nii.load(path_listInc[fileIndex])
    volumeMR = dataMR.get_fdata()

    # Convert the incidence data to a binary stroke volume.
    strokeVolume = thresholding(volumeMR, maskImg, thres,1)

    # Keep only parental atlas labels that overlap with the binary stroke volume.
    fValues_Anno = volumeAnno*strokeVolume

    # Save the labelled stroke overlap in subject space.
    scaledNiiData = nii.Nifti1Image(fValues_Anno, dataAnno.affine)
    hdrIn = scaledNiiData.header
    hdrIn.set_xyzt_units('mm')
    output_file =  os.path.join(outfile,os.path.basename(path_listMR[fileIndex]).split('.')[0]+ 'Anno_parmask.nii.gz')
    nii.save(scaledNiiData, output_file)

    # Extract the unique non-zero parental labels affected by the stroke.
    ref_Image = fValues_Anno
    fValues_Anno = np.unique(fValues_Anno)
    nullValues = np.argwhere(fValues_Anno<=0.0)
    fValues_Anno = np.delete(fValues_Anno, nullValues)

    # Calculate how much of each affected parental region is covered by the stroke.
    regionAffectPercent = np.zeros(np.size(fValues_Anno))
    for i in range(np.size(fValues_Anno)):
        # Existing parental behavior scales by 200 and caps values at 100.
        regionAffectPercent[i] = (np.sum(ref_Image == fValues_Anno[i]) / np.sum(volumeAnno == fValues_Anno[i])) * 200
        regionAffectPercent[regionAffectPercent > 100] = 100

    # Keep only label IDs that are actually present in the affected-label list.
    labCounterList = np.isin(labMat[:, 0], fValues_Anno)
    labMat = labMat[labCounterList,0]

    # Create an affected-region overlay in the reference atlas space.
    labCounterColor = np.isin(realAraImg, fValues_Anno)
    coloredAraLabels[labCounterColor] = realAraImg[labCounterColor]

    # Offset one hemisphere by 2000 to preserve left/right label separation.
    xdim = np.size(coloredAraLabels, 0)
    coloredAraLabels[int(xdim / 2):xdim, :, :] = coloredAraLabels[int(xdim / 2):xdim, :, :] + 2000
    coloredAraLabels[coloredAraLabels == 2000] = 0

    # Save the affected parental regions as a NIfTI overlay.
    scaledNiiData = nii.Nifti1Image(coloredAraLabels, araDataTemplate.affine)
    hdrIn = scaledNiiData.header
    hdrIn.set_xyzt_units('mm')
    output_file = os.path.join(outfile, 'affectedRegions_Parental.nii.gz')
    nii.save(scaledNiiData, output_file)

    # Stroke volume calculation
    betMask = nii.load(os.path.join(outfile,os.path.basename(path_listInc[fileIndex]).split('.')[0]+'_mask.nii.gz'))
    betMaskImg = betMask.get_fdata()
    oneValues = betMaskImg > 0.0
    betMaskImg[oneValues] = 1.0

    # Estimate volumes from voxel count times voxel volume.
    strokeVolumeInCubicMM = np.sum(maskImg * (dataMR.affine[0, 0] * dataMR.affine[1, 1] * dataMR.affine[2, 2]))
    brainVolumeInCubicMM = np.sum(betMaskImg * (dataMR.affine[0, 0] * dataMR.affine[1, 1] * dataMR.affine[2, 2]))

    # Load label names and write the text summary of affected parental regions.
    lines =open(os.path.abspath(os.path.join(os.getcwd(), os.pardir,os.pardir))+ '/lib/annoVolume.nii.txt').readlines()
    o=open(os.path.join(outfile, 'affectedRegions_Parental.txt'), 'w')
    o.write("Stroke: %0.2f %% - Stroke Volume: %0.2f mm^3\n"  % (((strokeVolumeInCubicMM/brainVolumeInCubicMM)*100),strokeVolumeInCubicMM,))
    matIndex = 0
    labelNamesAffected = ["" for x in range(np.size(fValues_Anno))]
    labelNames = ["" for x in range(np.size(lines))]
    for i in range(len(lines)):
        labelNames[i] = lines[i].split('\t')[1]
        if  np.isin(int(lines[i].split('\t')[0]),labMat):
            # Write label ID, label name, and affected percentage for each matched region.
            o.write(lines[i][:-1] + "\t %0.2f %%\n" % regionAffectPercent[matIndex])
            labelNamesAffected[matIndex] = lines[i].split('\t')[1]
            matIndex = matIndex + 1

            #o.write(str(int(lines[i].split('	')[0]) + 2000) + '	R_' + lines[i].split('	')[1])
    o.close()

    # Store the same region statistics in a MATLAB file for downstream workflows.
    rows = np.shape(labMat)[0]
    labMat = np.stack((labMat, regionAffectPercent[0:rows]))
    matFile['ABLAbelsIDsParental'] = labMat
    matFile['ABANamesPar'] = labelNamesAffected
    matFile['ABAlabels'] = labelNames
    matFile['volumePer'] = (strokeVolumeInCubicMM / brainVolumeInCubicMM) * 100
    matFile['volumeMM'] = strokeVolumeInCubicMM
    sc.savemat(os.path.join(outfile, 'labelCount_par.mat'), matFile)






def findIncData(path):
    # Find incidence maps in immediate subfolders of the input folder.
    regMR_list = []

    for filename in glob.iglob(path+'*/*IncidenceData.nii.gz', recursive=False):
        regMR_list.append(filename)

    return regMR_list

def findBETData(path):
    # Find skull-stripped anatomical images in immediate subfolders of the input folder.
    regMR_list = []

    for filename in glob.iglob(path+'*/*Bet.nii.gz', recursive=False):
        regMR_list.append(filename)

    return regMR_list



def findRegisteredData(path):
    # Find registered template-space anatomical images. This helper is currently unused.
    regMR_list = []

    for filename in glob.iglob(path+'*/*_Template.nii.gz', recursive=True):
        regMR_list.append(filename)

    return regMR_list

def findRegisteredAnno(path):
    # Find parental annotation images produced by registration_DTI/T2 preprocessing.
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

    # Use the input folder as both source folder and output folder.
    if args.inputFolder is not None:
        inputFolder = args.inputFolder
        outfile = args.inputFolder
    if not os.path.exists(inputFolder):
        sys.exit("Error: '%s' is not an existing directory." % (inputFolder,))


    if args.allenBrain_anno is not None:
        allenBrain_anno = args.allenBrain_anno
    if not os.path.exists(allenBrain_anno):
        sys.exit("Error: '%s' is not an existing directory." % (allenBrain_anno,))

    # Resolve static label/template resources from the repository lib folder.
    thres = args.threshold
    labels = os.path.abspath(os.path.join(os.getcwd(), os.pardir,os.pardir))+ '/lib/rsfMRILablelID.mat'
    araDataTemplate = os.path.abspath(os.path.join(os.getcwd(), os.pardir,os.pardir))+ '/lib/annoVolume.nii.gz'

    # The stroke mask defines the lesion/incidence region used for all overlap calculations.
    if len(glob.glob(inputFolder+'/*Stroke_mask.nii.gz')) > 0:
        incidenceMask = glob.glob(inputFolder+'/*Stroke_mask.nii.gz')[0]
    else:
        sys.exit("Error: '%s' has no affected or masked regions." % (inputFolder,))

    # Collect required subject files from the input folder.
    path = os.path.join(inputFolder)
    regMR_list = findBETData(path)
    regInc_list = findIncData(path)
    regANNO_list = findRegisteredAnno(path)


    print("'%i' folder will be proccessed..." % (len(regMR_list),))

    if not len(regANNO_list) == len(regMR_list):
        sys.exit("Error: For one or more annotations is no corresponding MR file defined in '%s'." % (inputFolder,))

    # Calculate parental-region lesion overlap and write NIfTI, TXT, and MAT outputs.
    incidenceMap(regMR_list,regInc_list,regANNO_list,araDataTemplate,incidenceMask,thres,outfile,labels)
