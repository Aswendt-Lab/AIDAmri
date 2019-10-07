"""
Created on 10/08/2017

@author: Niklas Pallast
Neuroimaging & Neuroengineering
Department of Neurology
University Hospital Cologne

"""

import sys,os
import glob
import shutil as sh


def regABA2rsfMRI(inputVolume, T2data, brain_template, brain_anno, splitedAnno, splitedAnno_rsfMRI, anno_rsfMRI,
                  bsplineMatrix, dref, outfile):
    outputT2w = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_T2w.nii.gz')
    outputAff = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + 'transMatrixAff.txt')

    if dref:
        pathT2 = glob.glob(os.path.dirname(outfile) + '*/DTI/*T2w.nii.gz', recursive=False)
        sh.copy(pathT2[0], outputT2w)
    else:
        os.system(
            'reg_aladin -ref ' + inputVolume + ' -flo ' + T2data + ' -res ' + outputT2w + ' -aff ' + outputAff + ' -rigOnly')  # + -rigOnly' -fmask ' +MPITemplateMask+ ' -rmask ' + find_mask(inputVolume))
        #  resample Annotation
        outputAnno = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_Anno.nii.gz')
        os.system(
        'reg_resample -ref ' + inputVolume + ' -flo ' + brain_anno +
        ' -cpp ' + outputAff + ' -inter 0 -res ' + outputAnno)

    # resample splited  Annotation
    outputAnnoSplit = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_AnnoSplit.nii.gz')
    if dref:
        pathT2 = glob.glob(os.path.dirname(outfile) + '*/DTI/*AnnoSplit.nii.gz', recursive=False)
        sh.copy(pathT2[0], outputAnnoSplit)
    else:
        os.system(
        'reg_resample -ref ' + brain_anno + ' -flo ' + splitedAnno +
        ' -trans ' + bsplineMatrix + ' -inter 0 -res ' + outputAnnoSplit)
        os.system(
        'reg_resample -ref ' + inputVolume + ' -flo ' + outputAnnoSplit +
        ' -trans ' + outputAff + ' -inter 0 -res ' + outputAnnoSplit)

    # resample splited rsfMRI Annotation
    outputAnnoSplit_rsfMRI = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_AnnoSplit_rsfMRI.nii.gz')
    if dref:
        pathT2 = glob.glob(os.path.dirname(outfile) + '*/DTI/*AnnoSplit_rsfMRI.nii.gz', recursive=False)
        sh.copy(pathT2[0], outputAnnoSplit_rsfMRI)
    else:
        os.system(
        'reg_resample -ref ' + brain_anno + ' -flo ' + splitedAnno_rsfMRI +
        ' -trans ' + bsplineMatrix + ' -inter 0 -res ' + outputAnnoSplit_rsfMRI)
        os.system(
        'reg_resample -ref ' + inputVolume + ' -flo ' + outputAnnoSplit_rsfMRI +
        ' -trans ' + outputAff + ' -inter 0 -res ' + outputAnnoSplit_rsfMRI)

    # resample rsfMRI Annotation
    outputAnno_rsfMRI = os.path.join(outfile,
                                          os.path.basename(inputVolume).split('.')[0] + '_Anno_rsfMRI.nii.gz')
    if dref:
        pathT2 = glob.glob(os.path.dirname(outfile) + '*/DTI/*Anno_rsfMRI.nii.gz', recursive=False)
        sh.copy(pathT2[0], outputAnno_rsfMRI)
    else:
        os.system(
        'reg_resample -ref ' + brain_anno + ' -flo ' + anno_rsfMRI +
        ' -trans ' + bsplineMatrix + ' -inter 0 -res ' + outputAnno_rsfMRI)
        os.system(
        'reg_resample -ref ' + inputVolume + ' -flo ' + outputAnno_rsfMRI +
        ' -trans ' + outputAff + ' -inter 0 -res ' + outputAnno_rsfMRI)
        # resample in-house developed tempalate
        outputTemplate = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_Template.nii.gz')
        os.system(
        'reg_resample -ref ' + inputVolume + ' -flo ' + brain_template +
        ' -cpp ' + outputAff + ' -res ' + outputTemplate)


    return outputAnnoSplit

def find_RefStroke(refStrokePath,inputVolume):
    path =  glob.glob(refStrokePath+'/' + os.path.basename(inputVolume)[0:9]+'*/T2w/*IncidenceData_mask.nii.gz', recursive=False)
    return path

def find_RefAff(inputVolume):
    path =  glob.glob(os.path.dirname(os.path.dirname(inputVolume))+'/T2w/*MatrixAff.txt', recursive=False)
    return path

def find_RefTemplate(inputVolume):
    path =  glob.glob(os.path.dirname(os.path.dirname(inputVolume))+'/T2w/*TemplateAff.nii.gz', recursive=False)
    return path


def find_relatedData(pathBase):
    pathT2 =  glob.glob(pathBase+'*/T2w/*Bet.nii.gz', recursive=False)
    pathStroke_mask = glob.glob(pathBase + '*/T2w/*Stroke_mask.nii.gz', recursive=False)
    pathAnno = glob.glob(pathBase + '*/T2w/*Anno.nii.gz', recursive=False)
    pathAllen = glob.glob(pathBase + '*/T2w/*Allen.nii.gz', recursive=False)
    bsplineMatrix =  glob.glob(pathBase + '*/T2w/*MatrixBspline.nii', recursive=False)
    return pathT2,pathStroke_mask,pathAnno,pathAllen,bsplineMatrix


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Registration of T2 dataset and Allen Brain Atlas to rsfMRI')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-i', '--inputVolume', help='file name of DTI data after preprocessing', required=True)
    parser.add_argument('-d', '--dtiasRef', action='store_true', help='use DTI as reference if data quality is low')
    parser.add_argument('-r', '--referenceDay', help='Refernce Stroke mask', nargs='?', type=str,
                        default=None)
    parser.add_argument('-s', '--splitedAnno', help='Splited annotations atlas', nargs='?', type=str,
                        default=os.path.abspath(os.path.join(os.getcwd(), os.pardir,os.pardir))+'/lib/ARA_annotationR+2000.nii.gz')
    parser.add_argument('-f', '--splitedAnno_rsfMRI', help='Splited annotations atlas for rsfMRI', nargs='?', type=str,
                        default=os.path.abspath(os.path.join(os.getcwd(), os.pardir,os.pardir))+'/lib/annoVolume+2000_rsfMRI.nii.gz')
    parser.add_argument('-a', '--anno_rsfMRI', help='Annotations atlas for rsfMRI', nargs='?', type=str,
                        default=os.path.abspath(os.path.join(os.getcwd(), os.pardir,os.pardir))+'/lib/annoVolume.nii.gz')



    args = parser.parse_args()



    stroke_mask = None
    inputVolume = None
    refStrokePath = None
    splitedAnno = None
    splitedAnno_rsfMRI = None
    anno_rsfMRI = None

    if args.inputVolume is not None:
        inputVolume = args.inputVolume
    if not os.path.exists(inputVolume):
        sys.exit("Error: '%s' is not an existing directory." % (inputVolume,))

    outfile = os.path.join(os.path.dirname(inputVolume))
    if not os.path.exists(outfile):
        os.makedirs(outfile)

    print("rsfMRI Registration  \33[5m...\33[0m (wait!)", end="\r")
    # generate log - file
    sys.stdout = open(os.path.join(os.path.dirname(inputVolume), 'registration.log'), 'w')


    # find related  data
    pathT2, pathStroke_mask, pathAnno, pathTemplate, bsplineMatrix = find_relatedData(os.path.dirname(outfile))
    if len(pathT2) is 0:
        T2data = []
        sys.exit("Error: %s' has no reference T2 template." % (os.path.basename(inputVolume),))
    else:
        T2data = pathT2[0]

    if len(pathStroke_mask) is 0:
        pathStroke_mask = []
        print("Error: %s' has no reference stroke mask and is treated as sham." % (os.path.basename(inputVolume),))
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
            print("'%s' has no reference stroke mask and is treated as sham." % (os.path.basename(inputVolume),))
        else:
            refStroke_mask = refStroke_mask[0]

    if args.splitedAnno is not None:
        splitedAnno = args.splitedAnno
    if not os.path.exists(splitedAnno):
        sys.exit("Error: '%s' is not an existing directory." % (splitedAnno,))

    if args.splitedAnno_rsfMRI is not None:
        splitedAnno_rsfMRI = args.splitedAnno_rsfMRI
    if not os.path.exists(splitedAnno_rsfMRI):
        sys.exit("Error: '%s' is not an existing directory." % (splitedAnno_rsfMRI,))

    if args.anno_rsfMRI is not None:
        anno_rsfMRI = args.anno_rsfMRI
    if not os.path.exists(anno_rsfMRI):
        sys.exit("Error: '%s' is not an existing directory." % (anno_rsfMRI,))

    output = regABA2rsfMRI(inputVolume, T2data, brain_template, brain_anno, splitedAnno, splitedAnno_rsfMRI,
                           anno_rsfMRI, bsplineMatrix, args.dtiasRef, outfile)
    print(output + '...DONE!')
    sys.stdout = sys.__stdout__
    print('rsfMRI Registration  \033[0;30;42m COMPLETED \33[0m')



