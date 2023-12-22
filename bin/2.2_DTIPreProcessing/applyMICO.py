"""
@Article{li2014multiplicative,
  author    = {Li, Chunming and Gore, John C and Davatzikos, Christos},
  title     = {Multiplicative intrinsic component optimization (MICO) for MRI bias field estimation and tissue segmentation},
  journal   = {Magnetic resonance imaging},
  year      = {2014},
  volume    = {32},
  number    = {7},
  pages     = {913--923},
  publisher = {Elsevier},
}


Created on 10/08/2017

@author: Niklas Pallast
Neuroimaging & Neuroengineering
Department of Neurology
University Hospital Cologne

"""

import numpy as np
import nibabel as nii
import sys,os
import MICO
import progressbar
import cv2
from tqdm import tqdm

def run_MICO(IMGdata,outputPath):
    data = nii.load(IMGdata)
    v = 8

    vol = data.get_fdata()
    biasCorrectedVol = np.zeros(vol.shape[0:3])

    ImgMe = np.mean(vol)
    
    if ImgMe > 10000:
        nCvalue = 1000
    elif ImgMe > 1000:
        nCvalue = 10
    else:
        nCvalue = 1

    progressbar = tqdm(total=vol.shape[2], desc='Biasfieldcorrection')
    for idx in range(vol.shape[2]):

        Img = vol[:,:,idx] / nCvalue
        kernel =np.ones((5,5),np.uint8)
        erosion = cv2.erode(Img,kernel,iterations = 1)

        iterNum = 100
        N_region = 1
        q = 1
        thres = 100

        A = 1
        Img_original = Img

        nrow = Img.shape[0]
        ncol = Img.shape[1]
        n = nrow * ncol
        ROIt = Img > thres
        ROI = np.zeros([nrow,ncol])
        ROI[ROIt] = 1

        Bas = getBasisOrder3(nrow, ncol)

        N_bas = Bas.shape[2]

        ImgG = np.zeros([nrow,ncol,10])
        GGT = np.zeros([nrow, ncol, 10,10])

        for ii in range(N_bas):
            ImgG[:,:,ii] = Img * Bas[:,:, ii]*ROI
            for jj in range(N_bas):
                GGT[:,:,ii, jj] = Bas[:,:, ii]*Bas[:,:, jj]*ROI
                GGT[:,:,jj, ii] = GGT[:,:,ii, jj]

        energy_MICO = np.zeros([3, iterNum])

        b = np.ones([nrow,ncol])

        for ini_num  in range(1):
            C = np.random.rand(3, 1)
            C = C * A
            M = np.random.rand(nrow, ncol, 3)
            a = np.sum(M, 2)
            for k in range(N_region):
                M[:,:, k]=M[:,:, k]/ a

            e_max = np.amax(M,2)
            N_max = np.argmax(M,2)

            M_old = M
            chg = 10000

            energy_MICO[ini_num, 0] = get_energy(Img, b, C, M, ROI, q)

            for n in range(1,iterNum):
                M, b, C = MICO.runMICO(Img, q, ROI, M, C, b, Bas, GGT, ImgG, 1, 1)
                energy_MICO[ini_num, n] = get_energy(Img, b, C, M, ROI, q)
                if np.mod(n, 1) == 0:
                    PC = np.zeros([nrow,ncol])
                    for k in range(N_region):
                        PC = PC + C[k] * M[:,:, k]
                        img_bc = Img /b # bias field corrected image
                        smV = img_bc < 0
                        img_bc[smV] = 0
                        smV = img_bc > 5000
                        img_bc[smV] = 0




        M, C = sortMemC(M, C)
        seg = np.zeros([nrow,ncol])
        for k in range(N_region):
            seg = seg + k * M[:,:, k] # label  the k-th region


        biasCorrectedVol[:, :, idx] = img_bc
        
        progressbar.update(1)

    progressbar.close()

    unscaledNiiData = nii.Nifti1Image(biasCorrectedVol, data.affine)
    hdrOut = unscaledNiiData.header
    hdrOut.set_xyzt_units('mm')


    outputData = os.path.join(outputPath,os.path.basename(IMGdata).split('.')[0]+'Mico.nii.gz')
    nii.save(unscaledNiiData,outputData)

    return outputData

def sortMemC(M, C):
    C_out =np.sort(C)
    IDX= np.argsort(C)

    if len(M.shape) == 4:
        M_out = np.zeros([M.shape[0], M.shape[1], M.shape[2], len(IDX)])
        for k in range(np.size(C)):
            M_out[:,:,:,k] = M[:,:,:,IDX[k]]

    elif len(M.shape) == 3:
        M_out = np.zeros([M.shape[0], M.shape[1], len(IDX)])
        for k in range(np.size(C)):
            M_out[:,:,k] = M[:,:,IDX[k]]

    else:
            sys.exit('Error: sortMemC: wrong dimension of the membership function')

    return M_out, C_out


def get_energy(Img,b,C,M,ROI,q):
    N = M.shape[2]
    energy = 0

    for k in range(N):
        C_k = C[k] * np.ones([Img.shape[0],Img.shape[1]])
        energy = energy + np.sum(np.sum((Img * ROI - b * C_k * ROI) ** 2 * M[:, :, k] ** q))

    return energy


def getBasisOrder3(Height,Wide):
    x = np.zeros([Height,Wide])
    y = np.zeros([Height,Wide])
    for i  in range(Height):
        x[i,:] =  np.linspace(-1,1,Wide)

    for i  in range(Wide):
        y[:,i] =  np.linspace(-1,1,Height)

    bais = np.zeros([Height,Wide,10])
    bais[:,:,0] = 1
    bais[:,:,1] = x
    bais[:,:,2] = (3*x*x - 1)/2
    bais[:,:,3] = (5*x*x*x - 3*x)/2
    bais[:,:,4] = y
    bais[:,:,5] = x*y
    bais[:,:,6] = y*(3*x*x -1)/2
    bais[:,:,7] = (3*y*y -1)/2
    bais[:,:,8] = (3*y*y -1)*x/2
    bais[:,:,9] = (5*y*y*y -3*y)/2
    B = bais
    for kk in range(10):
        A=bais[:,:,kk]**2
        r = np.sqrt(sum(sum(A)))
        B[:,:,kk]=bais[:,:,kk]/r

    return B


if __name__ == "__main__":
    import argparse


    parser = argparse.ArgumentParser(description='Bias Correction')

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-i','--input', help='Path to input file',required=True)

    args = parser.parse_args()


    if args.input is not None and args.input is not None:
        input = args.input
    if not os.path.exists(input):
        sys.exit("Error: '%s' is not an existing directory of file %s is not in directory." % (input, args.file,))

    result = run_MICO(input)
