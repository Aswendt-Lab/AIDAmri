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
import sys

def runMICO(Img,q,W,M,C,b,Bas,GGT,ImgG, Iter, iterCM):


    D = np.zeros(M.shape)
    for n in range(Iter):
        C = updateC(Img, W, b, M)
        for k in range (iterCM):
            N_class = M.shape[2]
            e = np.zeros(M.shape)
            for kk in range(N_class):
                D[:,:, kk] = (Img - C[kk] * b)** 2

            M = updateM(D, q)

    b_out = updateB(Img, q, C, M, Bas, GGT, ImgG)
    M_out = M
    C_out = C


    return M_out, b_out, C_out

def updateB(Img, q, C, M, Bas,GGT,ImgG):
    PC2 = np.zeros(Img.shape)
    PC = PC2

    N_class = M.shape[2]
    for kk in range(N_class):
        PC2 = PC2 + C[kk] ** 2 * M[:,:, kk]** q
        PC = PC + C[kk] * M[:,:, kk]** q
    N_bas = Bas.shape[2]
    V = np.zeros(N_bas)
    A = np.zeros([N_bas,N_bas])
    for ii in range(N_bas):
        ImgG_PC = ImgG[:,:,ii] * PC # Mask in ImgG
        V[ii] = np.sum(ImgG_PC) # inner product
        for jj in range(N_bas):
            B = GGT[:,:,ii, jj] * PC2 # Mask in GGT
            A[ii, jj] = np.sum(B) # inner product
            A[jj, ii] = A[ii, jj]

    w = np.dot(np.linalg.inv(A) , V)
    b = np.zeros(Img.shape)
    for kk in range (N_bas):
        b = b + np.dot(w[kk] , Bas[:,:, kk])

    return b





def updateC(Img, W,b, M):
    N_class=M.shape[2]
    C_new = np.zeros(N_class)
    for nn in range (N_class):
        N=b*Img*M[:,:,nn]
        D=(b**2) *M[:,:,nn]
        sN = np.sum(N*W)    # inner product
        sD = np.sum(D*W)   # inner product
        C_new[nn]=sN/(sD+(sD==0))

    return  C_new

def updateM(e, q):

    M = np.zeros(e.shape)
    N_class= e.shape[2]
    if q >1:
        epsilon=0.000000000001
        e=e+epsilon  # avoid division by zero
        p = 1/(q-1)
        f = 1/(e**p)
        f_sum = np.sum(f,2)
        for kk in range(N_class):
            M[:,:,kk] = f[:,:,kk]/f_sum

    elif q==1:
        e_min = np.amin(e,2)
        N_min = np.argmin(e,2)
        for kk in range (N_class):
            tempComp = (N_min == kk)
            M[:,:,kk] = tempComp

    else:
        sys.exit('Error: MICO: wrong fuzzifizer')

    return M