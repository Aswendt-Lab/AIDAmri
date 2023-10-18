"""
Created on 10/08/2017

@author: Niklas Pallast
Neuroimaging & Neuroengineering
Department of Neurology
University Hospital Cologne

"""
import os
from math import *
from lmfit import  Minimizer, Parameters
import matplotlib.pyplot as plt
import nibabel as nii
import numpy as np
import progressbar

from .ReferenceMethods import brummerSNR



plt.interactive(False)


def t2_monoexp3 (params, t, data):
    """
    # t2_monoexp3
    #
    # Define a mono-exponential decay equation with Y0 offset
    #------------------------------------------------------------------------------
    # params:   name (str, optional) – Name of the Parameter.
    #           value (float, optional) – Numerical Parameter value.
    #           vary (bool, optional) – Whether the Parameter is varied during a fit (default is True).
    #           min (float, optional) – Lower bound for value (default is -numpy.inf, no lower bound).
    #           max (float, optional) – Upper bound for value (default is numpy.inf, no upper bound).
    #           expr (str, optional) – Mathematical expression used to constrain the value during the fit.
    #           user_data (optional) – User-definable extra attribute used for a Parameter.
    # t: time axis in index units
    # data: scattert points
    """
    S0 = params['S0']
    T2 = params['T2']
    Y0 = params['Y0']
    model = S0 * np.exp(-t/T2)+Y0
    return model - data

def t2_monoexp2(params, t, data):
    """
    # t2_monoexp2
    #
    # Define a mono-exponential decay equation without Y0 offset
    # ------------------------------------------------------------------------------
    # params:   name (str, optional) – Name of the Parameter.
    #           value (float, optional) – Numerical Parameter value.
    #           vary (bool, optional) – Whether the Parameter is varied during a fit (default is True).
    #           min (float, optional) – Lower bound for value (default is -numpy.inf, no lower bound).
    #           max (float, optional) – Upper bound for value (default is numpy.inf, no upper bound).
    #           expr (str, optional) – Mathematical expression used to constrain the value during the fit.
    #           brute_step (float, optional) – Step size for grid points in the brute method.
    #           user_data (optional) – User-definable extra attribute used for a Parameter.
    # t: time axis in index units
    # data: scattert points
    """
    S0 = params['S0']
    T2 = params['T2']
    model = S0 * np.exp(-t/T2)
    return model - data

###############################################################################
# t2_fitmonoexp1
#
# Perform the data fitting using a mono-exponential model:
# S = Y0 + (S0 * EXP(-TE / T2))
###############################################################################


#def t2_fitmonoexp1 (T2, S0, Y0, T2bn, T2pe, img, snr, snrlim, nx, ny, slc, te, start, pinfo, FIXOFFSE):
#(slice, echoTime, model, curSnrMap, snrLim, slc)
def t2_fitmonoexp1(slice,te,snrMap,snrLim, model,uplim):

    dims = slice.shape
    nx = dims[1]
    ny = dims[0]

    T2 = np.zeros([nx,ny],dtype='int8') #Temporary store T2 map
    S0 = np.zeros([nx,ny],dtype='int8') #Temporary store S0 map


    # // FITTING PROCEDURE //
    bar = progressbar.ProgressBar()
    for i in bar(range(nx)):
        for j in range(ny):
            y = slice[i][j][:]
            if np.mean(snrMap[i, j]) >= snrLim:
                result=mpfitfun(y, te, model,uplim)
                T2[i, j] = result['T2'].value
                S0[i, j] = result['S0'].value
        bar.update(i)
    allResult = {'T2': T2, 'S0': S0, 'SNR': snrMap}
    #plt.imshow(T2, cmap='gray')
    return allResult

###############################################################################
# t2_fitmonoexp2
#
# Perform the data fitting using a mono-exponential model:
# S = S0 * EXP(-TE / T2)
#------------------------------------------------------------------------------
###############################################################################

#def t2_fitmonoexp2 (T2, S0, T2bn, T2pe, slice, snr, snrlim, nx, ny, slc, te, start, pinfo):
def t2_fitmonoexp2(slice,te,snrMap,snrLim, model,uplim):

    dims = slice.shape
    nx = dims[1]
    ny = dims[0]

    T2 = np.zeros([nx, ny],dtype='int8')  # Temporary store T2 map
    S0 = np.zeros([nx, ny],dtype='int8')  # Temporary store S0 map
    Y0 = np.zeros([nx, ny],dtype='int8')  # Temporary storeY0 map

    # // FITTING PROCEDURE //
    bar = progressbar.ProgressBar()
    for i in bar(range(nx)):
        for j in range(ny):
            y = slice[i][j][:]
            if np.mean(snrMap[i, j]) >= snrLim:
                result = mpfitfun(y, te, model, uplim)
                T2[i, j] = result['T2'].value
                S0[i, j] = result['S0'].value
                Y0[i, j] = result['Y0'].value
        bar.update(i)
    allResult = {'T2': T2, 'S0': S0, 'Y0': Y0, 'SNR': snrMap}
    #plt.imshow(T2, cmap='gray')
    return allResult


###############################################################################
# t2_mapping
#
# Main function for the T2 mapping project. Iterate over the slices
# T2 maps are generated from Bruker ParaVision data.
# Generate arrays to store data and call the fitting routines.
###############################################################################

def t2_mapping(data,echoTime, model, uplim, snrLim, SNRMethod):


    imgData = data.get_data()


    nx = imgData.shape[0] # Images size in x - direction
    ny = imgData.shape[1] # Images size in y - direction
    ns = imgData.shape[3] # Number of slices


    if 'T2_2p' in model:
         # Array to store the T2, S0 and Y0 maps
         pvMaps = np.zeros([nx, ny, ns, 3],dtype=data.get_data_dtype())


         #Loop to go through all slices
         for slc in range(ns):
            #   Print % of progress
            print('Slice: ' + str(slc + 1))

            # Temporal image containing all TE values for the selected slice
            slice = imgData[:, :, :, slc]


            # Temporal map containing the snr values for the selected slice
            if 'Chang' in SNRMethod:
                curSnrMap, estStdSijbers, estStdSijbersNorm = changSNR.calcSNR(slice, 0, 1)
            elif 'Brummer' in SNRMethod:
                curSnrMap, estStdSijbers, estStdSijbersNorm = brummerSNR.calcSNR(slice, 0, 1)
            elif 'Sijbers' in SNRMethod:
                curSnrMap, estStdSijbers, estStdSijbersNorm = sijbersSNR.calcSNR(slice, 0, 1)
            else:
                sys.exit("Error: No valid SNR model.")

            # Fit the data of the single slice (model 1)
            results = t2_fitmonoexp1(slice, echoTime, curSnrMap, snrLim, model, uplim)

            # Store data of slice in final image
            pvMaps[:, :, slc, 0] = results['T2']
            pvMaps[:, :, slc, 1] = results['S0']
            pvMaps[:, :, slc, 2] = results['SNR'][:, :, 0]

    elif 'T2_3p' in model:
        # Array to store the T2, S0 maps
        pvMaps = np.zeros([nx, ny, ns, 4],dtype=data.get_data_dtype())

        # Loop to go through all slices
        for slc in range(ns):
            #   Print % of progress
            print('Slice: ' + str(slc + 1) +'\n')

            # Temporal image containing all TE values for the selected slice
            slice = imgData[:, :, :, slc]

            # Temporal map containing the snr values for the selected slice
            if 'Chang' in SNRMethod:
                curSnrMap, estStdSijbers, estStdSijbersNorm = changSNR.calcSNR(slice, 0, 1)
            elif 'Brummer' in SNRMethod:
                curSnrMap, estStdSijbers, estStdSijbersNorm = brummerSNR.calcSNR(slice, 0, 1)
            elif 'Sijbers' in SNRMethod:
                curSnrMap, estStdSijbers, estStdSijbersNorm = sijbersSNR.calcSNR(slice, 0, 1)
            else:
                sys.exit("Error: No valid SNR model.")

            # Fit the data of the single slice (model 1)
            results = t2_fitmonoexp2(slice, echoTime, curSnrMap, snrLim, model, uplim)

            # Store data of slice in final image
            pvMaps[:, :, slc, 0] = results['T2']
            pvMaps[:, :, slc, 1] = results['S0']
            pvMaps[:, :, slc, 2] = results['Y0']
            pvMaps[:, :, slc, 3] = results['SNR'][:,:,0]

    else:
        sys.exit("Error: No valid model.")


    return pvMaps

def mpfitfun(data,te,model,uplim):

    y = data
    x = te
    params = Parameters()

    if 'T2_2p' in model:
        estT2 = (x[1]-x[0])/((np.log(y[0])/np.log(y[1])))
        est = [estT2, y[0]]
        params.add('T2', value=est[0], min=0, max=uplim)
        params.add('S0', value=est[1])
        minner = Minimizer(t2_monoexp2, params, fcn_args=(x, y))
        result = minner.minimize()
    elif 'T2_3p' in model:
        estT2 = (x[1]-x[0])/(np.log(y[0])/np.log(y[1]))
        est = [estT2, y[0], y[len(y)-1]]
        params.add('T2', value=est[0], min=0, max=70)
        params.add('S0', value=est[1])
        params.add('Y0', value=est[2])
        minner = Minimizer(t2_monoexp3, params, fcn_args=(x, y))
        result = minner.minimize()

    return result.params

def parsePV(filename):
    """
    Parser for Bruker ParaVision parameter files in JCAMP-DX format
    """

    # Read file 'filename' -> list 'lines'
    f = open(filename, 'r')
    lines = f.readlines()
    f.close()

    # Dictionary for parameters
    params = {}

    # Get STUDYNAME, EXPNO, and PROCNO
    #if filename[-9:] == 'visu_pars':
    if 'visu_pars' in filename:
        tmp = lines[6].split('/')
        params['studyname'] = [[], tmp[-5]]
        params['expno'] = [[], tmp[-4]]
        params['procno'] = [[], tmp[-2]]

    # Remove comment lines
    remove = [] # Index list
    for index, line in enumerate(lines): # Find lines
        if line[0:2] == '$$':
            remove.append(index)
    for offset, index in enumerate(remove): # Remove lines
        del lines[index-offset]

    # Create list of LDR (Labelled Data Record) elements
    lines = ''.join(lines).split('\n##') # Join lines and split into LDRs
    #lines = map(rstrip, lines) # Remove trailing whitespace from each LDR
    lines[0] = lines[0].lstrip('##') # Remove leading '##' from first LDR

    # Combine LDR lines
    for index, line in enumerate(lines):
        lines[index] = ''.join(line.split('\n'))

    # Fill parameter dictionary
    for line in lines:
        line = line.split('=', 1)
        if line[0][0] == '$':
            key = line[0].lstrip('$')
            dataset = line[1]
            params[key] = []
            pos = 0
            if (len(dataset) > 4) and (dataset[0:2] == '( '):
                pos = dataset.find(' )', 2)
                if pos > 2:
                    pardim = [int(dim) for dim in dataset[2:pos].split(',')]
                    params[key].append(pardim)
                    params[key].append(dataset[pos+2:])
            if pos <= 2:
                params[key].append([])
                params[key].append(dataset)

    # Remove specific elements from parameter dictionary
    if '$VisuCoreDataMin' in params: del params['$VisuCoreDataMin']
    if '$VisuCoreDataMax' in params: del params['$VisuCoreDataMax']
    if '$VisuCoreDataOffs' in params: del params['$VisuCoreDataOffs']
    if '$VisuCoreDataSlope' in params: del params['$VisuCoreDataSlope']
    if '$VisuAcqImagePhaseEncDir' in params: del params['$VisuAcqImagePhaseEncDir']

    for key in params.keys():
        pardim = params[key][0]
        parval = params[key][1]
        if (len(pardim) > 0) and (len(parval) > 0) and (parval[0] == '<'):
            params[key][1] = parval.replace('<', '"').replace('>', '"')
        elif (len(parval) > 0) and (parval[0] == '('):
            params[key][1] = parval.replace('<', '"').replace('>', '"')
        params[key] = params[key][1]

    return params


def getT2mapping(path,model,upLim,snrLim,SNRMethod,echoTime,output_path):

    data = nii.load(path)
    hdr = data.header
    raw = hdr.structarr
    if raw['dim'][3] < 2:
        sys.exit("Error: '%s' has wrong dimensions." % (path,))

    print('Start to  fit '+model+'-Map over TE %s ...' % (echoTime,) )

    map = t2_mapping(data, echoTime, model=model, uplim=upLim, snrLim=snrLim, SNRMethod=SNRMethod)
    pathT2Map = os.path.split(path)[0]
    map = map[:, :, :, 0] #delete this line if you want more outputdata
    mapNii =  nii.as_closest_canonical(nii.Nifti1Image(map, data.affine))
    hdr = mapNii.header
    hdr.set_xyzt_units('mm')
    study = os.path.split(path)[1].split('.')[0]
    nii.save(mapNii, output_path)



