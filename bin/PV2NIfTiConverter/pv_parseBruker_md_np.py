"""
Created on 10/08/2017

@author: Niklas Pallast
Neuroimaging & Neuroengineering
Department of Neurology
University Hospital Cologne

"""

from __future__ import print_function

import os,sys

import numpy as np

from dict2xml import createXML


# from string import split

def parsePV(filename):
    """
    Parser for Bruker ParaVision parameter files in JCAMP-DX format

    Prarmeters:
    ===========
    filename: 'acqp', 'method', 'd3proc', 'roi', 'visu_pars', etc.
    """
    if not os.path.exists(filename):
        return []

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
    if np.size(lines) == 1:
        sys.exit("Error: visu_pars is not readable")

    if 'subject' in filename:
        tmp = lines[32].split('#$Name,')
        params['coilname'] = tmp[1].split('#$Id')[0]
        return params

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

def getXML(filename, writeFile=False):
    """
    Writes header dictionary to xml format

    Parameters:
    ==========
    filename: Bruker ParaVision '2dseq' file
    writeFile: Boolean, if 'False' return string containing xml-header, else save to file
    """

    path = os.path.abspath(os.path.dirname(filename))

    # Parse all parameter files
    header_acqp = parsePV(os.path.join(path, '..', '..', 'acqp'))
    header_method = parsePV(os.path.join(path, '..', '..', 'method'))
    #header_d3proc = parsePV(os.path.join(path, 'd3proc'))   # removed for PV6
    header_visu = parsePV(os.path.join(path, 'visu_pars'))

    header = {'Scaninfo': {}}
    header['Scaninfo']['acqp'] = header_acqp
    header['Scaninfo']['method'] = header_method
    #header['Scaninfo']['d3proc'] = header_d3proc            # removed for PV6
    header['Scaninfo']['visu_pars'] = header_visu

    xml = createXML(header, '<?xml version="1.0"?>\n')

    if writeFile:
        f = open('scaninfo.xml', 'w')
        f.write(xml)
        f.close()
    else:
        return xml

def getNiftiHeader(params, sc=10):
    """
    Returns necessary header parameters for NIfTI generation ()

    Parameters:
    ===========
    filename: header returned from parser
    sc: scales pixel dimension (defaults to 10 for animal imaging)
    """
    # List of 'VisuCoreSize' parameter strings
    if params == []:
        return
    CoreSize = str.split(params['VisuCoreSize'])

    if params['VisuCoreDimDesc'] == 'spectroscopic':
        print("spectroscopic")
        #quit(42)
        return params['VisuStudyDate'], int(CoreSize[0]), 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 8

    # Dimensions
    nX = int(CoreSize[0])
    nY = int(CoreSize[1])
    nZ = 1
    nT = 1



    # FrameGroup dimensions
    if 'VisuFGOrderDescDim' in params:
        if int(params['VisuFGOrderDescDim']) > 0:

            FGOrderDesc = params['VisuFGOrderDesc'][1:-1].split(') (')
            #FGOrderDesc = map(lambda item: item.split(', '), FGOrderDesc)
            FGOrderDesc = [item.split(', ') for item in FGOrderDesc]
            #frameDims = map(lambda item: int(item[0]), FGOrderDesc)
            frameDims = [int(item[0]) for item in FGOrderDesc]
            # Number of slices
            nZ = frameDims[0]
            if int(params['VisuFGOrderDescDim']) > 1:
                nT = frameDims[1]
            if int(params['VisuFGOrderDescDim']) > 2:
                nT *= frameDims[2]

    # Voxel dimensions
    extent = params['VisuCoreExtent'].split()
    dX = sc * float(extent[0]) / nX
    dY = sc * float(extent[1]) / nY
    VisuCoreSlicePacksSliceDist = params.get('VisuCoreSlicePacksSliceDist')
    print("VisuCoreSlicePacksSliceDist",VisuCoreSlicePacksSliceDist)
    print("VisuCoreFrameThickness", params['VisuCoreFrameThickness'])
    if VisuCoreSlicePacksSliceDist is None:
        dZ = sc * float(params['VisuCoreFrameThickness'])
    else:
        # Slice thickness inclusive gap (PV6)
        VisuCoreSlicePacksSliceDist=VisuCoreSlicePacksSliceDist.split()[0]
        print("VisuCoreSlicePacksSliceDist",VisuCoreSlicePacksSliceDist)
        dZ = sc * float(VisuCoreSlicePacksSliceDist)

    if 'VisuAcqRepetitionTime' in params:
        if (nT > 1) and (float(params['VisuAcqRepetitionTime']) > 0 ):
            dT = float(params['VisuAcqRepetitionTime']) / 1000
        else:
            dT=0
    else:
        dT = 0

    if int(params['VisuCoreDim']) == 3:
        nZ = int(CoreSize[2])
        nT = 1
        frameDims = None
        if 'VisuFGOrderDescDim' in params:
            if int(params['VisuFGOrderDescDim']) > 0:
                nT = frameDims[0]
        dZ = sc * float(extent[2]) / nZ
        if (nT > 1) and (float(params['VisuAcqRepetitionTime']) > 1 ):
            dT = float(params['VisuAcqRepetitionTime']) / 1000
        else:
            dT = 0

    DT = 4
    if params['VisuCoreWordType'] == '_8BIT_UNSGN_INT': DT = 'int8'
    if params['VisuCoreWordType'] == '_16BIT_SGN_INT' : DT = 'int16'
    if params['VisuCoreWordType'] == '_32BIT_SGN_INT' : DT = 'int32'
    if params['VisuCoreWordType'] == '_32BIT_FLOAT'   : DT = 'float32'

    tmp = params['studyname'] + '.' + params['expno'] + '.' + params['procno']
    return tmp, nX, nY, nZ, nT, dX, dY, dZ, dT, 0, 0, 0, DT

''' def getNiftiHeader_md(params, sc=10):
    """
    Returns necessary header parameters for NIfTI generation ()

    Parameters:
    ===========
    filename: header returned from parser
    sc: scales pixel dimension (defaults to 10 for animal imaging)
    """
    # List of 'VisuCoreSize' parameter strings
    global frameDims
    CoreSize = str.split(params['VisuCoreSize'])
    if params['VisuCoreDimDesc'] == 'spectroscopic':
        return params['VisuStudyName'], int(CoreSize[0]), 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 8
    
    # Dimensions
    nX = int(CoreSize[0])
    nY = int(CoreSize[1])
    nZ = 1
    nT = 1

    # FrameGroup dimensions
    if int(params['VisuFGOrderDescDim']) > 0:
        FGOrderDesc = params['VisuFGOrderDesc'][1:-1].split(') (')
        #FGOrderDesc = map(lambda item: item.split(', '), FGOrderDesc)
        FGOrderDesc = [item.split(', ') for item in FGOrderDesc]
        #frameDims = map(lambda item: int(item[0]), FGOrderDesc)
        frameDims = [int(item[0]) for item in FGOrderDesc]
        # Number of slices
        nZ = frameDims[0]
        if int(params['VisuFGOrderDescDim']) > 1:
            nT = frameDims[1]

    # Voxel dimensions
    extent = params['VisuCoreExtent'].split()
    dX = sc * float(extent[0]) / nX
    dY = sc * float(extent[1]) / nY
    VisuCoreSlicePacksSliceDist = params.get('VisuCoreSlicePacksSliceDist')
    print("VisuCoreSlicePacksSliceDist",VisuCoreSlicePacksSliceDist)
    if VisuCoreSlicePacksSliceDist is None:
        dZ = sc * float(params['VisuCoreFrameThickness']) 
    else:
        # Slice thickness inclusive gap (PV6)
        dz = sc * float(VisuCoreSlicePacksSliceDist)
    print("dz",dz)
    if (nT > 1) and (float(params['VisuAcqRepetitionTime']) > 0 ):
        dT = float(params['VisuAcqRepetitionTime']) / 1000
    else:
        dT = 0

    if int(params['VisuCoreDim']) == 3:
        nZ = int(CoreSize[2])
        nT = 1
        if int(params['VisuFGOrderDescDim']) > 0:
            nT = frameDims[0]
        dZ = sc * float(extent[2]) / nZ
        if (nT > 1) and (float(params['VisuAcqRepetitionTime']) > 1 ):
            dT = float(params['VisuAcqRepetitionTime']) / 1000
        else:
            dT = 0

    CoreWordType = params['VisuCoreWordType']
    if CoreWordType   == '_8BIT_UNSGN_INT': DT = 'DT_UINT8'   #  2: 'int8'
    elif CoreWordType == '_16BIT_SGN_INT':  DT = 'DT_INT16'   #  4: 'int16'
    elif CoreWordType == '_32BIT_SGN_INT':  DT = 'DT_INT32'   #  8: 'int32'
    elif CoreWordType == '_32BIT_FLOAT':    DT = 'DT_FLOAT32' # 16: 'float32'
    else: DT = 4

    tmp = params['studyname'] + '.' + params['expno'] + '.' + params['procno']
    return (tmp, nX, nY, nZ, nT, dX, dY, dZ, dT, 0, 0, 0, DT)
'''
def getRotMatrix(filename):
    """
    Returns rotation matrix for image registration

    Parameters:
    ===========
    filename: visu_pars file to parse
    sc : scales pixel dimension (defaults to 10 for animal imaging)
    """
    params = parsePV(filename)
    orientation = map(float, params['VisuCoreOrientation'].split())
    if not 'VisuCorePosition' in params:
        return np.array([0.0, 0.0, 0.0, 0.0])
    position    = map(float, params['VisuCorePosition'].split())
    orientation = np.array(orientation[0:9]).reshape((3, 3))
    position    = np.array(position[0:3]).reshape((3, 1))
    rotMatrix   = np.append(orientation, position, axis=1)
    rotMatrix   = np.append(rotMatrix, np.array([0.0, 0.0, 0.0, 1.0]).reshape(1, 4), axis=0)
    return rotMatrix

def writeRotMatrix(rotMatrix, filename):
    fid = open(filename, 'w')
    np.savetxt(fid, rotMatrix, fmt='%-7.2f')
    fid.close()

"""
if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Parse Bruker parameter files')
    parser.add_argument('file', type=str, help='name of parameter file (2dseq for all)')
    parser.add_argument('-t', '--type', type=str, default = "xml", help='nifti, xml')
    args = parser.parse_args()

    if args.type == "nifti" and os.path.basename(args.file)=="visu_pars":
        print(getNiftiHeader(args.file))
    elif args.type == "mat" and os.path.basename(args.file)=="visu_pars":
        print(getRotMatrix(args.file))
    elif args.type == "xml":
        print(getXML(args.file))
    else:
        print("Hmmm, works not that way ;)")
"""
