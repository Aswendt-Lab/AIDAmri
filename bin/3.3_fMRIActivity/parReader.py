"""
Created on 10/08/2017

@author: Niklas Pallast
Neuroimaging & Neuroengineering
Department of Neurology
University Hospital Cologne

"""


import numpy as np
import os,sys
def getPar(filename):

    ## Open the text file.
    fileID = open(filename,'r')

    # Read columns of data according to the format.
    fileID.seek(0)
    lines = fileID.readlines()
    parData = np.zeros([len(lines),6])
    for row in range(len(lines)):
        for j in range(0, 12, 2):
            colume = int(j/2)
            parData[row][colume] = float(lines[row].split(' ')[j])
    fileID.close()


    # Create output variable
    return parData

if __name__ == "__main__":
    import argparse


    parser = argparse.ArgumentParser(description='par Reader')

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-i','--input', help='file name of data',required=True)

    args = parser.parse_args()


    if args.input is not None and args.input is not None:
        input = args.input
    if not os.path.exists(input):
        sys.exit("Error: '%s' is not an existing directory of file %s is not in directory." % (input, args.file,))

    result = getPar(input)
