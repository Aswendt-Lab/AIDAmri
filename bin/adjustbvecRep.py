"""
usage
python adjustbvecRep.py /aida/Data/output/sub-subname/ses-sesname
"""

import os
import sys
import glob
import numpy as np
import nibabel as nib


directory_use = sys.argv[1]

if os.path.isdir(directory_use):
    print(f"Processing directory: {directory_use}")
else:
    print(f"Error: {directory_use} is not a valid directory.")

os.chdir(directory_use)
os.chdir("dwi")

dwi_files = glob.glob('*dwi.nii.gz')

for dwifile in dwi_files:
    print(f"Processing File: {dwifile}")
    img = nib.load(dwifile)
    image_shape = img.shape
    dwibasename, _ = os.path.splitext(dwifile)
    dwibasename, _ = os.path.splitext(dwibasename)
    dwibvalname = dwibasename + ".bval"
    print(f"DWI bvalue name: {dwibvalname}")
    dwibvecname = dwibasename + ".bvec"
    print(f"DWI bvec name: {dwibvecname}")

    bval = np.loadtxt(dwibvalname, dtype=float)
    bvec = np.loadtxt(dwibvecname, dtype = float)
 
    if image_shape[3] == 1:
        print(f"Skipping {dwifile} as it only has 1 time point.")
    else: 
        bvaltile_fact = image_shape[3]/len(bval)
        print(f"Scaling up bval by a factor of {bvaltile_fact} repetitions.")

        #bvecsize = bvec.shape
        bvectile_fact = image_shape[3]/(bvec.shape[1])
        print(f"Scaling up bvec by a factor of {bvectile_fact} repetitions.")
    
        if bvaltile_fact == int(bvaltile_fact):
           bvaltile_fact = int(bvaltile_fact)
        else: 
           print("length of bval not compatible with length of diffusion data.")

        if bvectile_fact == int(bvectile_fact):
           bvectile_fact = int(bvectile_fact)
        else: 
           print("length of bvec not compatible with length of diffusion data.")

        bval = np.tile(bval,bvaltile_fact)
        bvec = np.tile(bvec,bvectile_fact)
        bval = bval.reshape(1,-1)
        bval2name = dwibvalname   # for debugging + "2"
        bvec2name = dwibvecname   # for debugging + "2"
    
        np.savetxt(bval2name, bval, fmt="%1.15f", delimiter=" ")
        np.savetxt(bvec2name, bvec, fmt="%1.17f", delimiter=" ")
  
   



