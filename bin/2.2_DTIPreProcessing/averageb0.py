"""
Averages all b0 images in a DTI dataset based on the b-values.

Given a 4D DWI NIfTI image file, this function locates the corresponding b-values file
(.bval or .btable), identifies all volumes with b-values less than 75 (b0 images), and
computes their mean to generate an averaged b0 image. If multiple b0 volumes are found,
each is saved as a separate NIfTI file, and all b0 volumes are aligned to the first using
FSL's MCFLIRT before averaging. The resulting mean b0 image is saved to disk.

Parameters
----------
input_file : str
    Path to the input 4D DWI NIfTI (.nii.gz) file. The corresponding .bval or .btable file
    must exist in the same directory.

use_mcflirt : bool, optional
    If True, uses FSL's MCFLIRT to align each b0 volume to the first b0 volume before averaging.
    Default is False, which averages the b0 volumes without alignment.

Returns
-------
output_file : str
    Path to the output NIfTI file containing the averaged b0 image.

Raises
------
SystemExit
    If the b-values file is not found, or if no b0 images are detected in the dataset.
"""

import os
import sys
import numpy as np
import nibabel as nii
import nipype.interfaces.fsl as fsl


def averageb0(input_file, b0_thresh=100, use_mcflirt=False):
    """
    Averages all b0 images in a DTI dataset, based on the bvals.
    Requires a 4D dwi image (input_file), with an existing bvals file in the same directory.
    """
    bvalsname = input_file.replace(".nii.gz", ".bval")
    if not os.path.exists(bvalsname):
        try:
            bvalsname = input_file.replace("Patch2SelfDenoised.nii.gz", ".bval")
            if not os.path.exists(bvalsname):
                bvalsname = input_file.replace("Patch2SelfDenoised.nii", ".bval")
            if not os.path.exists(bvalsname):
                bvalsname = input_file.replace(".nii.gz", ".btable")
                btable = np.loadtxt(bvalsname, dtype=float)
                bvalsname = os.path.splitext(bvalsname)[0] + ".bval"
                np.savetxt(bvalsname, btable[0, :], fmt='%.6f')
        except:
            sys.exit(f"Error: bvals file {bvalsname} not found.")
    bvals = np.loadtxt(bvalsname, dtype=float)
    if bvals.ndim > 1:
        bvals = bvals[0, :]
    # find b-values < b0_thresh
    b0_indices = np.where(bvals < b0_thresh)[0]
    if len(b0_indices) == 0:
        sys.exit(f'Error: No b0 images found (b-values < {str(b0_thresh)}).')
    data = nii.load(input_file)
    img = data.get_fdata()
    b0 = img[:, :, :, b0_indices]
    if b0.ndim < 4 or b0.shape[3] == 1:
        sys.exit(f'Error: No b0 images found (b-values < {str(b0_thresh)}).')
    if b0.shape[3] > 1:
        b0_filename = os.path.join(os.path.dirname(input_file),
                                   os.path.basename(input_file).split('.')[0] + '_b0.nii.gz')
        b0_nii = nii.Nifti1Image(b0, data.affine)
        b0_nii.header.set_xyzt_units('mm')
        nii.save(b0_nii, b0_filename)
        if use_mcflirt is True:
            # use FSL mcflirt to align each b0 volume to the first b0 volume
            mcflirt = fsl.MCFLIRT(in_file=b0_filename, out_file=b0_filename.replace(
                '.nii.gz', '_mcflirt.nii.gz'), ref_vol=0)
            mcflirt.inputs.save_rms = True
            mcflirt.inputs.save_plots = True
            mcflirt.inputs.cost = 'mutualinfo'
            mcflirt.run()
            b0_mc = nii.load(b0_filename.replace(
                '.nii.gz', '_mcflirt.nii.gz')).get_fdata()
            # average the b0 volumes
            b0mean = np.mean(b0_mc, axis=3)
        elif use_mcflirt is False:
            # average the b0 volumes without mcflirt alignment
            b0mean = np.mean(b0, axis=3)
    else:
        b0mean = np.mean(b0, axis=3)
    unscaledNiiData = nii.Nifti1Image(b0mean, data.affine)
    hdrOut = unscaledNiiData.header
    hdrOut.set_xyzt_units('mm')
    output_file = os.path.join(os.path.dirname(input_file),
                               os.path.basename(input_file).split('.')[0] + 'B0mean.nii.gz')
    nii.save(unscaledNiiData, output_file)
    return output_file


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Average b0 images in a DTI dataset')
    parser.add_argument('-i', '--input_file', help='Path to the input 4D DWI NIfTI file', required=True)
    parser.add_argument('-b', '--b0_thresh', default=100, help='B-value threshold under which volumes are treated as b0', required=False)
    parser.add_argument('-mcflirt', '--use_mcflirt', help='Use FSL MCFLIRT to align b0 volumes', action='store_true', required=False)
    args = parser.parse_args()

    if not os.path.exists(args.input_file):
        sys.exit(f"Error: Input file '{args.input_file}' does not exist.")
    
    output_file = averageb0(args.input_file, use_mcflirt=args.use_mcflirt)
    print(f"Averaged b0 image saved to: {output_file}")
    sys.exit(0)
