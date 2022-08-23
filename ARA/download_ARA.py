"""
Created on 17/03/2020

@author: Niklas Pallast
Neuroimaging & Neuroengineering
Department of Neurology
University Hospital Cologne

More information can be found here:
http://alleninstitute.github.io/AllenSDK/_modules/allensdk/api/queries/reference_space_api.html
"""

import os
import nrrd  # pip install pynrrd, if pynrrd is not already installed
import nibabel as nib  # pip install nibabel, if nibabel is not already installed
import numpy as np
from allensdk.api.queries.reference_space_api import ReferenceSpaceApi
from allensdk.config.manifest import Manifest

# the annotation download writes a file, so we will need somwhere to put it
annotation_dir = 'annotation'
Manifest.safe_mkdir(annotation_dir)

annotation_path = os.path.join(annotation_dir, 'annotation.nrrd')

# this is a string which contains the name of the latest ccf version
annotation_version = ReferenceSpaceApi.CCF_VERSION_DEFAULT

# download annotations
mcapi = ReferenceSpaceApi()
mcapi.download_annotation_volume(annotation_version, 50, annotation_path)
annotation = nrrd.read(annotation_path)

# read nrrd data and header
_nrrd = nrrd.read(annotation_path)
data = _nrrd[0]
header = _nrrd[1]

# create header and for RAS orientation
space_value = header['space directions']
affine = np.eye(4) * 0.001 * space_value[0, 0]
affine[3][3] = 1
# ensure RAS orientation
data = np.swapaxes(data, 2, 0)
data = np.flip(data, 2)

img = nib.Nifti1Image(data, affine)  #
img = nib.as_closest_canonical(img)
hdrIn = img.header
hdrIn.set_xyzt_units('mm')
# img = nib.Nifti1Image(data, img.affine,hdrIn)
scaledNiiData = nib.as_closest_canonical(img)
nib.save(scaledNiiData, os.path.join(annotation_dir, os.path.dirname(annotation_path) + '.nii.gz'))

# the template download writes a file, so we will need somwhere to put it
template_dir = 'template'
Manifest.safe_mkdir(template_dir)

template_path = os.path.join(template_dir, 'template.nrrd')

# download templates
mcapi.download_template_volume(50, template_path)
template = nrrd.read(template_path)

# read nrrd data and header
_nrrd = nrrd.read(template_path)
data = _nrrd[0]
header = _nrrd[1]

# create header and for RAS orientation
space_value = header['space directions']
affine = np.eye(4) * 0.001 * space_value[0, 0]
affine[3][3] = 1
# ensure RAS orientation
data = np.swapaxes(data, 2, 0)
data = np.flip(data, 2)

img = nib.Nifti1Image(data, affine)  #
img = nib.as_closest_canonical(img)
hdrIn = img.header
hdrIn.set_xyzt_units('mm')
# img = nib.Nifti1Image(data, img.affine,hdrIn)
scaledNiiData = nib.as_closest_canonical(img)
nib.save(scaledNiiData, os.path.join(template_dir, os.path.dirname(template_path) + '.nii.gz'))
