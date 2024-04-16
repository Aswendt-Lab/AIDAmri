'''
Created on 25.08.2020

Author:
Michael Diedenhofen
Max Planck Institute for Metabolism Research, Cologne

Description:
Helper tool to compare the number of voxels included in the peri-infarct region for each subject.
'''

from __future__ import print_function

try:
    zrange = xrange
except NameError:
    zrange = range

import os
import sys

import proc_tools as pt

def count_voxels(rois):
    values = []
    for k in zrange(rois.shape[3]):
        roi = rois[:, :, :, k]
        values.append(roi[roi>0].size)

    return values

def main():
    # output text file
    path_data = os.path.join(pt.proc_out_dir, 'ROIs_count_voxels.txt')

    # read label IDs
    _, labels = pt.read_labels(pt.path_labels_1)

    data = []
    for index_t, timepoint in enumerate(pt.timepoints):
        data.append([timepoint] + labels[0])
        for index_g, group in enumerate(pt.groups):
            for subject in pt.study[index_t][index_g]:
                if subject is not None:
                    in_dir = os.path.join(pt.proc_out_dir, timepoint, group, subject, 'fMRI')
                    if not os.path.isdir(in_dir):
                        sys.exit("Error: '%s' is not an existing directory." % (in_dir,))
                    # input ROIs file (NIfTI)
                    #path_in_rois = os.path.join(in_dir, subject + '_cortex_rois_2.nii.gz')
                    path_in_rois = os.path.join(in_dir, 'Seed_ROIs_peri.nii.gz')
                    if not os.path.isfile(path_in_rois):
                        sys.exit("Error: '%s' is not a regular file." % (path_in_rois,))
                    # read ROIs hyperstack (4D)
                    rois, _ = pt.read_data(path_in_rois)
                    values = count_voxels(rois)
                    data.append([subject] + values)

    pt.save_csv(path_data, data)

if __name__ == '__main__':
    main()
