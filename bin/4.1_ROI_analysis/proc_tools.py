'''
Created on 29.09.2020

Author:
Michael Diedenhofen
Max Planck Institute for Metabolism Research, Cologne
'''

from __future__ import print_function

import csv
import os
import sys

import nibabel as nib

from datetime import datetime

# directories
lib_in_dir = r'C:\Users\Public\Linux\shared_folder\AIDAmri\lib'
proc_in_dir = r'C:\Users\Public\Linux\shared_folder\proc_data'
#proc_out_dir = r'C:\Users\Public\Linux\shared_folder\proc_data'
proc_out_dir = r'C:\Users\Michael\Projects\Markus\Goeteborg\processed_data'
raw_in_dir = r'C:\Users\Public\Linux\shared_folder\raw_data'

# Enput labels text file with atlas index and seed regions (labels) in each line
# Atlas (1 or 2), Label 1, Label 2, ...
path_label_names_2000 = os.path.join(lib_in_dir, 'annoVolume+2000_rsfMRI.nii.txt')
path_labels = os.path.join(lib_in_dir, 'annotation_50CHANGEDanno_label_IDs+2000.txt')
path_labels_1 = r'C:\Users\Michael\Projects\Markus\Goeteborg\processed_data\cortex_labels_1.txt'
path_labels_2 = r'C:\Users\Michael\Projects\Markus\Goeteborg\processed_data\cortex_labels_2.txt'

# Enter all time points of the experiment
timepoints = ['Baseline', 'P7', 'P14', 'P28', 'P42', 'P56']

# Enter all experimental groups
groups = ['Treatment_C3a', 'Treatment_PBS']

# Enter subject name (same as in the Bruker folder structure) in the order of animals per time point for each group.
# For example: [[['subject1-group1-time point1', 'subject1-group1-time point1', ...]],
#               [['subject1-group2-time point1', 'subject1-group2-time point1', ...]],
#               [['subject1-group1-time point2', 'subject1-group1-time point2', ...]],
#               [['subject1-group2-time point2', 'subject1-group2-time point2', ...]]]
study = [[['GV_T3_12_1_1_3_20190808_093354', 'GV_T3_12_2_1_2_20190808_110831', 'GV_T3_12_3_1_2_20190808_120901', 'GV_T3_13_1_1_2_20190808_133527',
           'GV_T3_13_2_1_3_20190809_100309', 'GV_T3_13_3_1_2_20190808_144158', 'GV_T3_13_4_1_2_20190809_105931', 'GV_T3_16_1_1_1_20190903_123147',
           'GV_T3_16_2_1_1_20190903_131430', 'GV_T3_16_3_1_1_20190903_142337'],
          ['GV_T3_14_1_1_2_20190809_134721', 'GV_T3_14_2_1_1_20190809_145138', 'GV_T3_14_3_1_1_20190809_154126', 'GV_T3_14_4_1_1_20190812_103755',
           'GV_T3_16_4_1_1_20190904_080859', 'GV_T3_17_3_1_1_20190904_101208', 'GV_T3_17_4_1_1_20190904_105336']],
         [['GV_T3_12_1_1_4_20190820_090634', 'GV_T3_12_2_1_3_20190820_104019', 'GV_T3_12_3_1_3_20190820_112855', 'GV_T3_13_1_1_3_20190820_130848',
           'GV_T3_13_2_1_4_20190820_121817', 'GV_T3_13_3_1_3_20190820_141137', 'GV_T3_13_4_1_3_20190820_152016', 'GV_T3_16_1_1_2_20190923_111155',
           'GV_T3_16_2_1_2_20190923_144608', 'GV_T3_16_3_1_2_20190923_162928'],
          ['GV_T3_14_1_1_3_20190821_095620', 'GV_T3_14_2_1_2_20190821_105807', 'GV_T3_14_3_1_2_20190821_130012', 'GV_T3_14_4_1_2_20190821_134908',
           'GV_T3_16_4_1_2_20190924_112719', 'GV_T3_17_3_1_2_20190924_125208', 'GV_T3_17_4_1_2_20190924_133317']],
         [['GV_T3_12_1_1_5_20190829_091927', 'GV_T3_12_2_1_4_20190829_103718', 'GV_T3_12_3_1_4_20190829_112627', 'GV_T3_13_1_1_4_20190828_160508',
           'GV_T3_13_2_1_5_20190828_170236', 'GV_T3_13_3_1_4_20190828_174327', 'GV_T3_13_4_1_4_20190828_182452', 'GV_T3_16_1_1_3_20191001_091617',
           'GV_T3_16_2_1_3_20191001_100959', 'GV_T3_16_3_1_3_20191001_105158'],
          ['GV_T3_14_1_1_4_20190829_121905', 'GV_T3_14_2_1_3_20190829_130307', 'GV_T3_14_3_1_3_20190829_134612', 'GV_T3_14_4_1_3_20190829_143623',
           'GV_T3_16_4_1_3_20191002_083801', 'GV_T3_17_3_1_3_20191002_100801', 'GV_T3_17_4_1_3_20191002_105022']],
         [['GV_T3_12_1_1_6_20190910_084053', 'GV_T3_12_2_1_5_20190910_093051', 'GV_T3_12_3_1_5_20190910_101439', 'GV_T3_13_1_1_5_20190910_105531',
           'GV_T3_13_2_1_6_20190910_120953', 'GV_T3_13_3_1_5_20190910_125559', 'GV_T3_13_4_1_5_20190910_134736', 'GV_T3_16_1_1_4_20191015_092459',
           'GV_T3_16_2_1_5_20191015_151515', 'GV_T3_16_3_1_4_20191015_105651'],
          ['GV_T3_14_1_1_5_20190911_085005', 'GV_T3_14_2_1_4_20190911_094801', 'GV_T3_14_3_1_4_20190911_103054', 'GV_T3_14_4_1_4_20190911_111332',
           'GV_T3_16_4_1_4_20191015_114108', 'GV_T3_17_3_1_4_20191015_131043', 'GV_T3_17_4_1_4_20191015_135026']],
         [[None                            , 'GV_T3_12_2_1_6_20190924_163833', 'GV_T3_12_3_1_7_20190924_171919', 'GV_T3_13_1_1_6_20190925_082850',
           'GV_T3_13_2_1_7_20190925_091435', 'GV_T3_13_3_1_6_20190925_095841', 'GV_T3_13_4_1_6_20190925_104204', 'GV_T3_16_1_1_5_20191029_092504',
           'GV_T3_16_2_1_6_20191029_101500', 'GV_T3_16_3_1_5_20191029_110341'],
          ['GV_T3_14_1_1_6_20190925_113349', 'GV_T3_14_2_1_5_20190925_122019', 'GV_T3_14_3_1_5_20190925_132917', 'GV_T3_14_4_1_5_20190925_141548',
           'GV_T3_16_4_1_5_20191029_125150', 'GV_T3_17_3_1_5_20191029_133548', 'GV_T3_17_4_1_5_20191029_144222']],
         [['GV_T3_12_1_1_8_20191008_102322', 'GV_T3_12_2_1_7_20191008_125211', 'GV_T3_12_3_1_8_20191008_150900', 'GV_T3_13_1_1_7_20191008_161218',
           'GV_T3_13_2_1_8_20191009_084507', 'GV_T3_13_3_1_7_20191009_101526', 'GV_T3_13_4_1_7_20191009_112021', 'GV_T3_16_1_1_6_20191112_100410',
           'GV_T3_16_2_1_7_20191112_111819', 'GV_T3_16_3_1_6_20191112_121307'],
          ['GV_T3_14_1_1_8_20191009_121103', 'GV_T3_14_2_1_6_20191009_130701', 'GV_T3_14_3_1_6_20191009_140516', 'GV_T3_14_4_1_6_20191009_150737',
           'GV_T3_16_4_1_6_20191112_131231', 'GV_T3_17_3_1_6_20191112_141745', 'GV_T3_17_4_1_6_20191112_150712']]]

#study = [[['GV_T3_12_1_1_3_20190808_093354', 'GV_T3_13_3_1_2_20190808_144158', 'GV_T3_16_3_1_1_20190903_142337'], []], [['GV_T3_12_1_1_4_20190820_090634', 'GV_T3_13_3_1_3_20190820_141137', 'GV_T3_16_3_1_2_20190923_162928'], []], [[], []], [[], []], [[None, 'GV_T3_13_3_1_6_20190925_095841', 'GV_T3_16_3_1_5_20191029_110341'], []], [['GV_T3_12_1_1_8_20191008_102322', 'GV_T3_13_3_1_7_20191009_101526', 'GV_T3_16_3_1_6_20191112_121307'], []]]

# experiment number T2w
expno_T2w = [[[10, 8, 5, 5, 6, 6, 8, 6, 7, 10],
              [6, 6, 6, 6, 5, 6, 5]],
             [[7, 7, 6, 11, 5, 18, 6, 8, 5, 6],
              [11, 8, 6, 6, 5, 5, 10]],
             [[6, 5, 6, 5, 5, 5, 5, 6, 5, 8],
              [6, 5, 5, 6, 6, 5, 5]],
             [[6, 6, 5, 5, 8, 10, 5, 5, 5, 5],
              [10, 6, 5, 6, 5, 5, 5]],
             [[None, 6, 6, 7, 7, 6, 9, 7, 9, 7],
              [8, 10, 9, 6, 6, 15, 11]],
             [[6, 5, 10, 5, 10, 8, 5, 6, 6, 7],
              [5, 5, 6, 5, 10, 5, 6]]]

#expno_T2w = [[[10, 6, 10], []], [[7, 18, 6], []], [[], []], [[], []], [[None, 6, 7], []], [[6, 8, 7], []]]

# experiment number rsfMRI
expno_rsfMRI = [[[11, 9, 6, 6, 7, 7, 9, 7, 8, 11],
                 [7, 7, 9, 7, 6, 7, 6]],
                [[8, 10, 9, 12, 6, 19, 8, 9, 6, 7],
                 [12, 9, 7, 7, 6, 6, 11]],
                [[7, 6, 7, 6, 6, 6, 6, 7, 6, 9],
                 [7, 6, 6, 7, 7, 6, 6]],
                [[7, 7, 6, 6, 9, 11, 6, 6, 6, 6],
                 [11, 7, 6, 7, 6, 6, 6]],
                [[None, 7, 7, 8, 8, 7, 10, 8, 10, 8],
                 [9, 11, 10, 7, 7, 16, 12]],
                [[10, 6, 11, 6, 14, 9, 6, 7, 7, 10],
                 [6, 6, 7, 6, 11, 6, 7]]]

#expno_rsfMRI = [[[11, 7, 11], []], [[8, 19, 7], []], [[], []], [[], []], [[None, 7, 8], []], [[10, 9, 10], []]]

# experiment number DTI
expno_DTI = [[[12, 10, 7, 7, 8, 8, 10, 8, 9, 12],
              [8, 8, 8, 8, 7, 8, 7]],
             [[9, 9, 8, 13, 7, 20, 7, 10, 7, 8],
              [13, 10, 8, 8, 7, 7, 12]],
             [[8, 7, 8, 7, 7, 7, 7, 8, 7, 10],
              [8, 7, 7, 8, 8, 7, 7]],
             [[8, 8, 7, 7, 10, 12, 7, 8, 7, 7],
              [12, 8, 7, 8, 7, 7, 7]],
             [[None, 8, 8, 9, 9, 8, 11, 9, 11, 9],
              [10, 12, 11, 8, 8, 17, 13]],
             [[9, 8, 13, 8, 13, None, 8, 9, 9, 11],
              [8, 8, 9, 8, 13, 8, 9]]]

#expno_DTI = [[[12, 8, 12], []], [[9, 20, 8], []], [[], []], [[], []], [[None, 8, 9], []], [[9, None, 11], []]]

# processed images number
procno = 1

if not os.path.isdir(lib_in_dir):
    sys.exit("Error: '%s' is not an existing directory." % (lib_in_dir,))

if not os.path.isdir(proc_in_dir):
    sys.exit("Error: '%s' is not an existing directory." % (proc_in_dir,))

if not os.path.isdir(proc_out_dir):
    sys.exit("Error: '%s' is not an existing directory." % (proc_out_dir,))

if not os.path.isdir(raw_in_dir):
    sys.exit("Error: '%s' is not an existing directory." % (raw_in_dir,))

if not os.path.isfile(path_label_names_2000):
    sys.exit("Error: '%s' is not a regular file." % (path_label_names_2000,))

if not os.path.isfile(path_labels):
    sys.exit("Error: '%s' is not a regular file." % (path_labels,))

if not os.path.isfile(path_labels_1):
    sys.exit("Error: '%s' is not a regular file." % (path_labels_1,))

if not os.path.isfile(path_labels_2):
    sys.exit("Error: '%s' is not a regular file." % (path_labels_2,))

def get_date():
    now = datetime.now()
    pvDate = now.strftime("%a %d %b %Y")
    pvTime = now.strftime("%H:%M:%S")

    return pvDate + ' ' + pvTime

def read_csv(filename):
    if not os.path.isfile(filename):
        return None

    with open(filename, 'r') as fid:
        reader = csv.reader(fid, delimiter=',', dialect='excel', skipinitialspace=True)
        next(reader, None)
        data = list(reader)

    return data

def save_csv(filename, data):
    with open(filename, 'w') as fid:
        writer = csv.writer(fid, delimiter=';', dialect='excel',  lineterminator='\n')
        writer.writerows(data)

    print(filename)

def read_labels(filename):
    # read labels text file
    iatlas = []
    labels = []
    for row in read_csv(filename):
        iatlas.append(int(row[0].strip()))
        labels.append([int(label.strip()) for label in row[1:]])
    #print("iatlas:", iatlas)
    #print("labels:", labels)

    return (iatlas, labels)

def save_matrix(filename, matrix):
    lines = '\n'.join(('  '.join('%.12g' % (x,) for x in matrix[y]) + '  ') for y in range(matrix.shape[0]))

    # Open text file to write binary (Unix format)
    fid = open(filename, 'w')

    # Write text file
    for line in lines.splitlines():
        print(line, end=chr(10), file=fid)

    # Close text file
    fid.close()

    print(filename)

def read_text(filename):
    if not os.path.isfile(filename):
        return None

    # open file to read
    fid = open(filename, 'r')

    # read file -> list of lines
    lines = fid.readlines()

    # close file
    fid.close()

    return lines

def save_text(filename, lines):
    with open(filename, 'w') as fid:
        for line in lines:
            print(' '.join(line), end=chr(10), file=fid)

    print(filename)

def read_data(path_data):
    image = nib.load(path_data)
    data = image.get_data()

    header = image.get_header()
    voxel_dims = header.get_zooms()
    #print("header.get_data_shape():", header.get_data_shape())
    #print("header.get_data_dtype():", header.get_data_dtype())
    #print("header.get_zooms():", header.get_zooms())
    #print("header.get_data_offset():", header.get_data_offset())
    #print("header.get_xyzt_units():", header.get_xyzt_units())

    return (data, voxel_dims)

def save_data(data, voxel_dims, path_data, dtype='float32'):
    image = nib.Nifti1Image(data, None)

    header = image.get_header()
    if dtype is not None:
        header.set_data_dtype(dtype)
    if data.ndim == 3:
        header.set_zooms(voxel_dims)
    else:
        header.set_zooms(voxel_dims + (1.0,))
    header.set_xyzt_units(xyz='mm', t=None)

    image.to_filename(path_data)
    print(image.get_filename())

if __name__ == '__main__':
    pass
