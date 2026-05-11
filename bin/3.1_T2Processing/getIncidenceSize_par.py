
"""
Created on 10/08/2017

@author: Niklas Pallast
Neuroimaging & Neuroengineering
Department of Neurology
University Hospital Cologne


"""


import os,sys
import nibabel as nii
import glob
import numpy as np
import scipy.io as sc

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))


def voxel_volume_mm3(nifti_img):
    # Use the full affine so rotated/sheared registered images still get a valid voxel volume.
    voxel_volume = abs(np.linalg.det(nifti_img.affine[:3, :3]))
    if voxel_volume == 0:
        voxel_volume = np.prod(nifti_img.header.get_zooms()[:3])
    return voxel_volume


def ensure_same_grid(reference_img, moving_img, moving_description):
    if reference_img.shape != moving_img.shape:
        sys.exit("Error: Stroke mask and %s have different dimensions." % (moving_description,))
    if not np.allclose(reference_img.affine, moving_img.affine, atol=1e-4):
        sys.exit("Error: Stroke mask and %s have different affine matrices." % (moving_description,))


def nifti_name_without_extension(file_path):
    file_name = os.path.basename(file_path)
    if file_name.endswith('.nii.gz'):
        return file_name[:-len('.nii.gz')]
    if file_name.endswith('.nii'):
        return file_name[:-len('.nii')]
    sys.exit("Error: '%s' is not a NIfTI file." % (file_path,))


def find_brkraw_nifti(input_folder):
    brkraw_folder = os.path.join(input_folder, 'brkraw')
    matches = sorted(glob.glob(os.path.join(brkraw_folder, '*.nii')) +
                     glob.glob(os.path.join(brkraw_folder, '*.nii.gz')))
    if len(matches) == 0:
        matches = sorted(glob.glob(os.path.join(input_folder, '**', 'brkraw', '*.nii'), recursive=True) +
                         glob.glob(os.path.join(input_folder, '**', 'brkraw', '*.nii.gz'), recursive=True))
    if len(matches) == 0:
        sys.exit("Error: No NIfTI file found in a brkraw folder under '%s'." % (input_folder,))
    if len(matches) > 1:
        sys.exit("Error: Multiple NIfTI files found in brkraw folders under '%s': %s" %
                 (input_folder, ', '.join(matches),))
    return matches[0]


def calculate_parental_stroke_overlap(brain_file, parental_annotation_file, ara_template_file, stroke_mask_file,
                                      incidence_lesion_mask_file, output_folder, label_file):

    # Load the reference Allen/ARA template used for the affected-regions overlay.
    ara_template_img  = nii.load(ara_template_file)
    ara_template_labels = ara_template_img.get_fdata()
    affected_template_labels = np.zeros([np.size(ara_template_labels, 0), np.size(ara_template_labels, 1), np.size(ara_template_labels, 2)])

    # Load the parental label IDs from the MATLAB label file.
    label_mat = sc.loadmat(label_file)
    parental_label_ids = label_mat['ABLAbelsIDsParental']

    # Save the template-space IncidenceData lesion mask labelled by the parental rsfMRI atlas.
    incidence_lesion_mask_img = nii.load(incidence_lesion_mask_file)
    incidence_lesion_mask = incidence_lesion_mask_img.get_fdata()
    incidence_lesion_mask[incidence_lesion_mask > 0.0] = 1.0
    incidence_lesion_mask[incidence_lesion_mask <= 0.0] = 0.0

    parental_incidence_atlas_file = os.path.join(REPO_ROOT, 'lib', 'annoVolume+2000_rsfMRI.nii.gz')
    parental_incidence_atlas_img = nii.load(parental_incidence_atlas_file)
    ensure_same_grid(incidence_lesion_mask_img, parental_incidence_atlas_img, "parental incidence atlas")
    parental_incidence_atlas = np.round(parental_incidence_atlas_img.get_fdata())

    labelled_incidence_lesion = parental_incidence_atlas*incidence_lesion_mask
    labelled_incidence_lesion_img = nii.Nifti1Image(labelled_incidence_lesion, incidence_lesion_mask_img.affine)
    labelled_incidence_lesion_img.header.set_xyzt_units('mm')
    incidence_data_dir = os.path.dirname(incidence_lesion_mask_file)
    incidence_data_name = os.path.basename(incidence_lesion_mask_file)
    incidence_input_suffix = 'IncidenceData_Lesion_mask.nii.gz'
    incidence_output_suffix = 'IncidenceData_Anno_parental_lesion_mask.nii.gz'
    if not incidence_data_name.endswith(incidence_input_suffix):
        sys.exit("Error: Incidence lesion mask filename must end with '%s'." % (incidence_input_suffix,))
    incidence_name_prefix = incidence_data_name[:-len(incidence_input_suffix)]
    output_name = incidence_name_prefix + incidence_output_suffix
    output_file = os.path.join(incidence_data_dir, output_name)
    nii.save(labelled_incidence_lesion_img, output_file)

    # Load and binarize the externally provided stroke mask.
    stroke_mask_img = nii.load(stroke_mask_file)
    stroke_mask = stroke_mask_img.get_fdata()
    stroke_mask[stroke_mask > 0.0] = 1.0
    stroke_mask[stroke_mask <= 0.0] = 0.0

    # Load subject-space parental annotation and brain image.
    parental_annotation_img = nii.load(parental_annotation_file)
    parental_annotation = np.round(parental_annotation_img.get_fdata())
    brain_img = nii.load(brain_file)
    brain_volume = brain_img.get_fdata()

    ensure_same_grid(stroke_mask_img, parental_annotation_img, "parental annotation")
    ensure_same_grid(stroke_mask_img, brain_img, "brain image")

    # Keep only parental atlas labels that overlap with the subject-space stroke mask.
    labelled_stroke_overlap = parental_annotation * stroke_mask

    # Extract the unique non-zero parental labels affected by the stroke.
    affected_labels = np.unique(labelled_stroke_overlap)
    affected_labels = affected_labels[affected_labels > 0.0]

    # Calculate how much of each affected parental region is covered by the stroke.
    region_percent_by_label = {}
    for label_id in affected_labels:
        # Percentage of the parental region covered by stroke voxels.
        affected_voxels = np.sum(labelled_stroke_overlap == label_id)
        total_region_voxels = np.sum(parental_annotation == label_id)
        region_percent = (affected_voxels / total_region_voxels) * 100
        region_percent_by_label[int(label_id)] = min(region_percent, 100)

    # Keep only label IDs that are actually affected by stroke (in the affected-label list)
    affected_label_ids = np.array(sorted(region_percent_by_label))
    affected_label_mask = np.isin(parental_label_ids[:, 0], affected_label_ids)
    parental_label_ids = parental_label_ids[affected_label_mask,0]
    parental_label_id_set = set(int(label_id) for label_id in parental_label_ids)

    # Create an affected-region overlay in the reference atlas space.
    affected_template_mask = np.isin(ara_template_labels, affected_label_ids)
    affected_template_labels[affected_template_mask] = ara_template_labels[affected_template_mask]

    # Offset one hemisphere by 2000 to preserve left/right label separation.
    xdim = np.size(affected_template_labels, 0)
    affected_template_labels[int(xdim / 2):xdim, :, :] = affected_template_labels[int(xdim / 2):xdim, :, :] + 2000
    affected_template_labels[affected_template_labels == 2000] = 0

    # Save the affected parental regions as a NIfTI overlay.
    affected_regions_img = nii.Nifti1Image(affected_template_labels, ara_template_img.affine)
    affected_regions_img.header.set_xyzt_units('mm')
    brkraw_nifti_file = find_brkraw_nifti(output_folder)
    affected_regions_prefix = '%s_%s_' % (nifti_name_without_extension(brkraw_nifti_file),
                                          nifti_name_without_extension(ara_template_file))
    output_file = os.path.join(output_folder, affected_regions_prefix + 'affectedRegions_Parental.nii.gz')
    nii.save(affected_regions_img, output_file)

    # Stroke volume calculation
    brain_mask = brain_volume.copy()
    brain_mask[brain_mask > 0.0] = 1.0
    brain_mask[brain_mask <= 0.0] = 0.0

    # Estimate volumes from binary voxel count times each image's voxel volume.
    strokeVoxelVolumeMM3 = voxel_volume_mm3(stroke_mask_img)
    brainVoxelVolumeMM3 = voxel_volume_mm3(brain_img)
    strokeVolumeInCubicMM = np.sum(stroke_mask > 0) * strokeVoxelVolumeMM3
    brainVolumeInCubicMM = np.sum(brain_mask > 0) * brainVoxelVolumeMM3
    if brainVolumeInCubicMM == 0:
        sys.exit("Error: Brain mask volume is zero.")

    # Load label names and write the text summary of affected parental regions.
    lines =open(os.path.join(REPO_ROOT, 'lib', 'annoVolume.nii.txt')).readlines()
    o=open(os.path.join(output_folder, affected_regions_prefix + 'affectedRegions_Parental.txt'), 'w')
    o.write("Stroke: %0.2f %% - Stroke Volume: %0.4f mm^3\n"  % (((strokeVolumeInCubicMM/brainVolumeInCubicMM)*100),strokeVolumeInCubicMM,))
    affected_label_names_by_id = {}
    labelNames = ["" for x in range(np.size(lines))]
    for i in range(len(lines)):
        label_id = int(lines[i].split('\t')[0])
        label_name = lines[i].split('\t')[1]
        labelNames[i] = label_name
        if label_id in region_percent_by_label and label_id in parental_label_id_set:
            # Write label ID, label name, and affected percentage for each matched region.
            o.write(lines[i][:-1] + "\t %0.2f %%\n" % region_percent_by_label[label_id])
            affected_label_names_by_id[label_id] = label_name

            #o.write(str(int(lines[i].split('	')[0]) + 2000) + '	R_' + lines[i].split('	')[1])
    o.close()

    # Store the same region statistics in a MATLAB file for downstream workflows.
    regionAffectPercent = np.array([region_percent_by_label[int(label_id)] for label_id in parental_label_ids])
    affected_label_names = [affected_label_names_by_id.get(int(label_id), "") for label_id in parental_label_ids]
    parental_label_ids = np.stack((parental_label_ids, regionAffectPercent))
    label_mat['ABLAbelsIDsParental'] = parental_label_ids
    label_mat['ABANamesPar'] = affected_label_names
    label_mat['ABAlabels'] = labelNames
    label_mat['volumePer'] = (strokeVolumeInCubicMM / brainVolumeInCubicMM) * 100
    label_mat['volumeMM'] = strokeVolumeInCubicMM
    sc.savemat(os.path.join(output_folder, 'labelCount_par.mat'), label_mat)






def find_files(input_folder, pattern):
    # Search the input folder recursively and return deterministic file ordering.
    return sorted(glob.glob(os.path.join(input_folder, '**', pattern), recursive=True))


def find_single_file(input_folder, pattern, description):
    matches = find_files(input_folder, pattern)
    if len(matches) == 0:
        sys.exit("Error: No %s found in '%s'." % (description, input_folder,))
    if len(matches) > 1:
        sys.exit("Error: Multiple %s found in '%s': %s" % (description, input_folder, ', '.join(matches),))
    return matches[0]

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Calculate incidence sizes of parental regions. You do not need to enter single files, but the path to the .../T2w folder')
    requiredNamed = parser.add_argument_group('Required named arguments')
    requiredNamed.add_argument('-i', '--inputFolder', help='.../T2w', required=True)

    parser.add_argument('-a', '--allenBrain_anno', help='File: parental Allen/ARA annotation template', nargs='?', type=str,
                        default=os.path.join(REPO_ROOT, 'lib', 'annoVolume.nii.gz'))

    input_folder = None
    allen_template_file = None
    output_folder = None

    args = parser.parse_args()

    # Use the input folder as both source folder and output folder.
    if args.inputFolder is not None:
        input_folder = args.inputFolder
        output_folder = args.inputFolder
    if not os.path.exists(input_folder):
        sys.exit("Error: '%s' is not an existing directory." % (input_folder,))


    if args.allenBrain_anno is not None:
        allen_template_file = args.allenBrain_anno
    if not os.path.isfile(allen_template_file):
        sys.exit("Error: '%s' is not an existing file." % (allen_template_file,))

    # Resolve static label resources from the repository lib folder.
    label_file = os.path.join(REPO_ROOT, 'lib', 'rsfMRILablelID.mat')
    ara_template_file = allen_template_file

    # Collect exactly one required subject file from the input folder.
    stroke_mask_file = find_single_file(input_folder, '*Stroke_mask.nii.gz', 'stroke mask')
    brain_file = find_single_file(input_folder, '*Bet.nii.gz', 'BET image')
    parental_annotation_file = find_single_file(input_folder, '*_AnnoSplit_parental.nii.gz', 'parental annotation')
    incidence_lesion_mask_file = find_single_file(input_folder, '*IncidenceData_Lesion_mask.nii.gz', 'incidence lesion mask')

    print("1 folder will be processed...")

    # Calculate parental-region lesion overlap and write NIfTI, TXT, and MAT outputs.
    calculate_parental_stroke_overlap(brain_file, parental_annotation_file, ara_template_file, stroke_mask_file,
                                      incidence_lesion_mask_file, output_folder, label_file)
