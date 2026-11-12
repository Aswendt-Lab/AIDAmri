"""Shared filenames for per-region DTI CSV exports."""

import os


def nifti_stem(path):
    name = os.path.basename(path)
    for suffix in (".nii.gz", ".nii"):
        if name.lower().endswith(suffix):
            return name[:-len(suffix)]
    raise ValueError(f"Expected a .nii or .nii.gz file: '{path}'.")


def getOutfile(roi_file, img_file, exclude_stroke=False):
    """Append the ROI atlas marker and optional stroke-exclusion suffix."""
    image_stem = nifti_stem(img_file)
    roi_stem = nifti_stem(roi_file)
    atlas = roi_stem
    for marker in ("AnnoSplit_parental", "Anno_parental", "AnnoSplit", "Anno"):
        if marker.lower() in roi_stem.lower():
            atlas = marker
            break
    suffix = "_exclude-stroke" if exclude_stroke else ""
    return os.path.join(os.path.dirname(img_file), f"{image_stem}_{atlas}{suffix}.csv")
