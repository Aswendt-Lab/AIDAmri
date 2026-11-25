"""Extract median DTI parameter values per brain region into a CSV file.

Read a 3D NIfTI parameter map (e.g. FA, MD, RD or AD) and an aligned
3D ROI labelmap. Calculate region medians in translation-table order when
a table is supplied, otherwise in ascending order of nonzero ROI labels.

Optionally exclude voxels marked 1 in an aligned binary stroke mask with
--exclude-stroke-mask PATH. Retain fully excluded ROIs with an empty median
cell and a warning identifying the subject, session and metric.
Append _exclude-stroke to the CSV basename when used.
Input images are not modified; remaining metric values must still be finite.

Optionally translate ROI IDs into acronyms using a whitespace-separated
text file. For AnnoSplit and AnnoSplit_parental atlases, expand unsplit
translations into left (L_, original ID) and right (R_, ID + 2000)
hemisphere entries: all left entries first, then all right entries, each in
table order. Preserve the order of already hemisphere-labelled translations.
Keep table entries absent from the atlas with empty median cells and no warning.
Warn and skip atlas regions missing from a supplied translation table.
Always omit background ID 0.

Write roi_id, acronym and median columns to a CSV beside the input
image, using the image basename and atlas marker. Without a translation
table, leave the acronym column empty.

Validate image dimensions, voxel-grid alignment, integer ROI labels
and finite image values in exported regions before opening the CSV.
CLI messages are shown in the terminal and appended to the enclosing
dwi folder's analysis log, preserving stdout and stderr separately.
Use the batchProc.py timestamp format and Europe/Berlin time zone.
Outside a dwi folder, use an analysis log beside the image instead.
Create the analysis log if needed and append to it on subsequent runs.
Leave existing process logs unchanged.
For identifiable subject/dwi folders, record created and modified outputs
in the dwi manifest under stage extraction. Otherwise warn and skip tracking.
When a parent batch owns extraction tracking for the same dwi folder,
keep logging but leave manifest tracking to that parent.
Append Git information at the end of each standalone run, including failures.
During batch runs, leave the header and final Git block to the parent once per session.

Usage:
    python DTIdata_extract.py -i IMAGE -r ROIS [-t TRANSLATIONS] [--exclude-stroke-mask MASK]

Created on 06.04.2019 by Niklas Pallast.
Modified by J. Cramer 15-10-2026
latest Update: Julian Cramer
"""

import argparse
from contextlib import redirect_stderr, redirect_stdout
import csv
import os
from pathlib import Path
import sys

import nibabel as nii
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.artifact_manifest import manifest_filename, start_output_tracking
from common.dti_csv import getOutfile, nifti_stem as _nifti_stem
from common.script_logging import (
    _TeeStream, append_git_information, build_script_log_path, format_script_log_header,
)

PARENT_TRACKING_ENV = "AIDAMRI_DTI_EXTRACTION_TRACKED_BY_PARENT"
PARENT_GIT_LOG_ENV = "AIDAMRI_DTI_GIT_LOGGED_BY_PARENT"
PARENT_HEADER_LOG_ENV = "AIDAMRI_DTI_HEADER_LOGGED_BY_PARENT"


def _analysis_log_path(image_file):
    image = Path(image_file).expanduser().resolve()
    folder = next((parent for parent in image.parents if parent.name == "dwi"), image.parent)
    return Path(build_script_log_path(folder, "analysis.log"))


def _extraction_context(image_file):
    """Identify subject, session and metric from an optional input NIfTI path."""
    subject = session = metric = "unknown"
    if image_file is not None:
        image = Path(image_file).expanduser().resolve()
        dwi = next((parent for parent in image.parents if parent.name == "dwi"), None)
        if dwi is not None:
            subject = dwi.parent.parent.name or "unknown"
            session = dwi.parent.name or "unknown"
        try:
            suffix = _nifti_stem(image).rsplit(".", 1)[-1].lower()
        except ValueError:
            pass
        else:
            metric = {"fa": "FA", "fa0": "FA", "md": "MD", "rd": "RD", "ad": "AD"}.get(suffix, "unknown")
    return f"subject={subject}, session={session}, metric={metric}"


def _parent_owns_log_section(folder, env_name):
    parent_folder = os.environ.get(env_name, "")
    try:
        return bool(parent_folder) and Path(parent_folder).is_absolute() \
            and Path(parent_folder).resolve() == folder
    except (OSError, ValueError):
        return False


def _finalize_tracking(tracker):
    try:
        tracker.finalize()
    except (OSError, ValueError) as error:
        print(f"Error: Cannot finalize extraction manifest in '{tracker.folder}': {error}", file=sys.stderr)
        return False
    return True


def _split_mapping(mapping, txt_file):
    """Expand unsplit atlas IDs, or keep an already hemisphere-labelled table."""
    mapping = {roi_id: label for roi_id, label in mapping.items() if roi_id != 0}
    labelled = [label.startswith(("L_", "R_")) for label in mapping.values()]
    if any(labelled):
        if not all(labelled):
            raise ValueError(f"Mixed hemisphere-labelled and unlabelled entries in '{txt_file}'.")
        return mapping

    split = {}
    for offset, prefix in ((0, "L_"), (2000, "R_")):
        for roi_id, label in mapping.items():
            hemisphere_id = roi_id + offset
            if hemisphere_id in split:
                raise ValueError(
                    f"Hemisphere ID collision for ROI ID {hemisphere_id} in '{txt_file}'."
                )
            split[hemisphere_id] = prefix + label
    return split


def _read_mapping(txt_file, split_hemispheres=False):
    mapping = {}
    with open(txt_file, "r", encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, start=1):
            if not line.strip():
                continue
            parts = line.strip().split(maxsplit=1)
            location = f"'{txt_file}', line {line_number}"
            if len(parts) != 2:
                raise ValueError(f"Missing ROI ID or label in {location}.")
            try:
                roi_id = int(parts[0])
            except ValueError as error:
                raise ValueError(f"Invalid integer ROI ID '{parts[0]}' in {location}.") from error
            if roi_id in mapping and mapping[roi_id] != parts[1]:
                raise ValueError(f"Conflicting labels for ROI ID {roi_id} in {location}.")
            mapping[roi_id] = parts[1]
    mapping.pop(0, None)
    return _split_mapping(mapping, txt_file) if split_hemispheres else mapping


def extractDTIData(img, rois, outfile, txt_file, roi_file=None, stroke_mask=None, image_file=None):
    """Use table order when supplied; leave missing and fully excluded ROI medians blank."""
    if img.ndim != 3 or rois.ndim != 3:
        raise ValueError(f"Image and ROI data must be 3D; got {img.shape} and {rois.shape}.")
    if img.shape != rois.shape:
        raise ValueError(f"Image shape {img.shape} and ROI shape {rois.shape} do not match.")
    if not np.all(np.isfinite(rois)):
        raise ValueError("ROI labels must be finite integers; found NaN or infinity.")
    if not np.issubdtype(rois.dtype, np.integer) and not np.all(rois == np.floor(rois)):
        raise ValueError("ROI labels must be integers; found fractional labels.")
    if stroke_mask is not None:
        if stroke_mask.ndim != 3:
            raise ValueError(f"Stroke mask must be 3D; got {stroke_mask.shape}.")
        if stroke_mask.shape != img.shape:
            raise ValueError(f"Image shape {img.shape} and stroke mask shape {stroke_mask.shape} do not match.")
        if not np.all((stroke_mask == 0) | (stroke_mask == 1)):
            raise ValueError("Stroke mask must be binary with finite values 0 or 1.")

    split_hemispheres = roi_file is not None and "annosplit" in _nifti_stem(roi_file).lower()
    mapping = _read_mapping(txt_file, split_hemispheres) if txt_file is not None else None
    regions = np.unique(rois)
    regions = regions[regions != 0]
    present_ids = {int(region) for region in regions}
    if mapping is not None:
        for region in regions:
            if int(region) not in mapping:
                print(
                    f"Warning: No acronym for ROI ID {int(region)} in '{txt_file}'. Skipping region.",
                    file=sys.stderr,
                )
        regions = mapping
    rows = []
    for region in regions:
        roi_id = int(region)
        acronym = mapping[roi_id] if mapping is not None else ""
        if roi_id not in present_ids:
            rows.append((roi_id, acronym, None))
            continue
        selection = rois == region
        if stroke_mask is not None:
            selection &= stroke_mask == 0
        values = img[selection]
        if not np.all(np.isfinite(values)):
            raise ValueError(f"Image contains NaN or infinity in ROI ID {roi_id}.")
        if values.size == 0:
            print(
                f"Warning: [{_extraction_context(image_file)}] All voxels in ROI ID {roi_id} "
                "excluded by stroke mask. Leaving median empty.",
                file=sys.stderr,
            )
            median = None
        else:
            median = float(np.median(values))
        rows.append((roi_id, acronym, median))

    with open(outfile, "w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("roi_id", "acronym", "median"))
        writer.writerows(rows)
    return outfile


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Extract median DTI parameter values (FA, MD, RD, AD) per ROI into a CSV "
                    "beside the input image, using the image name and ROI atlas marker."
    )
    parser.add_argument("-i", "--image_file", required=True,
                        help="3D NIfTI DTI parameter map e.g. FA Nifti (.nii or .nii.gz)")
    parser.add_argument("-r", "--roi_file", required=True,
                        help="3D NIfTI ROI labelmap aligned with the parameter map e.g registered atlas")
    parser.add_argument(
        "-t", "--translatorTXT",
        help="Optional text file containing integer ROI IDs and labels separated by whitespace. "
             "CSV rows follow table order, omitting ID 0; regions absent from the atlas have "
             "empty medians. For AnnoSplit atlases, unsplit tables expand to all L_ entries "
             "then all R_ entries in table order, with a +2000 right-hemisphere offset. "
             "Already L_/R_-labelled tables keep their existing order.",
        type=str,
    )
    parser.add_argument(
        "--exclude-stroke-mask", metavar="PATH", type=Path,
        help="Optional aligned binary 3D NIfTI stroke mask (.nii or .nii.gz): exclude voxels "
             "marked 1 from ROI medians. Fully excluded ROIs have an empty median cell. "
             "Append _exclude-stroke to the output CSV basename.",
    )
    args = parser.parse_args(argv)

    tracker = None
    finalized = False
    try:
        log_path = _analysis_log_path(args.image_file)
        folder = log_path.parent
        tracking_warning = ""
        tracking_error = ""
        managed_by_parent = False
        try:
            if folder.name != "dwi":
                raise ValueError("No enclosing dwi folder")
            manifest_filename(folder, "dwi")
        except ValueError:
            tracking_warning = (
                f"Warning: Manifest tracking skipped: no recognizable dwi/subject structure in '{folder}'."
            )
        else:
            parent_folder = os.environ.get(PARENT_TRACKING_ENV, "")
            try:
                managed_by_parent = bool(parent_folder) and Path(parent_folder).is_absolute() \
                    and Path(parent_folder).resolve() == folder
            except (OSError, ValueError):
                managed_by_parent = False
            if not managed_by_parent:
                try:
                    # Snapshot before opening the log so a new log is also tracked.
                    tracker = start_output_tracking(folder, "dwi", "extraction")
                except (OSError, ValueError) as error:
                    tracking_error = f"Error: Cannot start extraction manifest tracking in '{folder}': {error}"
        with log_path.open("a", encoding="utf-8", buffering=1) as log:
            with redirect_stdout(_TeeStream(sys.stdout, log)), redirect_stderr(_TeeStream(sys.stderr, log)):
                if not _parent_owns_log_section(folder, PARENT_HEADER_LOG_ENV):
                    print(f"\n{format_script_log_header('DTIdata_extract.py')}")
                print(f"Image: {Path(args.image_file).expanduser().resolve()}")
                print(f"ROI: {Path(args.roi_file).expanduser().resolve()}")
                print(f"Translation: {Path(args.translatorTXT).expanduser().resolve() if args.translatorTXT else 'None'}")
                if args.exclude_stroke_mask is not None:
                    print(f"Exclude stroke mask: {args.exclude_stroke_mask.expanduser().resolve()}")
                if managed_by_parent:
                    print("Manifest tracking: parent batch")
                if tracking_warning:
                    print(tracking_warning, file=sys.stderr)
                status = 1
                try:
                    if tracking_error:
                        print(tracking_error, file=sys.stderr)
                    else:
                        status = _run_extraction(args)
                finally:
                    try:
                        if tracker is not None:
                            log.flush()
                            finalized = True
                            if not _finalize_tracking(tracker):
                                status = 1
                    finally:
                        # Include tracking errors before the final Git block.
                        if not _parent_owns_log_section(folder, PARENT_GIT_LOG_ENV):
                            try:
                                append_git_information(log, output_dir=folder)
                            except (OSError, UnicodeError) as error:
                                print(f"Error: Cannot append Git information to '{log_path}': {error}",
                                      file=sys.stderr)
                                status = 1
                return status
    except OSError as error:
        print(f"Error: Cannot write analysis log for '{args.image_file}': {error}", file=sys.stderr)
        return 1
    finally:
        # Also release the tracker's atexit hook if opening/writing the log failed.
        if tracker is not None and not finalized:
            _finalize_tracking(tracker)


def _run_extraction(args):
    try:
        outfile = getOutfile(args.roi_file, args.image_file,
                             exclude_stroke=args.exclude_stroke_mask is not None)
        image_data = nii.load(args.image_file)
        roi_data = nii.load(args.roi_file)
        if not np.allclose(image_data.affine, roi_data.affine, rtol=0, atol=1e-5):
            raise ValueError("Image and ROI NIfTI affines do not match; inputs must share the same voxel grid.")
        img = image_data.get_fdata(dtype=np.float32)
        rois = np.asanyarray(roi_data.dataobj)
        stroke_mask = None
        if args.exclude_stroke_mask is not None:
            stroke_data = nii.load(args.exclude_stroke_mask.expanduser())
            if not np.allclose(image_data.affine, stroke_data.affine, rtol=0, atol=1e-5):
                raise ValueError("Image and stroke mask NIfTI affines do not match; "
                                 "inputs must share the same voxel grid.")
            stroke_mask = np.asanyarray(stroke_data.dataobj)
        extractDTIData(img, rois, outfile, args.translatorTXT, roi_file=args.roi_file,
                       stroke_mask=stroke_mask, image_file=args.image_file)
    except (OSError, ValueError, nii.filebasedimages.ImageFileError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1

    print(f"Done: {outfile}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
