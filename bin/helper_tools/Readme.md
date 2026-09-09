# AIDAmri Helper Tools

Utility scripts for AIDAmri data preparation, file cleanup, quality control,
atlas-volume summaries, and small batch support tasks.

Most scripts are intended to be run from this directory or from the pipeline
environment so that relative paths to `lib/` resources resolve as expected.

## Overview

| Script | Purpose |
| --- | --- |
| `DistributeStrokeMasks.py` | Propagate existing stroke masks across subject timepoints with NiftyReg `reg_resample`. |
| `MRI_files_summarizer.py` | Create a CSV inventory of NIfTI files found below `**/brkraw/*.nii.gz`. |
| `ReorientBatch.py` | Reorient NIfTI files to a target orientation while mirroring the input folder tree. |
| `adjustbvecRep.py` | Repeat DWI `.bval` and `.bvec` sidecars to match the number of image volumes. |
| `batch_qc_reports.py` | Python helper module for project-level BET and registration HTML QC reports. |
| `crop_T2.py` | Crop T2-weighted images in x/y using FSL through Nipype and write quick-look PNGs. |
| `fieldmap_json_edit.py` | Populate BIDS fieldmap JSON `IntendedFor` entries for DWI and functional files. |
| `getAtlasRegionSize_BIDS.py` | Compute per-annotation atlas region volumes in BIDS-style folder trees. |
| `getAtlasRegionSize_noBIDS.py` | Compute aggregate atlas region volumes for non-BIDS T2w folders. |
| `plot_sourcedata_niftis.py` | Generate per-session NIfTI mosaic PNGs and an HTML QC report. |
| `remove_carets_spaces.sh` | Remove spaces and caret (`^`) characters from a plain-text file. |
| `Reset_proc_folder.py` | Safely reset processing stages using runtime-generated output manifests. |
| `reset_naming.py` | Clean Bruker `subject` files before PV-to-NIfTI conversion. |

## Dependencies

The scripts use a mix of Python packages and external neuroimaging tools:

- Python packages used across the tools: `nibabel`, `numpy`, `pandas`,
  `matplotlib`, `scipy`, and `nipype`.
- FSL is required by `crop_T2.py` (`fslroi` via Nipype and optional `slicer`
  PNG generation).
- NiftyReg is required by `DistributeStrokeMasks.py` (`reg_resample` on
  `PATH`).
- The atlas region scripts expect lookup resources under the repository-level
  `lib/` folder, resolved relative to the current working directory:
  `lib/ABALabelsIDchanged.mat`, `lib/annoVolume+2000_rsfMRI.nii.txt`, and,
  when present, `lib/ARA_annotationR+2000.nii.txt`.

## Atlas Region Size

### `getAtlasRegionSize_BIDS.py`

Computes region-wise voxel counts and volumes for each matching annotation file
found recursively under a BIDS-style input folder. Results are written to:

```text
<inputFolder>/output_region_size/
```

Supported annotation patterns:

- `**/*_Anno_parental.nii.gz`
- `**/*_AnnoSplit_parental.nii.gz`
- `**/*_AnnorsfMRI.nii.gz`
- `**/*_Anno.nii.gz`
- `**/*_AnnoSplit.nii.gz`

Mask files are inferred by replacing the annotation suffix with
`_mask.nii.gz`. If a mask is missing, the script falls back to annotation
foreground voxels for brain volume.

Usage:

```bash
python getAtlasRegionSize_BIDS.py -i /path/to/input_folder
```

Outputs:

- One `.txt` and one `.mat` file per annotation.
- Parental outputs receive an additional `_par` suffix before the extension.
- Volumes use the forced voxel size `0.068359 x 0.068359 x 0.5 mm`.

### `getAtlasRegionSize_noBIDS.py`

Computes aggregate region-wise voxel counts and volumes for a non-BIDS T2w
folder.

Usage:

```bash
python getAtlasRegionSize_noBIDS.py -i /path/to/T2w_folder
```

Inputs and outputs:

- Parental input pattern: `**/*_AnnorsfMRI.nii.gz`
- Regular input pattern: `**/*_Anno.nii.gz`
- Paired masks are expected as `_mask.nii.gz`.
- Writes `region_size_mm_par.txt/.mat` and/or `region_size_mm.txt/.mat`
  directly into the input folder.
- Volumes use the forced voxel size `0.068359 x 0.068359 x 0.5 mm`.

## Bruker Subject File Cleanup

### `reset_naming.py`

Prepares Bruker raw data before PV-to-NIfTI conversion. The script recursively
finds files named `subject` below the input folder and edits them in place.

Main operations:

- Removes the first underscore in values following `##$SUBJECT_id=` and
  `##$SUBJECT_study_name=`.
- Replaces study-name values matching `baseline...` with `PT0`
  case-insensitively.

Usage:

```bash
python reset_naming.py -i /path/to/raw_data
```

### `remove_carets_spaces.sh`

Removes spaces and caret (`^`) characters from each line of a plain-text input
file, typically a Bruker `subject` file.

Usage:

```bash
bash remove_carets_spaces.sh /path/to/subject
```

Outputs are written in the current working directory:

- `subject_clean.txt`: cleaned lines with spaces and caret characters removed.
- `subject_orig.txt`: copy of the original input file.

Note: despite the historical script comments, the current implementation does
not remove underscores and does not overwrite the original `subject` file.

## Manifest-based Processing Reset

### `Reset_proc_folder.py`

The anat, dwi, func, and t2map stage scripts take a recursive snapshot of their
modality folder before processing and update a dataset-specific manifest when
the process exits. For example, `sub-01/ses-02/anat` uses
`.sub-01_ses-02_anat_aidamri_manifest.json`. Tracking also works when a stage script is
launched directly instead of through `batchProc.py`.

The reset helper deletes only files registered as newly created by stages after
the requested target phase. Standard BIDS sidecars matching `*_T2w.json`, `*_EPI.json`, `*_dwi.json`, 
`*.bvec`, or `*.bval` and `*Stroke_mask.nii.gz` are not deleted.
Existing files that were modified are reported but
not deleted. Files without any manifest assignment are explicitly reported as
not managed and are preserved. Tracked directories are removed only when
empty. Before planning a reset, every existing folder of the selected mode must
contain its readable, correctly named manifest. For example,
`--mode anat` validates only `anat` folders and does not inspect manifests in
`dwi`, `func`, or `t2map`. A missing or invalid manifest in the selected mode
aborts the complete reset before anything is deleted.

Do not run a reset concurrently with preprocessing,
registration, processing, or `batchProc.py`.

At the end of every dry-run or applied reset, warnings and unmanaged files are
shown again in a consolidated summary grouped by subject, session, and mode.
Modified files and unmanaged files are listed with their full paths so the user
can inspect and handle them manually.

Use `--delete-unmanaged` only after reviewing a dry-run. It permanently deletes
files that are not referenced anywhere in the manifest. Standard BIDS sidecars matching 
`*_T2w.json`, `*_EPI.json`, `*_dwi.json`, `*.bvec`, or `*.bval` and `*Stroke_mask.nii.gz` 
are always protected from this deletion. Other manually supplied masks, notes, and user data can still be
deleted. Unmanaged directories are never removed recursively: their files are
deleted individually and the directories are removed only when empty. A
directory containing a protected file is preserved.

```bash
python Reset_proc_folder.py --input /path/to/project \
  --mode anat \
  --phase preprocessing \
  --delete-unmanaged \
  --dry-run
```

Always inspect a dry-run first:

```bash
python Reset_proc_folder.py --input /path/to/project \
  --mode anat \
  --phase preprocessing \
  --dry-run
```

Apply the displayed plan interactively:

```bash
python Reset_proc_folder.py --input /path/to/project \
  --mode anat \
  --phase preprocessing
```

Use `--yes` only for an already reviewed automated invocation. Projects that
were processed before manifest tracking was introduced cannot be reconstructed
retrospectively. By default, their untracked files remain untouched.

## Reorientation

### `ReorientBatch.py`

Batch-reorients `.nii` and `.nii.gz` files under an input root while mirroring
the folder structure to an output root.

Key behavior:

- Uses the active NIfTI transform, preferring `sform` over `qform`, and falls
  back to `img.affine`.
- Reorients NIfTI data to a target orientation.
- Default interactive target is AIDAmri's `LIP`
  (`Left-Inferior-Posterior`) orientation.
- Copies non-NIfTI files unchanged.
- Skips `.bval` and `.bvec` files during tree traversal because they are
  handled together with the matching NIfTI file.
- Copies `.bval` sidecars and reorients FSL-style `.bvec` sidecars when the
  image orientation changes.
- Writes a log file into the output root.
- If every NIfTI already matches the target orientation, the script logs this
  and aborts before copying the whole tree.
- Processes files in parallel. By default, it uses `50%` of available CPU cores.

Usage:

```bash
python ReorientBatch.py -i /path/to/input_root -o /path/to/output_root
```

Non-interactive usage:

```bash
python ReorientBatch.py \
  -i /path/to/input_root \
  -o /path/to/output_root \
  -t LIP \
  -n
```

Options:

- `-t ORI`: target orientation, for example `LIP`, `RAS`, or `LAS`.
- `-n`: non-interactive mode. Requires `-t`.
- `-l LOGFILE`: log filename written into the output root. Default:
  `reorient_log.txt`.
- `-p`, `--cpu-percent`: CPU percentage for parallel processing, for example
  `50` or `50%`. Default: `50`.

## Quality Control and Data Summaries

Generated report timestamps use German local time in the format
`DD Month YYYY HH:MM:SS CET/CEST`, for example `15 July 2026 15:14:12 CEST`.

### `adjustbvecRep.py`

Adjusts DWI sidecars in a session folder. The script changes into
`<session>/dwi`, loads files matching `*dwi.nii.gz`, and repeats `.bval` and
`.bvec` contents when the image volume count is an integer multiple of the
sidecar length.

Usage:

```bash
python adjustbvecRep.py /path/to/sub-001/ses-001
```

Special cases:

- Single-volume DWIs are moved with sidecars into `.single_volume/`.
- DWIs with fewer than 12 volumes are moved into `.low_volume_count/`.
- DWIs with exactly 23 volumes are moved into `.volumes_23_count/`.
- JSON sidecars are moved together with skipped DWI files when present.

This helper is called automatically by `conv2Nifti_auto.py` during conversion
when DWI data are present.

### `MRI_files_summarizer.py`

Creates a CSV overview of NIfTI files below `**/brkraw/*.nii.gz`. Filename
tokens are used to extract subject, session/timepoint, run, and modality
information.

Usage:

```bash
python MRI_files_summarizer.py -i /path/to/project -o /path/to/output_folder
```

Output:

```text
<output_folder>/MRI_files_overview.csv
```

CSV columns:

- `FileAddress`
- `Modality`
- `TimePoint`
- `SubjectID`
- `RunNumber`

### `plot_sourcedata_niftis.py`

Creates mosaic PNG files and an HTML report for NIfTI files in one
`sub-*/ses-*` folder.

Usage:

```bash
python plot_sourcedata_niftis.py /path/to/sub-001/ses-001 /path/to/qc_output
```

Optional argument:

- `--n_slices`: number of slices per orientation. Default: `10`.

Behavior:

- Scans `dwi`, `anat`, and `func` directories.
- Also scans subfolders inside `dwi`.
- Uses the first volume for 4D images.
- Adds orientation labels from the NIfTI affine.
- Writes `*_report.png` images and
  `sub-<subject>_ses-<session>_report.html` into the output directory.

### `batch_qc_reports.py`

Python helper module for project-level QC reports. 
Direct command-line use:

```bash
python bin/helper_tools/batch_qc_reports.py -i /path/to/proc_data
python bin/helper_tools/batch_qc_reports.py -i /path/to/proc_data --report bet --n-slices 7
python bin/helper_tools/batch_qc_reports.py -i /path/to/proc_data --report registration
python bin/helper_tools/batch_qc_reports.py -i /path/to/proc_data --custom-parameter t2-frac=0.1 --custom-parameter t2-bias-method=mico
```

`--report` accepts `all` (the default), `bet`, or `registration`. `--n-slices`
sets the number of slices per orientation and defaults to `10`. Repeat
`--custom-parameter NAME=VALUE` to record processing parameters in the custom
parameters section of the generated HTML reports. Parameter names without a
leading `--` are normalized automatically.

Import use:

Available functions:

```python
from batch_qc_reports import build_bet_qc_report, build_registration_qc_report

build_bet_qc_report(
    "/path/to/proc_data",
    n_slices=10,
    custom_parameters=[("--t2-frac", 0.1)],
)
build_registration_qc_report(
    "/path/to/proc_data",
    n_slices=10,
    custom_parameters=[("--t2-frac", 0.1)],
)
```

When reports are created by `batchProc.py`, options explicitly supplied on the
command line are listed in a **Custom parameters** section at the top of each
HTML report. The required `--input` option is omitted. Direct callers can pass
the optional `custom_parameters` sequence shown above.

BET report behavior:

- Searches `sub-*/ses-*/*/*Bet.nii.gz`.
- Skips files ending in `_mask.nii.gz`.
- Overlays a matching `*_mask.nii.gz` contour when present and shape-compatible.
- Writes PNGs and `bet_report.html` under `<project_dir>/Report/BET/`.

Registration report behavior:

- Searches `sub-*/ses-*/*/*_AnnoSplit_parental.nii.gz`.
- Expects the matching BET file to have the same filename without the
  `_AnnoSplit_parental` suffix.
- Overlays annotation labels on the BET image.
- Writes PNGs and `registration_report.html` under
  `<project_dir>/Report/Registration/`.

## T2 Cropping

### `crop_T2.py`

Crops a T2-weighted image with FSL `fslroi` through Nipype and writes quick-look
PNG images with FSL `slicer` when available.

Usage:

```bash
python crop_T2.py \
  -i /path/to/input_T2w.nii.gz \
  -x_max 120 \
  -y_max 120 \
  -z_max 40 \
  -o /path/to/cropped_T2w.nii.gz
```

Options:

- `-x_min`, default `0`
- `-x_max`, required
- `-y_min`, default `0`
- `-y_max`, required
- `-z_min`, default `0`
- `-z_max`, required by the CLI

Important current behavior:

- The script currently sets `crop_z = False`, so z cropping is disabled in the
  implementation. `z_min` is reset to `0` and `z_size` uses the full input
  z-dimension, even though `-z_max` is still required by the parser.
- The output image contains one timepoint (`t_min=0`, `t_size=1`).
- PNG quick-look files are written next to the input and output NIfTI files.

## Fieldmap JSON Editing

### `fieldmap_json_edit.py`

Updates BIDS fieldmap JSON files with `IntendedFor` entries for matching DWI and
functional files in a participant/session folder.

Usage:

```bash
python fieldmap_json_edit.py \
  --bids_root /path/to/proc_data \
  --participant 001 \
  --session PT0
```

Optional argument:

- `--overwrite`: replace existing non-empty `IntendedFor` fields.

Behavior:

- Looks for fieldmap JSON files in
  `<bids_root>/sub-<participant>/ses-<session>/fmap/*.json`.
- Adds relative paths from the BIDS root.
- DWI targets are matched in the session `dwi/` folder.
- Functional targets are matched as `*epi.nii*` in the session `func/` folder.

## Stroke Mask Distribution

### `DistributeStrokeMasks.py`

Finds existing stroke masks and propagates them to other timepoints of the same
subject by resampling with NiftyReg.

Search pattern:

```text
**/anat/*Stroke_mask.nii.gz
```

Workflow:

1. Resample each stroke mask into incidence space using
   `*IncidenceData.nii.gz` and `*MatrixInv.txt` from the source `anat` folder.
2. For each other timepoint of the same subject, resample from incidence space
   to that timepoint using `*MatrixBspline.nii` and a `*BiasBet.nii.gz`
   reference.
3. Skip target timepoints that already contain a `*Stroke_mask.nii.gz`.
4. Append missing incidence or affine-matrix information to
   `missing_files_log.txt` in the input root.

Usage:

```bash
python DistributeStrokeMasks.py -i /path/to/dataset_root
```

Outputs:

- Intermediate `*_StrokeM_IncidenceSpace.nii.gz` files in source `anat`
  folders.
- New `*Stroke_mask.nii.gz` files in target timepoint `anat` folders that do
  not already have a stroke mask.
- `missing_files_log.txt` in the input root when required files are missing.

Requirement:

- `reg_resample` must be available on `PATH`.
