# 🧠✨ AIDAmri helper tools (`bin/helper_tools`)

Small utility scripts for data preparation, QC, and batch operations used in AIDAmri workflows.

---

## 🧰 Scripts at a glance

- 🧠 **getAtlasRegionSize.py** — Compute region-wise volumes from registered annotation NIfTI files (parental + standard annotations).
- 🧬 **getAtlasRegionSize_BIDS.py** — BIDS variant: per-file outputs for all annotation variants (parental + standard).
- 🏷️ **reset_naming.py** — Fix Bruker `subject` files (naming cleanup) before PV-to-NIfTI conversion.
- 🧭 **ReorientBatch.py** — Batch-reorient NIfTI files to a target orientation (default: **LIP**) while preserving folder structure; handles FSL bvecs/bvals.
- 📋 **MRI_files_summarizer.py** — Create a CSV inventory of NIfTI files under `**/brkraw/*.nii.gz` with basic fields parsed from filenames.
- 🩹 **DistributeStrokeMasks.py** — Resample/propagate existing `*Stroke_mask.nii.gz` across timepoints using `reg_resample`.
- 🕵️ **searchmissingClyinder.py** — From a `files.txt` list, report study IDs missing expected time windows (BL/P3/P7/P14/P28/P56).

---

<details>
<summary>🧠 <strong>getAtlasRegionSize.py</strong> — Atlas region size outputs (parental + standard annotations)</summary>

> 📍 **Location tip:** Place both atlas scripts in `./bin/3.1_T2Processing/` so relative paths to `./lib/` resolve as expected.

`getAtlasRegionSize.py` computes **region-wise volumes** from registered annotation NIfTI files and writes both **voxel counts** and **mm³** per region.

### ✅ Supported annotation types and outputs

#### 1) 👪 Parental annotation (larger / coarser labels)

**Input pattern**
- `**/*_AnnorsfMRI.nii.gz`

**Mask pairing**
- `*_AnnorsfMRI.nii.gz` → `*_mask.nii.gz`

**Outputs**
- `region_size_mm_par.txt`
- `region_size_mm_par.mat`

These `_par` outputs correspond to the **parental atlas/annotation** (coarser / larger labels).

#### 2) 🧩 Standard annotation

**Input pattern**
- `**/*_Anno.nii.gz`

**Mask pairing**
- `*_Anno.nii.gz` → `*_mask.nii.gz`

**Outputs**
- `region_size_mm.txt`
- `region_size_mm.mat`

These outputs correspond to the **standard atlas/annotation**.

### 📐 Volume computation

- **Region voxel counts**: computed from annotation labels (occurrences of each `label_id`).
- **Brain voxel count**: computed from the paired mask (`mask > 0`).
  - If the mask is missing, the script falls back to annotation foreground (`label > 0`).

**Physical volumes** are derived from a **forced voxel size**:
- `0.068359 × 0.068359 × 0.5 mm`
- voxel volume: `0.0023364764 mm³`

This voxel size and voxel volume are written into the header of each output `.txt`.

### 🧾 Output format

Each `.txt` contains:
- a header line with **total brain volume** (vox + mm³)
- a table:

`ID  Name  Voxels  Unit  Volume  Unit`

</details>

---

<details>
<summary>🧬 <strong>getAtlasRegionSize_BIDS.py</strong> — BIDS region size outputs (all annotation variants, per-file outputs)</summary>

> 📍 **Location tip:** Place both atlas scripts in `./bin/3.1_T2Processing/` so relative paths to `./lib/` resolve as expected.

`getAtlasRegionSize_BIDS.py` is the BIDS-oriented variant. It searches **recursively** under the given input folder and produces **one output pair (.txt + .mat) per annotation file**, writing everything into an `output_region_size/` folder inside the input directory.

### 📂 Output folder

All outputs are written to:
- `<inputFolder>/output_region_size/`

### 🧾 Supported annotation variants (processed all)

#### 👪 Parental variants (larger / coarser labels)

**Input patterns**
- `**/*_Anno_parental.nii.gz`
- `**/*_AnnoSplit_parental.nii.gz`
- `**/*_AnnorsfMRI.nii.gz` (legacy parental naming)

**Mask pairing**
- `<annotation>.nii.gz` → `<annotation_basename>_mask.nii.gz`  
  (implemented by replacing the recognized suffix with `_mask.nii.gz`)

**Outputs (per input file)**
- `<basename>_par.txt`
- `<basename>_par.mat`

Parental outputs preserve **left/right hemisphere region IDs**, including the common **+2000 right-hemisphere convention**, by filtering region IDs using the **parental lookup table keys**.

#### 🧩 Standard variants

**Input patterns**
- `**/*_Anno.nii.gz`
- `**/*_AnnoSplit.nii.gz`

**Mask pairing**
- `<annotation>.nii.gz` → `<annotation_basename>_mask.nii.gz`

**Outputs (per input file)**
- `<basename>.txt`
- `<basename>.mat`

### 📐 Volume computation

- **Region voxel counts**: computed from annotation labels (occurrences of each `label_id`).
- **Brain voxel count**: computed from the paired mask (`mask > 0`).
  - If the mask is missing, the script falls back to annotation foreground (`label > 0`).

**Physical volumes** are derived from a **forced voxel size**:
- `0.068359 × 0.068359 × 0.5 mm`
- voxel volume: `0.0023364764 mm³`

This voxel size and voxel volume are written into the header of each output `.txt`.

### 🧾 Output format

Each per-file `.txt` contains:
- a header with **brain volume** (vox + mm³) for that file
- paths to the annotation and mask used
- a table:

`ID  Name  Voxels  Unit  Volume  Unit`

### ▶️ Usage

```bash
python getAtlasRegionSize_BIDS.py -i /path/to/sub-*/ses-*
```

</details>

---

<details>
<summary>🏷️ <strong>reset_naming.py</strong> — Fix Bruker <code>subject</code> files (pre PV-to-NIfTI)</summary>

Prepares Bruker raw data before running `1_PV2NIfTiConverter/pv_conv2Nifti.py`.

### 🔧 What it does (in-place)

1. Removes the **first underscore** `_` in `SUBJECT_id` and `SUBJECT_study_name` lines.
2. Replaces `baseline` (case-insensitive) with `PT0` in the study name.

### ▶️ Usage

```bash
python reset_naming.py -i /path/to/raw_data
```

### 📥 Inputs

- `-i / --input`: parent folder containing Bruker raw data. The script expects a structure like:
  - `projectfolder/subjects/ses/data` (and scans recursively for files named `subject`)

### 📤 Outputs

- Updated `subject` files written back to disk (in-place).

</details>

---

<details>
<summary>🧭 <strong>ReorientBatch.py</strong> — Batch reorient NIfTI files (default target: LIP)</summary>

Batch reorientation of `.nii` / `.nii.gz` under an input root while mirroring the directory structure to an output root.

### ✅ Key behavior

- Reorients NIfTI images to a target orientation (default intended for AIDAmri: **LIP**).
- Copies non-NIfTI files unchanged.
- If diffusion sidecars exist:
  - `.bval` files are copied
  - `.bvec` files are reoriented consistently with the image orientation change

### 🧩 Requirements

- Python packages: `nibabel`, `numpy`

### ▶️ Usage

```bash
python ReorientBatch.py -i /path/to/input_root -o /path/to/output_root
```

### ⚙️ Common options

- `-t / --target`: target orientation (3 letters from `{L,R,A,P,S,I}`), e.g. `LIP`
- `-n / --non_interactive`: run without prompts (requires `-t`)
- `-l / --logfile`: name of the log file written into the output root

### 📤 Outputs

- Reoriented dataset under `output_root` with the same relative structure as input.
- A log file summarizing processing.

</details>

---

<details>
<summary>📋 <strong>MRI_files_summarizer.py</strong> — Create a CSV overview of NIfTI files in <code>brkraw</code></summary>

Scans for:
- `**/brkraw/*.nii.gz`

### 🧾 Parsed fields (from filename tokens)

- `SubjectID` from `sub-...`
- `TimePoint` from `ses-...`
- `RunNumber` from `run-...`
- `Modality` from the final token before `.nii.gz`

### 🧩 Requirements

- Python package: `pandas`

### ▶️ Usage

```bash
python MRI_files_summarizer.py -i /path/to/project -o /path/to/output_folder
```

### 📤 Outputs

- `MRI_files_overview.csv` written to the output folder.

</details>

---

<details>
<summary>🩹 <strong>DistributeStrokeMasks.py</strong> — Distribute stroke masks across timepoints (via <code>reg_resample</code>)</summary>

Finds existing stroke masks and propagates them to other timepoints by resampling.

### 🔎 Search pattern

- `**/anat/*Stroke_mask.nii.gz`

### 🔁 Workflow (per found mask)

1. Resample mask into incidence space using:
   - `*IncidenceData.nii.gz`
   - `*MatrixInv.txt`
2. For each other timepoint folder, resample from incidence space into that timepoint using:
   - `*MatrixBspline.nii`
   - reference `*BiasBet.nii.gz`
3. Writes a log `missing_files_log.txt` when required files are missing.

### 🧩 Requirements

- `reg_resample` available on PATH (e.g., from NiftyReg)

### ▶️ Usage

```bash
python DistributeStrokeMasks.py -i /path/to/subject_root_or_dataset_root
```

### 📤 Outputs

- New `*Stroke_mask.nii.gz` files in timepoints that do not already have a stroke mask.
- `missing_files_log.txt` in the input root.

</details>

---

<details>
<summary>🕵️ <strong>searchmissingClyinder.py</strong> — Report missing time windows from a <code>files.txt</code> list</summary>

Reads a local `files.txt` (in the current working directory), extracts `studyID` and timepoint tokens from filenames, and reports which studies are missing expected time windows.

### 🗓️ Expected windows

- `BL`
- `P3` (accepts P2/P3/P4)
- `P7` (accepts P6/P7/P8)
- `P14` (accepts P13/P14/P15)
- `P28` (accepts P27–P32)
- `P56` (accepts P55/P56/P57)

### ▶️ Usage

```bash
python searchmissingClyinder.py
```

### 📥 Inputs

- `files.txt` in the current directory (one filename per line)

### 📤 Outputs

- Printed report to stdout, e.g.:
  - `GV_T3_...: missing BL, P7, ...`

</details>
