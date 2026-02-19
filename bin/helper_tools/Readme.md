## Atlas region size outputs (parental + standard annotations)

`getAtlasRegionSize.py` computes region-wise volumes from registered annotation NIfTI files and writes both **voxel counts** and **mm³** per region.

This script supports **two annotation types** and produces separate outputs:

### 1) Parental annotation (larger labels)
**Input pattern**
- `**/*_AnnorsfMRI.nii.gz`

**Mask pairing**
- `*_AnnorsfMRI.nii.gz` → `*_mask.nii.gz`

**Outputs**
- `region_size_mm_par.txt`
- `region_size_mm_par.mat`

These `_par` outputs correspond to the **parental** atlas/annotation (coarser / larger labels).

### 2) Standard annotation
**Input pattern**
- `**/*_Anno.nii.gz`

**Mask pairing**
- `*_Anno.nii.gz` → `*_mask.nii.gz`

**Outputs**
- `region_size_mm.txt`
- `region_size_mm.mat`

These outputs correspond to the **standard** atlas/annotation.

---

### Volume computation

- Region voxel counts are computed from the annotation labels (`label_id` occurrences).
- Brain voxel count is computed from the paired mask (`mask > 0`). If the mask is missing, the annotation foreground (`label > 0`) is used as a fallback.
- Physical volumes are derived from a **forced voxel size**:

`0.068359 × 0.068359 × 0.5 mm`  
(voxel volume `= 0.0023364764 mm³`)

This voxel size and voxel volume are written into the header of each output `.txt`.

---

### Output format

Each `.txt` contains:

- a header line with total brain volume (vox + mm³)
- a table:

`ID  Name  Voxels  Unit  Volume  Unit`

Example:
