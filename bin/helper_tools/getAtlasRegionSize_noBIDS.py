"""
Created on 10/08/2017
Modified on 04/23/2021
Updated on 02/19/2026

@author: Niklas Pallast and Markus Aswendt
Neuroimaging & Neuroengineering
Department of Neurology
University Hospital Cologne

#!/usr/bin/env python3

Two outputs:

1) Parental annotation (larger labels):
   - Input patterns: **/*_Anno_parental.nii.gz, **/*_AnnoSplit_parental.nii.gz,
     **/*_AnnorsfMRI.nii.gz
   - Mask pairing: known annotation suffix -> _mask.nii.gz
   - Output: region_size_mm_par.txt / region_size_mm_par.mat
   - Uses the SIGMA lookup-table keys for filtering IDs

2) Standard annotation:
   - Input patterns: **/*_Anno.nii.gz, **/*_AnnoSplit.nii.gz
   - Mask pairing: known annotation suffix -> _mask.nii.gz
   - Output: region_size_mm.txt / region_size_mm.mat

General:
- get_fdata() everywhere
- no stroke/incidence thresholding
- no -a option
- writes BOTH voxel counts and mm^3
- savemat warnings avoided (drops __header__/__version__/__globals__)
- voxel size forced to user-confirmed truth: 0.068359 × 0.068359 × 0.5 mm
latest Update: Julian Cramer
"""

import os
import sys
import glob
import argparse
import numpy as np
import nibabel as nib
import scipy.io as sc


# ---- USER-CONFIRMED TRUE VOXEL SIZE (mm) ----
TRUE_VOX_MM = (0.068359, 0.068359, 0.5)
TRUE_VOX_VOL_MM3 = float(TRUE_VOX_MM[0] * TRUE_VOX_MM[1] * TRUE_VOX_MM[2])


def find_files(root: str, pattern: str):
    """Recursive glob under root."""
    return sorted(glob.glob(os.path.join(root, "**", pattern), recursive=True))


def build_id_to_name_map(lookup_txt_path: str):
    """Map atlas ID -> name using tab-separated or ITK-SNAP lookup file."""
    id_to_name = {}
    with open(lookup_txt_path, "r") as f:
        for ln in f:
            line = ln.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split(None, 7)
            if len(parts) >= 8:
                try:
                    rid = int(parts[0])
                except ValueError:
                    continue
                id_to_name[rid] = parts[7].strip().strip('"')
                continue

            parts = line.split("\t")
            if len(parts) >= 2 and parts[0].strip().lstrip("-").isdigit():
                try:
                    rid = int(parts[0])
                except ValueError:
                    continue
                id_to_name[rid] = parts[1].strip()
    return id_to_name


def load_label_template(labels_mat_path: str, id_to_name: dict):
    """Load existing label MAT if present, otherwise derive IDs from lookup text."""
    if labels_mat_path and os.path.exists(labels_mat_path):
        mat_in = sc.loadmat(labels_mat_path)
        if "SIGLabelIDs" in mat_in:
            return mat_in, np.array(mat_in["SIGLabelIDs"][:, 0]).astype(np.int64)
        if "ABALabelIDs" in mat_in:
            return mat_in, np.array(mat_in["ABALabelIDs"][:, 0]).astype(np.int64)
        raise RuntimeError(f"'{labels_mat_path}' does not contain key 'SIGLabelIDs' or 'ABALabelIDs'")

    atlas_label_ids = np.array(sorted(rid for rid in id_to_name.keys() if rid > 0), dtype=np.int64)
    return {}, atlas_label_ids


def derive_mask_path(anno_path: str) -> str:
    suffixes = [
        "_Anno_parental.nii.gz",
        "_AnnoSplit_parental.nii.gz",
        "_AnnorsfMRI.nii.gz",
        "_Anno.nii.gz",
        "_AnnoSplit.nii.gz",
    ]
    for sfx in suffixes:
        if anno_path.endswith(sfx):
            return anno_path.replace(sfx, "_mask.nii.gz")

    if anno_path.endswith(".nii.gz"):
        return anno_path[:-7] + "_mask.nii.gz"
    if anno_path.endswith(".nii"):
        return anno_path[:-4] + "_mask.nii.gz"
    return anno_path + "_mask.nii.gz"


def compute_region_sizes_from_annotation(
    anno_paths,
    labels_mat_path: str,
    lookup_txt_path: str,
    outfile_dir: str,
    output_stem: str,
    mask_replace_from: str,
    mask_replace_to: str,
    use_lookup_filter: bool = False,
):
    """
    Computes region voxel counts and mm^3 from annotation label images.
    Brain volume is computed from paired mask voxel count (or annotation foreground fallback),
    then multiplied by forced TRUE_VOX_VOL_MM3.

    Filtering of output IDs:
      - If use_lookup_filter=True: keep only IDs that exist in lookup table keys
        (recommended for parental atlas files).
      - Else: keep only IDs that exist in SIGLabelIDs loaded from labels_mat_path.
    """

    id_to_name = build_id_to_name_map(lookup_txt_path)
    mat_in, atlas_label_ids = load_label_template(labels_mat_path, id_to_name)

    region_vox = {}   # rid -> vox
    region_mm3 = {}   # rid -> mm^3
    total_brain_vox = 0
    total_brain_mm3 = 0.0
    n_missing_masks = 0

    for anno_path in anno_paths:
        anno_img = nib.load(anno_path)
        anno = np.round(anno_img.get_fdata()).astype(np.int64)

        # Paired mask path
        if mask_replace_from:
            mask_path = anno_path.replace(mask_replace_from, mask_replace_to)
        else:
            mask_path = derive_mask_path(anno_path)

        if os.path.exists(mask_path):
            mask_img = nib.load(mask_path)
            mask = mask_img.get_fdata()
            brain_vox = int(np.count_nonzero(mask > 0))
        else:
            # Fallback: annotation foreground
            n_missing_masks += 1
            brain_vox = int(np.count_nonzero(anno > 0))

        brain_mm3 = brain_vox * TRUE_VOX_VOL_MM3
        total_brain_vox += brain_vox
        total_brain_mm3 += brain_mm3

        # Regions in this annotation
        rids = np.unique(anno)
        rids = rids[rids > 0]

        for rid in rids:
            vox = int(np.count_nonzero(anno == rid))
            mm3 = vox * TRUE_VOX_VOL_MM3
            region_vox[rid] = region_vox.get(rid, 0) + vox
            region_mm3[rid] = region_mm3.get(rid, 0.0) + mm3

    # ---- Filter what we output ----
    present_ids = np.array(sorted(region_vox.keys()), dtype=np.int64)

    if use_lookup_filter:
        # Recommended for parental atlas: lookup table defines the valid output IDs.
        allowed_ids = set(id_to_name.keys())
        present_ids = np.array([rid for rid in present_ids if int(rid) in allowed_ids], dtype=np.int64)
    else:
        present_ids = present_ids[np.isin(present_ids, atlas_label_ids)]

    # ---- TXT ----
    out_txt = os.path.join(outfile_dir, f"{output_stem}.txt")
    with open(out_txt, "w") as o:
        o.write(
            "Brain volume (sum across processed files): "
            f"{total_brain_vox}\tvox\t{total_brain_mm3:0.2f}\tmm^3\n"
        )
        o.write(f"Voxel size forced to (mm): {TRUE_VOX_MM[0]} {TRUE_VOX_MM[1]} {TRUE_VOX_MM[2]}\n")
        o.write(f"Voxel volume forced to (mm^3): {TRUE_VOX_VOL_MM3:.10f}\n")
        if n_missing_masks > 0:
            o.write(f"Warning: {n_missing_masks} file(s) missing mask; used annotation foreground as fallback.\n")
        o.write("ID\tName\tVoxels\tUnit\tVolume\tUnit\n")

        names, vox_list, mm3_list = [], [], []
        for rid in present_ids:
            rid_int = int(rid)
            name = id_to_name.get(rid_int, "NA")
            vox = int(region_vox[rid_int])
            mm3 = float(region_mm3[rid_int])

            o.write(f"{rid_int}\t{name}\t{vox}\tvox\t{mm3:0.2f}\tmm^3\n")

            names.append(name)
            vox_list.append(vox)
            mm3_list.append(mm3)

    # ---- MAT (clean dict to avoid scipy warnings) ----
    out_mat = os.path.join(outfile_dir, f"{output_stem}.mat")
    clean_mat = {k: v for k, v in mat_in.items() if not k.startswith("_")}

    ids_out = present_ids.astype(np.float64)
    mm3_out = np.array(mm3_list, dtype=np.float64)
    vox_out = np.array(vox_list, dtype=np.float64)

    # Keep legacy-ish field name but store mm^3; also store vox explicitly
    clean_mat["SIGLabelIDs"] = np.stack((ids_out, mm3_out))   # ids + mm^3
    clean_mat["SIGLabelVox"] = np.stack((ids_out, vox_out))   # ids + vox
    clean_mat["SIGNames"] = np.array(names, dtype=object)
    clean_mat["SIGlabels"] = np.array(names, dtype=object)
    clean_mat["volumeMM"] = float(total_brain_mm3)
    clean_mat["volumeVox"] = float(total_brain_vox)
    clean_mat["voxelSizeMM"] = np.array(TRUE_VOX_MM, dtype=np.float64)
    clean_mat["voxelVolMM3"] = float(TRUE_VOX_VOL_MM3)
    clean_mat["usedLookupFilter"] = bool(use_lookup_filter)

    sc.savemat(out_mat, clean_mat)

    return out_txt, out_mat, len(present_ids)


def main():
    parser = argparse.ArgumentParser(description="Calculate atlas region sizes from registered annotation volumes.")
    parser.add_argument("-i", "--inputFolder", required=True, help="Path to .../T2w folder")
    args = parser.parse_args()

    input_folder = args.inputFolder
    if not os.path.isdir(input_folder):
        sys.exit(f"Error: '{input_folder}' is not an existing directory.")

    # Resources relative to repository layout
    base = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
    labels_mat = os.path.join(base, "lib", "SIGLabelsIDchanged.mat")

    # Parental atlas lookup
    lookup_txt_par = os.path.join(base, "lib", "sigma", "SIGMA_InVivo_Anatomical_Brain_Atlas_Labels.txt")

    lookup_txt_regular = lookup_txt_par

    for p in (lookup_txt_par,):
        if not os.path.exists(p):
            sys.exit(f"Error: Missing required file: '{p}'")

    # -------- 1) Parental annotation (anno_par) --------
    anno_par_paths = []
    for pat in ("*_Anno_parental.nii.gz", "*_AnnoSplit_parental.nii.gz", "*_AnnorsfMRI.nii.gz"):
        anno_par_paths.extend(find_files(input_folder, pat))
    anno_par_paths = sorted(set(anno_par_paths))
    print(f"[PARENTAL] '{len(anno_par_paths)}' anno_par file(s) will be processed...")

    if len(anno_par_paths) > 0:
        out_txt_par, out_mat_par, nreg_par = compute_region_sizes_from_annotation(
            anno_paths=anno_par_paths,
            labels_mat_path=labels_mat,
            lookup_txt_path=lookup_txt_par,
            outfile_dir=input_folder,
            output_stem="region_size_mm_par",
            mask_replace_from=None,
            mask_replace_to=None,
            use_lookup_filter=True,
        )
        print(f"[PARENTAL] Wrote: {out_txt_par} ({nreg_par} regions)")
        print(f"[PARENTAL] Wrote: {out_mat_par}")
    else:
        print("[PARENTAL] No *_AnnorsfMRI.nii.gz found. Skipping parental outputs.")

    # -------- 2) Regular annotation (anno) --------
    anno_paths = []
    for pat in ("*_Anno.nii.gz", "*_AnnoSplit.nii.gz"):
        anno_paths.extend(find_files(input_folder, pat))
    anno_paths = sorted(set(
        p for p in anno_paths
        if not (p.endswith("_Anno_parental.nii.gz") or p.endswith("_AnnoSplit_parental.nii.gz"))
    ))
    print(f"[REGULAR]  '{len(anno_paths)}' anno file(s) will be processed...")

    if len(anno_paths) > 0:
        out_txt, out_mat, nreg = compute_region_sizes_from_annotation(
            anno_paths=anno_paths,
            labels_mat_path=labels_mat,
            lookup_txt_path=lookup_txt_regular,
            outfile_dir=input_folder,
            output_stem="region_size_mm",
            mask_replace_from=None,
            mask_replace_to=None,
            use_lookup_filter=False,
        )
        print(f"[REGULAR]  Wrote: {out_txt} ({nreg} regions)")
        print(f"[REGULAR]  Wrote: {out_mat}")
    else:
        print("[REGULAR]  No *_Anno.nii.gz found. Skipping regular outputs.")


if __name__ == "__main__":
    main()
