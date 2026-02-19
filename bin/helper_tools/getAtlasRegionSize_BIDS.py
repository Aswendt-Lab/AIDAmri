"""
Created on 10/08/2017
Modified on 04/23/2021
Updated on 02/19/2026

@author: Niklas Pallast and Markus Aswendt
Neuroimaging & Neuroengineering
Department of Neurology
University Hospital Cologne

#!/usr/bin/env python3

getAtlasRegionSize_BIDS.py

Process ALL annotation versions found recursively under the input folder and
write per-file outputs into: <inputFolder>/output_region_size/

Annotation variants searched (recursive):
  Parental:
    - *_Anno_parental.nii.gz
    - *_AnnoSplit_parental.nii.gz
    - *_AnnorsfMRI.nii.gz            (legacy parental)
  Regular:
    - *_Anno.nii.gz                  (excluding *_Anno_parental.nii.gz)
    - *_AnnoSplit.nii.gz             (excluding *_AnnoSplit_parental.nii.gz)

Mask pairing:
  - Derived by replacing the matched annotation suffix with "_mask.nii.gz"

Outputs:
  - One .txt and one .mat per annotation file.
  - Parental outputs have suffix "_par" before extension to distinguish:
      <basename>_par.txt / <basename>_par.mat
    Regular outputs:
      <basename>.txt / <basename>.mat

Volumes:
  - Forced voxel size (mm): 0.068359 × 0.068359 × 0.5
  - Forced voxel volume (mm^3): 0.0023364764
  - Brain volume: mask voxel count (>0) * forced voxel volume
    (fallback: annotation foreground >0 if mask missing)

Parental hemisphere IDs:
  - For parental outputs, filtering uses lookup-table keys (keeps +2000 R hemisphere IDs).
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
    """Map atlas ID -> name using tab-separated lookup file."""
    id_to_name = {}
    with open(lookup_txt_path, "r") as f:
        for ln in f:
            parts = ln.rstrip("\n").split("\t")
            if len(parts) >= 2:
                try:
                    rid = int(parts[0])
                except ValueError:
                    continue
                id_to_name[rid] = parts[1]
    return id_to_name


def strip_nii_ext(filename: str) -> str:
    """Remove .nii.gz or .nii from filename."""
    if filename.endswith(".nii.gz"):
        return filename[:-7]
    if filename.endswith(".nii"):
        return filename[:-4]
    return filename


def derive_mask_path(anno_path: str) -> str:
    """Derive mask path from known annotation suffixes."""
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

    # fallback heuristic
    if anno_path.endswith(".nii.gz"):
        return anno_path[:-7] + "_mask.nii.gz"
    if anno_path.endswith(".nii"):
        return anno_path[:-4] + "_mask.nii.gz"
    return anno_path + "_mask.nii.gz"


def compute_one_annotation(
    anno_path: str,
    labels_mat_path: str,
    lookup_txt_path: str,
    out_dir: str,
    is_parental: bool,
    use_lookup_filter: bool,
):
    """
    Compute region voxel counts and mm^3 for a single annotation file.
    Writes one .txt and one .mat in out_dir.
    """
    # Load lookup and labels
    id_to_name = build_id_to_name_map(lookup_txt_path)

    mat_in = sc.loadmat(labels_mat_path)
    if "ABALabelIDs" not in mat_in:
        raise RuntimeError(f"'{labels_mat_path}' does not contain key 'ABALabelIDs'")
    atlas_label_ids = np.array(mat_in["ABALabelIDs"][:, 0]).astype(np.int64)

    # Load annotation
    anno_img = nib.load(anno_path)
    anno = np.round(anno_img.get_fdata()).astype(np.int64)

    # Brain voxel count from mask (preferred)
    mask_path = derive_mask_path(anno_path)
    used_mask = False
    if os.path.exists(mask_path):
        mask = nib.load(mask_path).get_fdata()
        brain_vox = int(np.count_nonzero(mask > 0))
        used_mask = True
    else:
        brain_vox = int(np.count_nonzero(anno > 0))

    brain_mm3 = brain_vox * TRUE_VOX_VOL_MM3

    # Region voxel counts
    region_vox = {}
    rids = np.unique(anno)
    rids = rids[rids > 0]
    for rid in rids:
        rid_i = int(rid)
        region_vox[rid_i] = int(np.count_nonzero(anno == rid_i))

    # Filter output IDs
    present_ids = np.array(sorted(region_vox.keys()), dtype=np.int64)
    if use_lookup_filter:
        allowed_ids = set(id_to_name.keys())
        present_ids = np.array([rid for rid in present_ids if int(rid) in allowed_ids], dtype=np.int64)
    else:
        present_ids = present_ids[np.isin(present_ids, atlas_label_ids)]

    # Output names (per file)
    base = strip_nii_ext(os.path.basename(anno_path))
    out_stem = f"{base}_par" if is_parental else base
    out_txt = os.path.join(out_dir, f"{out_stem}.txt")
    out_mat = os.path.join(out_dir, f"{out_stem}.mat")

    # Write TXT
    with open(out_txt, "w") as o:
        o.write(f"Brain volume (single file): {brain_vox}\tvox\t{brain_mm3:0.2f}\tmm^3\n")
        o.write(f"Voxel size forced to (mm): {TRUE_VOX_MM[0]} {TRUE_VOX_MM[1]} {TRUE_VOX_MM[2]}\n")
        o.write(f"Voxel volume forced to (mm^3): {TRUE_VOX_VOL_MM3:.10f}\n")
        o.write(f"Annotation: {anno_path}\n")
        o.write(f"Mask: {mask_path} ({'found' if used_mask else 'missing -> used annotation foreground'})\n")
        o.write("ID\tName\tVoxels\tUnit\tVolume\tUnit\n")

        names, vox_list, mm3_list = [], [], []
        for rid in present_ids:
            rid_int = int(rid)
            name = id_to_name.get(rid_int, "NA")
            vox = int(region_vox[rid_int])
            mm3 = vox * TRUE_VOX_VOL_MM3
            o.write(f"{rid_int}\t{name}\t{vox}\tvox\t{mm3:0.2f}\tmm^3\n")
            names.append(name)
            vox_list.append(vox)
            mm3_list.append(mm3)

    # Write MAT (clean dict)
    clean_mat = {k: v for k, v in mat_in.items() if not k.startswith("_")}

    ids_out = present_ids.astype(np.float64)
    mm3_out = np.array(mm3_list, dtype=np.float64)
    vox_out = np.array(vox_list, dtype=np.float64)

    clean_mat["ABALabelIDs"] = np.stack((ids_out, mm3_out))  # ids + mm^3
    clean_mat["ABALabelVox"] = np.stack((ids_out, vox_out))  # ids + vox
    clean_mat["ABANames"] = np.array(names, dtype=object)
    clean_mat["ABAlabels"] = np.array(names, dtype=object)
    clean_mat["volumeMM"] = float(brain_mm3)
    clean_mat["volumeVox"] = float(brain_vox)
    clean_mat["voxelSizeMM"] = np.array(TRUE_VOX_MM, dtype=np.float64)
    clean_mat["voxelVolMM3"] = float(TRUE_VOX_VOL_MM3)
    clean_mat["usedLookupFilter"] = bool(use_lookup_filter)
    clean_mat["annotationPath"] = np.array([anno_path], dtype=object)
    clean_mat["maskPath"] = np.array([mask_path], dtype=object)

    sc.savemat(out_mat, clean_mat)

    return out_txt, out_mat, int(len(present_ids))


def main():
    parser = argparse.ArgumentParser(
        description="Compute atlas region sizes for ALL annotation versions in a BIDS folder."
    )
    parser.add_argument(
        "-i", "--inputFolder", required=True,
        help="Input folder to search recursively for annotation files."
    )
    args = parser.parse_args()

    input_folder = args.inputFolder
    if not os.path.isdir(input_folder):
        sys.exit(f"Error: '{input_folder}' is not an existing directory.")

    # Output folder inside input folder
    out_dir = os.path.join(input_folder, "output_region_size")
    os.makedirs(out_dir, exist_ok=True)

    # Resources relative to current working directory (pipeline layout)
    base = os.path.abspath(os.path.join(os.getcwd(), os.pardir, os.pardir))
    labels_mat = os.path.join(base, "lib", "ABALabelsIDchanged.mat")
    lookup_txt_par = os.path.join(base, "lib", "annoVolume+2000_rsfMRI.nii.txt")

    lookup_txt_regular = os.path.join(base, "lib", "ARA_annotationR+2000.nii.txt")
    if not os.path.exists(lookup_txt_regular):
        lookup_txt_regular = lookup_txt_par

    for p in (labels_mat, lookup_txt_par):
        if not os.path.exists(p):
            sys.exit(f"Error: Missing required file: '{p}'")

    # Collect all variants
    parental_patterns = [
        "*_Anno_parental.nii.gz",
        "*_AnnoSplit_parental.nii.gz",
        "*_AnnorsfMRI.nii.gz",  # legacy parental
    ]
    regular_patterns = [
        "*_Anno.nii.gz",
        "*_AnnoSplit.nii.gz",
    ]

    parental_files = []
    for pat in parental_patterns:
        parental_files.extend(find_files(input_folder, pat))
    parental_files = sorted(set(parental_files))

    regular_files = []
    for pat in regular_patterns:
        regular_files.extend(find_files(input_folder, pat))
    regular_files = sorted(set(regular_files))

    # Remove parental files from regular list
    regular_files = [
        p for p in regular_files
        if not (p.endswith("_Anno_parental.nii.gz") or p.endswith("_AnnoSplit_parental.nii.gz"))
    ]

    print(f"[INFO] Output folder: {out_dir}")
    print(f"[PARENTAL] Found {len(parental_files)} file(s).")
    print(f"[REGULAR]  Found {len(regular_files)} file(s).")

    # Process parental files
    for ap in parental_files:
        out_txt, out_mat, nreg = compute_one_annotation(
            anno_path=ap,
            labels_mat_path=labels_mat,
            lookup_txt_path=lookup_txt_par,
            out_dir=out_dir,
            is_parental=True,
            use_lookup_filter=True,   # keep +2000 hemisphere IDs
        )
        print(f"[PARENTAL] {os.path.basename(ap)} -> {os.path.basename(out_txt)} ({nreg} regions)")

    # Process regular files
    for ap in regular_files:
        out_txt, out_mat, nreg = compute_one_annotation(
            anno_path=ap,
            labels_mat_path=labels_mat,
            lookup_txt_path=lookup_txt_regular,
            out_dir=out_dir,
            is_parental=False,
            use_lookup_filter=False,
        )
        print(f"[REGULAR]  {os.path.basename(ap)} -> {os.path.basename(out_txt)} ({nreg} regions)")


if __name__ == "__main__":
    main()
