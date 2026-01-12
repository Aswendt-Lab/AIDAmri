#!/usr/bin/env python3
import os
import sys
import subprocess
import shutil
import traceback
import argparse

import numpy as np
import nibabel as nib
from nibabel import orientations as nio
from typing import Optional


"""
Batch reorientation of NIfTI files within a BIDS-like directory structure.

Usage:
    python ReorientBatch.py <INPUT_ROOT> <OUTPUT_ROOT>

- All .nii and .nii.gz files under INPUT_ROOT are processed.
- Non-NIfTI files are copied unchanged into the mirrored structure.
- The relative directory structure is preserved under OUTPUT_ROOT.
- Images are reoriented to a user-specified target orientation
  (default: LIP for AIDAmri).
- Only the orientation (affine and axis order) is changed; all other
  header information is preserved where possible.
"""

def strip_nii_ext(p: str) -> str:
    if p.endswith(".nii.gz"):
        return p[:-7]
    if p.endswith(".nii"):
        return p[:-4]
    return p


def reorient_bvecs_fsl(src_bvec: str, dst_bvec: str, ornt_trans: np.ndarray):
    """
    Reorient FSL-style bvecs (3xN) using nibabel orientation transform.

    ornt_trans is expected to map from NEW (target) axes to OLD (current) axes,
    i.e. the inverse transform relative to the data reorientation.
    """
    bvec = np.loadtxt(src_bvec)

    # Allow both 3xN (FSL) and Nx3
    if bvec.ndim != 2:
        raise ValueError(f"Unexpected bvec ndim: {bvec.ndim}")

    if bvec.shape[0] != 3:
        if bvec.shape[1] == 3:
            bvec = bvec.T
        else:
            raise ValueError(f"Unexpected bvec shape: {bvec.shape}")

    new = np.zeros_like(bvec)

    for new_ax in range(3):
        old_ax = int(ornt_trans[new_ax, 0])
        flip = float(ornt_trans[new_ax, 1])
        new[new_ax, :] = flip * bvec[old_ax, :]

    np.savetxt(dst_bvec, new, fmt="%.16f")


def copy_sidecars_if_present(src_base: str, dst_base: str, *, reorient: bool, ornt_trans=None, log=None):
    """
    Copy bval and (optionally reorient) bvec sidecars from src_base to dst_base.
    """
    src_bval = src_base + ".bval"
    src_bvec = src_base + ".bvec"
    dst_bval = dst_base + ".bval"
    dst_bvec = dst_base + ".bvec"

    has_bval = os.path.exists(src_bval)
    has_bvec = os.path.exists(src_bvec)

    if not (has_bval or has_bvec):
        return

    os.makedirs(os.path.dirname(dst_base), exist_ok=True)

    if has_bval:
        shutil.copy2(src_bval, dst_bval)
        if log:
            log("  Sidecar: copied .bval")

    if has_bvec:
        if reorient:
            if ornt_trans is None:
                raise ValueError("ornt_trans is required to reorient bvecs")
            reorient_bvecs_fsl(src_bvec, dst_bvec, ornt_trans)
            if log:
                log("  Sidecar: reoriented .bvec")
        else:
            shutil.copy2(src_bvec, dst_bvec)
            if log:
                log("  Sidecar: copied .bvec (no reorientation)")

def ask_target_orientation_with_default(non_interactive: bool, target_cli: Optional[str]) -> str:
    if non_interactive:
        if not target_cli:
            raise ValueError("-n was set but -t is missing.")
        return target_cli.upper()

    # if target_cli set, dont ask
    if target_cli:
        return target_cli.upper()

    while True:
        ans = input(
            "\nDo you want to reorient all images to the AIDAmri default orientation LIP "
            "(Left-Inferior-Posterior)? [Y/n]: "
        ).strip().lower()

        if ans in ("", "y", "yes"):
            return "LIP"
        elif ans in ("n", "no"):
            break
        else:
            print("Please answer 'y' (yes) or 'n' (no).")

    while True:
        ori = input(
            "Please enter the target orientation "
            "(three letters from {L, R, A, P, S, I}): "
        ).strip().upper()

        if len(ori) != 3 or not set(ori).issubset(set("LRAPSI")):
            print("Invalid orientation.")
            continue

        return ori



def get_current_orientation(img: nib.Nifti1Image):
    """
    Determine current orientation from the 'active' transform:
    prefer sform if sform_code > 0, else qform if qform_code > 0,
    else fall back to img.affine.

    Returns (ori_str, src_label, affine_used)
    """
    hdr = img.header
    s_code = int(hdr["sform_code"])
    q_code = int(hdr["qform_code"])

    if s_code > 0:
        A = img.get_sform()
        src = f"sform (code={s_code})"
    elif q_code > 0:
        A = img.get_qform()
        src = f"qform (code={q_code})"
    else:
        A = img.affine
        src = "img.affine (fallback)"

    ori = "".join(nio.aff2axcodes(A))
    return ori, src, A


def reorient_image(img: nib.Nifti1Image, target_ori: str, current_ori: str, base_affine: np.ndarray, log=None):
    """
    Returns: (img_new, did_reorient: bool, ornt_trans_for_bvecs_or_None)
    """
    data = img.get_fdata(dtype=np.float32)
    data = np.ascontiguousarray(data, dtype=np.float32)

    if log:
        log(f"  Target orientation: {target_ori}")

    if current_ori == target_ori:
        if log:
            log("  Current orientation already matches target. No reorientation is applied.")

        hdr_copy = img.header.copy()
        img_copy = nib.Nifti1Image(data, base_affine, header=hdr_copy)

        s_code = int(img.header["sform_code"])
        q_code = int(img.header["qform_code"])
        img_copy.set_sform(base_affine, code=(s_code if s_code > 0 else 1))
        img_copy.set_qform(base_affine, code=(q_code if q_code > 0 else 1))

        return img_copy, False, None

    curr_ornt = nio.axcodes2ornt(tuple(current_ori))
    target_ornt = nio.axcodes2ornt(tuple(target_ori))

    # data transform: current -> target
    ornt_trans = nio.ornt_transform(curr_ornt, target_ornt)

    # bvec transform: target -> current (inverse), as in your single-file script
    ornt_trans_bvec = nio.ornt_transform(target_ornt, curr_ornt)

    if log:
        log("  Applying orientation transform to image data...")

    data_reoriented = nio.apply_orientation(data, ornt_trans)

    if log:
        log("  Updating affine for new orientation...")

    inv_aff = nio.inv_ornt_aff(ornt_trans, img.shape[:3])
    new_affine = base_affine @ inv_aff

    hdr = img.header.copy()
    hdr.set_data_dtype(np.float32)

    # robustness after permutation/flip:
    hdr.set_data_shape(data_reoriented.shape)

    # optional but recommended: reset freq/phase/slice info (can become wrong after permutation)
    try:
        hdr["dim_info"] = 0
    except Exception:
        pass

    img_new = nib.Nifti1Image(data_reoriented, new_affine, header=hdr)
    img_new.set_sform(new_affine, code=1)
    img_new.set_qform(new_affine, code=1)

    if log:
        new_str = "".join(nio.aff2axcodes(new_affine))
        log(f"  New orientation (from affine): {new_str}")

    return img_new, True, ornt_trans_bvec


def reorient_single_image(src_path: str, dst_path: str, target_ori: str, log):
    log("")
    log("Processing file:")
    log(f"  Source:      {src_path}")
    log(f"  Destination: {dst_path}")

    img = nib.load(src_path)

    current_ori, ori_src, A_used = get_current_orientation(img)
    log(f"  Current orientation (from {ori_src}): {current_ori}")
    log(f"  Target orientation:                {target_ori}")

    img_out, did_reorient, bvec_ornt = reorient_image(
        img, target_ori, current_ori, A_used, log=log
    )

    nib.save(img_out, dst_path)
    log("  Saved NIfTI.")

    # Sidecars (.bval/.bvec) handled here
    src_base = strip_nii_ext(src_path)
    dst_base = strip_nii_ext(dst_path)

    copy_sidecars_if_present(
        src_base,
        dst_base,
        reorient=did_reorient,
        ornt_trans=bvec_ornt,
        log=log,
    )


def copy_non_nifti(src_path: str, dst_path: str, log):
    log("")
    log("Copying non-NIfTI file:")
    log(f"  Source:      {src_path}")
    log(f"  Destination: {dst_path}")
    shutil.copy2(src_path, dst_path)



def validate_target_ori(ori: str) -> str:
    ori = ori.strip().upper()
    if len(ori) != 3 or not set(ori).issubset(set("LRAPSI")):
        raise ValueError(f"Invalid orientation: {ori}")
    axes = set()
    for c in ori:
        if c in "LR": axes.add("x")
        elif c in "AP": axes.add("y")
        elif c in "SI": axes.add("z")
    if axes != {"x", "y", "z"}:
        raise ValueError(f"Invalid axis combination in orientation: {ori}")
    return ori

def parse_args():
    parser = argparse.ArgumentParser(
        description="Batch reorientation of NIfTI files (AIDAmri compatible)."
    )

    parser.add_argument(
        "-i",
        required=True,
        metavar="INPUT_ROOT",
        help="Input root directory (BIDS-like proc_data)"
    )

    parser.add_argument(
        "-o",
        required=True,
        metavar="OUTPUT_ROOT",
        help="Output root directory for reoriented data"
    )

    parser.add_argument(
        "-t",
        metavar="ORI",
        default=None,
        help="Target orientation (e.g. LIP). If omitted, ask interactively."
    )

    parser.add_argument(
        "-l",
        default="reorient_log.txt",
        metavar="LOGFILE",
        help="Log filename (written into output root)"
    )

    parser.add_argument(
        "-n",
        action="store_true",
        help="Non-interactive mode (requires -t)"
    )

    return parser.parse_args()

def main():
    args = parse_args()

    src_root = args.i
    dst_root = args.o
    LOG_FILENAME = args.l
    target_cli = args.t
    non_interactive = args.n

    if not os.path.isdir(src_root):
        print(f"Source root folder not found:\n{src_root}")
        sys.exit(1)

    os.makedirs(dst_root, exist_ok=True)
    log_path = os.path.join(dst_root, LOG_FILENAME)

    target_ori = ask_target_orientation_with_default(
        non_interactive=non_interactive,
        target_cli=target_cli
    )
    target_ori = validate_target_ori(target_ori)

    # --- Count total files for progress bar ---
    total_files = 0
    for root, dirs, files in os.walk(src_root):
        total_files += len(files)

    if total_files == 0:
        print("No files found in source root. Nothing to do.")
        return

    any_errors = False
    n_total_nifti = 0
    n_processed_nifti = 0
    n_copied_non_nifti = 0
    processed_files = 0
    bar_width = 40

    with open(log_path, "w") as log_fh:

        def log(msg: str):
            log_fh.write(msg + "\n")

        log(f"Target orientation for all images: {target_ori}")
        log(f"Source root: {src_root}")
        log(f"Destination root: {dst_root}")
        log(f"Total files (NIfTI + non-NIfTI): {total_files}")

        for root, dirs, files in os.walk(src_root):
            for fname in files:
                src_path = os.path.join(root, fname)
                rel_path = os.path.relpath(src_path, src_root)
                dst_path = os.path.join(dst_root, rel_path)

                os.makedirs(os.path.dirname(dst_path), exist_ok=True)

                is_nifti = fname.endswith(".nii") or fname.endswith(".nii.gz")
                is_sidecar = fname.endswith(".bvec") or fname.endswith(".bval")

                try:
                    if is_nifti:
                        n_total_nifti += 1
                        reorient_single_image(src_path, dst_path, target_ori, log)
                        n_processed_nifti += 1
                    elif is_sidecar:
                        # skip: handled with the image
                        pass
                    else:
                        copy_non_nifti(src_path, dst_path, log)
                        n_copied_non_nifti += 1
                except Exception as e:
                    any_errors = True
                    print(f"\nError processing file: {rel_path}: {e.__class__.__name__}: {e}", file=sys.stderr)
                    log("")
                    log("ERROR during processing file:")
                    log(f"  Source: {src_path}")
                    log(f"  Destination: {dst_path}")
                    log(f"  Exception: {e.__class__.__name__}: {e}")
                    log("  Traceback:")
                    log(traceback.format_exc())
                finally:
                    processed_files += 1
                    progress = processed_files / total_files
                    filled = int(bar_width * progress)
                    bar = "#" * filled + "-" * (bar_width - filled)
                    print(
                        f"\rProgress: [{bar}] {progress * 100:6.2f}% ({processed_files}/{total_files})",
                        end="",
                        flush=True,
                    )

    print()
    print("\nBatch processing completed.")
    print(f"Total NIfTI files found:     {n_total_nifti}")
    print(f"Total NIfTI files processed: {n_processed_nifti}")
    print(f"Non-NIfTI files copied:      {n_copied_non_nifti}")
    print(f"Reoriented data written to:  {dst_root}")
    print(f"Log file written to:         {log_path}")

    if any_errors:
        print("One or more errors occurred. See the log file for details.")

if __name__ == "__main__":
    main()
