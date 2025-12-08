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

# Log file name in destination root
LOG_FILENAME = "reorient_log.txt"


def ask_target_orientation_with_default() -> str:
    """
    Ask the user whether to reorient to the AIDAmri orientation (LIP).
    If not, ask for a custom 3-letter target orientation.

    Returns
    -------
    ori : str
        Valid orientation string of length 3 consisting of letters
        from {L, R, A, P, S, I}, e.g. 'LIP', 'RAS', 'LPI'.
    """
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

    # Custom target orientation
    while True:
        ori = input(
            "Please enter the target orientation for all images "
            "(three letters from {L, R, A, P, S, I}, e.g. 'RAS', 'LPI', 'LIP'): "
        ).strip().upper()

        if len(ori) != 3:
            print("Target orientation must consist of exactly 3 letters.")
            continue

        if not set(ori).issubset(set("LRAPSI")):
            print("Only letters from {L, R, A, P, S, I} are allowed.")
            continue

        # Check that each spatial axis (x, y, z) is represented exactly once
        axes = set()
        for c in ori:
            if c in "LR":
                axes.add("x")
            elif c in "AP":
                axes.add("y")
            elif c in "SI":
                axes.add("z")

        if axes != {"x", "y", "z"}:
            print(
                "Invalid combination: each of the axis pairs must appear exactly once:\n"
                "  - L or R for the x-axis\n"
                "  - A or P for the y-axis\n"
                "  - S or I for the z-axis"
            )
            continue

        return ori


def get_orientation_from_fsl(path: str, log) -> str:
    """
    Retrieve the image orientation using FSL (fslorient -getsform)
    and convert it into orientation codes (L/R, A/P, S/I).

    If fslorient does not return a valid 4x4 sform matrix, fall back to
    nibabel's affine.
    """
    cmd = ["fslorient", "-getsform", path]
    res = subprocess.run(cmd, capture_output=True, text=True)

    txt = res.stdout.strip()

    # If FSL failed or returned nothing, fall back to nibabel
    if res.returncode != 0 or not txt:
        log("  Warning: Could not retrieve SForm using fslorient. Falling back to nibabel affine.")
        img = nib.load(path)
        return "".join(nio.aff2axcodes(img.affine))

    # Parse lines containing numeric values
    lines = [ln.strip() for ln in txt.splitlines() if ln.strip()]
    rows = []
    for ln in lines:
        parts = ln.split()
        nums = []
        for p in parts:
            try:
                nums.append(float(p))
            except ValueError:
                # Ignore non-numeric tokens, if any
                continue
        if len(nums) >= 4:
            rows.append(nums[:4])

    # If we did not get any numeric rows, fall back
    if len(rows) == 0:
        log("  Warning: fslorient output could not be parsed as a matrix. Falling back to nibabel affine.")
        img = nib.load(path)
        return "".join(nio.aff2axcodes(img.affine))

    sform = np.array(rows, dtype=float)

    # Try to reshape to 4x4; if that fails, fall back
    if sform.size != 16:
        log("  Warning: fslorient did not return a 4x4 sform matrix. Falling back to nibabel affine.")
        img = nib.load(path)
        return "".join(nio.aff2axcodes(img.affine))

    sform = sform.reshape(4, 4)

    # Convert matrix -> axis codes
    axcodes = nio.aff2axcodes(sform)
    return "".join(axcodes)


def reorient_single_image(src_path: str, dst_path: str, target_ori: str, log):
    """
    Load a single NIfTI image, reorient it to the target orientation,
    and save it to dst_path, preserving header information where possible.
    All messages go to the log function.
    """
    log("")
    log("Processing file:")
    log(f"  Source:      {src_path}")
    log(f"  Destination: {dst_path}")

    img = nib.load(src_path)
    data = img.get_fdata(dtype=np.float32)
    data = np.ascontiguousarray(data, dtype=np.float32)

    # Determine current orientation
    curr_str = get_orientation_from_fsl(src_path, log)
    log(f"  Current orientation: {curr_str}")
    log(f"  Target orientation:  {target_ori}")

    # If the current orientation is already the target, just copy to dst
    if curr_str == target_ori:
        log("  Current orientation already matches target. Copying image with unchanged affine.")
        hdr_copy = img.header.copy()
        img_copy = nib.Nifti1Image(data, img.affine, header=hdr_copy)
        img_copy.set_sform(img.affine, code=img.get_sform(coded=True)[1])
        img_copy.set_qform(img.affine, code=img.get_qform(coded=True)[1])
        nib.save(img_copy, dst_path)
        return

    # Build orientation arrays
    curr_axcodes = tuple(curr_str)
    curr_ornt = nio.axcodes2ornt(curr_axcodes)
    target_axcodes = tuple(target_ori)
    target_ornt = nio.axcodes2ornt(target_axcodes)

    # Transformation from current orientation to target orientation
    ornt_trans = nio.ornt_transform(curr_ornt, target_ornt)

    # Reorient data
    log("  Applying orientation transform to image data...")
    data_reoriented = nio.apply_orientation(data, ornt_trans)

    # Compute new affine consistent with the reoriented data
    log("  Updating affine for new orientation...")
    inv_aff = nio.inv_ornt_aff(ornt_trans, img.shape[:3])
    new_affine = img.affine @ inv_aff

    # Build new image, preserving header except for affine/sform/qform
    hdr = img.header.copy()
    hdr.set_data_dtype(np.float32)

    img_new = nib.Nifti1Image(data_reoriented, new_affine, header=hdr)
    img_new.set_sform(new_affine, code=1)
    img_new.set_qform(new_affine, code=1)

    # Report new orientation
    new_axcodes = nio.aff2axcodes(new_affine)
    new_str = "".join(new_axcodes)
    log(f"  New orientation (from affine): {new_str}")

    nib.save(img_new, dst_path)


def copy_non_nifti(src_path: str, dst_path: str, log):
    """
    Copy a non-NIfTI file unchanged to the destination path.
    """
    log("")
    log("Copying non-NIfTI file:")
    log(f"  Source:      {src_path}")
    log(f"  Destination: {dst_path}")
    shutil.copy2(src_path, dst_path)


def batch_reorient(src_root: str, dst_root: str):
    """
    Batch-process all files under src_root and mirror them under dst_root.
    """
    # --- Sanity checks for root folders ---
    if not os.path.isdir(src_root):
        print(f"Source root folder not found:\n{src_root}")
        sys.exit(1)

    os.makedirs(dst_root, exist_ok=True)
    log_path = os.path.join(dst_root, LOG_FILENAME)

    # --- Ask once for target orientation (default option: LIP) ---
    target_ori = ask_target_orientation_with_default()

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

        # --- Walk through src_root and process all files ---
        for root, dirs, files in os.walk(src_root):
            for fname in files:
                src_path = os.path.join(root, fname)
                rel_path = os.path.relpath(src_path, src_root)
                dst_path = os.path.join(dst_root, rel_path)

                # Ensure target directory exists
                os.makedirs(os.path.dirname(dst_path), exist_ok=True)

                # NIfTI or not?
                is_nifti = fname.endswith(".nii") or fname.endswith(".nii.gz")

                try:
                    if is_nifti:
                        n_total_nifti += 1
                        reorient_single_image(src_path, dst_path, target_ori, log)
                        n_processed_nifti += 1
                    else:
                        copy_non_nifti(src_path, dst_path, log)
                        n_copied_non_nifti += 1
                except Exception as e:
                    any_errors = True
                    # Short message for terminal:
                    print(f"\nError processing file: {rel_path}: {e.__class__.__name__}: {e}", file=sys.stderr)
                    # Detailed message for logfile:
                    log("")
                    log("ERROR during processing file:")
                    log(f"  Source: {src_path}")
                    log(f"  Destination: {dst_path}")
                    log(f"  Exception: {e.__class__.__name__}: {e}")
                    log("  Traceback:")
                    log(traceback.format_exc())
                finally:
                    # Update progress bar
                    processed_files += 1
                    progress = processed_files / total_files
                    filled = int(bar_width * progress)
                    bar = "#" * filled + "-" * (bar_width - filled)
                    print(
                        f"\rProgress: [{bar}] {progress * 100:6.2f}% "
                        f"({processed_files}/{total_files})",
                        end="",
                        flush=True,
                    )

    # Zeile nach dem Fortschrittsbalken beenden
    print()

    # --- Final summary to terminal only ---
    print("\nBatch processing completed.")
    print(f"Total NIfTI files found:     {n_total_nifti}")
    print(f"Total NIfTI files processed: {n_processed_nifti}")
    print(f"Non-NIfTI files copied:      {n_copied_non_nifti}")
    print(f"Reoriented data written to:  {dst_root}")
    print(f"Log file written to:         {log_path}")

    if any_errors:
        print("One or more errors occurred. See the log file for details.")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Batch reorient NIfTI files in a BIDS-like directory tree."
    )
    parser.add_argument(
        "-i", "--input_root",
        help="Input root directory (e.g. proc_data).",
    )
    parser.add_argument(
        "-o", "--output_root",
        help="Output root directory for reoriented data.",
    )

    # also allow positional arguments for backward compatibility
    parser.add_argument(
        "positional_input",
        nargs="?",
        help="Positional input directory (optional)."
    )
    parser.add_argument(
        "positional_output",
        nargs="?",
        help="Positional output directory (optional)."
    )

    args = parser.parse_args()

    input_root = args.input_root or args.positional_input
    output_root = args.output_root or args.positional_output

    if not input_root or not output_root:
        parser.error("Must specify input and output directories via -i/-o or as positional arguments.")

    return input_root, output_root



if __name__ == "__main__":
    input_root, output_root = parse_args()
    batch_reorient(input_root, output_root)

