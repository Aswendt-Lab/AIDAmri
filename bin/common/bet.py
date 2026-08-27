import os
import subprocess

import nibabel as nib
import nipype.interfaces.fsl as fsl
import numpy as np





def set_default_xyzt_units_if_unknown(target):
    """
       Set NIfTI space/time units to mm/sec only if they are unknown.

       Accepts either:
       - a file path to a NIfTI file, or
       - a nibabel header object
       """
    if isinstance(target, str):
        img = nib.load(target)
        hdr = img.header

        space_unit, time_unit = hdr.get_xyzt_units()
        if not space_unit or space_unit.lower() == "unknown":
            space_unit = "mm"
        if not time_unit or time_unit.lower() == "unknown":
            time_unit = "sec"

        hdr.set_xyzt_units(space_unit, time_unit)
        nib.save(img, target)
        return target

    if hasattr(target, "get_xyzt_units") and hasattr(target, "set_xyzt_units"):
        space_unit, time_unit = target.get_xyzt_units()
        if not space_unit or space_unit.lower() == "unknown":
            space_unit = "mm"
        if not time_unit or time_unit.lower() == "unknown":
            time_unit = "sec"

        target.set_xyzt_units(space_unit, time_unit)
        return target

    raise TypeError("Expected a NIfTI file path or nibabel header object")



def estimate_center_intensity_based(nifti, percentile=60):
    """
    Estimate BET center (-c) using intensity-weighted center-of-gravity (fslstats -C),
    but excluding low-intensity voxels using a data-adaptive threshold (-l = P{percentile}).
    """
    # 1) Get intensity percentile
    p = subprocess.check_output(
        ["fslstats", nifti, "-P", str(percentile)]
    ).decode().strip()

    # 2) Compute center-of-gravity using only voxels > P{percentile}
    center = subprocess.check_output(
        ["fslstats", nifti, "-l", p, "-C"]
    ).decode().strip().split()

    cx, cy, cz = [float(v) for v in center]
    return [cx, cy, cz], float(p)


def transform_voxel_center_between_images(center, src_img_or_path, dst_img_or_path):
    """
    Map a voxel-space BET center from one NIfTI geometry to another while
    preserving the same physical point.
    """
    src_img = nib.load(src_img_or_path) if isinstance(src_img_or_path, str) else src_img_or_path
    dst_img = nib.load(dst_img_or_path) if isinstance(dst_img_or_path, str) else dst_img_or_path

    center = np.asarray(center, dtype=float)
    if center.shape != (3,):
        raise ValueError("BET center must contain exactly three coordinates: x y z")

    src_voxel = np.append(center, 1.0) # convert voxel coordinate to homogeneous coordinates
    world = src_img.affine @ src_voxel # convert voxel coordinates into world coordinates
    dst_voxel = np.linalg.inv(dst_img.affine) @ world # convert world coordinates into voxel coordinates of the destination image

    return dst_voxel[:3].tolist()


def skip_bet_function(input_file, return_mask=False):
    """
    Create BET-compatible outputs when BET is skipped.
    Reproduces key geometry/orientation steps so that downstream
    AIDAmri pipeline steps remain compatible.
    """

    print("Skipping BET")
    print("Creating BET-compatible outputs for pipeline compatibility")

    outputBET = os.path.join(
        os.path.dirname(input_file),
        os.path.basename(input_file).split('.')[0] + 'Bet.nii.gz'
    )

    src = nib.load(input_file)
    data = src.get_fdata(dtype=np.float32)
    aff = src.affine.copy()

    hdr = src.header.copy()
    hdr.set_data_dtype(np.float32)
    hdr["pixdim"][0] = 1
    hdr["pixdim"][4:8] = 1
    set_default_xyzt_units_if_unknown(hdr)

    final_img = nib.Nifti1Image(
        np.ascontiguousarray(data, dtype=np.float32),
        aff,
        header=hdr
    )

    final_img.set_qform(aff, code=1)
    final_img.set_sform(aff, code=1)

    nib.save(final_img, outputBET)

    print(f"BET skipped -> created compatibility image: {outputBET}")

    bet_mask_path = outputBET.replace(".nii.gz", "_mask.nii.gz")
    mask = (data > 0).astype(np.uint8)

    mask_hdr = hdr.copy()
    mask_hdr.set_data_dtype(np.uint8)
    mask_hdr["pixdim"][0] = 1
    mask_hdr["pixdim"][4:8] = 1
    set_default_xyzt_units_if_unknown(mask_hdr)

    mask_img = nib.Nifti1Image(
        np.ascontiguousarray(mask, dtype=np.uint8),
        aff,
        header=mask_hdr
    )

    mask_img.set_qform(aff, code=1)
    mask_img.set_sform(aff, code=1)

    nib.save(mask_img, bet_mask_path)

    print(f"BET mask created at {bet_mask_path}")
    return _format_result(outputBET, bet_mask_path, return_mask)


FSL_BET_WORLD_SWAPS = [(1, 2)]


def apply_world_ops(mat, swaps=()):
    out = mat.copy()

    # swaps first
    for a, b in swaps:
        out[[a, b], :] = out[[b, a], :]

    return out


def save_header_only_reoriented_copy(src_path, dst_path, swaps=()):
    img = nib.load(src_path)
    data = img.get_fdata(dtype=np.float32)
    aff = img.affine.copy()

    # header-only world-axis operation
    aff[:3, :] = apply_world_ops(aff[:3, :], swaps=swaps)

    hdr = img.header.copy()
    hdr["pixdim"][0] = 1
    hdr.set_data_dtype(np.float32)
    set_default_xyzt_units_if_unknown(hdr)

    out = nib.Nifti1Image(np.ascontiguousarray(data, np.float32), aff, header=hdr)
    out.set_qform(aff, code=1)
    out.set_sform(aff, code=1)
    nib.save(out, dst_path)

    return dst_path


def _format_result(output_file, mask_file, return_mask):
    '''
    Gives back imageBet.nii.gz for Preprocessing und mask for process_fMRI.py and regress.py
    '''
    if return_mask:
        return output_file, mask_file
    return output_file


def applyBET(
    input_file,
    frac,
    radius,
    horizontal_gradient=0.0,
    use_bet4animal=False,
    center=None,
    return_mask=False,
    output_path=None,
    verbose=False,
):
    """Apply the shared AIDAmri BET implementation."""
    if output_path is None:
        output_path = os.path.dirname(input_file)

    if use_bet4animal:
        print("Using BET for animal brains")
        print("Note: bet4animal requires that the AC-PC line of brain is parallel to Y-axis")
        w_value = 2  # smooth the surface (lissencephalic weighting)
        species_id = 6  #-> mouse # species_id = 5 -> rat
        output_file = os.path.join(
            output_path,
            os.path.basename(input_file).split(".")[0] + "AnimalBet.nii.gz",
        )

        tmp_bet = os.path.join(os.path.dirname(input_file), "bet4animal_tmp_out.nii.gz") # output name of bet4animnal
        tmp_mask = tmp_bet.replace(".nii.gz", "_mask.nii.gz") # bet4animal output mask
        mask_file = output_file.replace(".nii.gz", "_mask.nii.gz") #final mask

        # ----- fslreorient2std -----
        tmp_std = os.path.join(os.path.dirname(input_file), "bet4animal_reorient2std.nii.gz") #output of fslreorient2std/input in bet4animal

        cmd = ["fslreorient2std", input_file, tmp_std]
        subprocess.run(cmd, check=True)

        # OPTIONAL: sanity prints
        # print("tmp_hdr axcodes:", nib.aff2axcodes(nib.load(tmp_hdr).affine))
        # print("tmp_std axcodes:", nib.aff2axcodes(nib.load(tmp_std).affine))

        bet_in = tmp_std

        # print("Header-only reorientation saved:", tmp_hdr)
        # print("New axcodes:", nib.aff2axcodes(aff))
        if center is None:
            center, p = estimate_center_intensity_based(bet_in)
            #print(f"Using intensity-based center in reoriented image: {center}")
        else:
            #original_center = center
            center = transform_voxel_center_between_images(center, input_file, bet_in)
            #print(f"Using user-defined center transformed for reoriented image: {original_center} -> {center}")
            center = [int(round(c)) for c in center]
        cx, cy, cz = center

        cmd = [
            "/aida/bin/bet4animal",
            bet_in,
            tmp_bet,
            "-m",  # mask
            "-w", str(w_value),
            "-z", str(species_id),
            "-c", str(cx), str(cy), str(cz),
        ]
        subprocess.run(cmd, check=True)

        # ===== AFTER bet4animal =====
        # Nifti has to be reoriented to match the expected orientation and geometry of the AIDAmri pipeline (similar to real BET output) so that downstream steps remain compatible

        # ---------- (1) Reorient to original ----------
        input_img = nib.load(input_file)
        target_axcodes = nib.aff2axcodes(input_img.affine)

        img = nib.load(tmp_bet)
        data = img.get_fdata(dtype=np.float32)
        aff = img.affine

        ornt_cur = nib.orientations.io_orientation(aff)
        ornt_tgt = nib.orientations.axcodes2ornt(target_axcodes)
        transform = nib.orientations.ornt_transform(ornt_cur, ornt_tgt)

        data_lip = nib.orientations.apply_orientation(data, transform)
        aff_lip = aff @ nib.orientations.inv_ornt_aff(transform, img.shape)

        # ---------- (2) Set affine ----------
        hdr_final = img.header.copy()
        hdr_final.set_data_dtype(np.float32)
        hdr_final["pixdim"][0] = 1
        hdr_final["pixdim"][4:8] = 1
        set_default_xyzt_units_if_unknown(hdr_final)

        img_final = nib.Nifti1Image(
            np.ascontiguousarray(data_lip, dtype=np.float32),
            aff_lip,
            header=hdr_final
        )
        # To match FSL BET output
        img_final.set_qform(aff_lip, code=1)
        img_final.set_sform(aff_lip, code=1)

        nib.save(img_final, output_file)

        # print("Final orientation:", nib.aff2axcodes(aff_final))
        # print("Final offset:", aff_final[:3, 3])

        # ===== APPLY SAME POST-PROCESSING TO BET MASK =====
        if os.path.exists(tmp_mask):
            # (1) Reorient to LIP
            m_img = nib.load(tmp_mask)
            m_data = m_img.get_fdata(dtype=np.float32)
            m_aff = m_img.affine

            m_ornt_cur = nib.orientations.io_orientation(m_aff)
            m_transform = nib.orientations.ornt_transform(m_ornt_cur, ornt_tgt)

            m_data_lip = nib.orientations.apply_orientation(m_data, m_transform)
            m_aff_lip = m_aff @ nib.orientations.inv_ornt_aff(m_transform, m_img.shape)

            # Binarize mask + uint8
            m_bin = (m_data_lip > 0.5).astype(np.uint8)

            m_hdr_lip = m_img.header.copy()
            m_hdr_lip.set_data_dtype(np.uint8)
            m_hdr_lip["pixdim"][0] = 1
            m_hdr_lip["pixdim"][4:8] = 1
            set_default_xyzt_units_if_unknown(m_hdr_lip)

            m_out = nib.Nifti1Image(
                np.ascontiguousarray(m_bin, dtype=np.uint8),
                m_aff_lip,
                header=m_hdr_lip
            )
            # To match FSL BET output
            m_out.set_qform(m_aff_lip, code=1)
            m_out.set_sform(m_aff_lip, code=1)

            nib.save(m_out, mask_file)
            # print("Mask processed:", bet_mask_path)
        else:
            print("Warning: BET mask not found:", tmp_mask)

        # remove temp
        try:
            for tmp_file in [tmp_std, tmp_bet, tmp_mask]:
                if os.path.exists(tmp_file):
                    os.remove(tmp_file)
        except Exception:
            pass
    # FSL BET (human modified version)
    else:
        data = nib.load(input_file)
        imgTemp = data.get_fdata()
        # create 4x4 scaling matrix and scale by 10 to match human like brain size
        if verbose:
            print("Image dimensions before scaling:", data.header.get_zooms())

        scale = np.eye(4)
        scale[0, 0] = 10
        scale[1, 1] = 10
        scale[2, 2] = 10
        if verbose:
            print("Image dimensions after scaling:", (data.affine * scale)[:3, :3])

        # Create new Nifti image with scaled affine
        scaled_affine = data.affine @ scale
        scaledNiiData = nib.Nifti1Image(imgTemp, scaled_affine)
        hdrIn = scaledNiiData.header
        set_default_xyzt_units_if_unknown(hdrIn)

        fslPath = os.path.join(os.path.dirname(input_file), 'fslScaleTemp.nii.gz')
        nib.save(scaledNiiData, fslPath)

        # temporary BET input with header-only LIP -> LPI world swap
        # This insures that the vertical gradient works in horizontally (so snout to cerebellum). This is needed bc this is a issue in FSL BET used for mice
        # Any questions regarding this ask Julian and hope he still knows
        if verbose:
            print("Saved scaled image to:", fslPath)
            print("Image dimensions:", scaledNiiData.header.get_zooms())

        bet_input_tmp = os.path.join(os.path.dirname(input_file), 'fslScaleTemp_LPIhdr.nii.gz')
        save_header_only_reoriented_copy(
            fslPath,
            bet_input_tmp,
            swaps=FSL_BET_WORLD_SWAPS,
        )

        # final output path
        output_file = os.path.join(
            output_path,
            os.path.basename(input_file).split(".")[0] + "Bet.nii.gz",
        )

        # temporary BET output in LPI-header space
        bet_output_tmp = os.path.join(
            os.path.dirname(input_file),
            os.path.basename(input_file).split(".")[0] + 'Bet_LPIhdr_tmp.nii.gz',
        )

        if center is not None:
            # nipype BET requires ints
            center_int = [int(round(c)) for c in center]
            print(f"Using user-defined center (rounded): {center_int}")
            myBet = fsl.BET(
                in_file=bet_input_tmp,
                out_file=bet_output_tmp,
                frac=frac,
                radius=radius,
                vertical_gradient=horizontal_gradient,
                center=center_int,
                mask=True
            )
        else:
            print("Using robust center estimation (-R)")
            myBet = fsl.BET(
                in_file=bet_input_tmp,
                out_file=bet_output_tmp,
                frac=frac,
                radius=radius,
                vertical_gradient=horizontal_gradient,
                robust=True,
                mask=True
            )
        if verbose:
            print(myBet.cmdline)
        myBet.run()

        # backswap BET image: LPI -> LIP (same swap again, because swap is its own inverse)
        save_header_only_reoriented_copy(
            bet_output_tmp,
            output_file,
            swaps=FSL_BET_WORLD_SWAPS,
        )

        # backswap BET mask: LPI -> LIP
        mask_tmp_file = bet_output_tmp.replace(".nii.gz", "_mask.nii.gz")
        mask_file = output_file.replace(".nii.gz", "_mask.nii.gz")

        if os.path.exists(mask_tmp_file):
            save_header_only_reoriented_copy(
                mask_tmp_file,
                mask_file,
                swaps=FSL_BET_WORLD_SWAPS,
            )

        # unscale result data by factor 10ˆ(-1)
        dataOut = nib.load(output_file)
        imgOut = dataOut.get_fdata(dtype=np.float32)
        # rescale nifti
        inv_scale = np.eye(4)
        inv_scale[0, 0] = 0.1
        inv_scale[1, 1] = 0.1
        inv_scale[2, 2] = 0.1

        unscaled_affine = dataOut.affine @ inv_scale
        unscaledNiiData = nib.Nifti1Image(imgOut, unscaled_affine)
        unscaledNiiData.set_qform(unscaled_affine, code=1)
        unscaledNiiData.set_sform(unscaled_affine, code=1)
        hdrOut = unscaledNiiData.header
        set_default_xyzt_units_if_unknown(hdrOut)
        nib.save(unscaledNiiData, output_file)
        if verbose:
            print("Image dimensions after unscaling:", unscaledNiiData.header.get_zooms())

        # also unscale BET mask
        mask_file = output_file.replace('.nii.gz', '_mask.nii.gz')
        if os.path.exists(mask_file):
            mask_data = nib.load(mask_file)
            bet_ref = nib.load(output_file)
            # make binary mask and apply affine of BET NIFTI
            mask_img = (mask_data.get_fdata() > 0.5).astype(np.uint8)

            finalMask = nib.Nifti1Image(mask_img, bet_ref.affine)
            finalMask.set_qform(bet_ref.affine, code=1)
            finalMask.set_sform(bet_ref.affine, code=1)

            hdrMask = finalMask.header
            hdrMask.set_data_dtype(np.uint8)
            set_default_xyzt_units_if_unknown(hdrMask)

            nib.save(finalMask, mask_file)
        # delete temporary files
        for tmp_file in [fslPath, bet_input_tmp, bet_output_tmp, mask_tmp_file]:
            if os.path.exists(tmp_file):
                os.remove(tmp_file)

    print(f"Brain extraction completed, output saved to {output_file}")
    return _format_result(output_file, mask_file, return_mask)
