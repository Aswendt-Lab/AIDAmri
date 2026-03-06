"""
Created on 10/08/2017

@author: Niklas Pallast
Neuroimaging & Neuroengineering
Department of Neurology
University Hospital Cologne

"""


import nipype.interfaces.fsl as fsl
import os,sys
import nibabel as nib
import numpy as np
import applyMICO
import subprocess
import shutil
import nipype.interfaces.ants as ants

def reset_orientation(input_file):

    brkraw_dir = os.path.join(os.path.dirname(input_file), "brkraw")
    if os.path.exists(brkraw_dir):
        return 

    os.mkdir(brkraw_dir)
    dst_path = os.path.join(brkraw_dir, os.path.basename(input_file))

    shutil.copyfile(input_file, dst_path)

    data = nib.load(input_file)
    raw_img = data.dataobj.get_unscaled()

    raw_nii = nib.Nifti1Image(raw_img, data.affine)
    nib.save(raw_nii, input_file)

    delete_orient_command = f"fslorient -deleteorient {input_file}"
    subprocess.run(delete_orient_command, shell=True)

    # Befehl zum Festlegen der radiologischen Orientierung
    forceradiological_command = f"fslorient -forceradiological {input_file}"
    subprocess.run(forceradiological_command, shell=True)

def post_flip_and_canonicalize(in_path: str, out_path: str | None = None) -> str:
    """
    Apply the same processing convention as the FSL-BET branch:
    flip axis 2 (z) and convert to closest canonical orientation.
    If out_path is None, overwrite in_path.
    Returns output path.
    """
    if out_path is None:
        out_path = in_path

    img = nib.load(in_path)
    data = img.get_fdata()

    data = np.flip(data, 2)
    out_img = nib.Nifti1Image(data, img.affine, img.header)
    out_img = nib.as_closest_canonical(out_img)

    nib.save(out_img, out_path)
    return out_path

def n4biasfieldcorr(input_file):
    output_file = os.path.join(os.path.dirname(input_file), os.path.basename(input_file).split('.')[0] + 'Bias.nii.gz')
    # Note: shrink_factor is set to 4 to speed up the process, but can be adjusted
    myAnts = ants.N4BiasFieldCorrection(input_image=input_file, output_image=output_file,
                                        shrink_factor=2, bspline_fitting_distance=20,
                                        bspline_order=3, n_iterations=[50, 50, 50, 50, 0], dimension=3)
    myAnts.run()
    print("Biasfield correction completed")
    return output_file

def copy_xform(ref_file, dst_file):
    ref = nib.load(ref_file)
    dst = nib.load(dst_file)

    data = dst.get_fdata(dtype=np.float32)

    new = nib.Nifti1Image(data, ref.affine, header=dst.header)

    # Copy qform from reference
    qaff, qcode = ref.get_qform(coded=True)
    if qaff is not None:
        new.set_qform(qaff, int(qcode))

    # Set sform explicitly → Code = 2 (aligned anatomical)
    saff, _ = ref.get_sform(coded=True)
    if saff is None:
        saff = ref.affine

    new.set_sform(saff, code=2)

    nib.save(new, dst_file)

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


def applyBET(input_file,frac,radius,vertical_gradient,use_bet4animal=False, species='mouse', center= None):
    """Apply BET"""
    if use_bet4animal == True:
        # Use BET for animal brains
        print("Using BET for animal brains")
        print("Note: bet4animal requires that the AC-PC line of brain is parallel to Y-axis")
        w_value = 2 #smooth the surface (lissencephalic weighting)
        species_id = 6 if species == 'mouse' else 5
        output_file = os.path.join(os.path.dirname(input_file), os.path.basename(input_file).split('.')[0] + 'Bet.nii.gz')
        #----- Reorient-----#
        world_swaps = [(1, 2)]
        world_flips = [1, 2]
        # -----------------------------------------------

        tmp_hdr = os.path.join(os.path.dirname(input_file), "bet4animal_hdrtmp.nii.gz")

        img = nib.load(input_file)

        # keep data unchanged (header-only operation)
        try:
            data = img.dataobj.get_unscaled()
            data = np.asanyarray(data)
        except Exception:
            data = img.get_fdata(dtype=np.float32)

        aff = img.affine.copy()

        # swaps first
        for a, b in world_swaps:
            aff[[a, b], :] = aff[[b, a], :]

        # flips after
        for ax in world_flips:
            aff[ax, :] *= -1

        hdr = img.header.copy()
        hdr["pixdim"][0] = 1
        hdr["pixdim"][4:8] = 1
        hdr.set_data_dtype(np.float32)

        tmp_img = nib.Nifti1Image(np.ascontiguousarray(data, dtype=np.float32), aff, header=hdr)
        tmp_img.set_qform(aff, code=1)
        tmp_img.set_sform(aff, code=1)
        nib.save(tmp_img, tmp_hdr)

        # ----- fslreorient2std -----
        tmp_std = os.path.join(os.path.dirname(input_file), "bet4animal_reorient2std.nii.gz")

        cmd = ["fslreorient2std", tmp_hdr, tmp_std]
        subprocess.run(cmd, check=True)

        # OPTIONAL: sanity prints
        # print("tmp_hdr axcodes:", nib.aff2axcodes(nib.load(tmp_hdr).affine))
        # print("tmp_std axcodes:", nib.aff2axcodes(nib.load(tmp_std).affine))

        # ab jetzt tmp_std als Input für bet4animal verwenden:
        bet_in = tmp_std

        #print("Header-only reorientation saved:", tmp_hdr)
        #print("New axcodes:", nib.aff2axcodes(aff))
        if center is None:
            center, p = estimate_center_intensity_based(bet_in)
        cx, cy, cz = center

        command = (
            f"/aida/bin/bet4animal {bet_in} {output_file} "
            f"-f {frac} -m -w {w_value} -z {species_id} -c {cx} {cy} {cz}"
        )  # m = binary mask output
        subprocess.run(command, shell=True, check=True)

        # remove temp

        try:
            os.remove(tmp_hdr)
            os.remove(tmp_std)
        except Exception:
            pass
        


        # ===== AFTER bet4animal =====

        # ---------- (1) Reorient to LIP ----------
        target_axcodes = ('L', 'I', 'P')

        img = nib.load(output_file)
        data = img.get_fdata(dtype=np.float32)
        aff = img.affine

        ornt_cur = nib.orientations.io_orientation(aff)
        ornt_tgt = nib.orientations.axcodes2ornt(target_axcodes)
        transform = nib.orientations.ornt_transform(ornt_cur, ornt_tgt)

        data_lip = nib.orientations.apply_orientation(data, transform)
        aff_lip = aff @ nib.orientations.inv_ornt_aff(transform, img.shape)

        # ---------- (2) Header-only swaps/flips ----------
        world_swaps = [(1, 2)]
        world_flips = [1, 2]

        aff2 = aff_lip.copy()

        # swaps (safe)
        for a, b in world_swaps:
            tmp = aff2[[a, b], :].copy()
            aff2[[a, b], :] = tmp[::-1, :]

        # flips
        for ax in world_flips:
            aff2[ax, :] *= -1

        # ---------- (3) Flip data axis 2 ----------
        data_flip = np.flip(data_lip, 2)

        img_flip = nib.Nifti1Image(
            np.ascontiguousarray(data_flip, dtype=np.float32),
            aff2
        )

        # ---------- (4) Closest canonical ----------
        img_final = nib.as_closest_canonical(img_flip)

        # ---------- (5) Set affine offset to 0 ----------
        aff_final = img_final.affine.copy()
        aff_final[:3, 3] = 0

        hdr_final = img_final.header.copy()
        hdr_final.set_data_dtype(np.float32)
        hdr_final["pixdim"][0] = 1
        hdr_final["pixdim"][4:8] = 1
        hdr_final.set_xyzt_units('mm', 'sec')

        img_final2 = nib.Nifti1Image(
            np.ascontiguousarray(img_final.get_fdata(dtype=np.float32), dtype=np.float32),
            aff_final,
            header=hdr_final
        )

        img_final2.set_qform(aff_final, code=0)
        img_final2.set_sform(aff_final, code=2)

        nib.save(img_final2, output_file)

        #print("Final orientation:", nib.aff2axcodes(aff_final))
        #print("Final offset:", aff_final[:3, 3])

        # ===== APPLY SAME POST-PROCESSING TO BET MASK =====
        bet_mask_path = output_file.replace(".nii.gz", "_mask.nii.gz")
        if os.path.exists(bet_mask_path):
            # (1) Reorient to LIP
            m_img = nib.load(bet_mask_path)
            m_data = m_img.get_fdata(dtype=np.float32)
            m_aff = m_img.affine

            m_ornt_cur = nib.orientations.io_orientation(m_aff)
            m_ornt_tgt = nib.orientations.axcodes2ornt(target_axcodes)  # ('L','I','P')
            m_transform = nib.orientations.ornt_transform(m_ornt_cur, m_ornt_tgt)

            m_data_lip = nib.orientations.apply_orientation(m_data, m_transform)
            m_aff_lip = m_aff @ nib.orientations.inv_ornt_aff(m_transform, m_img.shape)

            # (2) Header-only swaps/flips
            m_aff2 = m_aff_lip.copy()
            for a, b in world_swaps:
                tmp = m_aff2[[a, b], :].copy()
                m_aff2[[a, b], :] = tmp[::-1, :]
            for ax in world_flips:
                m_aff2[ax, :] *= -1

            # (3) Flip data axis 2
            m_data_flip = np.flip(m_data_lip, 2)

            m_img_flip = nib.Nifti1Image(
                np.ascontiguousarray(m_data_flip, dtype=np.float32),
                m_aff2
            )

            # (4) Closest canonical
            m_img_final = nib.as_closest_canonical(m_img_flip)

            # (5) Set affine offset to 0
            m_aff_final = m_img_final.affine.copy()
            m_aff_final[:3, 3] = 0

            # Binarize mask + uint8
            m_bin = (m_img_final.get_fdata(dtype=np.float32) > 0.5).astype(np.uint8)

            m_hdr_final = m_img_final.header.copy()
            m_hdr_final.set_data_dtype(np.uint8)
            m_hdr_final["pixdim"][0] = 1
            m_hdr_final["pixdim"][4:8] = 1
            m_hdr_final.set_xyzt_units('mm', 'sec')

            m_out = nib.Nifti1Image(
                np.ascontiguousarray(m_bin, dtype=np.uint8),
                m_aff_final,
                header=m_hdr_final
            )
            m_out.set_qform(m_aff_final, code=0)
            m_out.set_sform(m_aff_final, code=2)

            nib.save(m_out, bet_mask_path)
            #print("Mask processed:", bet_mask_path)
        else:
            print("Warning: BET mask not found:", bet_mask_path)

    else:
        data = nib.load(input_file)
        imgTemp = data.get_fdata()
        # create 4x4 scaling matrix
        scale = np.eye(4) * 10
        #Set last element to 1 (important for affine matrix)
        scale[3][3] = 1
        imgTemp = np.flip(imgTemp, 2)

        #Create new Nifti image with flipped data and scaled affine
        scaledNiiData = nib.Nifti1Image(imgTemp, data.affine * scale)
        hdrIn = scaledNiiData.header
        hdrIn.set_xyzt_units('mm')
        scaledNiiData = nib.as_closest_canonical(scaledNiiData)

        fslPath = os.path.join(os.path.dirname(input_file), 'fslScaleTemp.nii.gz')
        nib.save(scaledNiiData, fslPath)

        # set output file path
        output_file = os.path.join(os.path.dirname(input_file), os.path.basename(input_file).split('.')[0] + 'Bet.nii.gz')
        if center is not None:
            # nipype BET requires ints
            center_int = [int(round(c)) for c in center]
            print(f"Using user-defined center (rounded): {center_int}")
            myBet = fsl.BET(
                in_file=fslPath,
                out_file=output_file,
                frac=frac,
                radius=radius,
                vertical_gradient=vertical_gradient,
                center=center_int,
                mask=True
            )
        else:
            print("Using robust center estimation (-R)")
            myBet = fsl.BET(
                in_file=fslPath,
                out_file=output_file,
                frac=frac,
                radius=radius,
                vertical_gradient=vertical_gradient,
                robust=True,
                mask=True
            )
        myBet.run()
        os.remove(fslPath) # remove temporary scaled file

        # unscale result data by factor 10ˆ(-1)
        dataOut = nib.load(output_file)
        imgOut = dataOut.get_fdata(dtype=np.float32)
        scale = np.eye(4) / 10
        scale[3][3] = 1
        #create unscaled Nifti image with unscaled affine and flip
        unscaledNiiData = nib.Nifti1Image(imgOut, dataOut.affine * scale)
        hdrOut = unscaledNiiData.header
        hdrOut.set_data_dtype(np.float32)
        hdrOut.set_xyzt_units('mm', 'sec')
        hdrOut["pixdim"][0] = 1
        hdrOut["pixdim"][4:8] = 1
        nib.save(unscaledNiiData, output_file)

        # also unscale BET mask
        mask_file = output_file.replace('.nii.gz', '_mask.nii.gz')
        if os.path.exists(mask_file):
            mask_data = nib.load(mask_file)
            mask_img = mask_data.get_fdata()
            #unscale mask affine
            scale = np.eye(4) / 10
            scale[3][3] = 1
            #make binary mask and apply unscaled affine
            unscaledMask = nib.Nifti1Image(
                (mask_img > 0.5).astype(np.uint8),
                mask_data.affine * scale
            )
            hdrMask = unscaledMask.header
            hdrMask.set_data_dtype(np.uint8)
            hdrMask.set_xyzt_units('mm', 'sec')
            hdrMask["pixdim"][0] = 1
            hdrMask["pixdim"][4:8] = 1
            #set offset to 0 (important for BET/mask overlay)
            aff_mask = unscaledMask.affine.copy()
            aff_mask[:3, 3] = 0

            finalMask = nib.Nifti1Image(
                (mask_img > 0.5).astype(np.uint8),
                aff_mask,
                header=hdrMask
            )
            #inline with Bet image
            finalMask.set_qform(aff_mask, code=0)
            finalMask.set_sform(aff_mask, code=2)

            nib.save(finalMask, mask_file)
    print(f"Brain extraction completed, output saved to {output_file}")
    return output_file

#%% Program

if __name__ == "__main__":
    import argparse


    parser = argparse.ArgumentParser(description='Preprocessing of T2 Data')

    requiredNamed = parser.add_argument_group('Required named arguments')
    requiredNamed.add_argument('-i','--input_file', help='path to input file',required=True)

    parser.add_argument(
        '-f',
        '--frac',
        help='Fractional intensity threshold - default=0.15  smaller values give larger brain outline estimates',
        nargs='?',
        type=float,
        default=0.15,
        )
    parser.add_argument(
        '-r', 
        '--radius',
        help='Head radius (mm not voxels) - default=45',
        nargs='?',
        type=int,
        default=45,
        )
    parser.add_argument(
        '-g',
        '--vertical_gradient',
        help='Vertical gradient in fractional intensity threshold - default=0.0   positive values give larger brain outlines at bottom and smaller brain outlines at top',
        nargs='?',
        type=float,
        default=0.0,
        )
    parser.add_argument(
        '-c', '--center',
        help='Brain center in voxel coordinates: x y z',
        nargs=3,
        type=float,
        default=None
    )

    parser.add_argument(
        '-b',
        '--bias_skip',
        help='Skip bias field correction',
        action='store_true'
    )

    parser.add_argument(
        '--bet_skip',
        help='Skip BET during T2 preprocessing (still creates *Bet.nii.gz as copy for pipeline compatibility)', #Output will stil be named as '*Bet.nii.gz' but will be identical to bias-corrected (or original) image
        action='store_true'
    )

    parser.add_argument(
        '--bias_method',
        help='Biasfield correction method - default="mico", other options are "mico" or "ants"',
        nargs='?',
        type=str,
        default="mico",
        )
    parser.add_argument(
        '--use_bet4animal',
        help='Use BET for animal brains',
        action='store_true'
    )

    args = parser.parse_args()

    # set Parameters
    input_file = None
    if args.input_file is not None and args.input_file is not None:
        input_file = args.input_file
    if not os.path.exists(input_file):
        sys.exit(f"Error: input file does not exist: {input_file}")


    frac = args.frac
    radius = args.radius
    vertical_gradient = args.vertical_gradient
    bias_skip = args.bias_skip
    bias_method = args.bias_method

    print(f"Frac: {frac} Radius: {radius} Gradient {vertical_gradient}")

    reset_orientation(input_file)
    print("Orientation resetted to RAS")

    #intensity correction using non parametric bias field correction algorithm
    print("Starting Biasfieldcorrection:")
    if not args.bias_skip:
        if bias_method == "mico":
            try:
                outputBiasCorr = applyMICO.run_MICO(input_file, os.path.dirname(input_file))
                print("Biasfield correction was successful")
            except Exception as e:
                print(f'Error in bias field correction\nError message: {str(e)}')
                raise
        elif bias_method == "ants":
            try:
                outputBiasCorr = n4biasfieldcorr(input_file=input_file)
                copy_xform(input_file, outputBiasCorr)
                print("Biasfield correction was successful")
            except Exception as e:
                print(f'Error in bias field correction\nError message: {str(e)}')
                raise
    else:
        outputBiasCorr = input_file
    
    print(os.path.exists(outputBiasCorr))

    use_bet4animal = args.use_bet4animal

    if args.bet_skip:
        print("Skipping BET")
        outputBET = os.path.join(
            os.path.dirname(outputBiasCorr),
            os.path.basename(outputBiasCorr).split('.')[0] + 'Bet.nii.gz'
        )
        # --- reproduce BET pre-steps: flip axis 2 + closest canonical ---
        src = nib.load(outputBiasCorr)
        data = src.get_fdata()

        data = np.flip(data, 2)  # match applyBET()
        bet_like = nib.Nifti1Image(data, src.affine)
        bet_like = nib.as_closest_canonical(bet_like)

        nib.save(bet_like, outputBET)
        print(f"BET skipped -> copied to {outputBET}")  # downstream "output" is the bias-corrected (or original) image

        # Dummy BET mask
        bet_mask_path = outputBET.replace('.nii.gz', '_mask.nii.gz')
        bet_data = bet_like.get_fdata()

        mask = (bet_data > 0).astype(np.uint8)
        mask_img = nib.Nifti1Image(mask, bet_like.affine)
        mask_img.set_data_dtype(np.uint8)
        nib.save(mask_img, bet_mask_path)

        print(f"BET mask created at {bet_mask_path}")
    else:
        # brain extraction
        print("Starting brain extraction")
        try:
            outputBET = applyBET(input_file=outputBiasCorr,frac=frac,radius=radius,vertical_gradient=vertical_gradient,use_bet4animal=use_bet4animal, center=args.center)
            print("Brain extraction was successful")
        except Exception as e:
            print(f'Error in brain extraction\nError messsage: {str(e)}')
            raise
    
    print("Preprocessing completed")
 










