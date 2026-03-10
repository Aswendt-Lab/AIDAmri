"""
Created on 10/08/2017

@author: Niklas Pallast
Neuroimaging & Neuroengineering
Department of Neurology
University Hospital Cologne

Edited by Paul Camacho 2025

"""


import nipype.interfaces.fsl as fsl
import os, sys
import nibabel as nib
import numpy as np
import applyMICO
import nipype.interfaces.ants as ants
import subprocess
import shutil
import averageb0
import dipy.denoise.patch2self as patch2self

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

def estimate_center_intensity_based(nifti, percentile=60):
    """
    Estimate BET center (-c) using intensity-weighted center-of-gravity (fslstats -C),
    excluding low-intensity voxels using a data-adaptive threshold (-l = P{percentile}).
    Returns center as floats (voxel coordinates).
    """
    p = subprocess.check_output(["fslstats", nifti, "-P", str(percentile)]).decode().strip()
    center = subprocess.check_output(["fslstats", nifti, "-l", p, "-C"]).decode().strip().split()
    cx, cy, cz = [float(v) for v in center]
    return [cx, cy, cz], float(p)

def skip_bet_function(input_file):
    """
    Create BET-like output and full-brain dummy mask for pipeline compatibility
    when BET is skipped.
    Reproduces the key geometry/orientation steps of applyBET(..., use_bet4animal=False).
    """
    output_file = os.path.join(
        os.path.dirname(input_file),
        os.path.basename(input_file).split('.')[0] + 'Bet.nii.gz'
    )

    src = nib.load(input_file)
    data = src.get_fdata(dtype=np.float32)

    # Reproduce main pre-BET manipulations from applyBET() non-bet4animal branch
    data = np.flip(data, 2)

    bet_like = nib.Nifti1Image(data, src.affine)
    bet_like.header.set_data_dtype(np.float32)
    bet_like.header.set_xyzt_units('mm', 'sec')

    bet_like = nib.as_closest_canonical(bet_like)

    # Match downstream expectations a bit more closely
    aff = bet_like.affine.copy()
    aff[:3, 3] = 0

    hdr = bet_like.header.copy()
    hdr.set_data_dtype(np.float32)
    hdr["pixdim"][0] = 1
    hdr["pixdim"][4:8] = 1
    hdr.set_xyzt_units('mm', 'sec')

    final_img = nib.Nifti1Image(
        np.ascontiguousarray(bet_like.get_fdata(dtype=np.float32), dtype=np.float32),
        aff,
        header=hdr
    )
    final_img.set_qform(aff, code=0)
    final_img.set_sform(aff, code=2)

    nib.save(final_img, output_file)

    # Create full mask for pipeline compatibility
    bet_mask_path = output_file.replace('.nii.gz', '_mask.nii.gz')
    img_data = final_img.get_fdata(dtype=np.float32)

    mask = (img_data > 0).astype(np.uint8)

    mask_hdr = final_img.header.copy()
    mask_hdr.set_data_dtype(np.uint8)
    mask_hdr["pixdim"][0] = 1
    mask_hdr["pixdim"][4:8] = 1
    mask_hdr.set_xyzt_units('mm', 'sec')

    mask_img = nib.Nifti1Image(
        np.ascontiguousarray(mask, dtype=np.uint8),
        aff,
        header=mask_hdr
    )
    mask_img.set_qform(aff, code=0)
    mask_img.set_sform(aff, code=2)

    nib.save(mask_img, bet_mask_path)

    print(f"BET skipped -> created compatibility image: {output_file}")
    print(f"BET skipped -> created compatibility mask: {bet_mask_path}")

    return output_file

def applyBET(input_file, frac=0.40, radius=6, vertical_gradient=0.0,
             use_bet4animal=False, species='mouse', verbose=True, center=None):
    """Apply BET"""
    if use_bet4animal:
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
        # scale Nifti data by factor 10
        data = nib.load(input_file)
        imgTemp = data.get_fdata()
        if verbose:
            print("Image dimensions before scaling:", data.header.get_zooms())
        scale = np.eye(4) * 10
        scale[3][3] = 1
        imgTemp = np.flip(imgTemp, 2)
        if verbose:
            print("Image dimensions after scaling:", (data.affine * scale)[:3,:3])

        scaledNiiData = nib.Nifti1Image(imgTemp, data.affine * scale)
        if verbose:
            print("Image dimensions after flipping:", scaledNiiData.header.get_zooms())
        hdrIn = scaledNiiData.header
        hdrIn.set_xyzt_units('mm')
        scaledNiiData = nib.as_closest_canonical(scaledNiiData)

        fslPath = os.path.join(os.path.dirname(input_file), 'fslScaleTemp.nii.gz')
        nib.save(scaledNiiData, fslPath)
        if verbose:
            print("Saved scaled image to:", fslPath)
            print("Image dimensions:", scaledNiiData.header.get_zooms())

        # extract brain
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
        os.remove(fslPath)

        # unscale result data by factor 10ˆ(-1)
        dataOut = nib.load(output_file)
        imgOut = dataOut.get_fdata()
        scale = np.eye(4) / 10
        scale[3][3] = 1

        unscaledNiiData = nib.Nifti1Image(imgOut, dataOut.affine * scale)
        hdrOut = unscaledNiiData.header
        hdrOut.set_xyzt_units('mm', 'sec')
        if verbose:
            print("Image dimensions after unscaling:", unscaledNiiData.header.get_zooms())
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


def denoise_patch2self(input_file, output_path, b0_thresh=100):
    """
    Denoises the input DTI image using Patch2Self from DIPY.
    Requires an appropriate input file (input_file) and the output path (output_path).
    """
    bvalsname = input_file.replace(".nii.gz", ".bval")
    if not os.path.exists(bvalsname):
        try:
            bvalsname = input_file.replace(".nii.gz", ".btable")
            btable = np.loadtxt(bvalsname, dtype=float)
            bvalsname = os.path.splitext(bvalsname)[0] + ".bval"
            np.savetxt(bvalsname, btable[0, :], fmt='%.6f')
        except:
            sys.exit(f"Error: bvals file {bvalsname} not found.")
    bvals = np.loadtxt(bvalsname, dtype=float)
    data = nib.load(input_file)
    img = data.get_fdata()
    affine = data.affine
    debug = True
    if debug:
        print("Debugging information:")
        print("Image header:", data.header)
        print("Affine matrix:", affine)
        print("Image sform:", data.header.get_sform())
    if img.ndim != 4:
        raise ValueError("Input image must be a 4D NIfTI file.")
    
    # Apply Patch2Self denoising
    denoised_img = patch2self.patch2self(img, bvals, b0_threshold=b0_thresh, model='ols', out_dtype=np.int16)
    
    # Save the denoised image
    output_file = os.path.join(output_path, os.path.basename(input_file).split('.')[0] + 'Patch2SelfDenoised.nii.gz')
    denoised_nii = nib.Nifti1Image(denoised_img, affine)
    denoised_nii.header.set_xyzt_units('mm')
    if debug:
        print("Denoised image header:", denoised_nii.header)
        print("Denoised affine matrix:", denoised_nii.affine)
        print("Denoised image sform:", denoised_nii.header.get_sform())
    nib.save(denoised_nii, output_file)

    # # Copy header from original image to denoised image using fslcpgeom
    # myFslCpGeom = fsl.utils.CopyGeom(dest_file=output_file, in_file=input_file)
    # myFslCpGeom.run()
    # print(f"Denoising completed, output saved to {output_file}")
    # if debug is True:
    #     output = nii.load(output_file)
    #     print("Final denoised image header after copying geometry:", output.header)
    #     print("Final denoised affine matrix after copying geometry:", output.affine)
    #     print("Final denoised image sform after copying geometry:", output.header.get_sform())
    return output_file



def dwibiasfieldcorr(input_file,outputPath):
    output_file = os.path.join(outputPath, os.path.basename(input_file).split('.')[0] + 'N4Bias.nii.gz')
    # Note: shrink_factor is set to 4 to speed up the process, but can be adjusted
    myAnts = ants.N4BiasFieldCorrection(input_image=input_file,output_image=output_file,
                                        shrink_factor=4,bspline_fitting_distance=20,
                                        bspline_order=3,n_iterations=[1000,0],dimension=3)
    myAnts.run()
    print("Biasfield correction completed")
    return output_file


def smoothIMG(input_file, output_path,skip_min=False):
    """
    Smoothes image via FSL. Only input and output has do be specified. Parameters are fixed to box shape and to the kernel size of 0.1 voxel.
    If skip_min is True, the median filter is not applied, and the image is directly smoothed (No "DN" files are produced).
    """
    data = nib.load(input_file)
    vol = data.get_fdata()
    if not skip_min:
        ImgSmooth = np.min(vol, 3)
        unscaledNiiData = nib.Nifti1Image(ImgSmooth, data.affine)
        hdrOut = unscaledNiiData.header
        hdrOut.set_xyzt_units('mm')
        output_file = os.path.join(os.path.dirname(input_file),
                                   os.path.basename(input_file).split('.')[0] + 'DN.nii.gz')
        nib.save(unscaledNiiData, output_file)
        input_file = output_file
    else:
        ImgSmooth = vol

    output_file = os.path.join(output_path, os.path.basename(input_file).split('.')[0] + 'Smooth.nii.gz')
    myGauss =  fsl.SpatialFilter(
        in_file = input_file,
        out_file = output_file, 
        operation = 'median',
        kernel_shape = 'box',
        kernel_size = 0.1
    )
    myGauss.run()
    return output_file

def thresh(input_file, output_path):
    #output_file = os.path.join(os.path.dirname(input_file),os.path.basename(input_file).split('.')[0]+ 'Thres.nii.gz')
    output_file = os.path.join(output_path, os.path.basename(input_file).split('.')[0] + 'Thres.nii.gz')
    myThres = fsl.Threshold(in_file=input_file,out_file=output_file,thresh=20)#,direction='above')
    myThres.run()
    return output_file

def cropToSmall(input_file,output_path):
    #output_file = os.path.join(os.path.dirname(input_file),os.path.basename(input_file).split('.')[0]  + 'Crop.nii.gz')
    output_file = os.path.join(output_path, os.path.basename(input_file).split('.')[0] + 'Crop.nii.gz')
    myCrop = fsl.ExtractROI(in_file=input_file,roi_file=output_file,x_min=40,x_size=130,y_min=50,y_size=110,z_min=0,z_size=12)
    myCrop.run()
    return  output_file


if __name__ == "__main__":
    import argparse


    parser = argparse.ArgumentParser(description='Preprocessing of DTI Data')

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '-i',
        '--input',
        help='Path to the raw NIfTI DTI file',
        required=True,
    )

    parser.add_argument(
        '-f',
        '--frac',
        help='Fractional intensity threshold - default=0.4, smaller values give larger brain outline estimates',
        type=float,
        default=0.4,
    )
    parser.add_argument(
        '-r',
        '--radius',
        help='Head radius (mm not voxels) - default=45',
        type=int,
        default=45,
    )
    parser.add_argument(
        '-g',
        '--vertical_gradient',
        help='Vertical gradient in fractional intensity threshold - default=0.0, positive values give larger brain outlines at bottom and smaller brain outlines at top',
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
        '-d',
        '--denoiser',
        help='Denoising method - default=None, other option is "patch2self"',
        choices = ["patch2self"],
        type=str.lower,
        default=None
    )
    parser.add_argument(
        '-b',
        '--bias_method',
        help='Biasfield correction method - default=None, other options are "mico" or "ants"',
        choices = ["mico", "ants"],
        type=str.lower,
        default=None,
    )
    parser.add_argument(
        '--use_bet4animal',
        help='Use BET for animal brains',
        action = 'store_true'
    )
    parser.add_argument(
        '--bet_skip',
        help='Skip BET during DTI preprocessing (still creates *Bet.nii.gz and *_mask.nii.gz for pipeline compatibility)',
        action='store_true'
    )
    parser.add_argument(
        '--average_b0',
        help='Average the b0 volumes',
        action='store_true'
    )
    parser.add_argument(
        '--skip_min',
        help='Skip the minimum filter before smoothing',
        action='store_true'
    )
    args = parser.parse_args()

    # set parameters
    input_file = args.input

    if not os.path.exists(input_file):
        sys.exit(f"Error: input file does not exist: {input_file}")
        
    frac = args.frac
    radius = args.radius
    vertical_gradient = args.vertical_gradient
    output_path = os.path.dirname(input_file)
    b0_thresh=100

    print(f"Frac: {frac} Radius: {radius} Gradient {vertical_gradient}")

    if reset_orientation_needed:
        reset_orientation(input_file)
        print("Orientation reset to RAS")
    
    if args.denoiser == "patch2self":
        # Denoising using Patch2Self
        denoised_image = denoise_patch2self(input_file, output_path, b0_thresh)
        print("Denoising completed, output saved to", denoised_image)
        # reset_orientation(denoised_image)
        input_file = denoised_image

    if args.average_b0:
        # Average b0 volumes
        b0image = averageb0.averageb0(input_file,b0_thresh)
        # # Copy header with fslcopygeom
        # myFslCpGeom = fsl.utils.CopyGeom(dest_file=b0image, in_file=input_file)
        # myFslCpGeom.run()
        input_file = b0image
        print("Averaging b0 volumes completed, output saved to", input_file)

    try:
        output_smooth = smoothIMG(input_file = input_file, output_path = output_path, skip_min=args.skip_min)
        print("Smoothing completed")
    except Exception as e:
        print(f'Fehler in der Biasfieldcorrecttion\nFehlermeldung: {str(e)}')
        raise
    
    if args.bias_method is None:
        print("No bias field correction applied")
        output_biascorr = output_smooth
    elif args.bias_method == "mico":
        # intensity correction using MICO
        try:
            output_biascorr = applyMICO.run_MICO(output_smooth, output_path)
            print("Biasfield correction was successful")
        except Exception as e:
            print(f'Error in bias field correction\nError message: {str(e)}')
            raise
    elif args.bias_method == "ants":
        # intensity correction using ANTs N4BiasFieldCorrection
        try:
            output_biascorr = dwibiasfieldcorr(input_file=output_smooth, outputPath=output_path)
            print("Biasfield correction was successful")
        except Exception as e:
            print(f'Error in bias field correction\nError message: {str(e)}')
            raise

    if args.bet_skip:
        print("Skipping brain extraction.")
        outputBET = skip_bet_function(output_biascorr)
    else:
        outputBET = applyBET(
            input_file=output_biascorr,
            frac=frac,
            radius=radius,
            vertical_gradient=vertical_gradient,
            use_bet4animal=args.use_bet4animal,
            center=args.center
        )
        print("Brain extraction was successful")