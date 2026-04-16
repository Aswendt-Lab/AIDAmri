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
import itertools
import sys
import threading
import time

FATAL_LIP_HEADER_EXIT_CODE = 86

def creat_brkraw_backup(input_file):

    brkraw_dir = os.path.join(os.path.dirname(input_file), "brkraw")
    if os.path.exists(brkraw_dir):
        return 

    os.mkdir(brkraw_dir)
    dst_path = os.path.join(brkraw_dir, os.path.basename(input_file))

    shutil.copyfile(input_file, dst_path)

    data = nib.load(input_file)
    raw_img = data.dataobj.get_unscaled().astype(np.float32)

    hdr = data.header.copy()
    hdr.set_data_dtype(np.float32)
    space_unit, time_unit = hdr.get_xyzt_units()

    if not space_unit or space_unit.lower() == "unknown":
        space_unit = "mm"
    if not time_unit or time_unit.lower() == "unknown":
        time_unit = "sec"

    hdr.set_xyzt_units(space_unit, time_unit)

    raw_nii = nib.Nifti1Image(raw_img, data.affine, header=hdr)
    raw_nii.set_qform(data.affine, code=1)
    raw_nii.set_sform(data.affine, code=1)
    nib.save(raw_nii, input_file)

def header_check(input_file):
    img = nib.load(input_file)
    axcodes = nib.aff2axcodes(img.affine)

    if axcodes != ("L", "I", "P"):
        print(
            f"Fatal header check failure: expected LIP orientation, found {axcodes} in {input_file}",
            file=sys.stderr,
        )
        sys.exit(FATAL_LIP_HEADER_EXIT_CODE)

    data = img.get_fdata(dtype=np.float32)

    out = nib.Nifti1Image(data, img.affine, header=img.header.copy())
    out.set_qform(img.affine, code=1)
    out.set_sform(img.affine, code=1)

    hdr = out.header
    hdr.set_data_dtype(np.float32)

    nib.save(out, input_file)
    return input_file

def set_default_xyzt_units_if_unknown(input_file):
    img = nib.load(input_file)
    hdr = img.header

    space_unit, time_unit = hdr.get_xyzt_units()
    if not space_unit or space_unit.lower() == "unknown":
        space_unit = "mm"
    if not time_unit or time_unit.lower() == "unknown":
        time_unit = "sec"

    hdr.set_xyzt_units(space_unit, time_unit)
    nib.save(img, input_file)
    return input_file

def n4biasfieldcorr(input_file):
    output_file = os.path.join(os.path.dirname(input_file), os.path.basename(input_file).split('.')[0] + 'AntsBias.nii.gz')
    # Note: shrink_factor is set to 4 to speed up the process, but can be adjusted
    myAnts = ants.N4BiasFieldCorrection(input_image=input_file,output_image=output_file,
                                        shrink_factor=4,bspline_fitting_distance=20,
                                        bspline_order=3,n_iterations=[1000,0],dimension=3)
    myAnts.run()
    print("Biasfield correction completed")
    return output_file

def spinner(stop_event, message="Working"):
    """
    Displays a simple terminal spinner while a long-running processing step is active.
    Does not report the actual progress of external tools such as FSL BET or ANTs.
    """
    for ch in itertools.cycle("|/-\\"):
        if stop_event.is_set():
            break
        sys.stdout.write(f"\r{message}... {ch}")
        sys.stdout.flush()
        time.sleep(0.1)
    sys.stdout.write(f"\r{message}... done\n")
    sys.stdout.flush()


def set_xform_codes_to_one(input_file):
    img = nib.load(input_file)
    img.set_qform(img.affine, code=1)
    img.set_sform(img.affine, code=1)
    nib.save(img, input_file)
    return input_file

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

def skip_bet_function(input_file):
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

    bet_like = nib.Nifti1Image(data, src.affine)

    # --- normalize header / affine like real BET output ---
    aff = bet_like.affine.copy()

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

    final_img.set_qform(aff, code=1)
    final_img.set_sform(aff, code=1)

    nib.save(final_img, outputBET)

    print(f"BET skipped -> created compatibility image: {outputBET}")

    # --- create dummy mask ---
    bet_mask_path = outputBET.replace('.nii.gz', '_mask.nii.gz')

    mask = (final_img.get_fdata(dtype=np.float32) > 0).astype(np.uint8)

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

    mask_img.set_qform(aff, code=1)
    mask_img.set_sform(aff, code=1)

    nib.save(mask_img, bet_mask_path)

    print(f"BET mask created at {bet_mask_path}")

    return outputBET

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

    out = nib.Nifti1Image(np.ascontiguousarray(data, np.float32), aff, header=hdr)
    out.set_qform(aff, code=1)
    out.set_sform(aff, code=1)
    nib.save(out, dst_path)

    return dst_path

def applyBET(input_file, frac, radius, horizontal_gradient,
             use_bet4animal=False, species='mouse', verbose=False, center=None):
    """Apply BET"""
    if use_bet4animal:
        # Use BET for animal brains
        print("Using BET for animal brains")
        print("Note: bet4animal requires that the AC-PC line of brain is parallel to Y-axis")
        w_value = 2 #smooth the surface (lissencephalic weighting)
        species_id = 6 if species == 'mouse' else 5
        output_file = os.path.join(os.path.dirname(input_file), os.path.basename(input_file).split('.')[0] + 'AnimalBet.nii.gz')

        tmp_bet = os.path.join(os.path.dirname(input_file), "bet4animal_tmp_out.nii.gz")
        tmp_mask = tmp_bet.replace(".nii.gz", "_mask.nii.gz")
        final_mask = output_file.replace(".nii.gz", "_mask.nii.gz")

        # ----- fslreorient2std -----
        tmp_std = os.path.join(os.path.dirname(input_file), "bet4animal_reorient2std.nii.gz")

        cmd = ["fslreorient2std", input_file, tmp_std]
        subprocess.run(cmd, check=True)

        # OPTIONAL: sanity prints
        # print("tmp_hdr axcodes:", nib.aff2axcodes(nib.load(tmp_hdr).affine))
        # print("tmp_std axcodes:", nib.aff2axcodes(nib.load(tmp_std).affine))

        bet_in = tmp_std

        #print("Header-only reorientation saved:", tmp_hdr)
        #print("New axcodes:", nib.aff2axcodes(aff))
        if center is None:
            center, p = estimate_center_intensity_based(bet_in)
        cx, cy, cz = center

        cmd = [
            "/aida/bin/bet4animal",
            bet_in,
            tmp_bet,
            "-f", str(frac),
            "-m", #mask
            "-w", str(w_value),
            "-z", str(species_id),
            "-c", str(cx), str(cy), str(cz),
        ]
        subprocess.run(cmd, check=True)

        # ===== AFTER bet4animal =====
        #Nifti has to be reoriented to match the expected orientation and geometry of the AIDAmri pipeline (similar to real BET output) so that downstream steps remain compatible

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
        hdr_final.set_xyzt_units('mm', 'sec')

        img_final = nib.Nifti1Image(
            np.ascontiguousarray(data_lip, dtype=np.float32),
            aff_lip,
            header=hdr_final
        )
        # To match FSL BET output
        img_final.set_qform(aff_lip, code=1)
        img_final.set_sform(aff_lip, code=1)

        nib.save(img_final, output_file)

        #print("Final orientation:", nib.aff2axcodes(aff_final))
        #print("Final offset:", aff_final[:3, 3])

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
            m_hdr_lip.set_xyzt_units('mm', 'sec')

            m_out = nib.Nifti1Image(
                np.ascontiguousarray(m_bin, dtype=np.uint8),
                m_aff_lip,
                header=m_hdr_lip
            )
            # To match FSL BET output
            m_out.set_qform(m_aff_lip, code=1)
            m_out.set_sform(m_aff_lip, code=1)

            nib.save(m_out, final_mask)
            #print("Mask processed:", bet_mask_path)
        else:
            print("Warning: BET mask not found:", bet_mask_path)

        # remove temp
        try:
            for tmp_file in [tmp_std, tmp_bet, tmp_mask]:
                if os.path.exists(tmp_file):
                    os.remove(tmp_file)
        except Exception:
            pass
    #FSL BET (human modified version)
    else:
        data = nib.load(input_file)
        imgTemp = data.get_fdata()
        if verbose:
            print("Image dimensions before scaling:", data.header.get_zooms())
        # create 4x4 scaling matrix and scale by 10 to match human like brain size
        scale = np.eye(4)
        scale[0, 0] = 10
        scale[1, 1] = 10
        scale[2, 2] = 10
        if verbose:
            print("Image dimensions after scaling:", (data.affine * scale)[:3,:3])

        #Create new Nifti image with scaled affine
        scaled_affine = data.affine @ scale
        scaledNiiData = nib.Nifti1Image(imgTemp, scaled_affine)
        hdrIn = scaledNiiData.header
        hdrIn.set_xyzt_units('mm')

        fslPath = os.path.join(os.path.dirname(input_file), 'fslScaleTemp.nii.gz')
        nib.save(scaledNiiData, fslPath)
        if verbose:
            print("Saved scaled image to:", fslPath)
            print("Image dimensions:", scaledNiiData.header.get_zooms())

        # temporary BET input with header-only LIP -> LPI world swap
        #This insures that the vertical gradient works in horizontally (so snout to cerebellum). This is needed bc this is a issue in FSL BET used for mice
        #Any questions regarding this ask Julian and hope he still knows
        bet_input_tmp = os.path.join(os.path.dirname(input_file), 'fslScaleTemp_LPIhdr.nii.gz')
        save_header_only_reoriented_copy(
            fslPath,
            bet_input_tmp,
            swaps=FSL_BET_WORLD_SWAPS,
        )

        # final output path
        output_file = os.path.join(
            os.path.dirname(input_file),
            os.path.basename(input_file).split('.')[0] + 'Bet.nii.gz'
        )

        # temporary BET output in LPI-header space
        bet_output_tmp = os.path.join(
            os.path.dirname(input_file),
            os.path.basename(input_file).split('.')[0] + 'Bet_LPIhdr_tmp.nii.gz'
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
        myBet.run()

        # backswap BET image: LPI -> LIP (same swap again, because swap is its own inverse)
        save_header_only_reoriented_copy(
            bet_output_tmp,
            output_file,
            swaps=FSL_BET_WORLD_SWAPS,
        )

        # backswap BET mask: LPI -> LIP
        mask_tmp_file = bet_output_tmp.replace('.nii.gz', '_mask.nii.gz')
        mask_file = output_file.replace('.nii.gz', '_mask.nii.gz')

        if os.path.exists(mask_tmp_file):
            save_header_only_reoriented_copy(
                mask_tmp_file,
                mask_file,
                swaps=FSL_BET_WORLD_SWAPS,
            )

        # unscale result data by factor 10ˆ(-1)
        dataOut = nib.load(output_file)
        imgOut = dataOut.get_fdata(dtype=np.float32)
        #rescale nifti
        inv_scale = np.eye(4)
        inv_scale[0, 0] = 0.1
        inv_scale[1, 1] = 0.1
        inv_scale[2, 2] = 0.1

        unscaled_affine = dataOut.affine @ inv_scale
        unscaledNiiData = nib.Nifti1Image(imgOut, unscaled_affine)
        unscaledNiiData.set_qform(unscaled_affine, code=1)
        unscaledNiiData.set_sform(unscaled_affine, code=1)
        hdrOut = unscaledNiiData.header
        hdrOut.set_xyzt_units('mm', 'sec')
        if verbose:
            print("Image dimensions after unscaling:", unscaledNiiData.header.get_zooms())
        nib.save(unscaledNiiData, output_file)

        # also unscale BET mask
        mask_file = output_file.replace('.nii.gz', '_mask.nii.gz')
        if os.path.exists(mask_file):
            mask_data = nib.load(mask_file)
            bet_ref = nib.load(output_file)
            #make binary mask and apply affine of BET NIFTI
            mask_img = (mask_data.get_fdata() > 0.5).astype(np.uint8)

            finalMask = nib.Nifti1Image(mask_img, bet_ref.affine)
            finalMask.set_qform(bet_ref.affine, code=1)
            finalMask.set_sform(bet_ref.affine, code=1)

            hdrMask = finalMask.header
            hdrMask.set_data_dtype(np.uint8)
            hdrMask.set_xyzt_units('mm', 'sec')

            nib.save(finalMask, mask_file)
        # delete temporary files
        for tmp_file in [fslPath, bet_input_tmp, bet_output_tmp, mask_tmp_file]:
            if os.path.exists(tmp_file):
                os.remove(tmp_file)
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
    debug = False
    if debug:
        print("Debugging information:")
        print("Image header:", data.header)
        print("Affine matrix:", affine)
        print("Image sform:", data.header.get_sform())
    if img.ndim != 4:
        raise ValueError("Input image must be a 4D NIfTI file.")
    
    # Apply Patch2Self denoising
    denoised_img = patch2self.patch2self(img, bvals, b0_threshold=b0_thresh, model='ols', out_dtype=np.float32)
    
    # Save the denoised image
    output_file = os.path.join(output_path, os.path.basename(input_file).split('.')[0] + 'Patch2SelfDenoised.nii.gz')
    denoised_nii = nib.Nifti1Image(denoised_img, affine)
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

def smoothIMG(input_file, output_path,skip_min=False):
    """
    Smoothes image via FSL. Only input and output has do be specified. Parameters are fixed to box shape and to the kernel size of 0.1 voxel.
    If skip_min is True, the median filter is not applied, and the image is directly smoothed (No "DN" files are produced).
    """
    data = nib.load(input_file)
    vol = data.get_fdata()
    if not skip_min:
        ImgSmooth = np.min(vol, 3).astype(np.float32)
        unscaledNiiData = nib.Nifti1Image(ImgSmooth, data.affine)
        unscaledNiiData.set_qform(data.affine, code=1)
        unscaledNiiData.set_sform(data.affine, code=1)

        hdrOut = unscaledNiiData.header
        hdrOut.set_data_dtype(np.float32)
        space_unit, time_unit = hdrOut.get_xyzt_units()

        if not space_unit or space_unit.lower() == "unknown":
            space_unit = "mm"
        if not time_unit or time_unit.lower() == "unknown":
            time_unit = "sec"
        hdrOut.set_xyzt_units(space_unit, time_unit)
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

    requiredNamed = parser.add_argument_group('Required named arguments')
    requiredNamed.add_argument(
        '-i',
        '--input_file',
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
        '--horizontal_gradient',
        help='Horizontal gradient in fractional intensity threshold - default=0.0. Not for bet4animals! Higher positive values make the BET stricter posterior and less stricter anterior (snout)',
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
        '--bet_skip',
        help='Skip BET during DTI preprocessing (still creates *Bet.nii.gz and *_mask.nii.gz for pipeline compatibility). '
             'If not set it uses FSL BET (modified human version)',
        action='store_true'
    )

    parser.add_argument(
        '-b',
        '--bias_method',
        help='Biasfield correction method - default="mico", other options are "ants" or "none"',
        choices = ["none", "mico", "ants"],
        type=str.lower,
        default=None,
    )

    parser.add_argument(
        '--use_bet4animal',
        help='Use BET for animal brains. '
             'If not set it use FSL (modified human version)',
        action='store_true'
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

    # set Parameters
    input_file = args.input_file
    if not os.path.exists(input_file):
        sys.exit(f"Error: input file does not exist: {input_file}")

    frac = args.frac
    radius = args.radius
    horizontal_gradient = args.horizontal_gradient
    bias_method = args.bias_method
    output_path = os.path.dirname(input_file)
    b0_thresh=100

    print(f"Frac: {frac} Radius: {radius} Gradient {horizontal_gradient}")

    creat_brkraw_backup(input_file)
    header_check(input_file)
    
    if args.denoiser == "patch2self":
        # Denoising using Patch2Self
        print("Starting denoising using patch2self")
        try:
            # start spinner
            stop_event = threading.Event()
            thread = threading.Thread(
                target=spinner,
                args=(stop_event, "Running Denoising with Patch2Self")
            )
            thread.start()

            try:
                denoised_image = denoise_patch2self(input_file, output_path, b0_thresh)
                set_xform_codes_to_one(denoised_image)
                set_default_xyzt_units_if_unknown(denoised_image)
            finally:
                stop_event.set()
                thread.join()

            print("Denoising completed, output saved to", denoised_image)
        except Exception as e:
            print(f'Error in Patch2Self denoising\nError message: {str(e)}')
            raise


        input_file = denoised_image

    if args.average_b0:
        # Average b0 volumes
        print("Starting averaging b0 volumes")
        try:
            # start spinner
            stop_event = threading.Event()
            thread = threading.Thread(
                target=spinner,
                args=(stop_event,)
            )
            thread.start()

            try:
                b0image = averageb0.averageb0(input_file,b0_thresh)
                set_xform_codes_to_one(b0image)
                set_default_xyzt_units_if_unknown(b0image)
            finally:
                stop_event.set()
                thread.join()
        # # Copy header with fslcopygeom
        # myFslCpGeom = fsl.utils.CopyGeom(dest_file=b0image, in_file=input_file)
        # myFslCpGeom.run()
            input_file = b0image
            print("Averaging b0 volumes completed, output saved to", input_file)
        except Exception as e:
            print(f'Error in averaging b0 volumes\nError message: {str(e)}')
            raise


    try:
        # start spinner
        stop_event = threading.Event()
        thread = threading.Thread(
            target=spinner,
            args=(stop_event, "Running smoothing")
        )
        thread.start()

        try:
            output_smooth = smoothIMG(input_file = input_file, output_path = output_path, skip_min=args.skip_min)
        finally:
            stop_event.set()
            thread.join()
        print(f"Smoothing completed, output saved to {output_path}")
    except Exception as e:
        print(f'Error in smoothing\nError message: {str(e)}')
        raise

    # intensity correction using non parametric bias field correction algorithm
    if bias_method == "none":
        print("No bias field correction applied")
        outputBiasCorr = output_smooth
    elif bias_method == "mico":
        print("Starting Biasfieldcorrection with MICO:")
        try:
            outputBiasCorr = applyMICO.run_MICO(output_smooth, output_path)
            set_xform_codes_to_one(outputBiasCorr)
            set_default_xyzt_units_if_unknown(outputBiasCorr)
            print("Biasfield correction was successful")
        except Exception as e:
            print(f'Error in bias field correction\nError message: {str(e)}')
            raise
    elif bias_method == "ants":
        # intensity correction using ANTs N4BiasFieldCorrection
        print("Starting Biasfieldcorrection with ANTS:")
        try:
            #start spinner
            stop_event = threading.Event()
            thread = threading.Thread(
                target=spinner,
                args=(stop_event, "Running N4 ANTS bias correction")
            )
            thread.start()

            try:
                outputBiasCorr = n4biasfieldcorr(input_file=output_smooth)
            finally:
                stop_event.set()
                thread.join()
            print("Biasfield correction was successful")
        except Exception as e:
            print(f'Error in bias field correction\nError message: {str(e)}')
            raise
    #print(os.path.exists(outputBiasCorr))

    use_bet4animal = args.use_bet4animal

    if args.bet_skip:
        print("Skipping brain extraction.")
        outputBET = skip_bet_function(outputBiasCorr)
    else:
        # brain extraction
        print("Starting brain extraction")
        try:
            stop_event = threading.Event()
            thread = threading.Thread(
                target=spinner,
                args=(stop_event, "Running Brain extraction")
            )
            thread.start()
            try:
                outputBET = applyBET(
                    input_file=outputBiasCorr,
                    frac=frac,
                    radius=radius,
                    horizontal_gradient=horizontal_gradient,
                    use_bet4animal=use_bet4animal,
                    center=args.center)
            finally:
                stop_event.set()
                thread.join()
            print("Brain extraction was successful")
        except Exception as e:
            print(f'Error in brain extraction\nError messsage: {str(e)}')
            raise

    print("Preprocessing completed")

