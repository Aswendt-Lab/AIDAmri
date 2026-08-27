import os
import sys
import glob
import numpy as np
import nibabel as nib

def hide_files(directory_use, dwifile, dwibasename, suffix):
    hidden_dir = os.path.join(directory_use, suffix)
    if not os.path.exists(hidden_dir):
        os.makedirs(hidden_dir)
    os.rename(dwifile, os.path.join(hidden_dir, dwifile))
    dwibvalname = dwibasename + ".bval"
    dwibvecname = dwibasename + ".bvec"
    if os.path.exists(dwibvalname):
        os.rename(dwibvalname, os.path.join(hidden_dir, dwibvalname))
    if os.path.exists(dwibvecname):
        os.rename(dwibvecname, os.path.join(hidden_dir, dwibvecname))
    json_name = dwibasename + ".json"
    if os.path.exists(json_name):
        os.rename(json_name, os.path.join(hidden_dir, json_name))

def adjust_bvec_rep(directory_use):
    if os.path.isdir(directory_use):
        print(f"Processing directory: {directory_use}")
    else:
        print(f"Error: {directory_use} is not a valid directory.")
        return

    os.chdir(directory_use)
    os.chdir("dwi")

    dwi_files = glob.glob('*dwi.nii.gz')
    dwi_file_length = []

    for dwifile in dwi_files:
        print(f"Processing File: {dwifile}")
        img = nib.load(dwifile)
        image_shape = img.shape
        dwibasename, _ = os.path.splitext(dwifile)
        dwibasename, _ = os.path.splitext(dwibasename)

        if len(image_shape) < 4 or image_shape[3] == 1:
            print(f"Skipping {dwifile} as it only has 1 volume.")
            hide_files(directory_use, dwifile, dwibasename, '.single_volume')
            continue
        elif image_shape[3] < 12:
            print(f"Skipping {dwifile} as it has less than 12 volumes.")
            hide_files(directory_use, dwifile, dwibasename, '.low_volume_count')
            continue
        elif image_shape[3] == 23:
            print(f"Skipping {dwifile} as it has 23 volumes.")
            hide_files(directory_use, dwifile, dwibasename, '.volumes_23_count')
            continue

        dwibvalname = dwibasename + ".bval"
        print(f"DWI bvalue name: {dwibvalname}")
        dwibvecname = dwibasename + ".bvec"
        print(f"DWI bvec name: {dwibvecname}")

        bval = np.loadtxt(dwibvalname, dtype=float)
        bvec = np.loadtxt(dwibvecname, dtype=float)

        bvaltile_fact = image_shape[3] / len(bval)
        print(f"Scaling up bval by a factor of {bvaltile_fact} repetitions.")

        bvectile_fact = image_shape[3] / (bvec.shape[1])
        print(f"Scaling up bvec by a factor of {bvectile_fact} repetitions.")

        if bvaltile_fact == int(bvaltile_fact):
            bvaltile_fact = int(bvaltile_fact)
        else:
            print("length of bval not compatible with length of diffusion data.")
            continue

        if bvectile_fact == int(bvectile_fact):
            bvectile_fact = int(bvectile_fact)
        else:
            print("length of bvec not compatible with length of diffusion data.")
            continue

        bval = np.tile(bval, bvaltile_fact)
        bvec = np.tile(bvec, bvectile_fact)
        bval = bval.reshape(1, -1)
        bval2name = dwibvalname
        bvec2name = dwibvecname

        np.savetxt(bval2name, bval, fmt="%1.15f", delimiter=" ")
        np.savetxt(bvec2name, bvec, fmt="%1.17f", delimiter=" ")
        print(f"Saved bval to {bval2name} and bvec to {bvec2name}.")
        dwi_file_length.append({'filename': dwifile, 'image_shape': image_shape})

    print("DWI files processed:")
    for item in dwi_file_length:
        print(f" - {item['filename']} with shape {item['image_shape']}")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        adjust_bvec_rep(sys.argv[1])
    else:
        print("Usage: python adjustbvecRep.py /path/to/directory")
