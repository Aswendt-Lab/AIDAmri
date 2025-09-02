import os
import json
import argparse
from glob import glob

def find_nifti_files(bids_root, participant, session):
    dwi_nifti_files = []
    # Diffusion files
    file_ending = 'dwi'
    #  If denoised data exists, change the file ending
    if os.path.exists(os.path.join(bids_root, f"sub-{participant}", f"ses-{session}", "dwi", f"sub-{participant}_ses-{session}_*_dwi*Denoised.nii.gz")):
        file_ending = 'Denoised'
    dwi_path = os.path.join(bids_root, f"sub-{participant}", f"ses-{session}", "dwi", "*"f"{file_ending}.nii*")
    dwi_nifti_files.extend([os.path.relpath(f, bids_root) for f in glob(dwi_path)])
    func_nifti_files = []
    # Functional files
    file_ending = 'epi'
    func_path = os.path.join(bids_root, f"sub-{participant}", f"ses-{session}", "func", "*"f"{file_ending}.nii*")
    func_nifti_files.extend([os.path.relpath(f, bids_root) for f in glob(func_path)])
    return dwi_nifti_files, func_nifti_files

# def update_fieldmap_jsons(fmap_folder, intended_for_list, bids_root):
#     json_files = glob(os.path.join(fmap_folder, "*.json"))
#     for json_file in json_files:
#         with open(json_file, "r") as f:
#             data = json.load(f)
#         data["IntendedFor"] = intended_for_list
#         with open(json_file, "w") as f:
#             json.dump(data, f, indent=4)
#         print(f"Updated {json_file}")

def update_fieldmap_jsons(fmap_folder, intended_for, overwrite=False):
    for json_file in glob(os.path.join(fmap_folder, "*.json")):
        print(f"Processing {json_file}")
        with open(json_file, "r") as f:
            data = json.load(f)
        if "IntendedFor" in data and data["IntendedFor"] and not overwrite:
            print(f"IntendedFor field is not empty in {json_file}, skipping update.")
            continue
        data["IntendedFor"] = intended_for
        with open(json_file, "w") as f:
            json.dump(data, f, indent=4)
        print(f"Updated {json_file} with IntendedFor field.")

# The following code block was removed because 'args' is only defined inside main().
# Please use the script via the command line or call main() directly.


def main():
    parser = argparse.ArgumentParser(description="Populate IntendedFor in fieldmap JSONs.")
    parser.add_argument("--bids_root", help="Path to BIDS root directory")
    parser.add_argument("--participant", help="Participant label (without 'sub-')")
    parser.add_argument("--session", help="Session label (without 'ses-')")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing IntendedFor fields")
    args = parser.parse_args()

    fmap_folder = os.path.join(args.bids_root, f"sub-{args.participant}", f"ses-{args.session}", "fmap")
    intended_for_dwi, intended_for_func = find_nifti_files(args.bids_root, args.participant, args.session)
    apply_to_all = True
    if apply_to_all is True:
        intended_for = intended_for_dwi + intended_for_func
    update_fieldmap_jsons(fmap_folder, intended_for, overwrite=args.overwrite)

if __name__ == "__main__":
    main()