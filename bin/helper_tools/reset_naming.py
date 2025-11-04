import os
import argparse
import glob

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='This script prepares raw data before running 1_PV2NIfTiConverter/pv_conv2Nifti.py. The raw data must follow the structure: projectfolder/days/subjects/data/. A customized groupMapping.csv is required, specifying the group name for each subject folder. Optional arguments allow selecting specific days (-d) and groups (-g). A new folder called proc_data will be created in the same directory as the raw_data folder.Example:python conv2Nifti_auto.py -i /Volumes/Desktop/MRI/raw_data -d Baseline P1 P7 P14 P28')
    parser.add_argument('-i', '--input', required=True,
                        help='Path to the parent project folder of the dataset, e.g. raw_data', type=str) 
    
    # read out parameters
    args = parser.parse_args()

    # get list of raw data in input folder
    list_of_raw = sorted([d for d in os.listdir(args.input) if os.path.isdir(os.path.join(args.input, d)) \
                               or (os.path.isfile(os.path.join(args.input, d)) and (('zip' in d) or ('PvDataset' in d)))])
    
    subject_files = glob.glob(os.path.join(args.input, "**", "subject"), recursive=True)
    print(subject_files)
    print(list_of_raw)

    for path in subject_files:
        subject_file = path
        subject_id = "##$SUBJECT_id="
        session_id = "##$SUBJECT_study_name="
        
        if os.path.exists(subject_file):
            with open(subject_file, 'r') as infile:
                lines = infile.readlines()
            
            for idx, line in enumerate(lines):
                if subject_id in line or session_id in line:
                    # Replace the first underscore with '' in the next line
                    lines[idx + 1] = lines[idx + 1].replace("_", "", 1)

            with open(subject_file, 'w') as outfile:
                outfile.writelines(lines)
                print(f"Modified: {subject_file}")
                print('success')

                        
                       
                

                        
