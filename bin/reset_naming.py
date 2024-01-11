import os
import argparse
import glob

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='This script automates the conversion from the raw bruker data format to the NIfTI format using 1_PV2NIfTiConverter/pv_conv2Nifti.py. The raw data needs to be in the following structure: projectfolder/days/subjects/data/. For this script to work, the groupMapping.csv needs to be adjusted, where the group name of every subject''s folder in the raw data structure needs to be specified. This script computes the converison either for all data in the raw project folder or for certain days and/or groups specified through the optional arguments -d and -g. During the processing a new folder called proc_data is being created in the same directory where the raw data folder is located. Example: python conv2Nifti_auto.py -f /Volumes/Desktop/MRI/raw_data -d Baseline P1 P7 P14 P28')
    parser.add_argument('-i', '--input', required=True,
                        help='Path to the parent project folder of the dataset, e.g. raw_data', type=str) 
						
	## read out parameters
    args = parser.parse_args()
    
    # get list of raw data in input folder
    list_of_raw = sorted([d for d in os.listdir(args.input) if os.path.isdir(os.path.join(args.input, d)) \
                               or (os.path.isfile(os.path.join(args.input, d)) and (('zip' in d) or ('PvDataset' in d)))])
    subject_files = glob.glob(os.path.join(args.input,"**","subject"),recursive=True)
    print(subject_files)
    print(list_of_raw)
    for path in subject_files:
    
        subject_file = path
        subject_id = "##$SUBJECT_id="
        session_id = "##$SUBJECT_study_name="
        if os.path.exists(subject_file):
            lines = []
            with open(subject_file, 'r') as infile:
                lines = infile.readlines() 
                for idx, line in enumerate(lines):
                    if subject_id in line:
                        lines[idx+1] = lines[idx+1].replace("_", "s", 1).replace("_", "c", 1).replace("_", "m", 1)
                    if session_id in line:
                        lines[idx+1] = lines[idx+1].replace("_","")
                
                with open(subject_file, 'w') as outfile:
                    outfile.writelines(lines)
                    print("success")
                        
                       
                

                        
