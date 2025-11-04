import os
import argparse
import glob
import pandas as pd

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='This script processes NIfTI files in a directory. It extracts relevant parts of the file name and creates a DataFrame.')
    parser.add_argument('-i', '--input', required=True, help='Path to the parent project folder of the dataset, e.g., raw_data', type=str)
    parser.add_argument('-o', '--output', required=True, help='Path where the output CSV file should be saved', type=str)
   
    args = parser.parse_args()

    # Search for all subject files in the input path
    temp_files = glob.glob(os.path.join(args.input, "**","brkraw", "*.nii.gz"), recursive=True)
    
    data = []

    for tt in temp_files:
        filename = os.path.basename(tt)
        list_split = filename.split("_")
        # Initialize the dictionary to collect file details
        file_info = {
            "FileAddress": tt,
            "Modality": None,
            "TimePoint": None,
            "SubjectID": None,
            "RunNumber": None
        }

        # Parse the split filename for specific identifiers
        for element in list_split:
            if "sub-" in element:
                file_info['SubjectID'] = element.replace("sub-", "")
            elif "ses-" in element:
                file_info['TimePoint'] = element.replace("ses-", "")
            elif "run-" in element:
                file_info['RunNumber'] = element.replace("run-", "")
            elif ".nii.gz" in element:
                file_info['Modality'] = element.replace(".nii.gz", "")

        data.append(file_info)
        
    # Create DataFrame
    df = pd.DataFrame(data)

    # Remove duplicate rows
    df = df.drop_duplicates()

    # Save the DataFrame to a CSV file
    output_file_path = os.path.join(args.output, "MRI_files_overview.csv")
    df.to_csv(output_file_path, index=False)
    
    print("Data processing complete. The DataFrame has been saved to", output_file_path)
