import os
import subprocess
import argparse


DATASET_PATH = "/mnt/data2/2024_Grandjean_Multiverse"  # Hardcoded dataset path

def get_subfolders(input_path):
    """Generate a list of all folders starting with 'sub-' in the given input path."""
    return sorted([
        os.path.join(input_path, f)
        for f in os.listdir(input_path)
        if os.path.isdir(os.path.join(input_path, f)) and f.startswith("sub-")
    ])


def execute_task(folder_path, script_path):
    """Execute the specified tasks on the given folder."""
    try:
        # Switch to the datalad dataset path
        os.chdir(DATASET_PATH)

        # Step 2: Datalad get and unlock the folder
        subprocess.run(["datalad", "get", folder_path], check=True)
        subprocess.run(["datalad", "unlock", folder_path], check=True)

        # Step 3: Execute the Python script
        subprocess.run(["python", script_path, "-i", folder_path], check=True)

        # Step 4: Save the folder
        subprocess.run(["datalad", "save", "-m", f"add the sub-folder {os.path.basename(folder_path)}", folder_path], check=True)

        # Step 5: Push changes to the remote (e.g., gin)
        subprocess.run(["datalad", "push", "--to", "gin"], check=True)

        # Step 6: Drop the folder
        subprocess.run(["datalad", "drop", folder_path], check=True)
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while processing {folder_path}: {e}")
    finally:
        # Return to the original working directory to avoid issues
        os.chdir(os.path.dirname(os.path.abspath(__file__)))


def main(input_path, script_path):
    while True:
        subfolders = get_subfolders(input_path)
        if not subfolders:
            print("No more folders to process.")
            break

        # Process the first folder
        first_folder = subfolders.pop(0)
        print(f"Processing folder: {first_folder}")
        execute_task(first_folder, script_path)

        if not subfolders:
            print("All folders have been processed.")
            break


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process folders starting with 'sub-'.")
    parser.add_argument("-i", "--input", required=True, help="Path to the input directory containing 'sub-' folders.")
    parser.add_argument("-s", "--script", required=True, help="Path to the Python script to execute.")

    args = parser.parse_args()

    main(args.input, args.script)
