1. Use this multiverse branch only. 
2. Download https://drive.google.com/drive/folders/1-5WheC8y0SxnamYIcaqo6noVtZOu_Cv0?usp=sharing and save the folder in the bin folder. 
3. In the provided dataset https://github.com/grandjeanlab/multiverse there are three different categories of image orientations, which needed correction for proper AIDAmri processing and generating Multiverse-specific output. 

4. Special steps needed for orient_3 files: 
anat
4.1.1. Create a corrected anat file first (once): correctAnatHeader.py 
4.1.2. Apply correct header information to sub-300301_T2w.nii.gz (the correct T2w file which will serve as a template for the orientation) fslreorient2std corrected_T2w.nii.gz sub-300301_T2w.nii.gz 
4.1.3. Apply from the template T2 file the corrected orientation to the actual T2w file: fslreorient2std corrected_orientation_sub-300301_T2w.nii.gz /Volumes/Publications/2024_Grandjean_Multiverse/input/mri/proc_data/version2/orient_3/sub-300302/ses-1/anat/sub-300302_T2w.nii.gz
4.1.4. Run preprocessing without bias field correction
4.1.5. Apply correct orientation to resulting be file: fslreorient2std sub-300301_T2wBet.nii.gz /Volumes/Publications/2024_Grandjean_Multiverse/input/mri/proc_data/version2/orient_3/sub-300302/ses-1/anat/sub-300302_T2wBet.nii.gz
4.1.6. Run registration python ../2.1_T2PreProcessing/registration_T2.py -i /Volumes/Publications/2024_Grandjean_Multiverse/input/mri/proc_data/version2/orient_3/sub-300302/ses-1/anat/sub-300302_T2wBet.nii.gz

func

4.2.1. Create a corrected function file first (once): 
4.2.2. Apply corrected header information to func file:  fslreorient2std corrected_task-rest_bold.nii.gz /Volumes/Publications/2024_Grandjean_Multiverse/input/mri/proc_data/version2/orient_3/sub-300302/ses-1/func/sub-300302_task-rest_bold.nii.gz   
4.2.3. Run preprocessing: python ../2.3_fMRIPreProcessing/preProcessing_fMRI.py -i /Volumes/Publications/2024_Grandjean_Multiverse/input/mri/proc_data/version2/orient_3/sub-300302/ses-1/func/sub-300302_task-rest_bold.nii.gz 
4.2.4. Apply correct orientation to resulting bet file:  fslreorient2std sub-300301_task-rest_boldSmoothBet.nii.gz /Volumes/Publications/2024_Grandjean_Multiverse/input/mri/proc_data/version2/orient_3/sub-300302/ses-1/func/sub-300302_task-rest_bold.nii.gz 
4.2.5. Run registration: python ../2.3_fMRIPreProcessing/registration_rsfMRI.py -i /Volumes/Publications/2024_Grandjean_Multiverse/input/mri/proc_data/version2/orient_3/sub-300302/ses-1/func/sub-300302_task-rest_boldSmoothBet.nii.gz
4.2.6. Run fmriproc: python ../3.3_fMRIActivity/process_fMRI.py -i /Volumes/Publications/2024_Grandjean_Multiverse/input/mri/proc_data/version2/orient_3/sub-300302/ses-1/func/sub-300302_task-rest_bold.nii.gz 
