Here you will find all the tools to create the AIDAmri parental atlas. Therefore two files are necessary: 

1. The original annotations of the Allen Brain Institute. These can be downloaded using AllenSDK using the download_ARA.py tool. This program can be executed without input values. The files annotation.nii.gz and template.nii.gz are generated automatically.

2. An xlsx file according to annotation_label_IDs_valid.xlsx. This file was created manually using AllenBrainAPI-master. The file has only been added for the sake of completeness and is not part of its own development. The implementation comes from https://github.com/SainsburyWellcomeCentre/AllenBrainAPI. Here are more details on how to use AllenBrianAPI-master.

The parental atlas is created using both files. All regions in annotation.nii.gz are merged according to the table annotation_label_IDs_valid.xlsx and saved as annotation_parent.nii.gz. Enter the following command in the command window of Matlab
getParentalARA ('./ annotation_label_IDs.xlsx', '. / annotation / annotation.nii.gz')

