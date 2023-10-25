#!/bin/bash
if [[ $# -eq 0 ]] ; then
    echo "Please enter a testfile"
    exit 0
fi

if [[ $# -ne 1 ]] ; then
    echo "Too many arguments"
    exit 0
fi

echo "You entered the file: $1"
if [[ $1 != *.nii ]] && [[ $1 != *.nii.gz ]] ; then
    echo "Wrong data type, please use Nifti (.nii) or packaged Nifti (.nii.gz) format."
    exit 0
fi

testrun_start=`date +%Y-%m-%d_%H-%M`
echo $testrun_start
logname=${testrun_start}_testrun_T2_preProc.log
ls -lh > $logname
echo "\n\n\nPreProcessing starts\n" >> $logname
python ../bin/2.1_T2PreProcessing/preProcessing_T2.py -i $1 >> $logname
testrun_end=`date +%Y-%m-%d_%H-%M`
echo "\n\n\nPreProcessing ends at ${testrun_end}\n" >> $logname
ls -lh >> $logname
