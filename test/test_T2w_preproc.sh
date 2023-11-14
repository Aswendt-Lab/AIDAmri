#!/bin/bash

readarray -t arr < contents_T2_preproc.txt
bindir="/aida/bin/2.1_T2PreProcessing"

testrun_start=`date +%Y-%m-%d_%H-%M`
echo "Pre-processing T2w test at ${testrun_start}"
cd /testData
logname=$(pwd)/${testrun_start}_testrun_T2_preProc.log
touch ${logname}
cd /testData/T2w
echo -e "Date: $testrun_start\nPresent files: \n" &>> ${logname}
ls -lh &>> ${logname}
echo -e "\n\n\nPreProcessing starts\n" &>> ${logname}

python $bindir/preProcessing_T2.py -i testData.5.1.nii.gz &>> ${logname}
ls -lh
for f in *; do
    if [[ ! " $arr " =~ " ${f}\n " ]]; then
        echo -e "ERROR: $f NOT FOUND; EXIT 1\n" >> ${logname}
        exit 1
    fi
    echo -e "Contains: $f\n" >> ${logname}
done
testrun_end=`date +%Y-%m-%d_%H-%M`
echo "\n\n\nPreProcessing ends at ${testrun_end}\n" >> ${logname}