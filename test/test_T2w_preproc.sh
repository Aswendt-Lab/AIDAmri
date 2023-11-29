#!/bin/bash

readarray -t ctable < output_contents.txt
while IFS= read -r line; do
    if [[ $line = T2preprocessing* ]]; then
        arr=($line)
    fi
done <<< $ctable
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

echo -e "\nCatching stream output:\n\n" &>> ${logname}
python $bindir/preProcessing_T2.py -i testData.5.1.nii.gz &>> ${logname}
count=`ls -1 *.log 2>/dev/null | wc -l`
if [ $count != 0 ]
then 
    echo -e "\n\nLog file generated. Check data folder output data directory.\n\n" &>> ${logname}
fi 
for i in "${arr[@]:1}"; do
    if ! [ -f ${i} ]; then
        echo -e "ERROR: $i NOT FOUND; EXIT 1\n" &>> ${logname}
        exit 1
    fi
    echo -e "Contains: $i\n" &>> ${logname}
done

testrun_end=`date +%Y-%m-%d_%H-%M`
echo -e "\n\n\nPreProcessing ends at ${testrun_end}\n" &>> ${logname}