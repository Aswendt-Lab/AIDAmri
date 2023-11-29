#!bin/bash

test_routine () {
    : '
    Arguments: Tool names to test (without .py ending).
    Arr[0]: Tool name
    Arr[1]: Binary location
    Arr[2]: Test data directory
    Arr[3]: Input file
    Arr[@]:4: Contents after processing (excluding input file)
    '
    readarray -t ctable < output_contents.txt
    for arg in "$@"; do
        while IFS="" read -r line || [ -n "$line" ]; do
            if [[ "$line" =~ ^"$arg" ]]; then
                arr=($line)
            fi
        done < output_contents.txt
        bindir=${arr[1]}

        testrun_start=`date +%Y-%m-%d_%H-%M`
        echo "${arg} testing starting at ${testrun_start}"
        cd /testData
        logname=$(pwd)/${testrun_start}_testrun_${arg}.log
        touch ${logname}
        cd ${arr[2]}
        echo -e "Date: $testrun_start\nPresent files: \n" &>> ${logname}
        ls -lh &>> ${logname}
        echo -e "\nCatching stream output:\n\n" &>> ${logname}
        python $bindir/$arg.py -i ${arr[3]} &>> ${logname}
        count=`ls -1 *.log 2>/dev/null | wc -l`
        if [ $count != 0 ]
        then 
            echo -e "\n\nLog file generated. Check data folder output data directory.\n\n" &>> ${logname}
        fi 
        for i in "${arr[@]:4}"; do
            if ! [ -f ${i} ]; then
                echo -e "ERROR: $i NOT FOUND; EXIT 1\n" &>> ${logname}
                exit 1
            fi
            echo -e "Contains: $i\n" &>> ${logname}
        done

        testrun_end=`date +%Y-%m-%d_%H-%M`
        echo -e "\n\n\nProcess ended at ${testrun_end}\n" &>> ${logname}
    done
}
test_routine preProcessing_T2