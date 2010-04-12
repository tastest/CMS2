#!/bin/bash
#!/bin/bash

while [ 1 ]; do
    # get time stamp
    DATE=`date +%Y%m%d%H%M%S`
    echo $DATE

    echo "Starting to check for new files to process..."
    find /tas03/disk03/slava77/reltestdata/CMSSW_3_5_6-cms2-data -maxdepth 1 -name "merged_ntuple_[0-9]*_[0-9]*.root" -printf "%p %s %C@\n" > AllRunsAvailable.txt
    touch RunsProcessed.txt
    diff RunsProcessed.txt AllRunsAvailable.txt | grep ">" | awk '{print substr($0,3,length($0))}' > RunsToProcess.txt
    # check how many new files are there, 
    # only run if there are actually new ones
    nFiles=`cat RunsToProcess.txt | wc -l`
    if [ $nFiles -lt 1 ];
    then
     echo "No new runs to process - stopping..."
    else
        while read line; do
            file=`echo $line | awk '{print $1}'`
            cmd="root -b -l -q 'UpdateSkim.C+(\"$file\")'"
            eval $cmd
            if [ $? -ne 0 ]; then
                echo "Error processing $file\n"
            else
                echo $line >>RunsProcessed.txt
            fi
        done <RunsToProcess.txt
    fi
sleep 600;
done
