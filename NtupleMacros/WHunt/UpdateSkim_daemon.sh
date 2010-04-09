#!/bin/bash

while [ 1 ]; do
    # get time stamp
    DATE=`date +%Y%m%d%H%M%S`
    echo $DATE

		echo "Starting to check for new files to process..."
		ls /store/disk03/slava77/reltestdata/CMSSW_3_5_6-cms2-data/*root > AllRunsAvailable.txt
		touch RunsProcessed.txt
		diff RunsProcessed.txt AllRunsAvailable.txt | grep ">" | awk '{print $2}' > RunsToProcess.txt
# check how many new files are there, 
# only run if there are actually new ones
		nFiles=`cat RunsToProcess.txt | wc -l`
		if [ $nFiles -lt 1 ];
				then
				echo "No new runs to process - stopping..."
		else
				echo "prepared list to process - starting to produce ntuple based on that list..."
				root -b -q UpdateSkim.C
		fi
    sleep 600
done


