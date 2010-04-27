#!/bin/bash

echo "Starting to check for new files to process..."
find /tas03/disk03/slava77/reltestdata/CMSSW_3_5_6-cms2-data -follow -maxdepth 1 -name "merged_ntuple_[0-9]*_[0-9]*.root" -printf "%p %s %C@\n" > AllRunsAvailable.txt
find /tas03/disk03/slava77/reltestdata/CMSSW_3_5_7-cms2-data -maxdepth 1 -name "merged_ntuple_[0-9]*_[0-9]*.root" -printf "%p %s %C@\n" >> AllRunsAvailable.txt
rm RunsToProcess.txt
touch RunsProcessed.txt

# loop over all runs available, check to find file name in RunsProcessed.txt,
# if not there, add to RunsToProcess.txt
while read line; do
    wasprocessed=`cat RunsProcessed.txt | grep "$line"`
    if [ `echo -n $wasprocessed | wc -c` -eq 0 ];
    then
      runnumber=`echo $line | awk -F "_" '{print $6}'`
  if [ $runnumber -gt 133222 ] && [ $runnumber -le 133250 ];
  then
    echo "Run with magnetic field off, skipping run: "$runnumber
  else 
        echo $line >> RunsToProcess.txt
  fi
    fi
done <AllRunsAvailable.txt

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

    # slava's merging and remerging may remove files in which case
    # we will have duplicates so rm parentless babies and skims
    echo "cleaning out parentless babies and skims"
    for baby in `ls baby/*.root`; do
        ident=`echo $baby | sed 's!^.*/dilepskim_baby_\(.*\).root.*$!\1!'`
        grep $ident AllRunsAvailable.txt >/dev/null 2>&1
        if [ $? -ne 0 ]; then
            echo "found parentless baby, removing $baby"
            rm $baby
        fi
    done
    for skim in `ls skim/*.root`; do
        ident=`echo $skim | sed 's!^.*/dilepskim_\(.*\).root.*$!\1!'`
        grep $ident AllRunsAvailable.txt >/dev/null 2>&1
        if [ $? -ne 0 ]; then
            echo "found parentless skim, removing $skim"
            rm $skim
        fi
    done

    # Run basic plotting and scanning script
    # in a separate process so that if it
    # stalls it does not affect this process
    #./makePlots.sh&
fi
