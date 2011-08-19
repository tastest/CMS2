#!/bin/bash

while [ 1 ]
do
    # get time stamp
    DATE=`date +%Y%m%d%H%M%S`
    echo $DATE

    echo "Checking for new files to process..."
    find /store/disk00/slava77/reltestdata/CMSSW_3_6_1-cms2-data -follow -maxdepth 1 -name "merged_ntuple_[0-9]*_[0-9]*.root" -printf "%p %s %C@\n" | awk -F'_' '{if ($6<138728) print $0}' >AllRunsAvailable.txt
    # The awk above is used because these two processings overlap
    # and we don't want to deal with the duplication
    # The second processing below is for runs >= 138728
    find /store/disk00/slava77/reltestdata/CMSSW_3_6_1-patch4-cms2-data -follow -maxdepth 1 -name "merged_ntuple_[0-9]*_[0-9]*.root" -printf "%p %s %C@\n" >>AllRunsAvailable.txt
    rm RunsToProcess.txt
    touch RunsProcessed.txt

    # loop over all AllRunsAvailable.txt, check RunsProcessed.txt,
    # if not there, add to RunsToProcess.txt
    while read line
    do
        grep "$line" RunsProcessed.txt >/dev/null 2>&1
        if [ $? -ne 0 ]
        then
            run=`echo $line | awk -F'_' '{print $6}'`
            if [ $run -gt 133222 ] && [ $run -le 133250 ];
            then
                echo "Magnetic field off, skipping run "$run
            else 
                echo $line >>RunsToProcess.txt
            fi
        fi
    done <AllRunsAvailable.txt

    # check if there are new files and process them gods willing
    nfiles=`cat RunsToProcess.txt | wc -l`
    if [ $nfiles -lt 1 ]
    then
        echo "No new runs to process..."
    else
        while read line
        do
            file=`echo $line | awk '{print $1}'`
            cmd="root -b -l -q 'UpdateSkims.C+(\"$file\")'"
            eval $cmd
            if [ $? -ne 0 ]
            then
                echo "Error processing $file"
            else
                echo $line >>RunsProcessed.txt
            fi
        done <RunsToProcess.txt

        # slava's [re-]merging may remove files in which case we
        # will have duplicates so rm parentless babies and skims
        echo "Cleaning out parentless babies and skims"
        # emu
        for baby in `ls emu_baby/*.root`
        do
            ident=`echo $baby | sed 's!^.*/emuskim_baby_\(.*\).root$!\1!'`
            grep $ident AllRunsAvailable.txt >/dev/null 2>&1
            if [ $? -ne 0 ]
            then
                echo "Found parentless baby, removing $baby"
                rm $baby
            fi
        done
        for skim in `ls emu_skim/*.root`
        do
            ident=`echo $skim | sed 's!^.*/emuskim_\(.*\).root$!\1!'`
            grep $ident AllRunsAvailable.txt >/dev/null 2>&1
            if [ $? -ne 0 ]
            then
                echo "Found parentless skim, removing $skim"
                rm $skim
            fi
        done
        # dilep
        for baby in `ls dilep_baby/*.root`
        do
            ident=`echo $baby | sed 's!^.*/dilepskim_baby_\(.*\).root$!\1!'`
            grep $ident AllRunsAvailable.txt >/dev/null 2>&1
            if [ $? -ne 0 ]
            then
                echo "Found parentless baby, removing $baby"
                rm $baby
            fi
        done
        for skim in `ls dilep_skim/*.root`
        do
            ident=`echo $skim | sed 's!^.*/dilepskim_\(.*\).root$!\1!'`
            grep $ident AllRunsAvailable.txt >/dev/null 2>&1
            if [ $? -ne 0 ]
            then
                echo "Found parentless skim, removing $skim"
                rm $skim
            fi
        done
        # trilep
        for baby in `ls trilep_baby/*.root`
        do
            ident=`echo $baby | sed 's!^.*/trilepskim_baby_\(.*\).root$!\1!'`
            grep $ident AllRunsAvailable.txt >/dev/null 2>&1
            if [ $? -ne 0 ]
            then
                echo "Found parentless baby, removing $baby"
                rm $baby
            fi
        done
        for skim in `ls trilep_skim/*.root`
        do
            ident=`echo $skim | sed 's!^.*/trilepskim_\(.*\).root$!\1!'`
            grep $ident AllRunsAvailable.txt >/dev/null 2>&1
            if [ $? -ne 0 ]
            then
                echo "Found parentless skim, removing $skim"
                rm $skim
            fi
        done
        # clean logs as well
        for log in `ls logs/*.txt`
        do
            ident=`echo $log | sed 's!^.*/log_\(.*\).txt$!\1!'`
            grep $ident AllRunsAvailable.txt >/dev/null 2>&1
            if [ $? -ne 0 ]
            then
                echo "Found parentless log, removing $log"
                rm $log
            fi
        done

        # Produce cands.txts and do some basic plotting
        # in a separate process so that if it stalls it
        # it does not affect this process
        ./makeCands.sh&

        # Last run processed
        cat RunsProcessed.txt | cut -d '_' -f 6 | sort -n -r | head -n 1 >lastrun.txt
        scp lastrun.txt uaf-4.t2.ucsd.edu:~/public_html/hunt
    fi

    sleep 600;
done
