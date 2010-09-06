#!/bin/bash

while [ 1 ]
do
    # get time stamp
    DATE=`date +%Y%m%d%H%M%S`
    echo $DATE

    echo "Checking for new files to process..."
    find /tas/cms2/MinimumBias_Commissioning10-SD_EG-Jun14thSkim_v1_RECO/V03-04-26-02/singleLepPt10Skim -follow -maxdepth 1 -name "skimmed_ntuple_[0-9]*_[0-9]*.root" -printf "%p %s %C@\n" | awk -F'_' '{if ($6<135446) print $0}' >>AllPromptRunsAvailable.txt
    find /tas/cms2/EG_Run2010A-Jun14thReReco_v1_RECO/V03-04-26-01/singleLepPt5Skim -follow -maxdepth 1 -name "skimmed_ntuple_[0-9]*_[0-9]*.root" -printf "%p %s %C@\n" | awk -F'_' '{(if ($6>135808) && if ($6<137437)) print $0}' >>AllPromptRunsAvailable.txt
    find /tas/cms2/EG_Run2010A-PromptReco-v4_RECO/V03-04-25/singleLepPt5Skim -follow -maxdepth 1 -name "skimmed_ntuple_[0-9]*_[0-9]*.root" -printf "%p %s %C@\n" | awk -F'_' '{((if ($6>137436) && if ($6<139779)) || (if ($6>140160) && if ($6<140401))) print $0}' >>AllPromptRunsAvailable.txt
    find /tas/cms2/EG_Run2010A-PromptReco-v4_RECO/V03-04-26-01/singleLepPt5Skim -follow -maxdepth 1 -name "skimmed_ntuple_[0-9]*_[0-9]*.root" -printf "%p %s %C@\n" | awk -F'_' '{((if ($6>137436) && if ($6<139779)) || (if ($6>140160) && if ($6<140401))) print $0}' >>AllPromptRunsAvailable.txt
    find /tas/cms2/EG_Run2010A-PromptReco-v4_RECO/V03-04-26-02/singleLepPt5Skim -follow -maxdepth 1 -name "skimmed_ntuple_[0-9]*_[0-9]*.root" -printf "%p %s %C@\n" | awk -F'_' '{((if ($6>137436) && if ($6<139779)) || (if ($6>140160) && if ($6<140401))) print $0}' >>AllPromptRunsAvailable.txt
    find /tas/slava77/cms2/EG_Run2010A-Jul16thReReco-v2_RECO/V03-04-26-03 -follow -maxdepth 1 -name "merged_ntuple_[0-9]*_[0-9]*.root" -printf "%p %s %C@\n" | awk -F'_' '{(if ($6>139778) && if ($6<140161)) print $0}' >>AllPromptRunsAvailable.txt

    find /tas/cms2/MinimumBias_Commissioning10-SD_Mu-Jun14thSkim_v1_RECO/V03-04-26-02/singleLepPt10Skim -follow -maxdepth 1 -name "skimmed_ntuple_[0-9]*_[0-9]*.root" -printf "%p %s %C@\n" | awk -F'_' '{if ($6<135446) print $0}' >>AllPromptRunsAvailable.txt
    find /tas/cms2/Mu_Run2010A-Jun14thReReco_v1_RECO/V03-04-26-01/singleLepPt5Ski/ -follow -maxdepth 1 -name "skimmed_ntuple_[0-9]*_[0-9]*.root" -printf "%p %s %C@\n" | awk -F'_' '{(if ($6>135808) && if ($6<137437)) print $0}' >>AllPromptRunsAvailable.txt
    find /tas/cms2/Mu_Run2010A-PromptReco-v4_RECO/V03-04-25/singleLepPt5Skim -follow -maxdepth 1 -name "skimmed_ntuple_[0-9]*_[0-9]*.root" -printf "%p %s %C@\n" | awk -F'_' '{((if ($6>137436) && if ($6<139779)) || (if ($6>140160) && if ($6<140401))) print $0}' >>AllPromptRunsAvailable.txt
    find /tas/cms2/Mu_Run2010A-PromptReco-v4_RECO/V03-04-26-01/singleLepPt5Skim -follow -maxdepth 1 -name "skimmed_ntuple_[0-9]*_[0-9]*.root" -printf "%p %s %C@\n" | awk -F'_' '{((if ($6>137436) && if ($6<139779)) || (if ($6>140160) && if ($6<140401))) print $0}' >>AllPromptRunsAvailable.txt
    find /tas/cms2/Mu_Run2010A-PromptReco-v4_RECO/V03-04-26-02/singleLepPt5Skim -follow -maxdepth 1 -name "skimmed_ntuple_[0-9]*_[0-9]*.root" -printf "%p %s %C@\n" | awk -F'_' '{((if ($6>137436) && if ($6<139779)) || (if ($6>140160) && if ($6<140401))) print $0}' >>AllPromptRunsAvailable.txt
    find /tas/slava77/cms2/Mu_Run2010A-Jul16thReReco-v2_RECO/V03-04-26-03 -follow -maxdepth 1 -name "merged_ntuple_[0-9]*_[0-9]*.root" -printf "%p %s %C@\n" | awk -F'_' '{(if ($6>139778) && if ($6<140161)) print $0}' >>AllPromptRunsAvailable.txt
    # The awk above is used because these two processings overlap
    # and we don't want to deal with the duplication
    # The second processing below is for runs >= 140401
    find /store/disk00/slava77/reltestdata/CMSSW_3_6_1-patch4-cms2-data -follow -maxdepth 1 -name "merged_ntuple_[0-9]*_[0-9]*.root" -printf "%p %s %C@\n" >>AllPromptRunsAvailable.txt
    rm PromptRunsToProcess.txt
    touch PromptRunsProcessed.txt

    # loop over all AllRunsAvailable.txt, check PromptRunsProcessed.txt,
    # if not there, add to PromptRunsToProcess.txt
    while read line
    do
        grep "$line" PromptRunsProcessed.txt >/dev/null 2>&1
        if [ $? -ne 0 ]
        then
            run=`echo $line | awk -F'_' '{print $6}'`
            if [ $run -gt 133222 ] && [ $run -le 133250 ];
            then
                echo "Magnetic field off, skipping run "$run
            else 
                echo $line >>PromptRunsToProcess.txt
            fi
        fi
    done <PromptAllRunsAvailable.txt

    # check if there are new files and process them gods willing
    nfiles=`cat PromptRunsToProcess.txt | wc -l`
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
                echo $line >>PromptRunsProcessed.txt
            fi
        done <PromptRunsToProcess.txt

        # slava's [re-]merging may remove files in which case we
        # will have duplicates so rm parentless babies and skims
        echo "Cleaning out parentless babies and skims"
        # emu
        for baby in `ls prompt_emu_baby/*.root`
        do
            ident=`echo $baby | sed 's!^.*/emuskim_baby_\(.*\).root$!\1!'`
            grep $ident PromptAllRunsAvailable.txt >/dev/null 2>&1
            if [ $? -ne 0 ]
            then
                echo "Found parentless baby, removing $baby"
                rm $baby
            fi
        done
        for skim in `ls prompt_emu_skim/*.root`
        do
            ident=`echo $skim | sed 's!^.*/emuskim_\(.*\).root$!\1!'`
            grep $ident PromptAllRunsAvailable.txt >/dev/null 2>&1
            if [ $? -ne 0 ]
            then
                echo "Found parentless skim, removing $skim"
                rm $skim
            fi
        done
        # dilep
        for baby in `ls prompt_dilep_baby/*.root`
        do
            ident=`echo $baby | sed 's!^.*/dilepskim_baby_\(.*\).root$!\1!'`
            grep $ident PromptAllRunsAvailable.txt >/dev/null 2>&1
            if [ $? -ne 0 ]
            then
                echo "Found parentless baby, removing $baby"
                rm $baby
            fi
        done
        for skim in `ls prompt_dilep_skim/*.root`
        do
            ident=`echo $skim | sed 's!^.*/dilepskim_\(.*\).root$!\1!'`
            grep $ident PromptAllRunsAvailable.txt >/dev/null 2>&1
            if [ $? -ne 0 ]
            then
                echo "Found parentless skim, removing $skim"
                rm $skim
            fi
        done
        # trilep
        for baby in `ls prompt_trilep_baby/*.root`
        do
            ident=`echo $baby | sed 's!^.*/trilepskim_baby_\(.*\).root$!\1!'`
            grep $ident PromptAllRunsAvailable.txt >/dev/null 2>&1
            if [ $? -ne 0 ]
            then
                echo "Found parentless baby, removing $baby"
                rm $baby
            fi
        done
        for skim in `ls prompt_trilep_skim/*.root`
        do
            ident=`echo $skim | sed 's!^.*/trilepskim_\(.*\).root$!\1!'`
            grep $ident PromptAllRunsAvailable.txt >/dev/null 2>&1
            if [ $? -ne 0 ]
            then
                echo "Found parentless skim, removing $skim"
                rm $skim
            fi
        done
        # clean logs as well
        for log in `ls prompt_logs/*.txt`
        do
            ident=`echo $log | sed 's!^.*/log_\(.*\).txt$!\1!'`
            grep $ident PromptAllRunsAvailable.txt >/dev/null 2>&1
            if [ $? -ne 0 ]
            then
                echo "Found parentless log, removing $log"
                rm $log
            fi
        done

        # Produce cands.txts and do some basic plotting
        # in a separate process so that if it stalls it
        # it does not affect this process
        #./makeCands.sh&

        # Last run processed
        #cat PromptRunsProcessed.txt | cut -d '_' -f 6 | sort -n -r | head -n 1 >lastrun.txt
        #scp lastrun.txt uaf-4.t2.ucsd.edu:~/public_html/hunt
    fi

    sleep 600;
done
