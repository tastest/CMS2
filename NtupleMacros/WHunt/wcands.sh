#!/bin/bash

# Turn output of TTree::Scan into something
# pickAnEvent.pl can process
# The scan output includes a dump of all
# baby branches, not just run, ls and evt

cd /afs/cern.ch/user/j/jribnik/scratch0/whunt
if [ ! -e wcands.txt ]; then
    echo wcands.txt does not exist, I have no reason to live
    exit 1
fi

tac wcands.txt | head -n-3 | tac | head -n-1 | sed 's/\*//g' | awk '{ORS=" ";for(i=2;i<=NF;i++) print $i;ORS="\n";print ""}' | sed 's/ $/ /g' >wcands_stripped.txt
touch wcands_stripped_processed.txt
while read line; do
    grep "$line" wcands_stripped_processed.txt >/dev/null 2>&1
    if [ $? -eq 0 ]; then
        continue
    fi

    echo Processing new candidate: $line

    # Choose the right dataset!!!
    run_=`echo $line | awk '{print $1}'`
    ls_=`echo $line | awk '{print $2}'`
    evt_=`echo $line | awk '{print $3}'`
    if [ $run_ -lt 132514 ]; then
        echo /ExpressPhysics/Commissioning10-Express-v7/FEVT >eventToPick.txt
    else
        echo /ExpressPhysics/Commissioning10-Express-v8/FEVT >eventToPick.txt
    fi

    # Need CMSSW environment for pickAnEvent.pl
    if [ ! -d CMSSW_3_5_6/src ]; then
        scramv1 p CMSSW CMSSW_3_5_6
    fi
    cd CMSSW_3_5_6/src
    eval `scramv1 ru -sh`
    cd -

    # Finally ready for pickAnEvent.pl
    echo "$run_ $ls_ $evt_" >>eventToPick.txt
    ./pickAnEvent.pl eventToPick.txt

    # Dump the event
    dumpFile="dump_${run_}_${ls_}_${evt_}.txt"
    echo pfmet `echo $line | awk '{print $4}'` >>$dumpFile
    echo njets `echo $line | awk '{print $5}'` >>$dumpFile
    echo jet1pt `echo $line | awk '{print $6}'` >>$dumpFile
    echo dphimetjet `echo $line | awk '{print $7}'` >>$dumpFile
    echo eormu `echo $line | awk '{print $8}'` >>$dumpFile
    echo type `echo $line | awk '{print $9}'` >>$dumpFile
    echo pt `echo $line | awk '{print $10}'` >>$dumpFile
    echo iso `echo $line | awk '{print $11}'` >>$dumpFile
    echo d0corr `echo $line | awk '{print $12}'` >>$dumpFile
    echo dphimet `echo $line | awk '{print $13}'` >>$dumpFile
    echo drjet `echo $line | awk '{print $14}'` >>$dumpFile
    echo mt `echo $line | awk '{print $15}'` >>$dumpFile
    echo mu_muonid `echo $line | awk '{print $16}'` >>$dumpFile
    echo mu_goodmask `echo $line | awk '{print $17}'` >>$dumpFile
    echo mu_gfitchi2 `echo $line | awk '{print $18}'` >>$dumpFile
    echo e_cand01 `echo $line | awk '{print $19}'` >>$dumpFile
    echo e_eopin `echo $line | awk '{print $20}'` >>$dumpFile
    echo e_hoe `echo $line | awk '{print $21}'` >>$dumpFile
    echo e_dphiin `echo $line | awk '{print $22}'` >>$dumpFile
    echo e_detain `echo $line | awk '{print $23}'` >>$dumpFile
    echo e_eMe55 `echo $line | awk '{print $24}'` >>$dumpFile
    echo e_nmHits `echo $line | awk '{print $25}'` >>$dumpFile

    echo $line >>wcands_stripped_processed.txt
done <wcands_stripped.txt

# Transfer back to cms-tas03
scp *.root cms-tas03:/tas03/disk01/whunt/picks
scp wcands_stripped_processed.txt cms-tas03:/tas03/disk01/whunt/
scp dump*.txt cms-tas03:/tas03/disk01/whunt/dumps

# Cleanup
rm *.root *.py dump*.txt
