#!/bin/bash

# Turn output of TTree::Scan into something
# pickAnEvent.pl can process
# The scan output includes a dump of all
# baby branches, not just run, ls and evt

cd /afs/cern.ch/user/j/jribnik/scratch0/whunt
if [ $# -ne 2 ]
then
    echo usage: processCands.sh cands.txt tas03loc
    echo
    echo e.g. processCands.sh wcands.txt /tas03/disk01/whunt
    exit 1
fi

candsfile=$1
candsfiletrunc=`echo $candsfile | awk -F'.' '{print $1}'`
tas03loc=$2

if [ ! -e $candsfile ]
then
    echo $candsfile does not exist, I have no reason to live
    exit 2
fi

fields=`tac $candsfile | head -n-1 | tac | head -n 1 | sed 's/\*//g' | awk '{ORS=" ";for(i=2;i<=NF;i++) print $i;ORS="\n";print ""}' | sed 's/ $//g'`
nf=`echo $fields | awk '{print NF}'`
tac $candsfile | head -n-3 | tac | head -n-1 | sed 's/\*//g' | awk '{ORS=" ";for(i=2;i<=NF;i++) print $i;ORS="\n";print ""}' | sed 's/ $//g' >${candsfiletrunc}_stripped.txt

touch ${candsfiletrunc}_stripped_processed.txt
while read line
do
    grep "$line" ${candsfiletrunc}_stripped_processed.txt >/dev/null 2>&1
    if [ $? -eq 0 ]
    then
        continue
    fi

    run_=`echo $line | awk '{print $1}'`
    ls_=`echo $line | awk '{print $2}'`
    evt_=`echo $line | awk '{print $3}'`
    echo Processing new candidate: $run_ $ls_ $evt_

    # Choose the right dataset!!!
    if [ $run_ -lt 132514 ]
    then
        echo /ExpressPhysics/Commissioning10-Express-v7/FEVT >eventToPick.txt
    else
        echo /ExpressPhysics/Commissioning10-Express-v8/FEVT >eventToPick.txt
    fi

    # Need CMSSW environment for pickAnEvent.pl
    if [ ! -d CMSSW_3_5_6/src ]
    then
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
    for i in `seq 4 $nf`; do
        fieldcmd="echo $fields | awk '{print \$$i}'"
        valuecmd="echo $line | awk '{print \$$i}'"
        field=`eval $fieldcmd`
        value=`eval $valuecmd`
        echo $field $value >>$dumpFile
    done

    echo $line >>${candsfiletrunc}_stripped_processed.txt
done <${candsfiletrunc}_stripped.txt

# Transfer back
scp *.root cms-tas03:${tas03loc}/picks
scp dump*.txt cms-tas03:${tas03loc}/dumps
scp ${candsfiletrunc}_stripped*.txt cms-tas03:${tas03loc}/

# Cleanup
rm *.root *.py dump*.txt
