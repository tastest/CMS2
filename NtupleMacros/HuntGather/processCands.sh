#!/bin/bash

# Turn output of TTree::Scan into something useful,
# i.e. a dump, and trigger pickAnEvent.pl
# E-mail whomever about existence of new candidate

#fireshotdir=/store/disk00/jribnik/fireShot
fireshotdir=/nfs-3/userdata/yanjuntu/hunt/fireShot
#huntdir=/store/disk00/jribnik/hunt
huntdir=/nfs-3/userdata/yanjuntu/hunt
cd $huntdir

if [ $# -ne 1 ]
then
    echo usage: processCands.sh cands
    exit 1
fi

candsfile=$1
candsfiletrunc=`echo $candsfile | awk -F'.' '{print $1}'`

if [ ! -e $candsfile ]
then
    echo $candsfile does not exist, I have no reason to live
    exit 2
fi

fields=`tac $candsfile | head -n-1 | tac | head -n 1 | sed 's/\*//g' | awk '{ORS=" ";for(i=2;i<=NF;i++) print $i;ORS="\n";print ""}' | sed 's/ $//g'`
nf=`echo $fields | awk '{print NF}'`
tac $candsfile | head -n-3 | tac | head -n-1 | sed 's/\*//g' | awk '{ORS=" ";for(i=2;i<=NF;i++) print $i;ORS="\n";print ""}' | sed 's/ $//g' >${candsfiletrunc}_stripped.cands

# Dump and email first, will process,
# i.e. pickEvent, afterward
touch ${candsfiletrunc}_stripped_processed.cands
while read line
do
    grep "$line" ${candsfiletrunc}_stripped_processed.cands >/dev/null 2>&1
    if [ $? -eq 0 ]
    then
        continue
    fi

    run_=`echo $line | awk '{print $2}'`
    ls_=`echo $line | awk '{print $3}'`
    evt_=`echo $line | awk '{print $4}'`
    echo Dumping new candidate in $run_ $ls_ $evt_

    # Dump the event
    # In case there are multiple cands
    # per event suffix with cand index
    candi=0
    dumpfile="${candsfiletrunc}_${run_}_${ls_}_${evt_}_${candi}.dump"
    while [ -e $dumpfile ]
    do
        candi=$(($candi+1))
        dumpfile="${candsfiletrunc}_${run_}_${ls_}_${evt_}_${candi}.dump"
    done

    for i in `seq 5 $nf`; do
        fieldcmd="echo $fields | awk '{print \$$i}'"
        valuecmd="echo $line | awk '{print \$$i}'"
        field=`eval $fieldcmd`
        value=`eval $valuecmd`
        echo $field $value >>$dumpfile
    done

    # No email if the given hunt is defunct
    isdefunct=0
    grep "${candsfiletrunc}" defunct.txt
    if [ $? -eq 0 ]
    then
        isdefunct=1
    fi

    # Email dump unless there exists the do
    # not email file or if hunt is defunct
    if [ ! -e .donotemail ] && [ $isdefunct -eq 0 ]
    then
        #echo "Check out: http://uaf-2.t2.ucsd.edu/~jribnik/hunt/index.php#${candsfiletrunc}_${run_}_${ls_}_${evt_}" >email.tmp
	echo "Check out: http://uaf-2.t2.ucsd.edu/~yanjuntu/hunt/index.php#${candsfiletrunc}_${run_}_${ls_}_${evt_}" >email.tmp
        echo               >>email.tmp
        echo "Here she is:">>email.tmp
        echo               >>email.tmp
        echo "run ${run_}" >>email.tmp
        echo "ls ${ls_}"   >>email.tmp
        echo "evt ${evt_}" >>email.tmp
        cat $dumpfile >>email.tmp
        #cat email.tmp | mail -s "[$candsfiletrunc ALERT] New candidate found!" jribnik@cern.ch
        #cat email.tmp | mail -s "[$candsfiletrunc ALERT] New candidate found!" fgolf@physics.ucsd.edu
        #cat email.tmp | mail -s "[$candsfiletrunc ALERT] New candidate found!" slava77@fnal.gov
        #cat email.tmp | mail -s "[$candsfiletrunc ALERT] New candidate found!" tas@fnal.gov
	cat email.tmp | mail -s "[$candsfiletrunc ALERT] New candidate found!" yanjuntu@cern.gov
    fi

    # Move dump
    mv $dumpfile $huntdir/dumps2xfer
done <${candsfiletrunc}_stripped.cands

# Process events, i.e. pickEvent
while read line
do
    grep "$line" ${candsfiletrunc}_stripped_processed.cands >/dev/null 2>&1
    if [ $? -eq 0 ]
    then
        continue
    fi

    dataset_=`echo $line | awk '{print $1}'`
    run_=`echo $line | awk '{print $2}'`
    ls_=`echo $line | awk '{print $3}'`
    evt_=`echo $line | awk '{print $4}'`
    echo Processing new candidate in $run_ $ls_ $evt_

    # Only pickAnEvent.pl events that have not
    # already been picked
    alreadypicked=0
    if [ -e $fireshotdir/picked.txt ]
    then
        grep "${run_} ${ls_} ${evt_}" $fireshotdir/picked.txt >/dev/null 2>&1
        if [ $? -eq 0 ]
        then
            alreadypicked=1
        fi
    fi

    # If we must pick, only do so for candi==0
    # as we have not yet transferred this pick
    # in order to know that it has been picked
    if [ $alreadypicked -eq 0 ] && [ $candi -eq 0 ]
    then
        # Make unique eventToPick.txt file
        eventtopick=`mktemp eventToPick.txt.XXXXXX`
        echo $dataset_ >$eventtopick

        # Need CMSSW environment for pickAnEvent.pl
        if [ ! -d CMSSW_3_6_1/src ]
        then
            scramv1 p CMSSW CMSSW_3_6_1
        fi
        cd CMSSW_3_6_1/src
        eval `scramv1 ru -sh`
        cd -

        # Finally ready for pickAnEvent.pl
        echo "$run_ $ls_ $evt_" >>$eventtopick
        ./pickAnEvent.pl $eventtopick

        if [ $? -eq 0 ]
        # Event is probably not in dbs yet
        then
            echo $line >>${candsfiletrunc}_stripped_processed.cands

            # Move pick
            tmp=`echo $dataset_ | sed 's!/!_!g'`
            pick="${tmp:1}_${run_}_${ls_}_${evt_}.root"
            mv $pick $fireshotdir/picks

            # Cleanup
            rm "${tmp:1}_${run_}_${ls_}_${evt_}.py" # py file
        fi

        # Cleanup
        rm $eventtopick
    else
        echo $line >>${candsfiletrunc}_stripped_processed.cands
    fi
done <${candsfiletrunc}_stripped.cands
