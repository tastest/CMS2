#! /bin/bash

MAXJOBS=5
if [ $1 != "" ]; then
    if [ $1 -le 1000 ] && [ $1 -ge 1 ]; then
	MAXJOBS=$1
    fi
fi


echo "Submitting $MAXJOBS test jobs."
echo "Use   condor_q $USER   to check your jobs" 
echo "When your jobs finish you can find their std out/err files in /data/tmp/${USER}/TestJob/std_logs"

JOBNUM=1
while : ; do
    ./lib_sh/submit.sh -e ./test/test.sh -a $JOBNUM -i ./test/test.txt -u "TestJob"
    if [ $JOBNUM = $MAXJOBS ]; then
	break
    fi
    JOBNUM=$((JOBNUM+1))
done
echo 
echo "Done submitting test jobs"