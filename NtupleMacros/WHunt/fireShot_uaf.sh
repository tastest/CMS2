#!/bin/bash

cd ~/devel/whunt

# Need CMSSW environment for cmsShow
if [ ! -d CMSSW_3_5_6/src ]; then
    scramv1 p CMSSW CMSSW_3_5_6
fi
cd CMSSW_3_5_6/src
eval `scramv1 ru -sh`
cd -

echo Waiting for the starting gun...
while [ 1 ]; do
    if [ -e .fireShot.bang ]; then
        echo Yeehaw time to ride!

        for f in `find ~/public_html/whunt/picks -name "*.root"`; do
            run_=`echo $f | sed 's/^.*FEVT_\(.*\).root$/\1/' | awk -F'_' '{print $1}'`
            ls_=`echo $f | sed 's/^.*FEVT_\(.*\).root$/\1/' | awk -F'_' '{print $2}'`
            evt_=`echo $f | sed 's/^.*FEVT_\(.*\).root$/\1/' | awk -F'_' '{print $3}'`

            # Already done?
            ls ~/public_html/whunt/shots/fireworks_${run_}_${ls_}_${evt_}.png >/dev/null 2>&1
            if [ $? -eq 0 ]; then
                continue;
            fi

            cmsShow -c whunt.fwc $f &
            #pid=`ps u | grep cmsShow.exe | awk '{print $2}'`
            #echo pid: $pid
            sleep 30

            ssh noms-1 "rm ~/.fireShot.done"
            ssh noms-1 "echo $run_ $ls_ $evt_ >~/.fireShot.prep"

            while [ 1 ]; do
                ssh noms-1 "ls ~/.fireShot.done >/dev/null 2>&1"
                if [ $? -eq 0 ]; then
                    #kill -9 $pid
                    killall cmsShow.exe
                    break
                fi

                sleep 10
            done

            if [ -e ~/devel/whunt/xwd/fireworks_${run_}_${ls_}_${evt_}.xwd ]; then
                echo Found dump, converting to png
                convert ~/devel/whunt/xwd/fireworks_${run_}_${ls_}_${evt_}.xwd ~/public_html/whunt/shots/fireworks_${run_}_${ls_}_${evt_}.png
            else
                echo Dump not found, not good
            fi
        done

        rm .fireShot.bang
    fi

    sleep 10
done
