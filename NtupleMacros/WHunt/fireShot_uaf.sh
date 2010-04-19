#!/bin/bash

cd ~/devel/fireShot

echo Waiting for the starting gun...
while [ 1 ]
do
    if [ -e .fireShot.bang ]
    then
        echo Yeehaw time to ride!

        for f in `find ~/devel/fireShot/picks -name "*.root"`
        do
            run_=`echo $f | sed 's/^.*FEVT_\(.*\).root$/\1/' | awk -F'_' '{print $1}'`
            ls_=`echo $f | sed 's/^.*FEVT_\(.*\).root$/\1/' | awk -F'_' '{print $2}'`
            evt_=`echo $f | sed 's/^.*FEVT_\(.*\).root$/\1/' | awk -F'_' '{print $3}'`

            # Already done?
            ls ~/devel/fireShot/shots/fireworks_${run_}_${ls_}_${evt_}.png >/dev/null 2>&1
            if [ $? -eq 0 ]
            then
                echo Skipping $run_ $ls_ $evt_
                continue
            fi

            # There _should_ be a cmsShow instance
            # ready and willing on localhost:$PORT
            echo "Psst, cmsShow, I've got the goods: $f"
            echo "$f" | nc -w 10 localhost 9999
            sleep 30

            ssh noms-1 "rm ~/.fireShot.done"
            ssh noms-1 "echo $run_ $ls_ $evt_ >~/.fireShot.prep"

            while [ 1 ]
            do
                ssh noms-1 "ls ~/.fireShot.done >/dev/null 2>&1"
                if [ $? -eq 0 ]
                then
                    break
                fi

                sleep 10
            done

            if [ -e ~/devel/fireShot/xwd/fireworks_${run_}_${ls_}_${evt_}.xwd ]
            then
                echo Found dump, converting to png
                convert ~/devel/fireShot/xwd/fireworks_${run_}_${ls_}_${evt_}.xwd ~/devel/fireShot/shots/fireworks_${run_}_${ls_}_${evt_}.png
            else
                echo Dump not found, not good
            fi
        done

        rm .fireShot.bang
        echo "That was fun, thanks!"
    fi

    sleep 10
done
