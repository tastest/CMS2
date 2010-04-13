#!/bin/bash

echo Waiting for launch codes on $DISPLAY...

while [ 1 ]; do
    if [ -e ~/.fireShot.prep ]; then
        echo Incoming!

        run_=`awk '{print $1}' ~/.fireShot.prep`
        ls_=`awk '{print $2}' ~/.fireShot.prep`
        evt_=`awk '{print $3}' ~/.fireShot.prep`

        xwd -display $DISPLAY -root -out fireworks_${run_}_${ls_}_${evt_}.xwd
        scp fireworks_${run_}_${ls_}_${evt_}.xwd uaf-6:~/devel/whunt/xwd
        rm fireworks_${run_}_${ls_}_${evt_}.xwd

        rm ~/.fireShot.prep
        touch ~/.fireShot.done
    fi

    sleep 10
done
