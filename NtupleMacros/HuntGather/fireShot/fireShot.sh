#!/bin/bash

echo "fireShot go on DISPLAY=${DISPLAY}"
fireShotdir="/store/disk00/jribnik/fireShot"
echo "fireShotdir=${fireShotdir}"

while [ 1 ]
do
    for f in `find ${fireShotdir}/picks -name "*.root"`
    do
        run_=`echo $f | sed 's/^.*FEVT_\(.*\).root$/\1/' | awk -F'_' '{print $1}'`
        ls_=`echo $f | sed 's/^.*FEVT_\(.*\).root$/\1/' | awk -F'_' '{print $2}'`
        evt_=`echo $f | sed 's/^.*FEVT_\(.*\).root$/\1/' | awk -F'_' '{print $3}'`

        # File name sans extension
        namesansext="fireworks_${run_}_${ls_}_${evt_}"

        # Already done?
        ls ${fireShotdir}/shots/${namesansext}.png >/dev/null 2>&1
        if [ $? -eq 0 ]
        then
            continue
        fi
        echo Yeehaw time to ride!

        # There _should_ be a cmsShow instance
        # ready and willing on localhost:$PORT
        echo "Psst, cmsShow, I've got the goods: $f"
        echo "$f" | nc -w 10 localhost 9999
        sleep 30

        echo Say CHEESE cmsShow!
        xwd -display $DISPLAY -root -out ${fireShotdir}/xwd/${namesansext}.xwd

        if [ -e ${fireShotdir}/xwd/${namesansext}.xwd ]
        then
            echo Found dump, converting to png
            convert ${fireShotdir}/xwd/${namesansext}.xwd ${fireShotdir}/shots/${namesansext}.png
        else
            echo Dump not found, not good
        fi

        # scp pick and shot to uaf
        if [ -e ${fireShotdir}/shots/${namesansext}.png ]
        then
            echo scp, as easy as 1, 2, 3...
            scp $f uaf-4.t2.ucsd.edu:~/public_html/hunt/picks/
            scp ${fireShotdir}/shots/${namesansext}.png uaf-4.t2.ucsd.edu:~/public_html/hunt/shots/
        else
            echo Shot not found, not good
        fi

        # update list of picked events so
        # we don't duplicate efforts
        ls picks | sed 's/\.root//' | awk -F'_' '{print $4,$5,$6}' >picked.txt

        echo "All done here, waiting again :("
    done

    sleep 10
done
