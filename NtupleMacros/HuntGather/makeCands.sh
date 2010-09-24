#!/bin/bash

# Produce cands.txts and do some basic plotting
# Process cands.txts with processCands.sh
# Transfer plots and dumps to uaf-4.t2.ucsd.edu

# Produce cands.txts
root -b -l -q makeCands.C
cd cands

for cands in `ls *.raw.cands`
do
    if [ ! -e $cands.lock ]
    # So that multiple instances of this
    # script don't interfere check for a
    # lock file. If no lock file then we
    # got here first so create one 
    then
        touch $cands.lock
        ./processCands.sh $cands
        rm $cands.lock
    fi
done

scp *.raw.cands uaf-4.t2.ucsd.edu:~/public_html/hunt/
cd ..
scp dumps2xfer/* uaf-4.t2.ucsd.edu:~/public_html/hunt/dumps
mv dumps2xfer/* dumps/