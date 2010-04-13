#!/bin/bash

# Run basic plotting and scanning script
# and transfer plots and dumps to uaf
# Note that because cms-tas0? is not slc5
# we are executing wcands.sh on lxplus5
# and transferring everything back over
root -b -l -q makePlots.C
scp wcands.txt lxplus303:~/scratch0/whunt
ssh lxplus303 /afs/cern.ch/user/j/jribnik/scratch0/whunt/wcands.sh
scp plots/* uaf-2.t2.ucsd.edu:~/public_html/whunt/plots
scp dumps/* uaf-2.t2.ucsd.edu:~/public_html/whunt/dumps
scp picks/* uaf-2.t2.ucsd.edu:~/public_html/whunt/picks
