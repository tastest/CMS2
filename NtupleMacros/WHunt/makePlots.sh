#!/bin/bash

# Run basic plotting and scanning scripts
# and transfer resulting plots and dumps
# to uaf. Note that because cms-tas0? is
# not slc5 we are executing processCands.sh
# on lxplus5 and transferring everything
# back over

# Produces wcands.txt
root -b -l -q makePlots.C

scp -o "StrictHostKeyChecking no" wcands.txt lxplus5:~/scratch0/whunt
ssh -o "StrictHostKeyChecking no" lxplus5 /afs/cern.ch/user/j/jribnik/scratch0/whunt/processCands.sh wcands.txt /tas03/disk01/whunt

scp plots/* uaf-4.t2.ucsd.edu:~/public_html/whunt/plots
scp dumps/* uaf-4.t2.ucsd.edu:~/public_html/whunt/dumps
