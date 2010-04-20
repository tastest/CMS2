#!/bin/bash

# Run basic plotting and scanning scripts
# and transfer resulting plots and dumps
# to uaf. Note that because cms-tas0? is
# not slc5 we are executing processCands.sh
# on lxplus5 and transferring everything
# back over

# Produces zcands.txt and topcands.txt
root -b -l -q makePlots.C

scp -o "StrictHostKeyChecking no" zcands.txt lxplus5:~/scratch0/dilephunt
ssh -o "StrictHostKeyChecking no" lxplus5 /afs/cern.ch/user/j/jribnik/scratch0/dilephunt/processCands.sh zcands.txt /tas03/disk01/dilephunt/zcands

scp -o "StrictHostKeyChecking no" topcands.txt lxplus5:~/scratch0/dilephunt
ssh -o "StrictHostKeyChecking no" lxplus5 /afs/cern.ch/user/j/jribnik/scratch0/dilephunt/processCands.sh topcands.txt /tas03/disk01/dilephunt/topcands

scp plots/* uaf-4.t2.ucsd.edu:~/public_html/dilephunt/plots
scp zcands/dumps/* uaf-4.t2.ucsd.edu:~/public_html/dilephunt/zcands/dumps
scp topcands/dumps/* uaf-4.t2.ucsd.edu:~/public_html/dilephunt/topcands/dumps
