#!/bin/bash

# Run basic plotting and scanning scripts
# and transfer resulting plots and dumps
# to uaf. Note that because cms-tas0? is
# not slc5 we are executing processCands.sh
# on lxplus5 and transferring everything
# back over

# Produces zcands.txt and topcands.txt
root -b -l -q makePlots.C

scp zcands.txt lxplus303:~/scratch0/dilephunt
ssh lxplus303 /afs/cern.ch/user/j/jribnik/scratch0/dilephunt/processCands.sh zcands.txt /tas03/disk01/dilephunt/zcands

scp topcands.txt lxplus303:~/scratch0/dilephunt
ssh lxplus303 /afs/cern.ch/user/j/jribnik/scratch0/dilephunt/processCands.sh topcands.txt /tas03/disk01/dilephunt/topcands

scp plots/* uaf-2.t2.ucsd.edu:~/public_html/dilephunt/plots

scp zcands/dumps/* uaf-2.t2.ucsd.edu:~/public_html/dilephunt/zcands/dumps
scp zcands/picks/* uaf-2.t2.ucsd.edu:~/public_html/dilephunt/zcands/picks

scp topcands/dumps/* uaf-2.t2.ucsd.edu:~/public_html/dilephunt/topcands/dumps
scp topcands/picks/* uaf-2.t2.ucsd.edu:~/public_html/dilephunt/topcands/picks

ssh uaf-6.t2.ucsd.edu "touch ~/devel/dilephunt/.fireShot.bang"
