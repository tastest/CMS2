#!/bin/bash

# get root
export ROOTSYS=/code/osgcode/UCSD_root/root_v5.24.00
export PATH=$PATH:$ROOTSYS/bin
export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH

# get crab
source /code/osgcode/ucsdt2/Crab/etc/crab.sh

# gather it
cd /home/users/dlevans/gathering/HuntGather
root -q -b gather_doAll.C

# copy plots to web dir
cp *.png ~jribnik/public_html/gather/

