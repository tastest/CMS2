# Create initial directory structure
mkdir cands dumps dumps2xfer logs emu_skim emu_baby dilep_skim dilep_baby trilep_skim trilep_baby

# Be  indepedent  of the  user  when
# when running root
mv rootrc .rootrc
# make sure you are using root 5.26 or later

# Check out the HEAD of CORE; if it
# doesn't work fix it
cvs co -rHEAD -d CORE UserCode/JRibnik/CMS2/NtupleMacros/CORE

# Turn off the paranoia!
sed 's!#define PARANOIA!//#define PARANOIA!' CORE/CMS2.h >CORE/CMS2.h.new
mv CORE/CMS2.h.new CORE/CMS2.h

# This script runs forever and  does
# all sorts of things.
# Read the script to understand what
# I mean by that
./UpdateSkims_daemon.sh

# FYI about some of these files:
#
# dead.txt: these hunts are dead which
# means that they should be  commented
# out of makeCands.C completely
#
# defunct.txt: these hunts are defunct
# which means that while they are  run
# email notifications are not sent for
# new candidates

# For the gathering...

# First, you need a goodruns JSON file.
cp /afs/cern.ch/user/s/slava77/public/jsons/august30/special/Cert_TopAug30_Merged_135059-144114_recover_noESDCS.txt ./json.txt
# You can also try
# /nfs-3/userdata/yanjuntu/hunt/runlists/
# for runlists.

# Put this in at the appropriate place
# at the top of gather_doAll.C

# Both gather_example_doAll.C, and
# gather_several_doAll.C have been
# deprecated. Instead, I recommend you
# run the full gathering code, which is
# found in the file gather_doAll.C, and
# if necessary, tailor it to your needs.

# Just do
# root gather_doAll.C

# This calls functions which are defined
# in gather.C. That in turn uses 
# BabySample.h, and other code.