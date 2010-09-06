# Create initial directory structure
mkdir cands dumps dumps2xfer logs emu_skim emu_baby dilep_skim dilep_baby trilep_skim trilep_baby

# Be  indepedent  of the  user  when
# when running root
mv rootrc .rootrc

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
