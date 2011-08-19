# Be  indepedent  of the  user  when
# when running root
mv rootrc .rootrc

# Check out the HEAD of CORE; if it
# doesn't work fix it
cvs co -rHEAD -d CORE UserCode/JRibnik/CMS2/NtupleMacros/CORE

# Turn off the paranoia!
sed 's!#define PARANOIA!//#define PARANOIA!' CORE/CMS2.h >CORE/CMS2.h.new
mv CORE/CMS2.h.new CORE/CMS2.h

# Now, for every MC sample for which
# you want to make skims and babies,
# create a MYSAMPLENAMEHERE.txt that
# contains a list of all merged*.root
# files, full paths please, then bah
./bah.sh
