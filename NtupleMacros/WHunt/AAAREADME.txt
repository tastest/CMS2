mkdir baby skim plots dumps
cvs co -rHEAD -d CORE UserCode/JRibnik/CMS2/NtupleMacros/CORE
# get Dave/Yanyan data MC comparison code
cvs up -r whunt_v4 HybridLooper/oldstyle
cd HybridLooper/oldstyle
make build -j 4
cd -
# fix local CMS2.h
sed 's!#define PARANOIA!//#define PARANOIA!' CORE/CMS2.h >CORE/CMS2.h.new; mv CORE/CMS2.h.new CORE/CMS2.h
