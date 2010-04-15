mkdir baby skim plots picks dumps
cvs co -rHEAD -d CORE UserCode/JRibnik/CMS2/NtupleMacros/CORE
# get Dave/Yanyna data MC comparison code
cd ..
cvs up -r whunt_v4 HybridLooper/oldstyle
cd HybridLooper/oldstyle
make build -j 4
cd ../../WHunt
# fix local CMS2.h
sed 's!#define PARANOIA!//#define PARANOIA!' CORE/CMS2.h >CORE/CMS2.h.new; mv CORE/CMS2.h.new CORE/CMS2.h
