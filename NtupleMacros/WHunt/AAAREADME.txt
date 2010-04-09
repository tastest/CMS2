mkdir baby skim
cvs co -rHEAD -d CORE UserCode/JRibnik/CMS2/NtupleMacros/CORE
sed 's!#define PARANOIA!//#define PARANOIA!' CORE/CMS2.h >CORE/CMS2.h.new; mv CORE/CMS2.h.new CORE/CMS2.h
