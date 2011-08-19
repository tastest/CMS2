#################################################################
                            SETUP
#################################################################
Set CMS2_NTUPLE_LOCATION environment variable that points to data 
location. 

export CMS2_NTUPLE_LOCATION=/data/tmp/cms2

You need the CMS2_LOCATION environment variable because lot's of filenames 
are searched for relative to it. It is expected that this envvar points to the CMS2
directory, i.e. the one that contains NtupleMacros directory.

Root version:
use 5.24.  This is needed to get the right roofit version.
UAF: source /home/users/spadhi/ROOT/slc4_ia32_gcc34/root/bin/thisroot.sh
cms-tas01: use /opt/root

Check out the looper:
cvs co -r  SameSignV00-00-04 -d SSDilep UserCode/JRibnik/CMS2/NtupleMacros/SSDilep
cvs co -d Tools UserCode/JRibnik/CMS2/NtupleMacros/Tools

Make the FWLite lib:
cd Tools/MiniFWLite
make

Run it:
root -b
root [0] gSystem->Load("../Tools/MiniFWLite/libMiniFWLite.so");
root [1] gSystem->SetAclicMode(TSystem::kDebug);     
root [2] .x processData.C


#################################################################
                     Make Default Histograms
#################################################################

1) Start ROOT:

  root -l
 
2) execute the processData.C (Any modification in the looper goes to
doAnalysis.C)

  .x processData.C
  .q

3) Show results
  .x doTable.C 

4) show the classification for ttbar
  .x doTTconstituents.C

#################################################################
                    CURRENT WORKING REFERENCE
#################################################################
Output for V01-03-01 samples with reference selection

|Same Sign | Total SM | DY | ttbar | Wjets | WW | WZ  | ZZ | TW | LM0 | LM1 | LM2 | LM3 | LM4 | LM5 | LM6 | LM7 | LM8 | LM9 |
| ee | 4.49627 | 0 | 4.37334 | 0 | 0 | 0 | 0.122926 | 0 | 99.3162|19.462|2.11564|14.4871|4.91576|1.65058|4.47621|2.21723|7.20882|4.97318|
| mumu | 1.67256 | 0 | 1.312 | 0 | 0 | 0.29909 | 0.0614628 | 0 | 119.939|24.2126|3.025|19.1106|6.25091|1.84277|4.54287|2.67476|8.78922|6.41888|
| emu | 4.80132 | 0 | 3.93601 | 0 | 0 | 0.448635 | 0.0614628 | 0.355207 | 225.225|45.5135|4.75092|30.9006|11.7735|3.43682|8.61909|4.15291|16.2201|12.6065|
| total | 10.9701 | 0 | 9.62136 | 0 | 0 | 0.747726 | 0.245851 | 0.355207 | 444.48|89.1881|9.89156|64.4982|22.9402|6.93018|17.6382|9.0449|32.2181|23.9985|

| Same Sign | Type I+II+III | Type-I | Type-II | Type-III | Type-II-SemiLep | Type-II-Fakes | 
| | 4.37334 | 0.874669 | 3.49868 | 0 | 1.74934 | 1.74934 | 
| | 1.312 | 0 | 1.312 | 0 | 1.312 | 0 | 
| | 3.93601 | 1.312 | 2.62401 | 0 | 2.62401 | 0 | 
| | 9.62136 | 2.18667 | 7.43469 | 0 | 5.68535 | 1.74934 | 


