#################################################################
                            SETUP
#################################################################
Set CMS2_NTUPLE_LOCATION environment variable that points to data 
location. 

  setenv CMS2_NTUPLE_LOCATION /data/tmp/
  setenv CMS2_LOCATION ../../../CMS2/

You need the CMS2_LOCATION environment variable because lot's of filenames 
are searched for relative to it. It is expected that this envvar points to the CMS2
directory, i.e. the one that contains NtupleMacros directory.

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

| | DY | ttbar | Wjets | WW | WZ  | ZZ | TW | LM0 | LM1 | LM2 | LM3 | LM4 | LM5 | LM6 | LM7 | LM8 | LM9 |
| ee | 0 | 7.87202 | 0 | 0 | 0.29909 | 0.122926 | 0 | 102.573|23.5996|2.61672|15.4888|5.58333|1.93321|5.98097|2.35801|7.41677|5.20449|
| mumu | 0 | 1.74934 | 0 | 0 | 0.747726 | 0.122926 | 0 | 124.824|29.8826|3.32193|20.2665|7.58605|2.19324|6.12383|2.70995|8.8724|7.05499|
| emu | 0 | 6.56002 | 0 | 0 | 1.49545 | 0 | 0.53281 | 223.597|57.4666|5.62316|32.5188|13.5335|4.08123|11.2286|4.46966|16.6775|12.7221|
| total | 0 | 16.1814 | 0 | 0 | 2.54227 | 0.245851 | 0.53281 | 450.993|110.949|11.5618|68.2741|26.7029|8.20768|23.3334|9.53762|32.9667|24.9816|


| Same Sign | Type I+II+III | Type-I | Type-II | Type-III | Type-II-SemiLep | Type-II-Fakes | 
| | 7.87202 | 2.18667 | 5.68535 | 0 | 4.37334 | 1.312 | 
| | 1.74934 | 0 | 1.74934 | 0 | 1.74934 | 0 | 
| | 6.56002 | 3.06134 | 3.49868 | 0 | 3.49868 | 0 | 
| | 16.1814 | 5.24801 | 10.9334 | 0 | 9.62136 | 1.312 | 


