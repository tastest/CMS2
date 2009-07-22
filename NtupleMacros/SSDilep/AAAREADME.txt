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

3) Show results
  .x .x doTable.C 


#################################################################
                    CURRENT WORKING REFERENCE
#################################################################
Output for V01-03-01 samples with reference selection

| Same Sign | DY | ttbar | Wjets | WW | WZ  | ZZ | TW | LM0 | LM1 | LM2 | LM3 | LM4 | LM5 | LM6 | LM7 | LM8 | LM9 |
| ee | 0 | 7.87202 | 0 | 0 | 0.29909 | 0.122926 | 0 | 103.658|23.4464|2.5796|15.4118|5.52265|1.94452|6.00002|2.35801|7.30587|5.20449|
| mumu | 0 | 1.74934 | 0 | 0 | 0.747726 | 0.122926 | 0 | 124.281|29.5761|3.32193|20.2665|7.64674|2.19324|6.12383|2.67476|8.81694|6.99716|
| emu | 0 | 6.56002 | 0 | 0 | 1.49545 | 0 | 0.53281 | 225.225|57.7731|5.58604|32.2106|13.4121|4.06992|11.2286|4.50485|16.5389|12.78|
| total | 0 | 16.1814 | 0 | 0 | 2.54227 | 0.245851 | 0.53281 | 453.164|110.796|11.4876|67.8888|26.5815|8.20768|23.3525|9.53762|32.6617|24.9816|

