                       FWLITE is _highly_ recommended  

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
 
2) load some needed libraries ( you may put this in your rootlogon.C )

  .x init.C

3) load and compile main looper. It produces default histograms. Output
is stored in a single ROOT file.

  .L doAnalysis.C++

4) process data. You may edit file before you run to adjust what samples 
to use 

  .x processData.C

4) Show results
  .x showResults.C


#################################################################
                    CURRENT WORKING REFERENCE
#################################################################
Output for V01-02-06 samples with reference selection

|      |  *DY ee*  | *DY mumu* |*DY tautau*|   *ttbar*  |   *Wjets*   |    *WZ*   |    *ZZ*   |     *WW*    |    *TW*    |
| ee   | 0.0+/-0.0 | 0.0+/-0.0 | 0.0+/-0.0 |  4.4+/-1.4 |  0.0+/-0.0  | 0.8+/-0.1 | 1.1+/-0.1 |  41.6+/-2.7 |  2.3+/-0.7 | 
| mumu | 0.0+/-0.0 | 3.3+/-3.3 | 0.0+/-0.0 |  7.5+/-1.8 |  0.0+/-0.0  | 2.3+/-0.2 | 1.7+/-0.1 |  71.5+/-3.5 |  4.0+/-0.9 | 
| em   | 3.3+/-3.3 | 0.0+/-0.0 | 6.7+/-4.7 | 20.3+/-3.0 | 57.0+/-17.2 | 6.3+/-0.3 | 0.1+/-0.0 | 259.7+/-6.7 | 12.3+/-1.6 | 
| total| 3.3+/-3.3 | 3.3+/-3.3 | 6.7+/-4.7 | 32.3+/-3.8 | 57.0+/-17.2 | 9.5+/-0.3 | 2.9+/-0.1 | 372.8+/-8.0 | 18.5+/-2.0 | 


CSA07 Expectations (scale to 10TeV)

|      |  *DY ee*  | *DY mumu* |*DY tautau* |  *ttbar*   |  *Wjets*   |    *WZ*   |    *ZZ*   |     *WW*    |    *TW*    |
| ee   | 1.4+/-1.0 | 0.0+/-0.0 |  0.0+/-0.0 |  7.5+/-1.2 |  2.8+/-1.4 | 1.5+/-0.4 | 1.4+/-0.3 |  40.1+/-2.0 |  2.9+/-0.4 |
| mumu | 0.0+/-0.0 | 7.6+/-2.8 |  0.6+/-0.6 | 14.3+/-1.7 |  0.7+/-0.7 | 2.2+/-0.4 | 1.1+/-0.3 |  58.0+/-2.4 |  5.0+/-0.6 |
| emu  | 0.0+/-0.0 | 1.5+/-1.1 | 14.9+/-3.8 | 30.2+/-2.4 | 27.5+/-6.3 | 5.1+/-0.7 | 0.2+/-0.1 | 225.8+/-4.7 | 13.5+/-0.9 |
| total| 1.4+/-1.0 | 9.1+/-3.0 | 15.5+/-3.8 | 51.9+/-3.2 | 31.0+/-6.5 | 8.9+/-0.9 | 2.7+/-0.4 | 323.9+/-5.7 | 21.4+/-1.2 |
