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
  .x showResults.C+


#################################################################
                    CURRENT WORKING REFERENCE
#################################################################
Output for V01-02-06 samples with reference selection


|      |  *DY ee*  | *DY mumu* |*DY tautau*|  *ttbar*  |  *Wjets*  |    *WZ*   |    *ZZ*   |    *WW*    |    *TW*   |
| ee   | 0.2+/-0.2 | 0.0+/-0.0 | 0.0+/-0.0 | 0.4+/-0.1 | 0.0+/-0.0 | 0.1+/-0.0 | 0.1+/-0.0 |  3.9+/-0.4 | 0.3+/-0.1 |
| mumu | 0.0+/-0.0 | 0.2+/-0.2 | 0.0+/-0.0 | 0.8+/-0.2 | 0.0+/-0.0 | 0.2+/-0.1 | 0.1+/-0.0 |  6.8+/-0.5 | 0.4+/-0.1 |
| em   | 0.0+/-0.0 | 0.5+/-0.3 | 0.7+/-0.4 | 2.1+/-0.3 | 4.7+/-1.6 | 0.7+/-0.1 | 0.1+/-0.0 | 24.4+/-0.9 | 1.3+/-0.2 |
| total| 0.2+/-0.2 | 0.9+/-0.4 | 0.7+/-0.4 | 3.3+/-0.4 | 4.7+/-1.6 | 1.0+/-0.1 | 0.3+/-0.0 | 35.0+/-1.1 | 1.9+/-0.2 |

CSA07 Expectations (scale to 10TeV)

|      |   *DY ee*   |  *DY mumu*  | *DY tautau* |   *ttbar*   |   *Wjets*   |     *WZ*    |     *ZZ*    |      *WW*    |      *TW*    |
| ee   | 0.14+/-0.10 |  0.0+/-0.0  |  0.0+/-0.0  | 0.75+/-0.12 | 0.28+/-0.14 | 0.15+/-0.04 | 0.14+/-0.03 |  4.01+/-0.20 |  0.29+/-0.04 |
| mumu |  0.0+/-0.0  | 0.76+/-0.28 | 0.06+/-0.06 | 1.43+/-0.17 | 0.07+/-0.07 | 0.22+/-0.04 | 0.11+/-0.03 |  5.80+/-0.24 |  0.50+/-0.06 |
| emu  |  0.0+/-0.0  | 0.15+/-0.11 | 1.49+/-0.38 | 3.02+/-0.24 | 2.75+/-0.63 | 0.51+/-0.07 | 0.02+/-0.01 | 22.58+/-0.47 |  1.35+/-0.09 |
| total| 0.14+/-0.10 | 0.91+/-0.30 | 1.55+/-0.38 | 5.19+/-0.32 | 3.10+/-0.65 | 0.89+/-0.09 | 0.27+/-0.04 | 32.39+/-0.57 |  2.14+/-0.12 |
