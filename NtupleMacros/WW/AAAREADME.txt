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


#################################################################
                     Produce Default Tables
#################################################################

root -l init.C
.x doTables.C

#################################################################
                    CURRENT WORKING REFERENCE
#################################################################
Output for V01-02-01 samples:
  * CSA07 WW selection
  * tcMET->MET
  * el d0 is in, but data is bad, so electron efficiency is 20-30% lower.
  * d0->d0corr
  * isolation from PAT

|      |   *DY ee*   |  *DY mumu*  | *DY tautau* |    *ttbar*   |    *Wjets*    |    *WZ*     |    *ZZ*     |      *WW*     |     *TW*     | 
| emu  | 0.0 +/- 0.0 | 0.0 +/- 0.0 | 0.0 +/- 0.0 | 26.6 +/- 3.6 | 30.5 +/- 12.5 | 3.2 +/- 0.2 | 0.1 +/- 0.0 | 129.8 +/- 3.1 | 12.8 +/- 1.5 |  
| mumu | 0.0 +/- 0.0 | 4.2 +/- 4.2 | 0.0 +/- 0.0 |  6.3 +/- 1.7 |  0.0 +/- 0.0  | 1.7 +/- 0.1 | 1.4 +/- 0.1 |  59.2 +/- 2.1 |  5.0 +/- 0.9 |  
| ee   | 0.0 +/- 0.0 | 0.0 +/- 0.0 | 0.0 +/- 0.0 |  6.3 +/- 1.7 |  0.0 +/- 0.0  | 0.5 +/- 0.1 | 0.6 +/- 0.1 |  19.9 +/- 1.2 |  2.2 +/- 0.6 |  
| total| 0.0 +/- 0.0 | 4.2 +/- 4.2 | 0.0 +/- 0.0 | 39.2 +/- 4.4 | 30.5 +/- 12.5 | 5.4 +/- 0.2 | 2.1 +/- 0.1 | 208.9 +/- 4.0 | 20.0 +/- 1.9 | 

#################################################################
         EXPECTED RESULT BASED ON CSA07 SCALED to 10TeV
#################################################################

	   DY ee        DY mumu       DY tautau        ttbar         Wjets            WZ           ZZ              WW             TW
emu -    0.0 ± 0.0     1.5 ± 1.1      14.9 ± 3.8    30.2 ± 2.4     27.5 ± 6.3      5.1 ± 0.7    0.2 ± 0.1    225.8 ± 4.7      13.5 ± 0.9
mumu -   0.0 ± 0.0     7.6 ± 2.8       0.6 ± 0.6    14.3 ± 1.7      0.7 ± 0.7      2.2 ± 0.4    1.1 ± 0.3     58.0 ± 2.4       5.0 ± 0.6
ee -     1.4 ± 1.0     0.0 ± 0.0       0.0 ± 0.0     7.5 ± 1.2      2.8 ± 1.4      1.5 ± 0.4    1.4 ± 0.3     40.1 ± 2.0       2.9 ± 0.4
