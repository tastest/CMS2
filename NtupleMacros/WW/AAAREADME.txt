                       FWLITE is required 

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

root -l doTables.C

Output as of 11/18/2008:
  * feb selection
  * nTrk>2
  * MET->tcMET
  * el with calo iso - 0.92
  * trkJet veto
  * extra muon veto

 Njet = 0
| | DY ee| DY mumu | DY tautau | ttbar | Wjets | WZ | ZZ | WW | TW | 
| emu | 0 | 2.23398 | 21.9215 | 67.1341 | 40.4379 | 7.85087 | 0.337496 | 347.383 | 30.0962 | 
| mumu | 0 | 11.1401 | 0.906839 | 31.6739 | 1.03004 | 3.44336 | 1.68748 | 89.2269 | 11.1625 | 
| ee | 1.99738 | 0 | 0 | 16.629 | 4.06634 | 2.34149 | 2.13747 | 61.7374 | 6.35838 | 
| total | 1.99738 | 13.3741 | 22.8283 | 115.437 | 45.5343 | 13.6357 | 4.16245 | 498.348 | 47.6171 | 
  
 Njet = 1
| | DY ee| DY mumu | DY tautau | ttbar | Wjets | WZ | ZZ | WW | TW | 
| emu | 0 | 0 | 1.11699 | 93.4652 | 8.08419 | 4.68298 | 0.224997 | 39.9296 | 26.5638 | 
| mumu | 0 | 9.76207 | 0.601962 | 39.2934 | 1.05197 | 0.688673 | 0.78749 | 13.0539 | 9.46692 | 
| ee | 2.62537 | 0 | 0 | 25.1244 | 0 | 1.10188 | 0 | 7.67875 | 5.93449 | 
| total | 2.62537 | 9.76207 | 1.71895 | 157.883 | 9.13616 | 6.47353 | 1.01249 | 60.6623 | 41.9652 | 
  
 Njet = 2
| | DY ee| DY mumu | DY tautau | ttbar | Wjets | WZ | ZZ | WW | TW | 
| emu | 0 | 0 | 0.601962 | 44.9553 | 0 | 0.964142 | 0 | 3.6858 | 7.77135 | 
| mumu | 0 | 0 | 0 | 13.7726 | 0 | 0.413204 | 0.337496 | 1.38218 | 1.83687 | 
| ee | 0 | 0 | 0 | 6.49742 | 0 | 0.275469 | 0 | 0.614301 | 1.69557 | 
| total | 0 | 0 | 0.601962 | 65.2253 | 0 | 1.65281 | 0.337496 | 5.68228 | 11.3038 | 
