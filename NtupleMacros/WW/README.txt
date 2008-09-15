This is main CVS area for all common tools for WW and related 
analyses.

How to make default histograms
==========================================

0) make dataset.txt file which contains prefix of the full path or name
of all data. Example: /home/users/dmytro/ntuples/V00-04-01/
All files are assumed to be in this directory (dcache is allowed as well).

1) Start ROOT:

root -l
 
2) load some needed libraries ( you may put this in your rootlogon.C )

.x init.C

3) load and compile main looper. It produces default histograms. Output
is stored in a single ROOT file.

.L fkwLoopingFunctionFast.C++

4) process data

.x processData.C


How to produce default tables
=======================================
root -l doTables.C

Output for V00-04-01:
 Njet = 0
| | DY ee| DY mumu | DY tautau | ttbar | Wjets | WZ | ZZ | WW | TW |
| emu | 0 | 18.1821 | 29.9597 | 237.354 | 190.519 | 0 | 0 | 443.054 | 66.5507 |
| mumu | 0 | 61.8788 | 1.62716 | 75.0155 | 1.03004 | 0 | 0 | 102.434 | 21.1946 |
| ee | 50.772 | 0 | 0.92972 | 72.237 | 53.6584 | 0 | 0 | 89.9947 | 17.8035 |
| total | 50.772 | 80.0609 | 32.5165 | 384.606 | 245.207 | 0 | 0 | 635.483 | 105.549 |

 Njet = 1
| | DY ee| DY mumu | DY tautau | ttbar | Wjets | WZ | ZZ | WW | TW |
| emu | 0 | 18.6932 | 22.3281 | 1581.01 | 75.6349 | 0 | 0 | 149.889 | 278.213 |
| mumu | 0 | 121.746 | 3.42805 | 474.38 | 2.42419 | 0 | 0 | 41.619 | 84.0715 |
| ee | 77.2118 | 0 | 3.64942 | 466.77 | 18.3212 | 0 | 0 | 37.0117 | 67.3985 |
| total | 77.2118 | 140.44 | 29.4056 | 2522.17 | 96.3803 | 0 | 0 | 228.52 | 429.683 |

 Njet = 2
| | DY ee| DY mumu | DY tautau | ttbar | Wjets | WZ | ZZ | WW | TW |
| emu | 0 | 2.00568 | 15.5862 | 3149.71 | 18.6777 | 0 | 0 | 43.3083 | 169.697 |
| mumu | 0 | 103.874 | 3.8032 | 969.05 | 0 | 0 | 0 | 14.8968 | 54.5405 |
| ee | 57.9894 | 0 | 3.84994 | 817.379 | 5.80704 | 0 | 0 | 8.44663 | 40.1283 |
| total | 57.9894 | 105.88 | 23.2394 | 4936.14 | 24.4848 | 0 | 0 | 66.6517 | 264.366 |

Note: WZ and ZZ was not processed, so ignore their values.






############################################################################
                          Old Information
############################################################################


==========================================

As of 8/30/08 the state of affairs is as follows:

root -l
.x init.C
.L fkwLoopingFunctionFast.C++
.x doWWFast_diet.C
.q

root -l
.x doTable.C

This combination will get you the table:

| | DY ee| DY mumu | DY tautau | ttbar | Wjets | WZ | ZZ | WW | TW | 
| emu | 0 | 20.3318 | 28.0355 | 221.59 | 171.629 | 13.3602 | 0.78749 | 443.054 | 66.5507 | 
| mumu | 0 | 93.7201 | 3.20012 | 71.0864 | 1.03094 | 8.81501 | 2.24997 | 102.434 | 21.1946 | 
| ee | 46.6652 | 0 | 0.944629 | 70.355 | 49.8651 | 4.82071 | 2.47497 | 89.9947 | 17.8035 | 
| total | 46.6652 | 114.052 | 32.1802 | 363.032 | 222.525 | 26.996 | 5.51243 | 635.483 | 105.549 | 

To reproduce this yourself, you need to:

1.) move the data files from uaf-3.t2.ucsd.edu:/data/tmp/dietcms2 
    to wherever you work, and change the doWWFast_diet.C accordingly.
2.) change doTable.C to pick up whatever file you wrote your histograms into.
