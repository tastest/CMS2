This directory is used for pushing the WW analysis as of February
forwards into the soup area, and a public presentation in CMS
that is more or less complete.

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
