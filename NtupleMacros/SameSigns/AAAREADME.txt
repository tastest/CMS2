Made a skeleton dilepton analysis in SameSigns

For a quick start, do the following:
   cd SameSigns
   make -j 8 Results.tbl

Cut names are defined in Looper.h.
Cuts are defined in Looper.cc.
Histograms are booked and filled in Looper.cc.
ConfigAndRun.h defines samples and cuts to use for tables.

To make a table, add a rule in ConfigAndRun.h and run:
   make -j 8 <your table name>.tbl
