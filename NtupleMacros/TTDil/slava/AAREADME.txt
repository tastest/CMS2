To run the looper
root
// load FWLite stuff
.L doAll.C
doAll(cutBits)
// doAll(cutBits, skipFWLite=true) -- is possible to turn off loading FWlite --
// this only seems to work in compiled doAll case:
root -l -b -q "doAll.C+(bits, true)"; //need to do it twice though (the first time it confuses itself with libs)
   The default value for the bitmask is : 2**10 + 2**13 + 2**14 + 2**15 + 2**18 + 2**19 + 2**20  
   See doAll.C for info on what those bits mean
// the option to avoid this reqmnt to "compile" and then "load and run" is to pre-load all relevant libs ... lazy
//this will make myHist_cutBits.root
-------------------------

To make hisotgrams:

root -l -b -q "makeAllPlots.C(\"myHist_cutBits.root\")" # will make linear-scale plots in out/ directory
root -l -b -q "makeAllPlots.C(\"myHist_cutBits.root\", true)" # will make log-scale plots in out/directory

#there are more flags to makeAllPlots ....

---------------------------

