# GENERAL COMMENT HOW TO DEAL WITH THE CODE HERE
# checkout the code either based on the prescribed tag

#########
cd ../../
cvs -nq update  -rPUTATAGHERE CORE TTDil NtupleTools Tools Templates |\
  grep -v "root$\|log$\|^\? [^$]*file\|^\?[^$]*eventsCMS2\|\.d$\|\_o$\|_b$\|\.list$\|\.tmp$\|^\? [^$]*out$" 
then commit what's still missing or update wrt whateve someone else inserted in these packs
cvs tag [-b]TAGNAME CORE TTDil NtupleTools Tools Templates

######################################

To run the looper
root
// load FWLite stuff
.x setup.C
.L doAll.C++
doAllCombined(cutBits)
// this only seems to work in compiled doAll case:
echo -e ".x setup.C \n .L doAll.C++ \n doAllCombined(1957888)" root -l -b
   The baseline value for the bitmask (as in TOP-09-002) is : 1957888 = 2**10 + 2**13 + 2**14 + 2**15 + 2**18 + 2**19 + 2**20  
   See doAll.C for info on what those bits mean
//this will make myHist_cutBits_aHumanReadableString.root
-------------------------

To make hisotgrams:

root -l -b -q "histscripts/makeAllPlots.C(\"myHist_cutBits.root\")" # will make linear-scale plots in out/ directory
root -l -b -q "histscripts/makeAllPlots.C(\"myHist_cutBits.root\", true)" # will make log-scale plots in out/directory
root -l -b -q "histscripts/makeAllPlots.C(\"myHist_cutBits.root\", 0,0,0,1,2,1,0,0)" # will order and color hists as in TOP-09-002
#there are more flags to makeAllPlots ....

---------------------------

