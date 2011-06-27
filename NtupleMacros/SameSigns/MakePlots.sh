#!/bin/bash
root -b -q MakePlots.C\(\"SUMET_10.root\"\)
cd SUMET_10.root_out
../webify_pngfiles.sh .
cd ..

root -b -q MakePlots.C\(\"SUMET_1.root\"\)
cd SUMET_1.root_out
../webify_pngfiles.sh .
cd ..

root -b -q MakePlots.C\(\"MET_10.root\"\)
cd MET_10.root_out
../webify_pngfiles.sh .
cd ..

root -b -q MakePlots.C\(\"MET_1.root\"\)
cd MET_1.root_out
../webify_pngfiles.sh .
cd ..

root -b -q MakePlots.C\(\"BaseLine.root\"\)
cd BaseLine.root_out
../webify_pngfiles.sh .
cd ..
