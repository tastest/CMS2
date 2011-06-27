#!/bin/bash
if [ "$1" = "long" ]
then
rootfiles="\"OSSUSY.root\" \"SSSUSY.root\" \"TtbarBaseline_SS_FlipObserved.root\" \"TtbarBaseline_SS_FlipPredicted.root\" \"TtbarLoose_SS_FlipObserved.root\" \"TtbarLoose_SS_FlipPredicted.root\" \"Ttbar_SS_Fakerate.root\" \"Ttbar_SS_FOs_Not_MuonNumerator.root\" \"Ttbar_SS_FOs_Not_Numerator.root\" \"Ttbar_SS_MuonFakerate.root\" \"Ttbar_SS_MuonNumerator.root\" \"Ttbar_SS_Numerator.root\" \"Wjets_Fakerate.root\" \"Wjets_FOs_Not_MuonNumerator.root\" \"Wjets_FOs_Not_Numerator.root\" \"Wjets_MuonFakerate.root\" \"Wjets_MuonNumerator.root\" \"Wjets_Numerator.root\" \"Wjets_SS_Fakerate.root\" \"Wjets_SS_FOs_Not_MuonNumerator.root\" \"Wjets_SS_FOs_Not_Numerator.root\" \"Wjets_SS_MuonFakerate.root\" \"Wjets_SS_MuonNumerator.root\" \"Wjets_SS_Numerator.root\""
echo ' LONG mode!'
else
rootfiles="\"SSSUSY.root\""
fi

for rootfile in $rootfiles
do
 echo $rootfile
 echo "Processing "$rootfile" in "$PWD
 cd $PWD
 root -b -q /nfs-1/users/ibloch/MakePlots.C\($rootfile\)
 rootfile_clean=`echo $rootfile | sed s+"\""+""+g`
 dir=$rootfile_clean"_out"
# echo $dir
 cd $dir
 /nfs-1/users/ibloch/webify_pngfiles.sh .
 cd ..
done
