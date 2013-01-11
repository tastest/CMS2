#! /bin/bash
# preamble to set up cms things
# don't know if all of these are necessary
source /code/osgcode/cmssoft/cmsset_default.sh  > /dev/null 2>&1
export SCRAM_ARCH=slc5_amd64_gcc462

CMSSWRelease="CMSSW_5_3_2_patch4"

#checkout cmssw and setup release
scram project -n $CMSSWRelease CMSSW $CMSSWRelease
if [ $? != 0 ]; then
       echo "ERROR: Failed to check out CMSSW release $CMSSWRelease.
Exiting job without running."
               scram list
               exit 1
fi

cd $CMSSWRelease
eval `scramv1 ru -sh`
cd ..



ntupleI=$1
ntupleO=$2
rootfile=${ntupleI##*/}
rootskim=${rootfile%.*}_skim.root
skim_C=$3
libminifwlite=$4
isData=$5
echo "Skimming"
echo "In  $ntupleI"
echo "Out $ntupleO"
echo "Root $rootfile"
echo "skim $skim_C"
echo "libminifwlite $libminifwlite"
echo `ls`
echo "Host $HOSTNAME"
echo `env`
if [ ! -e $ntupleI ]; then
	echo "Error Skimming $rootfile. Job exit code 1. Cannot find file $ntupleI on hadoop."
	exit 1
fi
tar -xvzf CORE.tgz
root -b -q -l "makeSkim.C(\"$ntupleI\",\"$rootskim\",\"$isData\",\"$skim_C\",\"$libminifwlite\")"
error=$?
echo "Done Skimming $1."


## now copy everything back to hadoop
if [ $error = 0 ]; then
	lcg-cp -b -D srmv2 --vo cms -t 2400 --verbose file:`pwd`/${rootskim} srm://bsrm-1.t2.ucsd.edu:8443/srm/v2/server?SFN=$ntupleO/$rootskim
    stageout_error=$?
	if [ $stageout_error != 0 ]; then
		echo "Error skimming $rootfile. Job exit code $stageout_error. Stageout with lcg-cp failed."
	fi
else
    echo "Error skimming $rootfile. Job exit code $error. Error occurred while running makeSkime.C."
fi

if [ ! -e $rootskim ]; then
	echo "Cannot find final skimmed root file. Listing the directory for debug purposes."
	echo `ls`
	echo "Error skimming $rootfile. Job exit code 1. Cannot find output skimmed file $rootskim."
fi


## clean up the work area a little
rm *.d *.so $rootskim .root_hist



