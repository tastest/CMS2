
#!/bin/bash

# ucsd or cern
GATHER_SITE=$1
# tag name
GATHER_TAG=$2
# list of samples in input location to process
GATHER_SAMPLE=$3

#
# check input exists
#

if [ ! $# -eq 3 ]; then
    echo "USAGE: ./makeGatherMC.sh GATHER_SITE GATHER_SAMPLE
    GATHER_SITE - UCSD or CERN: e.g. CERN
    GATHER_TAG - name of subdir for babies
    GATHER_SAMPLE - sample: e.g. WZtoAnything_TuneZ2_7TeV-pythia6-tauola_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/"
    exit 1
fi

export SCRAM_ARCH=slc5_amd64_gcc434

GATHER_INPUT=""
GATHER_OUTPUT=""
if [ "$GATHER_SITE" == "HADOOP" ]; then
    GATHER_INPUT=/hadoop/cms/store/user/yanjuntu/
    GATHER_OUTPUT=/nfs-3/userdata/cms2/gather/$GATHER_TAG/
    export ROOTSYS=/code/osgcode/UCSD_root/root_v5.24.00
    source /code/osgcode/cmssoft/cms/cmsset_default.sh
elif [ "$GATHER_SITE" == "UCSD" ]; then
    GATHER_INPUT=/nfs-3/userdata/cms2/
    GATHER_OUTPUT=/nfs-3/userdata/cms2/gather/$GATHER_TAG/
    export ROOTSYS=/code/osgcode/UCSD_root/root_v5.24.00
    source /code/osgcode/cmssoft/cms/cmsset_default.sh
elif [ "$GATHER_SITE" == "CERN" ]; then
    GATHER_INPUT=/tas/cms2/
    GATHER_OUTPUT=/tas/cms2/gather/$GATHER_TAG/
    #export ROOTSYS=/afs/cern.ch/sw/lcg/app/releases/ROOT/5.26.00/x86_64-slc5-gcc34-opt/root/
    source /afs/cern.ch/cms/cmsset_default.sh

else
    echo "ERROR: GATHER_SITE either UCSD or CERN"
    exit 1
fi

    eval `scramv1 runtime -sh`

#export PATH=$ROOTSYS/bin:$PATH
#export LD_LIBRARY_PATH=$ROOTSYS/lib


#
# done with checks, now make the baby
#

# make a directory to hold the babies for this samples
mkdir -p $GATHER_OUTPUT/$GATHER_SAMPLE

# construct the input and output file name for making these babies
FILE_IN=$GATHER_INPUT/$GATHER_SAMPLE/*.root
FILE_OUT=$GATHER_OUTPUT/$GATHER_SAMPLE/baby_gather.root

# run root to make the baby
CMD="root -b -l -q 'makeGatherBaby.C(\"$FILE_IN\",\"$FILE_OUT\")'"
eval $CMD

exit 0

