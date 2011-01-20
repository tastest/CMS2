
#!/bin/bash

# base location of input ntuples
GATHER_MC_INPUT="/tas/cms2"

# base location of output babies
GATHER_MC_OUTPUT="/tas/cms2/gather/mc"

# list of samples in input location to process
GATHER_MC_SAMPLES="
WZtoAnything_TuneZ2_7TeV-pythia6-tauola_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/"
#ZZtoAnything_TuneZ2_7TeV-pythia6-tauola_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/"

# loop on input samples
for SAMPLE in $GATHER_MC_SAMPLES;
do

    # make a directory to hold the babies for this samples
    mkdir -p $GATHER_MC_OUTPUT/$SAMPLE

    # construct the input and output file name for making these babies
    FILE_IN=$GATHER_MC_INPUT/$SAMPLE/*.root
    FILE_OUT=$GATHER_MC_OUTPUT/$SAMPLE/baby_gather.root

    # run root to make the baby
    CMD="root -b -l -q 'makeGatherBaby.C(\"$FILE_IN\",\"$FILE_OUT\")'"
    eval $CMD

done


