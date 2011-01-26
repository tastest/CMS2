
#!/bin/bash

# base location of input ntuples
GATHER_DATA_INPUT="/tas/cms2"

# base location of output babies
GATHER_DATA_OUTPUT="/tas/cms2/gather/mc"

# list of samples in input location to process
GATHER_DATA_SAMPLES="

"

# loop on input samples
for SAMPLE in $GATHER_DATA_SAMPLES;
do

    # make a directory to hold the babies for this samples
    mkdir -p $GATHER_DATA_OUTPUT/$SAMPLE

    # construct the input and output file name for making these babies
    FILE_IN=$GATHER_DATA_INPUT/$SAMPLE/*.root
    FILE_OUT=$GATHER_DATA_OUTPUT/$SAMPLE/baby_gather.root

    # run root to make the baby
    CMD="root -b -l -q 'makeGatherBaby.C(\"$FILE_IN\",\"$FILE_OUT\")'"
    eval $CMD

done

