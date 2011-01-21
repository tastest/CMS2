#!/bin/bash

#
# set up locations needed by the script
#

# location of the gathering scripts
TOOL_DIR=$1
# name of the dataset to process
DATASET=$2
# CMS2 tag of the dataset to process
CMS2_TAG=$3

# location of the ntuples
BASE=/tas/cms2

# full location of the input ntuples
# and the output babies
INPUT_DATA_DIR=$BASE/$DATASET/$CMS2_TAG
OUTPUT_DATA_DIR=$BASE/gather/data/$DATASET/$CMS2_TAG

# location of ROOT
export ROOTSYS=/afs/cern.ch/sw/lcg/app/releases/ROOT/5.26.00/x86_64-slc5-gcc34-opt/root/
export PATH=$ROOTSYS/bin:$PATH
export LD_LIBRARY_PATH=$ROOTSYS/lib

#
# this script is now going to run forever
#

while [ 1 ]
do

    # check that the baby output directory exists.  If not create it
    [ ! -d "$OUTPUT_DATA_DIR" ] && echo Create $OUTPUT_DATA_DIR && mkdir -p $OUTPUT_DATA_DIR

    # all the runs/ntuples available for this dataset will be re listed
    # so start by clearing it
    rm      $OUTPUT_DATA_DIR/AllRunsAvailable.txt
    touch   $OUTPUT_DATA_DIR/AllRunsAvailable.txt

    # likewise the list of what needs to be processed because it is new
    rm      $OUTPUT_DATA_DIR/RunsToProcess.txt
    touch   $OUTPUT_DATA_DIR/RunsProcessed.txt

    # get time stamp
    DATE=`date +%Y%m%d%H%M%S`
    echo $DATE

    # get the list of all available files and write them to AllRunsAvailable.txt
    echo "Checking for new files to process in $INPUT_DATA_DIR"
    find $INPUT_DATA_DIR -follow -maxdepth 1 -name "merged_ntuple_[0-9]*.root" -printf "%p %s %C@\n"  >>  $OUTPUT_DATA_DIR/AllRunsAvailable.txt

    # loop over all AllRunsAvailable.txt, check RunsProcessed.txt,
    # if not there, add to RunsToProcess.txt
    while read line
    do
        grep "$line" $OUTPUT_DATA_DIR/RunsProcessed.txt >/dev/null 2>&1
        if [ $? -ne 0 ]
        then
            run=`echo $line | awk -F'_' '{print $5}'`
            if [ $run -gt 133222 ] && [ $run -le 133250 ];
            then
                echo "Magnetic field off, skipping run "$run
                else 
                    echo $line >> $OUTPUT_DATA_DIR/RunsToProcess.txt
            fi
        fi
    done < $OUTPUT_DATA_DIR/AllRunsAvailable.txt

    #
    # check if there are new files and process them gods willing
    #

    N_FILES=`cat $OUTPUT_DATA_DIR/RunsToProcess.txt | wc -l`
    if [ $N_FILES -lt 1 ]
    then
        echo "No new runs to process..."
    else

        # there are new files to process
        # so loop on the list of runs to process
        while read line
        do

            # input file
            INPUT_FILE=`echo $line | awk '{print $1}'`

            # output file
            NTUPLE_NAME=`echo $INPUT_FILE | awk -F'/' '{print $NF}'`
            OUTPUT_FILE=$OUTPUT_DATA_DIR/baby_gather_$NTUPLE_NAME

            # run the baby maker for this ntuple
            CMD="root -b -l -q '$TOOL_DIR/makeGatherBaby.C(\"$INPUT_FILE\", \"$OUTPUT_FILE\")'"
            eval $CMD

            # check the baby maker ran correctly
            if [ $? -ne 0 ]
            then
                echo "Error processing $FILE"
                exit 55

            # if it ran correctly then write this run in the list of runs 
            # that have been successfully processed
            else
                echo $line >> $OUTPUT_DATA_DIR/RunsProcessed.txt
            fi

        # finish loop on the list of runs to process
        done < $OUTPUT_DATA_DIR/RunsToProcess.txt

        echo cat RunsProcessed.txt | cut -d '_' -f 3 | sort -n -r | head -n 1 > $OUTPUT_DATA_DIR/lastrun.txt
    fi

    # wait an hour and then start again
    echo sleep 3600;

# repeat while loop
done

