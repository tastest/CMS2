#!/bin/bash

#
# set up locations needed by the script
#

# ucsd or cern
GATHER_SITE=$1
# location of the gathering tools
GATHER_TOOL=$2
# list of samples in input location to process
GATHER_SAMPLE=$3
# input ntuple prefix
NTUPLE_PFX=$4

#
# check input exists
#

if [ ! $# -eq 4 ]; then
    echo "USAGE: ./makeGatherPromptData.sh GATHER_SITE GATHER_TOOL GATHER_SAMPLE 
    GATHER_SITE - UCSD or CERN: e.g. CERN
    GATHER_TOOL - location of gathering scripts: e.g. /home/username/gathering/CMS2/HuntGather2011/babymaker/
    GATHER_SAMPLE - sample: e.g. Electron_Run2010B-Nov4ReReco_v1_RECO/V03-06-16/diLepPt1020Skim/
    NTUPLE_PFX - ntuple prexif: e.g. if ntuples called merged_ntuple_run_n.root then 'merged_ntuple'"
    exit 1
fi

INPUT_DATA_DIR=""
OUTPUT_DATA_DIR=""

export SCRAM_ARCH=slc5_amd64_gcc434
    source /code/osgcode/cmssoft/cms/cmsset_default.sh
    eval `scramv1 runtime -sh`

if [ "$GATHER_SITE" == "UCSD" ]; then
    #INPUT_DATA_DIR=/nfs-3/userdata/cms2/$GATHER_SAMPLE
    INPUT_DATA_DIR=/hadoop/cms/store/user/yanjuntu/$GATHER_SAMPLE
    OUTPUT_DATA_DIR=/nfs-3/userdata/cms2/gather/data/$GATHER_SAMPLE
    #export ROOTSYS=/code/osgcode/UCSD_root/root_v5.24.00


elif [ "$GATHER_SITE" == "CERN" ]; then
    INPUT_DATA_DIR=/tas/cms2/$GATHER_SAMPLE
    OUTPUT_DATA_DIR=/tas/cms2/gather/data/$GATHER_SAMPLE
else 
    echo "ERROR: GATHER_SITE either UCSD or CERN"
    exit 1
fi

#export PATH=$ROOTSYS/bin:$PATH
#export LD_LIBRARY_PATH=$ROOTSYS/lib

#
# this script is now going to run forever
#

while [ 1 ]
do

    # check that the baby output directory exists.  If not create it
    [ ! -d "$OUTPUT_DATA_DIR" ] && echo Create $OUTPUT_DATA_DIR && mkdir -p $OUTPUT_DATA_DIR

    # all the runs/ntuples available for this dataset will be re listed
    # so start by clearing it
    touch   $OUTPUT_DATA_DIR/AllRunsAvailable.txt
    rm      $OUTPUT_DATA_DIR/AllRunsAvailable.txt
    touch   $OUTPUT_DATA_DIR/AllRunsAvailable.txt

    # likewise the list of what needs to be processed because it is new
    touch   $OUTPUT_DATA_DIR/RunsToProcess.txt
    rm      $OUTPUT_DATA_DIR/RunsToProcess.txt
    touch   $OUTPUT_DATA_DIR/RunsProcessed.txt

    # get time stamp
    DATE=`date +%Y%m%d%H%M%S`
    echo $DATE

    # get the list of all available files and write them to AllRunsAvailable.txt
    echo "Checking for new files to process in $INPUT_DATA_DIR"
    find $INPUT_DATA_DIR -follow -maxdepth 1 -name ""$NTUPLE_PFX"_[0-9]*.root" -printf "%p %s %C@\n" >>  $OUTPUT_DATA_DIR/AllRunsAvailable.txt

    # loop over all AllRunsAvailable.txt, check RunsProcessed.txt,
    # if not there, add to RunsToProcess.txt
    while read line
    do
        grep "$line" $OUTPUT_DATA_DIR/RunsProcessed.txt >/dev/null 2>&1
        if [ $? -ne 0 ]
        then
            #run=`echo $line | awk -F'_' '{print $5}'`
            run=`echo $line | sed 's/.*ntuple_\(.*\)_.*\.root.*/\1/'`
            echo $run
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
        sleep 3600

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
            CMD="root -b -l -q '$GATHER_TOOL/makeGatherBaby.C(\"$INPUT_FILE\", \"$OUTPUT_FILE\")'"
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

