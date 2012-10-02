#!/bin/bash

#
# Write a crab config and a wrapper script
# to submit a looper job
#


if [ ! $# -eq 4 ]; then
    echo "
USAGE: ./batchProcessWithCrab.sh TASK OUTFILE INDEX ARGS
    TASK    - Unique name for this task
    OUTFILE - Name of outfile 
    INDEX   - Input file index card
    ARGS    - Argument values for wrapper script (comma separated)
"
    exit 1
fi

TASK=$1
OUTFILE=$2
INDEX=$3
ARGS=$4
CRABCFG=crab_${TASK}.cfg
WRAPPER=wrapper_${TASK}.sh
LOOPER=looper.tar.gz

#
# This is the crab configuration
# to sumbit the looper job
#

NJOBS=`cat ${INDEX} | grep root | wc -l`
cat > ${CRABCFG} << EOF
[CRAB]
jobtype   = cmssw
scheduler = condor

[CMSSW]
datasetpath             = None
pset                    = None
output_file             = ${OUTFILE}
### not used but must be specified
events_per_job          = 1414 
### this is used and determined from index...
number_of_jobs          = ${NJOBS}

[USER]
script_exe              = ${WRAPPER}
script_arguments        = ${ARGS}
return_data             = 1
ui_working_dir          = ${TASK}
additional_input_files  = ${LOOPER},${INDEX}
EOF

#
# make the tar file to run the looper
# include everything needed to run looper
#

tar -chzf ${LOOPER} files/ *.so processData.exe

#
# This is the wrapper that will run
# the looper on the remote WN
# - the first argument is always the job index
# - the latter arguments are as provided in the crab cfg
#

cat > ${WRAPPER} << EOF
#!/bin/bash

tar -zxf ${LOOPER}

#
# args
#

JobIndex=\$1
DataType=\$2
File=\`awk 'FNR=='\${JobIndex} ${INDEX}\`

echo "[wrapper] JobIndex    = " \${JobIndex}
echo "[wrapper] DataType    = " \${DataType}
echo "[wrapper] File        = " \${File}

#
# run it
#

./processData.exe \${DataType} \${File}

# outfile must be in working directory
# if not already
mv smurf_0_999999/${OUTFILE} ${OUTFILE}
EOF
chmod +x ${WRAPPER}

