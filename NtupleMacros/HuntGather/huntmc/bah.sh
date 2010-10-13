#!/bin/bash

for txt in `ls *.txt`;
do
    sample=`echo $txt | cut -d'.' -f 1`
    mkdir -p $sample/emu_skim $sample/emu_baby $sample/dilep_skim $sample/dilep_baby $sample/trilep_skim $sample/trilep_baby $sample/logs

    while read file
    do
        cmd="root -b -l -q 'UpdateSkims.C+(\"$sample\",\"$file\")'"
        eval $cmd
    done <$sample.txt
done
