#!/bin/bash
# how to run: source allSkim.sh &
# files need to be changed: allSkim.sh and ntupleFilter.cc 
# files don't need change: makeSkim.C and checkAndRunSkim.sh
 
# check the existance of CORE
[ ! -d "CORE" ] && echo CORE is missing ! Please check out from CVS

# the following are the input datasets, you can add what you need to skim for
dataset_names="`echo DoubleElectron_Run2011A-05Aug2011-v1_AOD/V04-02-30/DoubleElectronTriggerSkim`
    `echo DoubleMu_Run2011A-05Aug2011-v1_AOD/V04-02-30/DoubleMuTriggerSkim`
    `echo DoubleElectron_Run2011A-PromptReco-v6_AOD/V04-02-30/DoubleElectronTriggerSkim`
    `echo DoubleMu_Run2011A-PromptReco-v6_AOD/V04-02-30/DoubleMuTriggerSkim`
"
# you need to specify input_dir (where the datasets locate) and output_dir (where the analysis skims go)
# The skimming is run in a looper, which is awoke every 7200 seconds. This deals with the rolling-in data files. You can only stop the looper when 
# you are sure that the dataset is complete and no more files are rolling in. 

while [ 1 ]
  do
  for dataset_name in $dataset_names; do
    input_dir=`echo /nfs-6/userdata/cms2/$dataset_name `
    output_dir=`echo /nfs-6/userdata/yanjuntu/TPrimeSkim/$dataset_name`
    [ ! -d "${output_dir}" ] && echo Create ${output_dir} && mkdir ${output_dir}
    echo skimming $dataset_name from $input_dir to $output_dir
    source checkAndRunSkim.sh $dataset_name $input_dir $output_dir 
  done
  sleep 7200 
done


