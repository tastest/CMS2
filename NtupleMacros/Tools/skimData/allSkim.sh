#!/bin/bash
# how to run: source allSkim.sh &
# files need to be changed: allSkim.sh and ntupleFilter.cc 
# files don't need change: makeSkim.C and checkAndRunSkim.sh

# you need to specify input_dir (where the datasets locate), output_dir (where the analysis skims go), pattern of input file name
miniFWlib=`echo ~/MiniFWlib/libMiniFWLite_5.27.06b-cms10.so`
skim_C=`echo /home/users/yanjuntu/CMS2011/NtupleMacros/skimData/ntupleFilter.cc`
input_dir=`echo /nfs-6/userdata/cms2`
output_dir=`echo /nfs-7/userdata/yanjuntu/test` 
input_pattern=`echo skimmed\*.root`

# check the existance of CORE and miniFWlib
[ ! -d "CORE" ] && echo CORE is missing ! Please check out from CVS && exit 
[ ! -f "CORE/CMS2.h" ] && echo CORE/CMS2.h is missing ! && exit 
[ ! -f "CORE/CMS2.cc" ] && echo CORE/CMS2.C is missing ! && exit 
[ ! -f "${miniFWlib}" ] && miniFWlib is missing! && exit
[ ! -f "${skim_C}" ] && skim Macro is missing! && exit

# the following are the input datasets, you can add what you need to skim for
#dataset_names="`echo DoubleElectron_Run2011A-05Aug2011-v1_AOD/V04-02-30/DoubleElectronTriggerSkim`
#   `echo DoubleMu_Run2011A-05Aug2011-v1_AOD/V04-02-30/DoubleMuTriggerSkim`
#    `echo DoubleElectron_Run2011A-PromptReco-v6_AOD/V04-02-30/DoubleElectronTriggerSkim`
#    `echo DoubleMu_Run2011A-PromptReco-v6_AOD/V04-02-30/DoubleMuTriggerSkim`
#"
dataset_names="`echo DoubleElectron_Run2011A-05Aug2011-v1_AOD/V04-02-30/DoubleElectronTriggerSkim` "


# The skimming is run in a looper, which is awoke every 7200 seconds. This deals with the rolling-in data files. You can only stop the looper when 
# you are sure that the dataset is complete and no more files are rolling in. 

#while [ 1 ]
 # do
  for dataset_name in $dataset_names; do
    input_dir_sub=`echo ${input_dir}/${dataset_name} `
    output_dir_sub=`echo ${output_dir}/${dataset_name}`
    [ ! -d "${output_dir_sub}" ] && echo Create ${output_dir_sub} && mkdir -p ${output_dir_sub}
    echo skimming $dataset_name from $input_dir to $output_dir
    source checkAndRunSkim.sh $dataset_name $input_dir_sub $output_dir_sub $input_pattern $miniFWlib $skim_C 
  done
  #sleep 7200 
#done


