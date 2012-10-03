#!/bin/bash

#
# YJ
#

DATASETS="DoubleMu_Run2012C-PromptReco-v2_AOD \
DoubleElectron_Run2012C-PromptReco-v2_AOD \
MuEG_Run2012C-PromptReco-v2_AOD \
DoubleMu_Run2012C-24Aug2012-v1_AOD \
MuEG_Run2012C-24Aug2012-v1_AOD \
DoubleMu_Run2012B-13Jul2012-v4_AOD \
DoubleElectron_Run2012B-13Jul2012-v1_AOD \
MuEG_Run2012B-13Jul2012-v1_AOD \
DoubleMu_Run2012A-13Jul2012-v1_AOD \
DoubleElectron_Run2012A-13Jul2012-v1_AOD \
MuEG_Run2012A-13Jul2012-v1_AOD \
DoubleElectron_Run2012A-recover-06Aug2012-v1_AOD \
DoubleMu_Run2012A-recover-06Aug2012-v1_AOD \
MuEG_Run2012A-recover-06Aug2012-v1_AOD"

for DATASET in ${DATASETS}; do
    echo ./batchProcessWithCrab.sh ${DATASET} data.root index_${DATASET}.txt \"0\"
    find  /hadoop/cms/store/user/yanjuntu/CMSSW_5_3_2_patch4_V05-03-13/${DATASET}/merged  -type f -print \
        | sed 's/\/hadoop\/cms/root:\/\/xrootd.unl.edu\//' > index_${DATASET}.txt
done

#
# Jae
#

DATASETS="SingleMu_Run2012C-PromptReco-v2_AOD \
SingleMu_Run2012C-24Aug2012-v1_AOD \
SingleMu_Run2012B-13Jul2012-v1_AOD \
SingleMu_Run2012A-13Jul2012-v1_AOD \
SingleMu_Run2012A-recover-06Aug2012-v1_AOD"

for DATASET in ${DATASETS}; do
    echo ./batchProcessWithCrab.sh ${DATASET} data.root index_${DATASET}.txt \"0\"
    find  /hadoop/cms/store/user/jaehyeok/CMSSW_5_3_2_patch4_V05-03-13/${DATASET}/merged  -type f -print \
        | sed 's/\/hadoop\/cms/root:\/\/xrootd.unl.edu\//' > index_${DATASET}.txt
done

#
# Vince
#

DATASETS="SingleElectron_Run2012C-24Aug2012-v1_AOD \
SingleElectron_Run2012C-PromptReco-v2_AOD \
SingleElectron_Run2012B-13Jul2012-v1_AOD \
SingleElectron_Run2012A-13Jul2012-v1_AOD \
SingleElectron_Run2012A-recover-06Aug2012-v1_AOD"

for DATASET in ${DATASETS}; do
    echo ./batchProcessWithCrab.sh ${DATASET} data.root index_${DATASET}.txt \"0\"
    find  /hadoop/cms/store/user/cwelke/CMSSW_5_3_2_patch4_V05-03-13/${DATASET}/merged  -type f -print \
        | sed 's/\/hadoop\/cms/root:\/\/xrootd.unl.edu\//' > index_${DATASET}.txt
done


