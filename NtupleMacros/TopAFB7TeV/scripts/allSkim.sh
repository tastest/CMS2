#!/bin/bash 
export VDT_LOCATION=/data/vdt
export EDG_WL_LOCATION=$VDT_LOCATION/edg
source /data/vdt/setup.sh 

runOption=`echo skimNew`
if [ "$runOption" == "skimMC" ]; then 
    dataset_names=" `echo DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29`
   `echo DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Summer11-PU_S4_START42_V11-v1/V04-02-29`
   `echo DYToEE_M-10To20_CT10_TuneZ2_7TeV-powheg-pythia_Summer11-PU_S4_START42_V11-v1/V04-02-29`   
   `echo DYToMuMu_M-10To20_CT10_TuneZ2_7TeV-powheg-pythia_Summer11-PU_S4_START42_V11-v1/V04-02-29`
   `echo DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Summer11-PU_S4_START42_V11-v1/V04-02-29`
   `echo DYToTauTau_M-10To20_TuneZ2_7TeV-pythia6-tauola_Summer11-PU_S3_START42_V11-v2/V04-02-29`
   `echo DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29`   
   `echo WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29_singleLepton`
   `echo T_TuneZ2_s-channel_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29`
   `echo Tbar_TuneZ2_s-channel_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29`
   `echo T_TuneZ2_t-channel_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29`
   `echo Tbar_TuneZ2_t-channel_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29`
   `echo T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29`
   `echo Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29`
   `echo WWJetsTo2L2Nu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29`
   `echo WZJetsTo2L2Q_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29`
   `echo WZJetsTo3LNu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29`
   `echo ZZJetsTo2L2Nu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29`
   `echo ZZJetsTo2L2Q_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29` 
   `echo ZZJetsTo4L_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29`
"
 # dataset_names=" 
 #  `echo WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29_singleLepton` 
 #  `echo ZZ_TuneZ2_7TeV_pythia6_tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29`                                                                                                                  
 #  `echo WWTo2L2Nu_TuneZ2_7TeV_pythia6_tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29`                                                                                                           
 #   `echo WZ_TuneZ2_7TeV_pythia6_tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29`   
 #   `echo WWJetsTo2L2Nu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29`                                                                                                     
 #   `echo WZJetsTo2L2Q_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29`                                                                                                      
 #   `echo ZZJetsTo2L2Nu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29`                                                                                                     
 #   `echo ZZJetsTo2L2Q_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29`  
 #"
 #backup MC samples
 #   `echo TTJets_TuneZ2_7TeV-madgraph-tauola_Fall11-PU_S6_START42_V14B-v2/V04-02-29`   
 #  `echo TTJets_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29` 
 #  `echo TTTo2L2Nu2B_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1`
 #  `echo ZZ_TuneZ2_7TeV_pythia6_tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29`
 #  `echo WWTo2L2Nu_TuneZ2_7TeV_pythia6_tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29`
 #   `echo WZ_TuneZ2_7TeV_pythia6_tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29`
elif [ "$runOption" == "skimData" ]; then
    dataset_names="`echo DoubleElectron_Run2011A-05Aug2011-v1_AOD/V04-02-30/DoubleElectronTriggerSkim`
    `echo DoubleMu_Run2011A-05Aug2011-v1_AOD/V04-02-30/DoubleMuTriggerSkim`
    `echo DoubleElectron_Run2011A-PromptReco-v6_AOD/V04-02-30/DoubleElectronTriggerSkim`
    `echo DoubleMu_Run2011A-PromptReco-v6_AOD/V04-02-30/DoubleMuTriggerSkim`
"
elif [ "$runOption" == "skimHadoopData" ]; then
    dataset_names="`echo MuEG_Run2011A-05Aug2011-v1_AOD ` 
     `echo MuEG_Run2011A-PromptReco-v6_AOD `
"
elif [ "$runOption" == "skimHadoopMC" ]; then
    #dataset_names="`echo SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1 ` "
    #dataset_names="`echo WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1`"
    dataset_names="`echo TT_TuneZ2_7TeV-mcatnlo_Fall11-PU_S6_START42_V14B-v1`"
elif [ "$runOption" == "skimMuEG" ]; then
    dataset_names="`echo MuEG_Run2011A-May10ReReco-v1_AOD/V04-02-15 ` "
elif [ "$runOption" == "skimNew" ]; then
    dataset_names=" `echo  ZZJetsTo2L2Q_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLepton ` "
#    dataset_names="`echo ZZJetsTo4L_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29`

#  `echo WZJetsTo3LNu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29`
   # dataset_names="`echo TT_TuneZ2_7TeV-mcatnlo_Fall11-PU_S6_START42_V14B-v1/V04-02-29_fix_dilepton`"
   # dataset_names="`echo DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Summer11-PU_S4_START42_V11-v1/V04-02-29`"
 #"
    #dataset_names="`echo TT_TuneZ2_7TeV-mcatnlo_Fall11-PU_S6_START42_V14B-v1_genfix/V04-02-29`"
else 
    echo failed to define runOption && exit 25
fi



#while [ 1 ]                                                                                                                                  
 # do 

for dataset_name in $dataset_names; do
    if [ "$runOption" == "skimMC" ]; then 
	#input_dir=`echo /nfs-7/userdata/cms2/$dataset_name`
	input_dir=`echo /hadoop/cms/store/group/snt/papers2011/Summer11MC/$dataset_name`
	pattern=`echo merged\*.root `
    elif [ "$runOption" == "skimData" ]; then
	input_dir=`echo /nfs-6/userdata/cms2/$dataset_name `
	pattern=`echo skimmed\*.root ` 
    elif [ "$runOption" == "skimHadoopData" ]; then
	input_dir=`echo /hadoop/cms/store/user/yanjuntu/CMSSW_4_2_7_patch1_V04-02-30/$dataset_name/CMSSW_4_2_7_patch1_V04-02-30_merged/V04-02-30 `
	pattern=`echo merged\*.root`
    elif [ "$runOption" == "skimHadoopMC" ]; then
        #input_dir=`echo /hadoop/cms/store/user/benhoob/CMS2_V04-02-20-04/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1 `
	#input_dir=`echo /hadoop/cms/store/user/imacneill/Summer11MC/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/DileptonHyp`
	input_dir=`echo /nfs-6/userdata/cms2/TT_TuneZ2_7TeV-mcatnlo_Fall11-PU_S6_START42_V14B-v1/V04-02-29/preprocessing`
        pattern=`echo ntuple\*.root`
    elif [ "$runOption" == "skimMuEG" ]; then
        input_dir=`echo /nfs-4/userdata/cms2/$dataset_name `
        pattern=`echo merged\*.root`
    elif [ "$runOption" == "skimNew" ]; then
	#input_dir=`echo /nfs-7/userdata/cms2/$dataset_name`
	input_dir=`echo /hadoop/cms/store/group/snt/papers2011/Summer11MC/$dataset_name`
	pattern=`echo merged\*.root `
    fi
    
    output_dir=`echo /nfs-6/userdata/yanjuntu/AfbSkimv2/$dataset_name`
    echo skimming $dataset_name from $input_dir to $output_dir
    source checkAndRunSkim.sh $dataset_name $input_dir $output_dir $pattern
  done
 # sleep 7201 
#done


