{
gROOT->LoadMacro("myBabyMaker.C+");

TChain* chain1 = new TChain("Events");
chain1->Add("/tas07/disk00/cms2/MinimumBias_Commissioning10-May6thPDSkim2_SD_JetMETTau-v1/V03-04-09/skim/filtered_ntuple.root");
chain1->Add("/tas07/disk00/yanjuntu/Prompt_Data_CMSSW_3_6_0_V03-04-09-01/MinimumBias_Commissioning10-SD_JetMETTau-v9_RECO/singleLepPt5Skim/*.root");
chain1->Add("/tas07/disk00/yanjuntu/Prompt_Data_CMSSW_3_6_0_V03-04-09-01/JetMETTau_Run2010A-PromptReco-v1_RECO/singleLepPt5Skim/*.root");
chain1->Add("/tas07/disk00/yanjuntu/Prompt_Data_CMSSW_3_6_1_V03-04-10-01/JetMETTau_Run2010A-PromptReco-v1_RECO/singleLepPt5Skim/*.root");
chain1->Add("/tas07/disk00/yanjuntu/Prompt_Data_CMSSW_3_6_1_V03-04-10-01/JetMETTau_Run2010A-PromptReco-v2_RECO/singleLepPt5Skim/*.root");

myBabyMaker * baby1 = new myBabyMaker();
baby1->ScanChain(chain1, "JMT.root", -1);

TChain* chain2 = new TChain("Events");
chain2->Add("/tas07/disk00/cms2/MinimumBias_Commissioning10-May6thPDSkim2_SD_JetMETTauMonitor-v1/V03-04-09/skim/filtered_ntuple.root");
chain2->Add("/tas07/disk00/yanjuntu/Prompt_Data_CMSSW_3_6_0_V03-04-09-01/MinimumBias_Commissioning10-SD_JetMETTauMonitor-v9_RECO/singleLepPt5Skim/*.root");
chain2->Add("/tas07/disk00/yanjuntu/Prompt_Data_CMSSW_3_6_0_V03-04-09-01/JetMETTauMonitor_Run2010A-PromptReco-v1_RECO/singleLepPt5Skim/*.root");
chain2->Add("/tas07/disk00/yanjuntu/Prompt_Data_CMSSW_3_6_1_V03-04-10-01/JetMETTauMonitor_Run2010A-PromptReco-v1_RECO/singleLepPt5Skim/*.root");
chain2->Add("/tas07/disk00/yanjuntu/Prompt_Data_CMSSW_3_6_1_V03-04-10-01/JetMETTauMonitor_Run2010A-PromptReco-v2_RECO/singleLepPt5Skim/*.root");

myBabyMaker * baby2 = new myBabyMaker();
baby2->ScanChain(chain2, "JMTMonitor.root", -1);

}
