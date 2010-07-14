#include "ChainFromText.cc"
void doit(){

gROOT->LoadMacro("myBabyMaker.C++");

//TChain* chain1 = new TChain("Events");
//chain1->Add("/tas07/disk00/cms2/MinimumBias_Commissioning10-May6thPDSkim2_SD_JetMETTau-v1/V03-04-09/skim/filtered_ntuple.root");
//chain1->Add("/tas07/disk00/yanjuntu/Prompt_Data_CMSSW_3_6_0_V03-04-09-01/MinimumBias_Commissioning10-SD_JetMETTau-v9_RECO/singleLepPt5Skim/*.root");
//chain1->Add("/tas07/disk00/yanjuntu/Prompt_Data_CMSSW_3_6_0_V03-04-09-01/JetMETTau_Run2010A-PromptReco-v1_RECO/singleLepPt5Skim/*.root");
//chain1->Add("/tas07/disk00/yanjuntu/Prompt_Data_CMSSW_3_6_1_V03-04-10-01/JetMETTau_Run2010A-PromptReco-v1_RECO/singleLepPt5Skim/*.root");
//chain1->Add("/tas07/disk00/yanjuntu/Prompt_Data_CMSSW_3_6_1_V03-04-10-01/JetMETTau_Run2010A-PromptReco-v2_RECO/singleLepPt5Skim/*.root");

// old
//chain1->Add("/nfs-3/userdata/cms2/MinimumBias_Commissioning10-May6thPDSkim2_SD_JetMETTau-v1/V03-04-09/skim/filtered_ntuple.root");
//chain1->Add("/hadoop/cms/store/group/snt/MinimumBias_Commissioning10-SD_JetMETTau-v9_RECO/singleLepPt5Skim/*.root");
//chain1->Add("/hadoop/cms/store/group/snt/yanjuntu/Prompt_Data_CMSSW_3_6_0_V03-04-09-01/JetMETTau_Run2010A-PromptReco-v1_RECO/singleLepPt5Skim/*.root");
//chain1->Add("/hadoop/cms/store/group/snt/yanjuntu/Prompt_Data_CMSSW_3_6_1_V03-04-10-01/JetMETTau_Run2010A-PromptReco-v1_RECO/singleLepPt5Skim/*.root");
//chain1->Add("/hadoop/cms/store/group/snt/yanjuntu/Prompt_Data_CMSSW_3_6_1_V03-04-10-01/JetMETTau_Run2010A-PromptReco-v2_RECO/singleLepPt5Skim/*.root");

// new
//chain1->Add("/hadoop/cms/store/user/spadhi/CMS2_V03-04-26-02/Commissioning10-SD_JetMETTau-Jun14thSkim_v1/*.root");
//chain1->Add("/nfs-3/userdata/cms2/JetMETTau_Run2010A-Jun14thReReco_v2_RECO/V03-04-26-01/singleLepPt5Skim/*.root");
//chain1->Add("/nfs-3/userdata/cms2/JetMETTau_Run2010A-PromptReco-v4_RECO/V03-04-25/singleLepPt5Skim/*.root");
//chain1->Add("/nfs-3/userdata/cms2/JetMETTau_Run2010A-PromptReco-v4_RECO/V03-04-26-01/singleLepPt5Skim/*.root");
//chain1->Add("/nfs-3/userdata/cms2/JetMETTau_Run2010A-PromptReco-v4_RECO/V03-04-26-02/singleLepPt10Skim/*.root");



 // 
//TChain *chain1 = ChainFromText("jmt_july6.txt");
//cout << chain1->GetEntries() << " total entries." << endl;

//These are the EG 
//TChain *chain1 = new TChain("Events");
 // chain1->Add("/tas/cms2/EG_Run2010A-Jun14thReReco_v1_RECO/V03-04-26-01/singleLepPt5Skim/*.root");
 //chain1->Add("/tas/cms2/EG_Run2010A-PromptReco-v4_RECO/V03-04-25/singleLepPt5Skim/*.root");
 //chain1->Add("/tas/cms2/EG_Run2010A-PromptReco-v4_RECO/V03-04-26-01/singleLepPt5Skim/*.root");
 //chain1->Add("/tas/cms2/EG_Run2010A-PromptReco-v4_RECO/V03-04-26-02/singleLepPt10Skim/*.root");
 //chain1->Add("/tas/cms2/MinimumBias_Commissioning10-SD_EG-Jun14thSkim_v1_RECO/V03-04-26-02/singleLepPt10Skim/*.root");

 TChain *chain1 = new TChain("Events");
 chain1->Add("/tas/cms2/JetMETTau_Run2010A-Jun14thReReco_v2_RECO/V03-04-26-01/singleLepPt5Skim/*.root");
 chain1->Add("/tas/cms2/JetMETTau_Run2010A-PromptReco-v4_RECO/V03-04-25/singleLepPt5Skim/*.root");
 chain1->Add("/tas/cms2/JetMETTau_Run2010A-PromptReco-v4_RECO/V03-04-26-01/singleLepPt5Skim/*.root");
 chain1->Add("/tas/cms2/JetMETTau_Run2010A-PromptReco-v4_RECO/V03-04-26-02/singleLepPt10Skim/*.root");
 chain1->Add("/tas/cms2/MinimumBias_Commissioning10-SD_JetMETTau-Jun14thSkim_v1_RECO/V03-04-26-02/singleLepPt10Skim/*.root");


myBabyMaker * baby1 = new myBabyMaker();
baby1->ScanChain(chain1, "JMT.root", -1);
//baby1->ScanChain(chain1, "EG.root", -1);

//TChain* chain2 = new TChain("Events");
//chain2->Add("/tas07/disk00/cms2/MinimumBias_Commissioning10-May6thPDSkim2_SD_JetMETTauMonitor-v1/V03-04-09/skim/filtered_ntuple.root");
//chain2->Add("/tas07/disk00/yanjuntu/Prompt_Data_CMSSW_3_6_0_V03-04-09-01/MinimumBias_Commissioning10-SD_JetMETTauMonitor-v9_RECO/singleLepPt5Skim/*.root");
//chain2->Add("/tas07/disk00/yanjuntu/Prompt_Data_CMSSW_3_6_0_V03-04-09-01/JetMETTauMonitor_Run2010A-PromptReco-v1_RECO/singleLepPt5Skim/*.root");
//chain2->Add("/tas07/disk00/yanjuntu/Prompt_Data_CMSSW_3_6_1_V03-04-10-01/JetMETTauMonitor_Run2010A-PromptReco-v1_RECO/singleLepPt5Skim/*.root");
//chain2->Add("/tas07/disk00/yanjuntu/Prompt_Data_CMSSW_3_6_1_V03-04-10-01/JetMETTauMonitor_Run2010A-PromptReco-v2_RECO/singleLepPt5Skim/*.root");

//old
//chain1->Add("/nfs-3/userdata/cms2/MinimumBias_Commissioning10-May6thPDSkim2_SD_JetMETTauMonitor-v1/V03-04-09/skim/filtered_ntuple.root");
//chain1->Add("/hadoop/cms/store/group/snt/MinimumBias_Commissioning10-SD_JetMETTauMonitor-v9_RECO/singleLepPt5Skim/*.root");
//chain1->Add("/hadoop/cms/store/group/snt/yanjuntu/Prompt_Data_CMSSW_3_6_0_V03-04-09-01/JetMETTauMonitor_Run2010A-PromptReco-v1_RECO/singleLepPt5Skim/*.root");
//chain1->Add("/hadoop/cms/store/group/snt/yanjuntu/Prompt_Data_CMSSW_3_6_1_V03-04-10-01/JetMETTauMonitor_Run2010A-PromptReco-v1_RECO/singleLepPt5Skim/*.root");
//chain1->Add("/hadoop/cms/store/group/snt/yanjuntu/Prompt_Data_CMSSW_3_6_1_V03-04-10-01/JetMETTauMonitor_Run2010A-PromptReco-v2_RECO/singleLepPt5Skim/*.root");

// new
//chain2->Add("/hadoop/cms/store/user/spadhi/CMS2_V03-04-26-02/Commissioning10-SD_JetMETTauMonitor-Jun14thSkim_v1/*.root");
//chain2->Add("/nfs-3/userdata/cms2/JetMETTauMonitor_Run2010A-Jun14thReReco_v1_RECO/V03-04-26-01/singleLepPt5Skim/*.root");
//chain2->Add("/nfs-3/userdata/cms2/JetMETTauMonitor_Run2010A-PromptReco-v4_RECO/V03-04-25/singleLepPt5Skim/*.root");
//chain2->Add("/nfs-3/userdata/cms2/JetMETTauMonitor_Run2010A-PromptReco-v4_RECO/V03-04-26-01/singleLepPt5Skim/*.root");
//chain2->Add("/nfs-3/userdata/cms2/JetMETTauMonitor_Run2010A-PromptReco-v4_RECO/V03-04-26-02/singleLepPt10Skim/*.root");

//TChain *chain2 = ChainFromText("jmtm_july6.txt");
//cout << chain2->GetEntries() << " total entries." << endl;

// These are the muons
//TChain *chain2 = new TChain("Events");
//chain2->Add("/tas/cms2/Mu_Run2010A-Jun14thReReco_v1_RECO/V03-04-26-01/singleLepPt5Skim/*.root");
//chain2->Add("/tas/cms2/Mu_Run2010A-PromptReco-v4_RECO/V03-04-25/singleLepPt5Skim/*.root");
//chain2->Add("/tas/cms2/Mu_Run2010A-PromptReco-v4_RECO/V03-04-26-01/singleLepPt5Skim/*.root");
//chain2->Add("/tas/cms2/Mu_Run2010A-PromptReco-v4_RECO/V03-04-26-02/singleLepPt10Skim/*.root");
//chain2->Add("/tas/cms2/MinimumBias_Commissioning10-SD_Mu-Jun14thSkim_v1_RECO/V03-04-26-02/singleLepPt10Skim/*.root");


 TChain *chain2 = new TChain("Events");
 chain2->Add("/tas/cms2/JetMETTauMonitor_Run2010A-Jun14thReReco_v1_RECO/V03-04-26-01/singleLepPt5Skim/*.root");
 chain2->Add("/tas/cms2/JetMETTauMonitor_Run2010A-PromptReco-v4_RECO/V03-04-25/singleLepPt5Skim/*.root");
 chain2->Add("/tas/cms2/JetMETTauMonitor_Run2010A-PromptReco-v4_RECO/V03-04-26-01/singleLepPt5Skim/*.root");
 chain2->Add("/tas/cms2/JetMETTauMonitor_Run2010A-PromptReco-v4_RECO/V03-04-26-02/singleLepPt10Skim/*.root");
 chain2->Add("/tas/cms2/MinimumBias_Commissioning10-SD_JetMETTauMonitor-Jun14thSkim_v1_RECO/V03-04-26-02/singleLepPt10Skim/*.root");


myBabyMaker * baby2 = new myBabyMaker();
baby2->ScanChain(chain2, "JMTMonitor.root", -1);
//baby2->ScanChain(chain2, "Mu.root", -1);

}
