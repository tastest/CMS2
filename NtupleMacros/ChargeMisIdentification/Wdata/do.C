{

  using namespace std;

  // Load various tools  
  gROOT->ProcessLine(".x setup.C(true)");
  gROOT->LoadMacro("ScanChain.C+");

  TChain *chain = new TChain("Events");
  // data 

  chain->Add("/nfs-3/userdata/cms2/EG_Run2010A-Jun14thReReco_v1_RECO/V03-04-26-01/singleLepPt5Skim/skimmed_ntuple_*.root");
/*
  chain->Add("/nfs-3/userdata/cms2/EG_Run2010A-PromptReco-v4_RECO/V03-04-25/singleLepPt5Skim/skimmed_*.root");
  chain->Add("/nfs-3/userdata/cms2/EG_Run2010A-PromptReco-v4_RECO/V03-04-26-01/singleLepPt5Skim/skimmed_*.root");
  chain->Add("/nfs-3/userdata/cms2/EG_Run2010A-PromptReco-v4_RECO/V03-04-26-02/singleLepPt10Skim/skimmed_ntuple_139020_0.root");
  chain->Add("/nfs-3/userdata/cms2/EG_Run2010A-PromptReco-v4_RECO/V03-04-26-02/singleLepPt10Skim/skimmed_ntuple_139085_0.root");
  chain->Add("/nfs-3/userdata/cms2/EG_Run2010A-PromptReco-v4_RECO/V03-04-26-02/singleLepPt10Skim/skimmed_ntuple_139094_0.root");
  chain->Add("/nfs-3/userdata/cms2/EG_Run2010A-PromptReco-v4_RECO/V03-04-26-02/singleLepPt10Skim/skimmed_ntuple_139096_0.root");
  chain->Add("/nfs-3/userdata/cms2/EG_Run2010A-PromptReco-v4_RECO/V03-04-26-02/singleLepPt10Skim/skimmed_ntuple_139100_0.root");
  chain->Add("/nfs-3/userdata/cms2/EG_Run2010A-PromptReco-v4_RECO/V03-04-26-02/singleLepPt10Skim/skimmed_ntuple_139096_1.root");
  chain->Add("/nfs-3/userdata/cms2/EG_Run2010A-PromptReco-v4_RECO/V03-04-26-02/singleLepPt10Skim/skimmed_ntuple_139098_0.root");
  chain->Add("/nfs-3/userdata/cms2/EG_Run2010A-PromptReco-v4_RECO/V03-04-26-02/singleLepPt10Skim/skimmed_ntuple_139100_1.root");
  chain->Add("/nfs-3/userdata/cms2/EG_Run2010A-PromptReco-v4_RECO/V03-04-26-02/singleLepPt10Skim/skimmed_ntuple_139100_2.root");
  chain->Add("/nfs-3/userdata/cms2/EG_Run2010A-PromptReco-v4_RECO/V03-04-26-02/singleLepPt10Skim/skimmed_ntuple_139102_0.root");
  chain->Add("/nfs-3/userdata/cms2/EG_Run2010A-PromptReco-v4_RECO/V03-04-26-02/singleLepPt10Skim/skimmed_ntuple_139103_0.root");
*/
  // MC
  //chain->Add("/nfs-3/userdata/cms2/Wenu_Spring10-START3X_V26_S09-v1/V03-04-08-01/merged_ntuple.root");
  //chain->Add("/nfs-3/userdata/cms2/Wenu_Spring10-START3X_V26_S09-v1/V03-04-08-01/merged_*.root");
  ScanChain(chain);

  //save all the histograms
    
  const char* outFile = "myHist.root";
  hist::saveHist(outFile);
  hist::deleteHistos();
}
