int doAll() {

  gROOT->ProcessLine(".L TBitSet.cc+");
  gROOT->ProcessLine(".L ScanChain.C+");
  gROOT->ProcessLine(".L ../Tools/MiniFWLite/libMiniFWLite.so");


  double Lumi = 1000.0;
  
  //
  // Flags for files to run over 
  // (0 and 1 are easier to modify)
  //
  bool runHWW130 = 0;
  bool runHWW160 = 1;
  bool runHWW200 = 0;
  bool runWW    = 1;
  bool runWZ    = 0;
  bool runZZ    = 0;
  bool runWjets = 0;
  bool runWgamma = 0;
  bool runDYee  = 0;
  bool runDYmm  = 0;
  bool runDYtt  = 0;
  bool runttbar = 0;
  bool runtW    = 0;
  bool runQCD   = 0; 
  bool runData  = 0;
  
  if (runHWW160)
    ProcessSample("ggH160", "data/GluGluToHToWWTo2L2Nu_M-160_7TeV-powheg-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-18/*.root", -1 ,Lumi, 0.8664429*4.0/9.0, -1); 
  
  if (runHWW130)
    ProcessSample("ggH130","data/GluGluToHToWWTo2L2Nu_M-130_7TeV-powheg-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-18/*.root", -1 ,Lumi, 0.4520899*4.0/9.0, -1); 
  
  if (runHWW200) 
    ProcessSample("ggH200","data/GluGluToHToWWTo2L2Nu_M-200_7TeV-powheg-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-18/*.root", -1 ,Lumi, 0.4083049*4.0/9.0, -1); 


  if (runWW) {
    //ch->Add("data/WWTo2L2Nu_TuneZ2_7TeV-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-14/merged*.root");
    ProcessSample("WW", "data/VVJetsTo4L_TuneD6T_7TeV-madgraph-tauola_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/diLepPt2020/merged*.root", -1, Lumi, 4.5*0.919, 963356*(682015./963356.), "", false, true);
  }
  
  if (runttbar) {
    ch->Add("data/TTJets_TuneZ2_7TeV-madgraph-tauola_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v3/V03-06-18/diLepPt2020/*.root");
    ProcessSample("ttbar", ch, -1 ,Lumi, 157.5, -1); 
  } 
  if  (runWjets) {
    ch->Add("data//WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-18/diLepPt2020/merged_ntuple*root");
    ProcessSample("Wjets", ch, NprocEvents ,Lumi, 31314.0, -1); 
  }

  
  if (runData) {
    ch->Add("data/Mu_Run2010B-PromptReco-v2_RECO/V03-06-14/diLepPt1020Skim/*.root");
    ch->Add("data/Electron_Run2010B-PromptReco-v2_RECO/V03-06-14/diLepPt1020Skim/*.root");
    ProcessSample("data", ch, NprocEvents ,Lumi, -1, -1); 
  }
  
  if (runDYmm) {
    ch->Add("data/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/diLepPt2020/*.root");
    ProcessSample("Zmm", ch, -1 ,Lumi, 1666, -1); 
  }

  if (runDYee) {
    ch->Add("data/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/diLepPt2020/*.root");
    ProcessSample("Zee", ch, -1 ,Lumi, 1666, -1); 
  }

  if (runDYtt) 
    ProcessSample("Ztt", "data/DYToTauTau_M-20_TuneZ2_7TeV-pythia6-tauola_Fall10-START38_V12-v1/V03-06-14/*.root",  -1 ,Lumi, 1666, -1); 
  
  if (runWZ) 
    ProcessSample("WZ", "data/WZtoAnything_TuneZ2_7TeV-pythia6-tauola_Fall10-START38_V12-v1/V03-06-14/*.root",  -1 ,Lumi, 18.2, -1); 
  
  
  if (runZZ) 
    ProcessSample("ZZ","data/ZZtoAnything_TuneZ2_7TeV-pythia6-tauola_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/diLepPt2020/*.root", -1 ,Lumi, 5.9, -1); 
  
  if (runWgamma)
    ProcessSample("Vg", "data/PhotonVJets_7TeV-madgraph_Fall10-START38_V12-v1/CMS2_V03-06-17/*.root",  -1 ,Lumi, 165.0, -1); 
	      
}
