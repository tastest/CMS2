int doAll(std::string process){

  gROOT->ProcessLine(".L TBitSet.cc+");
  gROOT->ProcessLine(".L ScanChain.C+");


  //double Lumi=1000.0;
  double Lumi = 35.5;
  int NprocEvents=-1;

  double Xsec;

  TChain *ch = new TChain("Events"); 
  
  if (process=="WW") {
    //ch->Add("data/VVJetsTo4L_TuneD6T_7TeV-madgraph-tauola_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/diLepPt2020/*.root");
    // ch->Add("data/GluGluToWWTo4L_TuneZ2_7TeV-gg2ww-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-18/*.root");
    ch->Add("data/WWTo2L2Nu_TuneZ2_7TeV-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-14/merged*.root");
    Xsec= 4.5;
  }
  else if (process=="ttbar") {
    //    ch->Add("data/TT_TuneZ2_7TeV-pythia6-tauola_Fall10-START38_V12-v1_GEN-SIM-RECO/V03-06-14/*.root");
    //    Xsec=165.0;
    // ch->Add("data/TTTo2L2Nu2B_7TeV-powheg-pythia6_Fall10-START38_V12-v1/V03-06-14/*.root");
    ch->Add("data/TTJets_TuneZ2_7TeV-madgraph-tauola_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v3/V03-06-18/diLepPt2020/*.root");
    Xsec=157.5;
  } 
  else if (process=="Wjets") {
    ch->Add("data//WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-18/diLepPt2020/merged_ntuple*root");
    Xsec=31314;
  }
  else if (process=="data") {
    ch->Add("data/Mu_Run2010B-PromptReco-v2_RECO/V03-06-14/diLepPt1020Skim/*.root");
    ch->Add("data/Electron_Run2010B-PromptReco-v2_RECO/V03-06-14/diLepPt1020Skim/*.root");
    Xsec=-1;
  }
  else if (process=="Zmm") {
    ch->Add("data/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/diLepPt2020/*.root");
    // ch->Add("data/DYToMuMu_M-20_TuneZ2_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-14/*.root");
    Xsec=1666;
  }
  else if (process=="Zee") {
    //ch->Add("data/DYToEE_M-20_TuneZ2_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-14/*.root");
    ch->Add("data/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/diLepPt2020/*.root");
    Xsec=1666;
  }
  else if (process=="Ztt") {
    ch->Add("data/DYToTauTau_M-20_TuneZ2_7TeV-pythia6-tauola_Fall10-START38_V12-v1/V03-06-14/*.root");
    Xsec=1666;
  }
  else if (process=="WZ") {
    ch->Add("data/WZtoAnything_TuneZ2_7TeV-pythia6-tauola_Fall10-START38_V12-v1/V03-06-14/*.root");
    Xsec=18.2;
  }
  else if (process=="ZZ") {
    ch->Add("data/ZZtoAnything_TuneZ2_7TeV-pythia6-tauola_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/diLepPt2020/*.root");
    Xsec=5.9;
  }
  else if (process=="Vg") {
    ch->Add("data/PhotonVJets_7TeV-madgraph_Fall10-START38_V12-v1/CMS2_V03-06-17/*.root");
    Xsec=165.0;
  }
  else if (process=="ggH160") {
    ch->Add("data/GluGluToHToWWTo2L2Nu_M-160_7TeV-powheg-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-18/*.root");
    Xsec=0.8664429*4.0/9.0;
  }

  else if (process=="ggH130") {
    ch->Add("data/GluGluToHToWWTo2L2Nu_M-130_7TeV-powheg-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-18/*.root");
    Xsec=0.4520899*4.0/9.0;
  }

else if (process=="ggH200") {
    ch->Add("data/GluGluToHToWWTo2L2Nu_M-200_7TeV-powheg-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-18/*.root");
    Xsec=0.4083049*4.0/9.0;
  }
  
  cout<<"Total Events Before Selection "<<ch->GetEntries()<<"\n";

   ScanChain(process, ch, NprocEvents ,Lumi, Xsec, -1); 
}
