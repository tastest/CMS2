#include "TChain.h"
#include "looper.C"

void runLooper(char* prefix){

  TChain* ch = new TChain("Events");
  bool isData = false;

  //const char* path = "/nfs-3/userdata/cms2";
  const char* path = "/tas/cms2";
  

  if( TString(prefix).Contains("WToLNu") ){
    cout << "Running pythia wjets sample" << endl;
    ch->Add(Form("%s/WToENu_TuneZ2_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-14/merged_ntuple*.root",path));
    ch->Add(Form("%s/WToMuNu_TuneZ2_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-14/merged_ntuple*.root",path));
    ch->Add(Form("%s/WToTauNu_TuneZ2_7TeV-pythia6-tauola_Fall10-START38_V12-v1/V03-06-14/merged_ntuple*.root",path));
  }

  else if( TString(prefix).Contains("WJetsToLNu_PU") ){
    cout << "Running madgraph wjets sample" << endl;
    ch->Add(Form("%s/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-18/merged_ntuple*.root",path));
  }

  else if( strcmp( prefix , "WJetsToLNu" ) == 0 ){
    cout << "Running madgraph wjets sample, no pileup" << endl;
    ch->Add(Form("%s/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Fall10-START38_V12-v1/V03-06-18/merged_ntuple*.root",path));
  }

  else if( strcmp( prefix , "WWTo2L2Nu" ) == 0 ){
    ch->Add(Form("%s/WWTo2L2Nu_TuneZ2_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-14/merged_ntuple*.root",path));
  }

  else if( strcmp( prefix , "WWTo2L2Nu_PU" ) == 0 ){
    ch->Add(Form("%s/WWTo2L2Nu_TuneZ2_7TeV-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-14/merged_ntuple*.root",path));
  }

  else if( strcmp( prefix , "GluGluToWWTo4L_PU" ) == 0 ){
    ch->Add(Form("%s/GluGluToWWTo4L_TuneZ2_7TeV-gg2ww-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-18/merged_ntuple*.root",path));
  }

  else if( strcmp( prefix , "VVJetsTo4L_PU" ) == 0 ){
    ch->Add(Form("%s/VVJetsTo4L_TuneD6T_7TeV-madgraph-tauola_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/merged_ntuple*.root",path));
  }

  else if( strcmp( prefix , "TTJets" ) == 0 ){
    ch->Add(Form("%s/TTJets_TuneZ2_7TeV-madgraph-tauola_Fall10-START38_V12-v3/V03-06-18/merged_ntuple*.root",path));
  }

  else if( strcmp( prefix , "TTJets_PU" ) == 0 ){
    ch->Add(Form("%s/TTJets_TuneZ2_7TeV-madgraph-tauola_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v3/V03-06-18/merged_ntuple*.root",path));
  }

  else if( strcmp( prefix , "tW_PU" ) == 0 ){
    ch->Add(Form("%s/TToBLNu_TuneZ2_tW-channel_7TeV-madgraph_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/merged_ntuple*.root",path));
  }

  else if( strcmp( prefix , "ZZ_PU" ) == 0 ){
    ch->Add(Form("%s/ZZtoAnything_TuneZ2_7TeV-pythia6-tauola_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/merged_ntuple*.root",path));
  }

  else if( strcmp( prefix , "WZ_PU" ) == 0 ){
    ch->Add(Form("%s/WZtoAnything_TuneZ2_7TeV-pythia6-tauola_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/merged_ntuple*.root",path));
  }

  else if( strcmp( prefix , "DYToEEM10To20_PU" ) == 0 ){
    ch->Add(Form("%s/DYToEE_M-10To20_TuneZ2_7TeV-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/dilepPt2010Skim/skim*.root",path));
  }

  else if( strcmp( prefix , "DYToEEM20_PU" ) == 0 ){
    ch->Add(Form("%s/DYToEE_M-20_TuneZ2_7TeV-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/merged_ntuple*.root",path));
  }

  else if( strcmp( prefix , "DYToMuMuM10To20_PU" ) == 0 ){
    ch->Add(Form("%s/DYToMuMu_M-10To20_TuneZ2_7TeV-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/dilepPt2010Skim/skim*.root",path));
  }

  else if( strcmp( prefix , "DYToMuMuM20_PU" ) == 0 ){
    ch->Add(Form("%s/DYToMuMu_M-20_TuneZ2_7TeV-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/merged_ntuple*.root",path));
  }

  else if( strcmp( prefix , "DYToTauTauM10To20_PU" ) == 0 ){
    ch->Add(Form("%s/DYToTauTau_M-10To20_TuneZ2_7TeV-pythia6-tauola_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/dilepPt2010Skim/skim*.root",path));
  }

  else if( strcmp( prefix , "DYToTauTauM20_PU" ) == 0 ){
    ch->Add(Form("%s/DYToTauTau_M-20_TuneZ2_7TeV-pythia6-tauola_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/merged_ntuple*.root",path));
  }

  else if( strcmp( prefix , "HToWWTo2L2NuM130" ) == 0 ){
    ch->Add(Form("%s/GluGluToHToWWTo2L2Nu_M-130_7TeV-powheg-pythia6_Fall10-START38_V12-v1/V03-06-18/merged_ntuple*.root",path));
  }

  else if( strcmp( prefix , "HToWWTo2L2NuM130_PU" ) == 0 ){
    ch->Add(Form("%s/GluGluToHToWWTo2L2Nu_M-130_7TeV-powheg-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-18/merged_ntuple*.root",path));
  }

  else if( strcmp( prefix , "HToWWTo2Tau2NuM130_PU" ) == 0 ){
    ch->Add(Form("%s/GluGluToHToWWTo2Tau2Nu_M-130_7TeV-powheg-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-18/merged_ntuple*.root",path));
  }

  else if( strcmp( prefix , "HToWWToLNuTauNuM130_PU" ) == 0 ){
    ch->Add(Form("%s/GluGluToHToWWToLNuTauNu_M-130_7TeV-powheg-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-18/merged_ntuple*.root",path));
  }

  else if( strcmp( prefix , "HToWWTo2L2NuM160" ) == 0 ){
    ch->Add(Form("%s/GluGluToHToWWTo2L2Nu_M-160_7TeV-powheg-pythia6_Fall10-START38_V12-v1/V03-06-18/merged_ntuple*.root",path));
  }

  else if( strcmp( prefix , "HToWWTo2L2NuM160_PU" ) == 0 ){
    ch->Add(Form("%s/GluGluToHToWWTo2L2Nu_M-160_7TeV-powheg-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-18/merged_ntuple*.root",path));
  }

  else if( strcmp( prefix , "HToWWTo2Tau2NuM160_PU" ) == 0 ){
    ch->Add(Form("%s/GluGluToHToWWTo2Tau2Nu_M-160_7TeV-powheg-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-18/merged_ntuple*.root",path));
  }

  else if( strcmp( prefix , "HToWWToLNuTauNuM160_PU" ) == 0 ){
    ch->Add(Form("%s/GluGluToHToWWToLNuTauNu_M-160_7TeV-powheg-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-18/merged_ntuple*.root",path));
  }

  else if( strcmp( prefix , "HToWWTo2L2NuM200" ) == 0 ){
    ch->Add(Form("%s/GluGluToHToWWTo2L2Nu_M-200_7TeV-powheg-pythia6_Fall10-START38_V12-v1/V03-06-18/merged_ntuple*.root",path));
  }

  else if( strcmp( prefix , "HToWWTo2L2NuM200_PU" ) == 0 ){
    ch->Add(Form("%s/GluGluToHToWWTo2L2Nu_M-200_7TeV-powheg-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-18/merged_ntuple*.root",path));
  }

  else if( strcmp( prefix , "HToWWTo2Tau2NuM200_PU" ) == 0 ){
    ch->Add(Form("%s/GluGluToHToWWTo2Tau2Nu_M-200_7TeV-powheg-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-18/merged_ntuple*.root",path));
  }

  else if( strcmp( prefix , "HToWWToLNuTauNuM200_PU" ) == 0 ){
    ch->Add(Form("%s/GluGluToHToWWToLNuTauNu_M-200_7TeV-powheg-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-18/merged_ntuple*.root",path));
  }

  else if( strcmp( prefix , "MinBias" ) == 0 ){
    ch->Add("/tas/cerati/PostProcMinBias/merged*root");
    isData = true;
  }

  else{
    cout << "UNRECOGNIZED SAMPLE " << prefix << ", QUITTING" << endl;
    exit(0);
  }
  
  looper* mylooper = new looper();
  
  cout << "Running on sample " << prefix << endl;
  mylooper->ScanChain(ch, prefix, isData);
  
}
