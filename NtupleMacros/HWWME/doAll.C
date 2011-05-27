#include "TFile.h"

/*

The naming of the samples follow this convention in DataType in 
http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/Smurf/Core/SmurfTree.h

*/


int doAll() {

  gROOT->ProcessLine(".L TBitSet.cc+");
  gROOT->ProcessLine(".L ScanChain.C+");
  
  using namespace std;
 
  double Lumi = 1000.0;
  
  TFile *utilFile = new TFile("Util.root","RECREATE");
  
  //
  // Flags for files to run over 
  // (0 and 1 are easier to modify)
  //
  bool runHWW120 = 0;
  bool runHWW130 = 0;
  bool runHWW140 = 0;
  bool runHWW150 = 0;
  bool runHWW160 = 0;
  bool runHWW170 = 0;
  bool runHWW180 = 0;
  bool runHWW190 = 0;
  bool runHWW200 = 0;
  bool runHWW210 = 0;
  bool runHWW220 = 0;
  bool runHWW230 = 0;
  bool runHWW250 = 0;
  bool runHWW300 = 0;
  bool runWW    = 0;
  bool runWZ    = 0;
  bool runZZ    = 0;
  bool runWjets = 1;
  bool runWgamma = 0;
  bool runDYee  = 0;
  bool runDYmm  = 0;
  bool runDYtt  = 0;
  bool runttbar = 0;
  bool runtW    = 0;
  bool runQCD   = 0; 
  bool runData  = 0;
  
  if (runHWW120)
    ProcessSample("hww120", "data/GluGluToHToWWTo2L2Nu_M-120_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v2/V04-01-01/*.root", utilFile,  -1 ,Lumi, 0.250*4.0/9.0, -1); 
  
  if (runHWW130)
    // ProcessSample("hww130","data/GluGluToHToWWTo2L2Nu_M-130_7TeV-powheg-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-18/*.root", utilFile,  -1 ,Lumi, 0.4520899*4.0/9.0, -1); 
    ProcessSample("hww130","data/GluGluToHToWWTo2L2Nu_M-130_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/*.root", utilFile,  -1 ,Lumi, 0.4520899*4.0/9.0, -1); 

  if (runHWW140)
    ProcessSample("hww140","data/GluGluToHToWWTo2L2Nu_M-140_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/*.root", utilFile,  -1 ,Lumi, 0.642*4.0/9.0, -1); 

  if (runHWW150)
    ProcessSample("hww150","data/GluGluToHToWWTo2L2Nu_M-150_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/*.root", utilFile,  -1 ,Lumi, 0.770*4.0/9.0, -1); 
  
  if (runHWW160)
    //ProcessSample("hww160", "data/GluGluToHToWWTo2L2Nu_M-160_7TeV-powheg-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-18/*.root", utilFile,  -1 ,Lumi, 0.8664429*4.0/9.0, -1); 
    ProcessSample("hww160", "data/GluGluToHToWWTo2L2Nu_M-160_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/*.root", utilFile,  -1 ,Lumi, 0.8664429*4.0/9.0, -1); 

  if (runHWW170)
    ProcessSample("hww170","data/GluGluToHToWWTo2L2Nu_M-170_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/*.root", utilFile,  -1 ,Lumi, 0.783*4.0/9.0, -1); 
  
  if (runHWW180)
    ProcessSample("hww180","data/GluGluToHToWWTo2L2Nu_M-180_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/*.root", utilFile,  -1 ,Lumi, 0.659*4.0/9.0, -1); 
  

  if (runHWW190)
    ProcessSample("hww190","data/GluGluToHToWWTo2L2Nu_M-190_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/*.root", utilFile,  -1 ,Lumi, 0.486*4.0/9.0, -1); 
  
  if (runHWW200) 
    // ProcessSample("hww200","data/GluGluToHToWWTo2L2Nu_M-200_7TeV-powheg-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-18/*.root", utilFile,  -1 ,Lumi, 0.4083049*4.0/9.0, -1); 
    ProcessSample("hww200","data/GluGluToHToWWTo2L2Nu_M-200_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/*.root", utilFile,  -1 ,Lumi, 0.4083049*4.0/9.0, -1); 

  if (runHWW210)
    ProcessSample("hww210","data/GluGluToHToWWTo2L2Nu_M-210_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/*.root", utilFile,  -1 ,Lumi,0.358*4.0/9.0, -1); 
  
  if (runHWW220)
    ProcessSample("hww220","data/GluGluToHToWWTo2L2Nu_M-220_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/*.root", utilFile,  -1 ,Lumi,0.321*4.0/9.0, -1); 

  if (runHWW230)
    ProcessSample("hww230","data/GluGluToHToWWTo2L2Nu_M-230_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/*.root", utilFile,  -1 ,Lumi,0.290*4.0/9.0, -1); 
  
  if (runHWW250) 
    // ProcessSample("hww250","data/GluGluToHToWWTo2L2Nu_M-250_7TeV-powheg-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-18/*.root", utilFile,  -1 ,Lumi, 0.244*4.0/9.0, -1); 
    ProcessSample("hww250","data/GluGluToHToWWTo2L2Nu_M-250_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-00/*.root", utilFile,  -1 ,Lumi, 0.244*4.0/9.0, -1); 

  if (runHWW300)
    ProcessSample("hww300","data/GluGluToHToWWTo2L2Nu_M-300_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/*.root", utilFile,  -1 ,Lumi,0.176*4.0/9.0, -1); 

  if (runWW) 
    // ProcessSample("ww", "data/VVJetsTo4L_TuneD6T_7TeV-madgraph-tauola_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/*.root", utilFile,  -1, Lumi, 4.5*0.919, 963356*(682015./963356.), "", false, true);
    ProcessSample("ww", "data/WWTo2L2Nu_TuneZ2_7TeV-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/*.root", utilFile,  -1, Lumi, 4.5, -1);
  
  if (runttbar) 
    ProcessSample("ttbar", "data/TTJets_TuneZ2_7TeV-madgraph-tauola_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v3/V03-06-18/*.root", -1 ,Lumi, 157.5, -1); 
  
  if  (runWjets) {
    // ProcessSample("wjets","data/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Fall10-START38_V12-v1/V03-06-18/merged_ntuple*root", utilFile,  -1, Lumi, 31314.0, -1); // single lepton 10 filtered
    // ProcessSample("wjets","data/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/SingleLepton/merged_ntuple*root", utilFile,  -1, Lumi, 31314.0, -1); // single lepton 10 filtered
    ProcessSample("wjets","data/G_Pt_15to30_TuneZ2_7TeV_pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-00/merged_ntuple*root", utilFile,  -1, Lumi, 31314.0, -1); // single lepton 10 filtered
  }
  if (runDYmm)
    ProcessSample("dymm", "data/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/*.root", -1 ,Lumi, 1666, -1); 
  
  if (runDYee) 
    ProcessSample("dyee", "data/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/*.root", ch, -1 ,Lumi, 1666, -1); 
  
  if (runDYtt) 
    ProcessSample("dytt", "data/DYToTauTau_M-20_TuneZ2_7TeV-pythia6-tauola_Fall10-START38_V12-v1/V03-06-14/*.root", utilFile,   -1 ,Lumi, 1666, -1); 
  
  if (runWZ) 
    ProcessSample("wz", "data/WZtoAnything_TuneZ2_7TeV-pythia6-tauola_Fall10-START38_V12-v1/V03-06-14/*.root", utilFile,   -1 ,Lumi, 18.2, -1); 
  
  
  if (runZZ) 
    ProcessSample("zz","data/ZZtoAnything_TuneZ2_7TeV-pythia6-tauola_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/*.root", utilFile,  -1 ,Lumi, 5.9, -1); 
  
  if (runWgamma)
    ProcessSample("wgamma", "data/PhotonVJets_7TeV-madgraph_Fall10-START38_V12-v1/CMS2_V03-06-17/*.root", utilFile,   -1 ,Lumi, 165.0, -1); 

 if (runData) {
  std::vector<string> dataSamples;
  dataSamples.push_back(dataset+"/Mu_Run2010A-Sep17ReReco_v2_RECO/V03-06-14/diLepPt1020Skim/"+version+"/*.root");
  dataSamples.push_back(dataset+"/Mu_Run2010B-PromptReco-v2_RECO/V03-06-14-00/diLepPt1020Skim/"+version+"/*.root");
  dataSamples.push_back(dataset+"/Mu_Run2010B-PromptReco-v2_RECO/V03-06-14/diLepPt1020Skim/"+version+"/*.root");
  dataSamples.push_back(dataset+"/EG_Run2010A-Sep17ReReco_v2_RECO/V03-06-14/diLepPt1020Skim/"+version+"/*.root");
  dataSamples.push_back(dataset+"/Electron_Run2010B-PromptReco-v2_RECO/V03-06-14-00/diLepPt1020Skim/"+version+"/*.root");
  dataSamples.push_back(dataset+"/Electron_Run2010B-PromptReco-v2_RECO/V03-06-14/diLepPt1020Skim/"+version+"/*.root");
  ProcessSample("data", dataSamples, NprocEvents ,Lumi, -1, -1); 
 }
 
 utilFile->Close();	      
  
}
