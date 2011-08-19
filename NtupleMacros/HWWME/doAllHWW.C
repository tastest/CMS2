/*

The naming of the samples follow this convention in DataType in 
http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/Smurf/Core/SmurfTree.h

*/
#if defined(__CINT__) && !defined(__MAKECINT__)
{
  gSystem->Load("libCMS2NtupleMacrosCORE.so");
  gSystem->Load("libCMS2NtupleMacrosLooper.so");
  gSystem->CompileMacro("doAllHWW.C","k");
  doAllHWW();
}
#endif 


#include "TChain.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TRegexp.h"
#include "TFile.h"
#include "TBitSet.hh"

#ifndef __CINT__
#include "doAnalysis.h"
#endif

using namespace std;

void doAllHWW() {

  // load various tools
  gSystem->Load("../Tools/MiniFWLite/libMiniFWLite");
  gROOT->ProcessLine(".L TBitSet.cc++");
  
  double Lumi = 1000.0;
  
  TFile *utilFile = new TFile("Util_HWW.root","RECREATE");
  
  //
  // Flags for files to run over 
  // (0 and 1 are easier to modify)
  //
  bool runWW    = 0;
  bool runWjets = 1;
  bool runWZ    = 0;
  bool runZZ    = 0;
  bool runHWW115 = 0;
  bool runHWW120 = 0;
  bool runHWW130 = 0;
  bool runHWW140 = 0;
  bool runHWW150 = 0;
  bool runHWW160 = 0;
  bool runHWW170 = 0;
  bool runHWW180 = 0;
  bool runHWW190 = 0;
  bool runHWW200 = 0;
  bool runHWW250 = 0;
  bool runHWW300 = 0;
  bool runHWW350 = 0;
  bool runHWW400 = 0;
  bool runHWW450 = 0;
  bool runHWW500 = 0;
  bool runHWW550 = 0;
  bool runHWW600 = 0;
  
  if (runWW) 
    ProcessSample("ww", "data/VVJetsTo4L_TuneD6T_7TeV-madgraph-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-01-12/wwfilter/merged_ntuple.root", utilFile,  -1, Lumi, 4.5, -1, "HWW","",false,true);
  
  // single lepton 10 filtered
  if  (runWjets) 
    ProcessSample("wjets","data/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Spring11-PU_S1_START311_V1G1-v1_SingleLepton/V04-01-01/PartonSkim/merged_ntuple.root", utilFile,  -1, Lumi, 31314.0, -1); 
  
  
  if (runWZ) 
    ProcessSample("wz", "data//WZtoAnything_TuneZ2_7TeV-pythia6-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/merged_ntuple*.root", utilFile,   -1 ,Lumi, 18.2, -1); 
    
  if (runZZ) 
    ProcessSample("zz","data/ZZtoAnything_TuneZ2_7TeV-pythia6-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/merged_ntuple*.root", utilFile,  -1 ,Lumi, 7.41, -1);

  if (runHWW115)
    ProcessSample("hww115", "data/GluGluToHToWWTo2L2Nu_M-115_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/*.root", utilFile,  -1 ,Lumi, 0.165*4.0/9.0, -1); 
  
  if (runHWW120)
    ProcessSample("hww120", "data/GluGluToHToWWTo2L2Nu_M-120_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v2/V04-01-01/*.root", utilFile,  -1 ,Lumi, 0.250*4.0/9.0, -1); 
  
  if (runHWW130)
    ProcessSample("hww130","data/GluGluToHToWWTo2L2Nu_M-130_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/*.root", utilFile,  -1 ,Lumi, 0.4520899*4.0/9.0, -1); 

  if (runHWW140)
    ProcessSample("hww140","data/GluGluToHToWWTo2L2Nu_M-140_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/*.root", utilFile,  -1 ,Lumi, 0.642*4.0/9.0, -1); 

  if (runHWW150)
    ProcessSample("hww150","data/GluGluToHToWWTo2L2Nu_M-150_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/*.root", utilFile,  -1 ,Lumi, 0.770*4.0/9.0, -1); 
  
  if (runHWW160)
    ProcessSample("hww160", "data/GluGluToHToWWTo2L2Nu_M-160_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/*.root", utilFile,  -1 ,Lumi, 0.8664429*4.0/9.0, -1); 
  if (runHWW170)
    ProcessSample("hww170","data/GluGluToHToWWTo2L2Nu_M-170_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/*.root", utilFile,  -1 ,Lumi, 0.783*4.0/9.0, -1); 
  
  if (runHWW180)
    ProcessSample("hww180","data/GluGluToHToWWTo2L2Nu_M-180_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/*.root", utilFile,  -1 ,Lumi, 0.659*4.0/9.0, -1); 
  
  if (runHWW190)
    ProcessSample("hww190","data/GluGluToHToWWTo2L2Nu_M-190_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/*.root", utilFile,  -1 ,Lumi, 0.486*4.0/9.0, -1); 
  
  if (runHWW200) 
    ProcessSample("hww200","data/GluGluToHToWWTo2L2Nu_M-200_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/*.root", utilFile,  -1 ,Lumi, 0.4083049*4.0/9.0, -1); 

  if (runHWW250)
    ProcessSample("hww250","data/GluGluToHToWWTo2L2Nu_M-250_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/*.root", utilFile,  -1 ,Lumi, 0.244*4.0/9.0, -1); 

  if (runHWW300)
    ProcessSample("hww300","data/GluGluToHToWWTo2L2Nu_M-300_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/*.root", utilFile,  -1 ,Lumi,0.176*4.0/9.0, -1); 

  if (runHWW350)
    ProcessSample("hww350","data/GluGluToHToWWTo2L2Nu_M-350_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/*.root", utilFile,  -1 ,Lumi,0.160*4.0/9.0, -1); 

  if (runHWW400)
    ProcessSample("hww400","data/GluGluToHToWWTo2L2Nu_M-400_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/*.root", utilFile,  -1 ,Lumi,0.124330*4.0/9.0, -1); 

  if (runHWW450)
    ProcessSample("hww450","data/GluGluToHToWWTo2L2Nu_M-450_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/*.root", utilFile,  -1 ,Lumi,0.078433*4.0/9.0, -1); 

  if (runHWW500)
    ProcessSample("hww500","data/GluGluToHToWWTo2L2Nu_M-500_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/*.root", utilFile,  -1 ,Lumi,0.048702*4.0/9.0, -1); 

  if (runHWW550)
    ProcessSample("hww550","data/GluGluToHToWWTo2L2Nu_M-550_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/*.root", utilFile,  -1 ,Lumi,0.030364*4.0/9.0, -1); 

  if (runHWW600)
    ProcessSample("hww600","data/GluGluToHToWWTo2L2Nu_M-600_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/*.root", utilFile,  -1 ,Lumi,0.019184*4.0/9.0, -1); 
  
  utilFile->Close();	      
  
}
