/*

The naming of the samples follow this convention in DataType in 
http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/Smurf/Core/SmurfTree.h

*/
#if defined(__CINT__) && !defined(__MAKECINT__)
{
  gSystem->Load("libCMS2NtupleMacrosCORE.so");
  gSystem->Load("libCMS2NtupleMacrosLooper.so");
  gSystem->CompileMacro("doAllHZZ.C","k");
  doAllHZZ();
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

void doAllHZZ() {

  // load various tools
  gSystem->Load("../Tools/MiniFWLite/libMiniFWLite");
  gROOT->ProcessLine(".L TBitSet.cc++");
  
  double Lumi = 1000.0;
  
  TFile *utilFile = new TFile("Util_HZZ.root","RECREATE");
  
  //
  // Flags for files to run over 
  // (0 and 1 are easier to modify)
  //
  bool runWW    = 1;
  bool runWZ    = 1;
  bool runZZ    = 1;
  bool runHZZ250 = 1;
  bool runHZZ300 = 1;
  bool runHZZ350 = 1;
  bool runHZZ400 = 1;
  bool runHZZ500 = 1;
  bool runHZZ600 = 1;
  
  if (runWW) 
    ProcessSample("ww", "data/VVJetsTo4L_TuneD6T_7TeV-madgraph-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-01-12/wwfilter/merged_ntuple*.root", utilFile,  -1, Lumi, 4.5, -1, "HZZ","",false,true);
  
  if (runWZ) 
    ProcessSample("wz", "data//WZtoAnything_TuneZ2_7TeV-pythia6-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/merged_ntuple*.root", utilFile,   -1 ,Lumi, 18.2, -1, "HZZ"); 
    
  if (runZZ) 
    ProcessSample("zz","data/ZZtoAnything_TuneZ2_7TeV-pythia6-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/merged_ntuple*.root", utilFile,  -1 ,Lumi, 7.41, -1, "HZZ");

  if (runHZZ250)
    ProcessSample("hzz250","data/GluGluToHToZZTo2L2Nu_M-250_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-12/wwfilter/merged_ntuple*.root", utilFile,  -1 ,Lumi, 0.03974, -1, "HZZ"); 
  
  if (runHZZ300)
    ProcessSample("hzz300","data/GluGluToHToZZTo2L2Nu_M-300_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-12/wwfilter/merged_ntuple*.root", utilFile,  -1 ,Lumi, 0.0300396, -1, "HZZ"); 
  
 if (runHZZ350)
    ProcessSample("hzz350","data/GluGluToHToZZTo2L2Nu_M-350_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-12/wwfilter/merged_ntuple*.root", utilFile,  -1 ,Lumi, 0.0286009, -1, "HZZ"); 

 if (runHZZ400)
   ProcessSample("hzz400","data/GluGluToHToZZTo2L2Nu_M-400_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-12/wwfilter/merged_ntuple*.root", utilFile,  -1 ,Lumi, 0.0220830, -1, "HZZ"); 

 if (runHZZ500)
   ProcessSample("hzz500","data/GluGluToHToZZTo2L2Nu_M-500_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-12/wwfilter/merged_ntuple*.root", utilFile,  -1 ,Lumi, 0.0089522, -1, "HZZ"); 

 if (runHZZ600)
   ProcessSample("hzz600","data/GluGluToHToZZTo2L2Nu_M-600_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-12/wwfilter/merged_ntuple*.root", utilFile,  -1 ,Lumi, 0.0035900, -1, "HZZ"); 

  
  utilFile->Close();	      
  

}
