#include "ChainFromText.cc"

#include "TChain.h"

void runFR_os7June2011_SingleMu(){
  
  gROOT->LoadMacro("myBabyMaker.C++");
  // Single Muon
  TChain *chain3 = new TChain("Events");
  chain3->Add("/hadoop/cms/store/user/jaehyeok/CMSSW_4_2_3_patch1_V04-02-15/SingleMu_Run2011A-May10ReReco-v1_AOD/CMSSW_4_2_3_patch1_V04-02-15_merged/V04-02-15/merged*root");

  myBabyMaker* baby3 = new myBabyMaker();
  baby3->ScanChain(chain3, "SingleMu_os7June2011.root", true, -1);

}
