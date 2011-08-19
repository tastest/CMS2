#include "ChainFromText.cc"

#include "TChain.h"

void runFR_os7June2011_DoubleMu(){
  
  gROOT->LoadMacro("myBabyMaker.C++");
  // Double Muon
  
  TChain *chain2 = new TChain("Events");
  chain2->Add("/hadoop/cms/store/user/yanjuntu/CMSSW_4_2_3_patch1_V04-02-15/DoubleMu_Run2011A-May10ReReco-v1_AOD/CMSSW_4_2_3_patch1_V04-02-15_merged/V04-02-15/merged*root");

  myBabyMaker* baby2 = new myBabyMaker();
  baby2->ScanChain(chain2, "DoubleMu_os7June2011.root", true, -1);

}
