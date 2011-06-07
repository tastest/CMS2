#include "ChainFromText.cc"

#include "TChain.h"

void runFR_os7June2011_DoubleElectron(){
  
  gROOT->LoadMacro("myBabyMaker.C++");

  // Double Electron
  TChain *chain1 = new TChain("Events");
  chain1->Add("/hadoop/cms/store/user/yanjuntu/CMSSW_4_2_3_patch1_V04-02-15/DoubleElectron_Run2011A-May10ReReco-v1_AOD/CMSSW_4_2_3_patch1_V04-02-15_merged/V04-02-15/merged*root");

  myBabyMaker* baby1 = new myBabyMaker();
  baby1->ScanChain(chain1, "DoubleElectron_os7June2011.root", true, -1);

}
