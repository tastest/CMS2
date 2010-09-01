#include "ChainFromText.cc"
void doData(){

gROOT->LoadMacro("myBabyMaker.C++");

// Data

  // EG
  TChain *chain1 = ChainFromText( "input_data/eg_ucsd.txt" );
  myBabyMaker * baby1 = new myBabyMaker();
  baby1->ScanChain(chain1, "EG.root", true, -1);

  // Mu
  TChain *chain2 = ChainFromText( "input_data/mu_uaf.txt" );
  myBabyMaker * baby2 = new myBabyMaker();
  baby2->ScanChain(chain2, "Mu.root", true, -1);

  // JMT, JM
  TChain *chain3 = ChainFromText( "input_data/jmt_uaf.txt" );
  myBabyMaker * baby3 = new myBabyMaker();
  baby3->ScanChain(chain3, "JMT.root", true, -1);

// Monte Carlo

  // QCD 15
  TChain *chain4 = ChainFromText( "input_data/qcd15_uaf.txt" );
  myBabyMaker * baby4 = new myBabyMaker();
  baby4->ScanChain(chain4, "qcd15.root", false, -1);

  // QCD 30
  TChain *chain5 = ChainFromText( "input_data/qcd30_uaf.txt" );
  myBabyMaker * baby5 = new myBabyMaker();
  baby5->ScanChain(chain5, "qcd30.root", false, -1);

  // Wmunu
  TChain *chain6 = ChainFromText( "input_data/wmunu_uaf.txt" );
  myBabyMaker * baby6 = new myBabyMaker();
  baby6->ScanChain(chain6, "Wmunu.root", false, -1);

  // Wenu
  TChain *chain7 = ChainFromText( "input_data/wenu_uaf.txt" );
  myBabyMaker * baby7 = new myBabyMaker();
  baby7->ScanChain(chain7, "Wenu.root", false, -1);

  // InclusiveMu 15
  TChain *chain8 = ChainFromText( "input_data/inclMu_uaf.txt" );
  myBabyMaker * baby8 = new myBabyMaker();
  baby8->ScanChain(chain8, "inclMu.root", false, -1);

}
