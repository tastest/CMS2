#include "ChainFromText.cc"
void doData(){

gROOT->LoadMacro("myBabyMaker.C++");

  TChain *chain1 = ChainFromText( "eg.txt" );
  myBabyMaker * baby1 = new myBabyMaker();
  baby1->ScanChain(chain1, "EG.root", true, -1);

  TChain *chain2 = ChainFromText( "mu.txt" );
  myBabyMaker * baby2 = new myBabyMaker();
  baby2->ScanChain(chain2, "Mu.root", true, -1);

  TChain *chain3 = ChainFromText( "jmt.txt" );
  myBabyMaker * baby3 = new myBabyMaker();
  baby3->ScanChain(chain3, "JMT.root", true, -1);

  TChain *chain4 = ChainFromText( "qcd30.txt" );
  myBabyMaker * baby4 = new myBabyMaker();
  baby4->ScanChain(chain4, "qcd30.root", false, -1);


//    TChain *chain4 = ChainFromText( "qcd15.txt" );
//    myBabyMaker * baby4 = new myBabyMaker();
//    baby4->ScanChain(chain4, "qcd15.root", false, -1);
//
//
//  TChain *chain5 = ChainFromText( "inclMu.txt" );
//  myBabyMaker * baby5 = new myBabyMaker();
//  baby5->ScanChain(chain5, "inclMu.root", false, -1);
//
//  TChain *chain6 = ChainFromText( "wmunu.txt" );
//  myBabyMaker * baby6 = new myBabyMaker();
//  baby6->ScanChain(chain6, "Wmunu.root", false, -1);
//
//  TChain *chain7 = ChainFromText( "wenu.txt" );
//  myBabyMaker * baby7 = new myBabyMaker();
//  baby7->ScanChain(chain7, "Wenu.root", false, -1);

}
