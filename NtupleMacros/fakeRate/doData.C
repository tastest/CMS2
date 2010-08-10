#include "ChainFromText.cc"
void doData(){

gROOT->LoadMacro("myBabyMaker.C++");

//  TChain *chain1 = ChainFromText( "eg.txt" );
//  myBabyMaker * baby1 = new myBabyMaker();
//  baby1->ScanChain(chain1, "EG.root", true, -1);
//
//  TChain *chain2 = ChainFromText( "mu.txt" );
//  myBabyMaker * baby2 = new myBabyMaker();
//  baby2->ScanChain(chain2, "Mu.root", true, -1);

  TChain *chain3 = ChainFromText( "jmt.txt" );
  myBabyMaker * baby3 = new myBabyMaker();
  baby3->ScanChain(chain3, "JMT.root", true, -1);

//  TChain *chain4 = ChainFromText( "jmtm_july6.txt" );
//  myBabyMaker * baby4 = new myBabyMaker();
//  baby4->ScanChain(chain4, "JMTMonitor.root", true, -1);

}
