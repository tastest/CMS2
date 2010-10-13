#include <iostream>
#include <fstream>
#include <string>

#include "myBabyMaker.C"
#include "ChainFromText.cc"

using namespace std;

void doit() {
  //gSystem->Load("../Tools/MiniFWLite/libMiniFWLite.so"); //put this in your root logon, foo
  //gROOT->LoadMacro("myBabyMaker.C++");

  //TChain* chain1 = new TChain("Events"); //replace with cool script from Derek
  //chain1->Add("/tas/cms2/Mu_Run2010A-Jun9thReReco_v1_RECO/V03-05-01/singleLepPt5Skim/*.root");
  //chain1->Add("/tas/cms2/Mu_Run2010A-Jun9thReReco_v1_RECO/V03-05-01/singleLepPt5Skim/skimmed_ntuple_34.root");

  const bool doeg = true;
  const bool domu = true;

  const bool dobefore = true;
  const bool doafter  = true;

  const string gdruntxt = "goodruns.txt"; // was : jsonlist_132440_138751.txt

  if( dobefore ) {
	if( doeg ) {
	  myBabyMaker * baby_Egbefore = new myBabyMaker();
	  TChain* chain_Egbefore = ChainFromText("EG_before.txt");
	  //bools are electrons, muons
	  baby_Egbefore->ScanChain(chain_Egbefore, "validate_els_before.root", gdruntxt.c_str(), true, false);
	}
	if( domu ) {
	  myBabyMaker * baby_Mubefore = new myBabyMaker();
	  TChain* chain_Mubefore = ChainFromText("Mu_before.txt");
	  baby_Mubefore->ScanChain(chain_Mubefore, "validate_mus_before.root", gdruntxt.c_str(), false, true);
	}
  }

  if( doafter ) {
	if( doeg ) {
	  myBabyMaker * baby_Egafter = new myBabyMaker();
	  TChain* chain_Egafter = ChainFromText("EG_after.txt");
	  //bools are electrons, muons
	  baby_Egafter->ScanChain(chain_Egafter, "validate_els_after.root", gdruntxt.c_str(), true, false);
	}
	if( domu ) {
	  myBabyMaker * baby_Muafter = new myBabyMaker();
	  TChain* chain_Muafter = ChainFromText("Mu_after.txt");
	  baby_Muafter->ScanChain(chain_Muafter, "validate_mus_after.root", gdruntxt.c_str(), false, true);
	}
  }

}
