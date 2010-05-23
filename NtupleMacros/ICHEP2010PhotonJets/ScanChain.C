#include <iostream>
#include <vector>
#include <algorithm>
#include <math.h>

#include "TChain.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TH2F.h"
#include "CMS2.h"
#include "TProfile.h"
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "/home/users/wandrews/CMSSW_3_3_6/src/CMS2/NtupleMacros/Tools/goodrun.cc"
#include "CMS2.cc"
#include "ScanChainUtilities.cc"
//#include "../CORE/metSelections.h"
//#include "../CORE/jetSelections.cc"
//#include "../CORE/metSelections.cc"
//#include "../CORE/trackSelections.cc"


using namespace tas;
//typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;


TString ScanChain( TChain* chain, bool isGEN, bool requireTrackCuts = true, bool is2tev=false, int nEvents = -1) {

  //cout << "starting" << endl;
  i_permille_old = 0; //clear tty counter once per ScanChain
  TObjArray *listOfFiles = chain->GetListOfFiles();

  unsigned int nEventsChain=0;
  if(nEvents==-1) 
    nEvents = chain->GetEntries();
  nEventsChain = nEvents;
  unsigned int nEventsTotal = 0;
  //TDirectory *rootdir = gDirectory->GetDirectory("Rint:"); //wasn't used

  // file loop
  TIter fileIter(listOfFiles);

  vector<string> v_prefix;
  vector<string> v_title;
  //talk to Warren about how to deal with this map stuff. I want to count the number of evts that pass any of our favorite triggers
  map<string,int> v_numEvtsPerTrigger;

  if(isGEN) {
    v_prefix.push_back("mcpt_");
    v_prefix.push_back("mcft_"); 
    v_prefix.push_back("mcall_");
    v_title.push_back(" for MC events passing the trigger Requirements");
    v_title.push_back(" for MC events failing the trigger Requirements");
    v_title.push_back(" for all MC events");
  } else {
    v_prefix.push_back("HLT_L1Jet6U_");
    v_prefix.push_back("HLT_L1Jet10U_");
    v_prefix.push_back("HLT_Jet15U_");
    v_prefix.push_back("HLT_Jet30U_");
    v_prefix.push_back("HLT_Photon10_");
    v_prefix.push_back("L1_SingleEG5_");
    v_prefix.push_back("HLT_Photon15_");
    v_title.push_back(" for events passing L1Jet6U");
    v_title.push_back(" for events passing L1Jet10U");
    v_title.push_back(" for events passing Jet15U");
    v_title.push_back(" for events passing Jet30U");
    v_title.push_back(" for events passing Photon10");
    v_title.push_back(" for events passing L1SingleEG5");
    v_title.push_back(" for events passing Photon15");
  }

  if(v_prefix.size() != v_title.size() ) {
    cout << "The vector of prefixes and the vector of title are not the same size!!! Exiting!" << endl;
    return "";
  }

  unsigned int aSize = v_prefix.size();
  const double pi = TMath::Pi();

  const double crystalunit = 0.017453292519943295;
  const double barrelend = 85.5*crystalunit;
  
  //1d hist TH1F *h_bla1d[aSize];       
  //2d hist TH2F *h_bla2d[aSize];

  TH1F* h_prescales[aSize];

  const float metmaxw = 2500.;
  const int metbinsw = 250;

  for(unsigned int i = 0; i < v_prefix.size(); i++) {
    //1d hist NewHist(h_bla1d[i],(v_prefix.at(i)+"bla1d").c_str(), ";bla 1d", metbins, 0.0, metmax, isGEN);
    //2d hist: NewHist(h_bla2d[i],(v_prefix.at(i)+"bla2d").c_str(), ";bla x ;bla y", 100, -pi, pi, d0bins, -d0max, d0max, isGEN);

    NewHist(h_prescales[i],(v_prefix.at(i)+"prescales").c_str(), ";log of product of prescales",10, -0.01, 9.99, isGEN);

  }

  TDirectory* histdir = gDirectory;
  //pass fail counters
  int nGoodEvents = 0;
  int npassgoodrun = 0;

  set<int> v_goodRuns;
  map<int, int> nGoodEventsPerRun;
  //initialize all
  v_goodRuns.insert(0); //just do this to align with the maps below
  nGoodEventsPerRun[0] = 0;

  while ( TFile *currentFile = (TFile*)fileIter.Next() ) {
	//cout << "starting file loop" << endl;
    TFile f(currentFile->GetTitle());
    TTree *tree = (TTree*)f.Get("Events");
    cms2.Init(tree);
    
    //Event Loop
    unsigned int nEvents = tree->GetEntries();

    for( unsigned int event = 0; event < nEvents; ++event) {
      cms2.GetEntry(event);
      ++nEventsTotal;

	  progressbar( nEventsTotal, nEventsChain );

	  //THIS IS TO LOOK FOR ABSURD MET
	  //if( evt_met() < 100. )
	  //continue;

	  //insert jmu's goodrun fn
	  if( isGEN || goodrun( evt_run(), evt_lumiBlock() ) ){
	    npassgoodrun++;
	    v_goodRuns.insert(evt_run());
	  }
      
	  if( !(passesTrackCuts() && requireTrackCuts)) cout << "failed tracking cuts - this should never happen" <<  endl;

	  nGoodEvents++;

	  int pL1Jet6U = 1;
	  int pL1Jet10U = 1;
	  
	  
      //cout << "end of event loop" << endl << endl;
    }//event loop
  }//file loop

  cout << "\n\n********************SUMMARY********************" << endl;
  cout << "Total number of events: " << nEventsTotal << endl;
  cout << "Total number of events that pass good run: " << npassgoodrun 
	   << " (" << 100*(double)npassgoodrun/nEventsTotal << "% of total)" << endl;

  cout << "Total number of events that pass run and trigger cuts: " << nGoodEvents
	   << " (" << 100*(double)nGoodEvents/npassgoodrun << "%)" << endl;

  if ( nEventsChain != nEventsTotal ) {
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  }

  return "";
}

