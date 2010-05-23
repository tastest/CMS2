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

  vector<string> v_trignames;
  map<string,int> v_numEvtsPerTrigger;

  v_trignames.push_back("HLT_L1Jet6U");
  v_trignames.push_back("HLT_L1Jet10U");
  v_trignames.push_back("HLT_Jet15U");
  v_trignames.push_back("HLT_Jet30U");
  v_trignames.push_back("HLT_Photon10_L1R");
  v_trignames.push_back("L1_SingleEG5");
  v_trignames.push_back("HLT_Photon15_L1R");

  vector<string> v_l1assoc;
  v_l1assoc.push_back("L1_SingleJet6");
  v_l1assoc.push_back("L1_SingleJet10");
  v_l1assoc.push_back("L1_SingleJet6");
  v_l1assoc.push_back("L1_SingleJet20");
  v_l1assoc.push_back("L1_SingleEG5");
  v_l1assoc.push_back("NA");
  v_l1assoc.push_back("L1_SingleEG8");

  map<string, TH1F*> m_prescales;
  map<string, TH1F*> m_maxpfjetpt;
  for( unsigned int i=0;i<v_trignames.size(); i++ ) {
    m_prescales[v_trignames[i]] = MakeHist((v_trignames[i]+"prescales").c_str(),  ("prescales for "+v_trignames[i]+";log of product of prescales").c_str(), 10, -0.01, 9.99);
    m_maxpfjetpt[v_trignames[i]] = MakeHist((v_trignames[i]+"maxpfjetpt").c_str(),  ("max pfjet pT "+v_trignames[i]+";pT").c_str(), 500, 0.0, 500.0);
  }

  const double pi = TMath::Pi();

  const double crystalunit = 0.017453292519943295;
  const double barrelend = 85.5*crystalunit;

  const float metmaxw = 2500.;
  const int metbinsw = 250;

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
	  else continue;
	  
	  //if( !(passesTrackCuts() && requireTrackCuts)) //cout << "failed tracking cuts - this should never happen" <<  endl;
	  //  continue;

	  nGoodEvents++;

	  int pL1Jet6U = 1;
	  int pL1Jet10U = 1;

	  map<string, int> m_l1hlt;
	  for( unsigned int i=0; i < v_trignames.size(); i++ ) { //for all hlt trigs we want
		bool foundl1 = false;
		for( unsigned int j=0; j < l1_trigNames().size(); j++ ) {
		  //cout << l1_trigNames()[j] << endl;
		  if( v_l1assoc[i] == "NA" ) {
			foundl1 = true;
			continue;
		  }
		  else if( v_l1assoc[i] == l1_trigNames()[j] ) { //find associated l1 trig
			m_l1hlt[v_trignames[i]] = j;
			foundl1 = true;
			break;
		  }
		}
		if( !foundl1 )
		  cout << "Missed l1 assoc for " << v_trignames[i] << endl;
	  }
	  
	  //hlt
	  for( unsigned int i=0; i<hlt_trigNames().size(); i++ ) {
		for( unsigned int j=0; j<v_trignames.size(); j++ ) {
		  if( hlt_trigNames()[i] == v_trignames[j] && passHLTTrigger(hlt_trigNames()[i]) ){
		    //here's where we fill hists fr hlt triggers
		    m_prescales[v_trignames[j]]->Fill( log10( hlt_prescales()[i] * l1_prescales()[ m_l1hlt[(string)hlt_trigNames()[i]] ] ) );
		    float maxpfjetpt = 0;
		    for ( unsigned int k=0; k<pfjets_p4().size(); k++){
		      if( maxpfjetpt < pfjets_p4()[k].pt() ) maxpfjetpt = pfjets_p4()[k].pt();
		    }
		    m_maxpfjetpt[v_trignames[j]]->Fill(maxpfjetpt);
		  }//end of filling hists for hlt triggers
		}
	  }

	  //l1
	  for( unsigned int i=0; i<l1_trigNames().size(); i++ ) {
		for( unsigned int j=0; j<v_trignames.size(); j++ ) {
		  if( l1_trigNames()[i] == v_trignames[j] && passL1Trigger( l1_trigNames()[i] ) ){
		    //here's where we fill hists for l1 triggers
		    m_prescales[v_trignames[j]]->Fill( log10(l1_prescales()[j]) );
		    float maxpfjetpt = 0;
		    for ( unsigned int k=0; k<pfjets_p4().size(); k++){
		      if( maxpfjetpt < pfjets_p4()[k].pt() ) maxpfjetpt = pfjets_p4()[k].pt();
		    }
		    m_maxpfjetpt[v_trignames[j]]->Fill(maxpfjetpt);

		  }//end of filling hists for l1 triggers
		}
	  }
	  
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

