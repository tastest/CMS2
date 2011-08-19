/* Usage:
   root [0] .L ScanChain.C++
   root [1] TFile *_file0 = TFile::Open("merged_ntuple.root")
   root [2] TChain *chain = new TChain("Events")
   root [3] chain->Add("merged_ntuple.root")

   There are several places where one may create CMS2 cms2
   It can be done here (in a doAll.C script), i.e.:

   root [4] CMS2 cms2 

   It can be done in the source as is done below, or it can be
   ascertained by including CORE/CMS2.cc as is commented out
   below.  They are all the same, and everything will work so
   long as it is created somewhere globally.

   root [5] ScanChain(chain)
*/
#include <iostream>
#include <vector>

#include "TChain.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TROOT.h"

#include "CMS2.h"
#include "branches.h"
CMS2 cms2;

//#include "../CORE/CMS2.cc"
#include "../CORE/selections.cc"
#include "../CORE/utilities.cc"

using namespace tas;

int ScanChain( TChain* chain, int nEvents = -1) {

  TObjArray *listOfFiles = chain->GetListOfFiles();

  unsigned int nEventsChain=0;
  if(nEvents==-1) 
    nEvents = chain->GetEntries();
  nEventsChain = nEvents;
  unsigned int nEventsTotal = 0;
  InitSkimmedTree();
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");

  // array to count events
  unsigned int event_count[6] = {0};

  enum { E, M, EE, EM, MM, OTHER}; // enum of lepton final states

  const int mu_shift = 0x8000000;

  std::cout << "\n\n" << std::endl;

  // file loop
  TIter fileIter(listOfFiles);
  TFile *currentFile = 0;
  while ( currentFile = (TFile*)fileIter.Next() ) {
    TFile f(currentFile->GetTitle());
    TTree *tree = (TTree*)f.Get("Events");
    cms2.Init(tree);

    //Event Loop
    unsigned int nEvents = tree->GetEntries();
    for( unsigned int event = 0; event < nEvents; ++event) {
      cms2.GetEntry(event);
      ++nEventsTotal;

      if( !(nEventsTotal % 100000) ) std::cout << "Processed " << nEventsTotal << " events..." << std::endl;

      std::vector<int> lep_idx;

      for (unsigned int mus = 0; mus < mus_p4().size(); mus++) {

	if( mus_p4()[mus].pt() < 20 ) continue;

	if( !goodMuonIsolated(mus) ) continue;

	lep_idx.push_back( mus + mu_shift);

      }

      for (unsigned int els = 0; els < els_p4().size(); els++) {

	if( els_p4()[els].pt() < 20 ) continue;

	if( !goodElectronIsolated(els, true) ) continue;

	lep_idx.push_back( els );

      }

      if( lep_idx.size() == 0 ) continue;

      if( lep_idx.size() > 2 ) continue;

      if( lep_idx.size() == 2 ) {

	switch( ( lep_idx[0] + lep_idx[1] ) / mu_shift ) {

	case 0:
	  if( els_charge()[ lep_idx[0] ] * els_charge()[ lep_idx[1] ] < 0 ) {
	    /*
	    const LorentzVector zp4_ = els_p4()[ lep_idx[0] ] + els_p4()[ lep_idx[1] ];
	    
	    if( zp4_.M() > 76 && zp4_.M() < 106 ) {

	      if( evt_tcmet() < 30 ) ++event_count[EE];
	      }*/

	    ++event_count[EE];
	  }

	  break;

	case 1:
	  ++event_count[EM];
	  break;

	case 2:
	  if( mus_charge()[ lep_idx[0] % mu_shift ] * mus_charge()[ lep_idx[1] % mu_shift ] < 0 ) {
	    /*
	      const LorentzVector zp4_ = mus_p4()[ lep_idx[0] % mu_shift ] + mus_p4()[ lep_idx[1] % mu_shift ];
	      
	      if( zp4_.M() > 76 && zp4_.M() < 106 ) {

		if( evt_tcmet() < 30 ) ++event_count[MM];
		}*/

	      ++event_count[MM];
	  }
	  
	  break;	  
	}
      }

      else if( lep_idx.size() == 1 ) {
      
	switch( lep_idx[0] / mu_shift ) {

	case 0:
	  /*{
	    
	    bool isZcand = false;

	    const LorentzVector lepp4_ = els_p4()[ lep_idx[0] ];

	    for( unsigned int trks = 0; trks < trks_trk_p4().size(); trks++ ) {

	      LorentzVector trkp4_ = trks_trk_p4()[trks];

	      LorentzVector wp4_ = lepp4_ + trkp4_;

	      if( wp4_.M() > 76 && wp4_.M() < 106 ) {
		isZcand = true;
		break;
	      }
	    }

	    if( !isZcand ) {

	      if( evt_tcmet() > 20 ) {
	      
		const LorentzVector metp4_ (evt_tcmet() * cos(evt_tcmetPhi()), evt_tcmet() * sin(evt_tcmetPhi()), 0, evt_tcmet());

		wmt_ = TMath::Sqrt( 2. * lepp4_.pt() * metp4_.pt() * (1 - cos(metp4_.phi() - lepp4_.phi())));

		//		if( wmt_ > 40 && wmt_ < 100 ) ++event_count[E];
	      }
	    }
	  }
	    */
	  ++event_count[E];

	  break;

	case 1:
	  /*{
	    
	    bool isZcand = false;

	    const LorentzVector lepp4_ = mus_p4()[ lep_idx[0] % mu_shift ];

	    for( unsigned int trks = 0; trks < trks_trk_p4().size(); trks++ ) {

	      LorentzVector trkp4_ = trks_trk_p4()[trks];

	      LorentzVector wp4_ = lepp4_ + trkp4_;

	      if( wp4_.M() > 76 && wp4_.M() < 106 ) {
		isZcand = true;
		break;
	      }
	    }

	    if( !isZcand ) {

	      if( evt_tcmet() > 20 ) {
	      
		const LorentzVector metp4_ (evt_tcmet() * cos(evt_tcmetPhi()), evt_tcmet() * sin(evt_tcmetPhi()), 0, evt_tcmet());

		wmt_ = TMath::Sqrt( 2. * lepp4_.pt() * metp4_.pt() * (1 - cos(metp4_.phi() - lepp4_.phi())));

		//		if( wmt_ > 40 && wmt_ < 100 ) ++event_count[M];
	      }
	    }
	  }
	    */
	  ++event_count[M];

	  break;
	}
      }

      nJPTS_ = 0;

      for( unsigned int jpts = 0; jpts < jpts_p4().size(); jpts++ ) {

	if( jpts_p4()[jpts].Et() < 20 ) continue;

	if( fabs( jpts_p4()[jpts].eta() ) > 2.4 ) continue;

	++nJPTS_;
      }
    }
  }

  if ( nEventsChain != nEventsTotal ) {
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  }

  std::cout << "\n\n" << std::endl;
  std::cout << "Processed " << nEventsTotal << " total W+jets events." << std::endl;

  std::cout << "W->enu candidates passing selection: " << event_count[E]     << std::endl;
  std::cout << "W->mnu candidates passing selection: " << event_count[M]     << std::endl;
  std::cout << "Z->ee candidates passing selection : " << event_count[EE]    << std::endl;
  std::cout << "Z->mm candidates passing selection : " << event_count[MM]    << std::endl;

  std::cout << "\n\n" << std::endl;

  outFile_->cd();
  outTree_->Write();
  outFile_->Close();
  return 0;
}
