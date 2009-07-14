#include <iostream>
#include <vector>

#include "TMath.h"
#include "TChain.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TROOT.h"

#include "CMS2.h"
CMS2 cms2;
//#include "CORE/CMS2.cc"
#include "../CORE/selections.cc"
#include "../CORE/utilities.cc"
#include "../Tools/tools.cc"

using namespace tas;

int ScanChain( TChain* chain) {

  TObjArray *listOfFiles = chain->GetListOfFiles();

  unsigned int nEventsChain=chain->GetEntries();;
  unsigned int nEventsTotal = 0;

  // book histograms
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  unsigned int nBinsPt 	= 55;
  float lowBinsPt 	= 0.;
  float highBinsPt 	= 110.;

  unsigned int nBinsEta = 52;
  float lowBinsEta      = -2.6;
  float highBinsEta     =  2.6;
  
  TH1F *els_pt_sim 			= book1DHist("els_pt_sim", 
						     "true Electron: p_{T}^{true}",
						     nBinsPt,
						     lowBinsPt,
						     highBinsPt,
						     "p_{T}^{true} [GeV]",
						     "Electrons",2);
  TH1F *els_eta_sim 			= book1DHist("els_eta_sim", 
						     "true Electron: #eta^{true}",
						     nBinsEta,
						     lowBinsEta,
						     highBinsEta,
						     "#eta^{true}",
						     "Electrons",2);
  TH1F *els_pt_recosim 			= book1DHist("els_pt_recosim", 
						     "true Electron: p_{T}^{true}",
						     nBinsPt,
						     lowBinsPt,
						     highBinsPt,
						     "p_{T}^{true} [GeV]",
						     "Electrons",2);
  TH1F *els_eta_recosim 		= book1DHist("els_eta_recosim", 
						     "true Electron: #eta^{true}",
						     nBinsEta,
						     lowBinsEta,
						     highBinsEta,
						     "#eta^{true}",
						     "Electrons",2);
  TH1F *els_pt_reco 			= book1DHist("els_pt_reco", 
						     "reco Electron: p_{T}^{true}",
						     nBinsPt,
						     lowBinsPt,
						     highBinsPt,
						     "p_{T}^{true} [GeV]",
						     "Electrons",2);
  TH1F *els_eta_reco 			= book1DHist("els_eta_reco", 
						     "reco Electron: #eta^{true}",
						     nBinsEta,
						     lowBinsEta,
						     highBinsEta,
						     "#eta^{true}",
						     "Electrons",2);
  TH1F *els_pt_recosim_incorCharge 	= book1DHist("els_pt_recosim_incorCharge", 
						     "true Electron with incorrect reconstructed Charge: p_{T}^{true}",
						     nBinsPt,
						     lowBinsPt,
						     highBinsPt,
						     "p_{T}^{true} [GeV]",
						     "Electrons",2);
  TH1F *els_eta_recosim_incorCharge 	= book1DHist("els_eta_recosim_incorCharge", 
						     "true Electron with incorrect reconstructed Charge: #eta^{true}",
						     nBinsEta,
						     lowBinsEta,
						     highBinsEta,
						     "#eta^{true}",
						     "Electrons",2);
  TH1F *els_pt_recosim_corCharge 	= book1DHist("els_pt_recosim_corCharge", 
						     "true Electron with correct reconstructed Charge: p_{T}^{true}",
						     nBinsPt,
						     lowBinsPt,
						     highBinsPt,
						     "p_{T}^{true} [GeV]",
						     "Electrons",2);
  TH1F *els_eta_recosim_corCharge 	= book1DHist("els_eta_recosim_corCharge", 
						     "true Electron with correct reconstructed Charge: #eta^{true}",
						     nBinsEta,
						     lowBinsEta,
						     highBinsEta,
						     "#eta^{true}",
						     "Electrons",2);
  TH1F *els_pt_reco_corCharge 		= book1DHist("els_pt_reco_corCharge", 
						     "reco Electron with correct reconstructed Charge: p_{T}^{true}",
						     nBinsPt,
						     lowBinsPt,
						     highBinsPt,
						     "p_{T}^{true} [GeV]",
						     "Electrons",2);
  TH1F *els_eta_reco_corCharge 		= book1DHist("els_eta_reco_corCharge", 
						     "reco Electron with correct reconstructed Charge: #eta^{true}",
						     nBinsEta,
						     lowBinsEta,
						     highBinsEta,
						     "#eta^{true}",
						     "Electrons",2);

  TH1F *els_trkId         		= book1DHist("els_trkId", 
						     "Track ID of associated to reco Electron",
						     30,
						     -1050.,
						     450.,
						     "Track ID",
						     "Electrons",2);


  // file loop
  TIter fileIter(listOfFiles);
  TFile *currentFile = 0;
  while ( currentFile = (TFile*)fileIter.Next() ) {
    TFile f(currentFile->GetTitle());
    TTree *tree = (TTree*)f.Get("Events");
    cms2.Init(tree);
    
    //Event Loop
    unsigned int nEvents = tree->GetEntries();
    // nEvents = 100;
    for( unsigned int event = 0; event < nEvents; ++event) {
      cms2.GetEntry(event);
      ++nEventsTotal;
      
      if ( nEventsTotal%10000 == 0 ) {
	std::cout << "Event: " << nEventsTotal << endl;
      }

      // loop over true electrons
      for ( unsigned int els = 0;
	    els < genps_p4().size();
	    ++els ) {

	// check that electron is final state electron
	if ( genps_status()[els] != 1 ) continue;
	
	// check for true electron
	if ( TMath::Abs(genps_id()[els]) != 11 ) continue;

	// fill true histrograms
	els_pt_sim->Fill(genps_p4()[els].pt());
	els_eta_sim->Fill(genps_p4()[els].eta());

      }

      // loop over reco electrons
      for (unsigned int els = 0; 
           els < els_p4().size(); 
	   ++els) {

	// cuts
	//if ( !goodElectronWithoutIsolation(els) ) continue;
	if ( !goodElectronIsolated(els) ) continue;
	if ( conversionElectron(els) ) continue;
	// Yanjun's conversion removal
	if ( conversionElectron_PIXHIT(els) ) continue;

	// check how many electrons don't have an associated track
	els_trkId->Fill(els_trkidx().at(els));

	// tmp charge variable
	double charge = els_charge().at(els);

	// if electron has associated track and track charge is not equal to electron charge, veto
 	int trk = els_trkidx().at(els);
	if ( (trk >= 0) && (charge != trks_charge().at(trk)) ) {
	  continue;
	}

	// fill reco
	els_pt_reco->Fill(els_p4().at(els).Pt());
	els_eta_reco->Fill(els_p4().at(els).eta());

	// fill reco_corCharge
	if ( (charge == -1 && els_mc_id().at(els) == 11) || (charge == 1 && els_mc_id().at(els) == -11) ) {
	  els_pt_reco_corCharge->Fill(els_p4().at(els).Pt());
	  els_eta_reco_corCharge->Fill(els_p4().at(els).eta());
	}

	// exclude reco which has no true electron match
	if ( els_mc_id()[els] == -999 ) continue;

	// fill recosim
	els_pt_recosim->Fill(els_mc_p4().at(els).Pt());
	els_eta_recosim->Fill(els_mc_p4().at(els).eta());

	// correct charge identified
	if ( (charge == -1 && els_mc_id().at(els) == 11) || (charge == 1 && els_mc_id().at(els) == -11) ) {

	  els_pt_recosim_corCharge->Fill(els_mc_p4().at(els).Pt());
	  els_eta_recosim_corCharge->Fill(els_mc_p4().at(els).eta());

	  // incorrect charge identified
	} else {

	  els_pt_recosim_incorCharge->Fill(els_mc_p4().at(els).Pt());
	  els_eta_recosim_incorCharge->Fill(els_mc_p4().at(els).eta());
	}

      }
      
    }
    
  }

  if ( nEventsChain != nEventsTotal ) {
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  }

  return 0;
}
