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

  TH1F *els_nValidPixelHits_corCharge	= book1DHist("els_nValidPixelHits_corCharge", 
						     "Number of valid pixel hits, correct electron charge",
						     8,
						     0.,
						     8.,
						     "# valid pixel hits",
						     "Electrons",2);
  TH1F *els_nValidPixelHits_incorCharge	= book1DHist("els_nValidPixelHits_incorCharge", 
						     "Number of valid pixel hits, incorrect electron charge",
						     8,
						     0.,
						     8.,
						     "# valid pixel hits",
						     "Electrons",2);
  TH1F *els_nValidPixelHits_corCharge_forward	= book1DHist("els_nValidPixelHits_corCharge_forward", 
						     "Number of valid pixel hits, correct electron charge, forward",
						     8,
						     0.,
						     8.,
						     "# valid pixel hits",
						     "Electrons",2);
  TH1F *els_nValidPixelHits_incorCharge_forward	= book1DHist("els_nValidPixelHits_incorCharge_forward", 
						     "Number of valid pixel hits, incorrect electron charge, forward",
						     8,
						     0.,
						     8.,
						     "# valid pixel hits",
						     "Electrons",2);

  TH1F *els_nValidPixelHits_corCharge_barrel	= book1DHist("els_nValidPixelHits_corCharge_barrel", 
						     "Number of valid pixel hits, correct electron charge, barrel",
						     8,
						     0.,
						     8.,
						     "# valid pixel hits",
						     "Electrons",2);
  TH1F *els_nValidPixelHits_incorCharge_barrel	= book1DHist("els_nValidPixelHits_incorCharge_barrel", 
						     "Number of valid pixel hits, incorrect electron charge, barrel",
						     8,
						     0.,
						     8.,
						     "# valid pixel hits",
						     "Electrons",2);

  TH1F *els_layerFirstPixelHit_corCharge	= book1DHist("els_layerFirstPixelHit_corCharge", 
						     "layer of first valid pixel hit, correct electron charge",
						     8,
						     0.,
						     8.,
						     "layer",
						     "Electrons",2);
  TH1F *els_layerFirstPixelHit_incorCharge	= book1DHist("els_layerFirstPixelHit_incorCharge", 
						     "layer of first valid pixel hit, incorrect electron charge",
						     8,
						     0.,
						     8.,
						     "layer",
						     "Electrons",2);
  TH1F *els_layerFirstPixelHit_corCharge_forward	= book1DHist("els_layerFirstPixelHit_corCharge_forward", 
						     "layer of first valid pixel hit, correct electron charge, forward",
						     8,
						     0.,
						     8.,
						     "layer",
						     "Electrons",2);
  TH1F *els_layerFirstPixelHit_incorCharge_forward	= book1DHist("els_layerFirstPixelHit_incorCharge_forward", 
						     "layer of first valid pixel hit, incorrect electron charge, forward",
						     8,
						     0.,
						     8.,
						     "layer",
						     "Electrons",2);

  TH1F *els_layerFirstPixelHit_corCharge_barrel	= book1DHist("els_layerFirstPixelHit_corCharge_barrel", 
						     "layer of first valid pixel hit, correct electron charge, barrel",
						     8,
						     0.,
						     8.,
						     "layer",
						     "Electrons",2);
  TH1F *els_layerFirstPixelHit_incorCharge_barrel	= book1DHist("els_layerFirstPixelHit_incorCharge_barrel", 
						     "layer of first valid pixel hit, incorrect electron charge, barrel",
						     8,
						     0.,
						     8.,
						     "layer",
						     "Electrons",2);

  TH1F *els_chargeFirstPixelHit_corCharge	= book1DHist("els_chargeFirstPixelHit_corCharge", 
							     "Charge of first valid pixel hit, correct electron charge",
							     200,0.,200000.,
						     "charge",
						     "Electrons",2);
  TH1F *els_chargeFirstPixelHit_incorCharge	= book1DHist("els_chargeFirstPixelHit_incorCharge", 
						     "Charge of first valid pixel hit, incorrect electron charge",
							     200,0.,200000.,
						     "charge",
						     "Electrons",2);
  TH1F *els_chargeFirstPixelHit_corCharge_forward	= book1DHist("els_chargeFirstPixelHit_corCharge_forward", 
						     "Charge of first valid pixel hit, correct electron charge, forward",
								     200,0.,200000.,
						     "charge",
						     "Electrons",2);
  TH1F *els_chargeFirstPixelHit_incorCharge_forward	= book1DHist("els_chargeFirstPixelHit_incorCharge_forward", 
						     "Charge of first valid pixel hit, incorrect electron charge, forward",
								     200,0.,200000.,
						     "charge",
						     "Electrons",2);

  TH1F *els_chargeFirstPixelHit_corCharge_barrel	= book1DHist("els_chargeFirstPixelHit_corCharge_barrel", 
						     "Charge of first valid pixel hit, correct electron charge, barrel",
								     200,0.,200000.,
						     "charge",
						     "Electrons",2);
  TH1F *els_chargeFirstPixelHit_incorCharge_barrel	= book1DHist("els_chargeFirstPixelHit_incorCharge_barrel", 
						     "Charge of first valid pixel hit, incorrect electron charge, barrel",
								     200,0.,200000.,
						     "charge",
						     "Electrons",2);
  TH1F *els_detIdFirstPixelHit_corCharge                = book1DHist("els_detIdFirstPixelHit_corCharge", 
						     "DetId of first pixel hit, correct electron charge",
						     5,
						     0.,
						     5.,
						     "DetId",
						     "Electrons",2);
  TH1F *els_detIdFirstPixelHit_incorCharge              = book1DHist("els_detIdFirstPixelHit_incorCharge", 
						     "DetId of first pixel hit, incorrect electron charge",
						     5,
						     0.,
						     5.,
						     "DetId",
						     "Electrons",2);

  TH1F *els_detIdFirstPixelHit_corCharge_barrel         = book1DHist("els_detIdFirstPixelHit_corCharge_barrel", 
						     "DetId of first pixel hit, correct electron charge, barrel",
						     5,
						     0.,
						     5.,
						     "DetId",
						     "Electrons",2);
  TH1F *els_detIdFirstPixelHit_incorCharge_barrel       = book1DHist("els_detIdFirstPixelHit_incorCharge_barrel", 
						     "DetId of first pixel hit, incorrect electron charge, barrel",
						     5,
						     0.,
						     5.,
						     "DetId",
						     "Electrons",2);

  TH1F *els_detIdFirstPixelHit_corCharge_forward         = book1DHist("els_detIdFirstPixelHit_corCharge_forward", 
						     "DetId of first pixel hit, correct electron charge, forward",
						     5,
						     0.,
						     5.,
						     "DetId",
						     "Electrons",2);
  TH1F *els_detIdFirstPixelHit_incorCharge_forward       = book1DHist("els_detIdFirstPixelHit_incorCharge_forward", 
						     "DetId of first pixel hit, incorrect electron charge, forward",
						     5,
						     0.,
						     5.,
						     "DetId",
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
// 	  if ( (charge == -1 && els_mc_id().at(els) == 11) || (charge == 1 && els_mc_id().at(els) == -11) ) {
// 	    els_nValidPixelHits_corCharge->Fill(els_valid_pixelhits().at(els));
// 	    els_layerFirstPixelHit_corCharge->Fill(els_layer1_layer().at(els));
// 	    els_chargeFirstPixelHit_corCharge->Fill(els_layer1_charge().at(els));
// 	    if ( TMath::Abs(els_mc_p4().at(els).eta()) <= 1.479 ) {
// 	      els_nValidPixelHits_corCharge_barrel->Fill(els_valid_pixelhits().at(els));
// 	      els_layerFirstPixelHit_corCharge_barrel->Fill(els_layer1_layer().at(els));
// 	      els_chargeFirstPixelHit_corCharge_barrel->Fill(els_layer1_charge().at(els));
// 	    } else {
// 	      els_nValidPixelHits_corCharge_forward->Fill(els_valid_pixelhits().at(els));
// 	      els_layerFirstPixelHit_corCharge_forward->Fill(els_layer1_layer().at(els));
// 	      els_chargeFirstPixelHit_corCharge_forward->Fill(els_layer1_charge().at(els));
// 	    }
// 	    // incorrect charge identified
// 	  } else {
// 	    els_nValidPixelHits_incorCharge->Fill(els_valid_pixelhits().at(els));
// 	    els_layerFirstPixelHit_incorCharge->Fill(els_layer1_layer().at(els));
// 	    els_chargeFirstPixelHit_incorCharge->Fill(els_layer1_charge().at(els));
// 	    if ( TMath::Abs(els_mc_p4().at(els).eta()) <= 1.479 ) {
// 	      els_nValidPixelHits_incorCharge_barrel->Fill(els_valid_pixelhits().at(els));
// 	      els_layerFirstPixelHit_incorCharge_barrel->Fill(els_layer1_layer().at(els));
// 	      els_chargeFirstPixelHit_incorCharge_barrel->Fill(els_layer1_charge().at(els));
// 	    } else {
// 	      els_nValidPixelHits_incorCharge_forward->Fill(els_valid_pixelhits().at(els));
// 	      els_layerFirstPixelHit_incorCharge_forward->Fill(els_layer1_layer().at(els));
// 	      els_chargeFirstPixelHit_incorCharge_forward->Fill(els_layer1_charge().at(els));
// 	    }
// 	  }

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

	  els_nValidPixelHits_corCharge->Fill(els_valid_pixelhits().at(els));
	  els_detIdFirstPixelHit_corCharge->Fill(els_layer1_det().at(els));
	  els_layerFirstPixelHit_corCharge->Fill(els_layer1_layer().at(els));
	  els_chargeFirstPixelHit_corCharge->Fill(els_layer1_charge().at(els));
	  if ( TMath::Abs(els_mc_p4().at(els).eta()) <= 1.479 ) {
	    els_nValidPixelHits_corCharge_barrel->Fill(els_valid_pixelhits().at(els));
	    els_detIdFirstPixelHit_corCharge_barrel->Fill(els_layer1_det().at(els));
	    els_layerFirstPixelHit_corCharge_barrel->Fill(els_layer1_layer().at(els));
	    els_chargeFirstPixelHit_corCharge_barrel->Fill(els_layer1_charge().at(els));
	  } else {

	    els_nValidPixelHits_corCharge_forward->Fill(els_valid_pixelhits().at(els));
	    // veto forward entries with 0 valid pixelhits
	    if ( els_valid_pixelhits().at(els) == 0 ) continue;

	    // make following plots for 1,2 valid pixel hits
	    if ( els_valid_pixelhits().at(els) == 1 || els_valid_pixelhits().at(els) == 2 ) {
	      els_detIdFirstPixelHit_corCharge_forward->Fill(els_layer1_det().at(els));
	      if ( els_layer1_det().at(els) == 2 ) {
		els_layerFirstPixelHit_corCharge_forward->Fill(els_layer1_layer().at(els));
		els_chargeFirstPixelHit_corCharge_forward->Fill(els_layer1_charge().at(els));
	      }
	    }
	  }

	  els_pt_recosim_corCharge->Fill(els_mc_p4().at(els).Pt());
	  els_eta_recosim_corCharge->Fill(els_mc_p4().at(els).eta());

	  // incorrect charge identified
	} else {

	  els_nValidPixelHits_incorCharge->Fill(els_valid_pixelhits().at(els));
	  els_detIdFirstPixelHit_incorCharge->Fill(els_layer1_det().at(els));
	  els_layerFirstPixelHit_incorCharge->Fill(els_layer1_layer().at(els));
	  els_chargeFirstPixelHit_incorCharge->Fill(els_layer1_charge().at(els));
	  if ( TMath::Abs(els_mc_p4().at(els).eta()) <= 1.479 ) {
	    els_nValidPixelHits_incorCharge_barrel->Fill(els_valid_pixelhits().at(els));
	    els_detIdFirstPixelHit_incorCharge_barrel->Fill(els_layer1_det().at(els));
	    els_layerFirstPixelHit_incorCharge_barrel->Fill(els_layer1_layer().at(els));
	    els_chargeFirstPixelHit_incorCharge_barrel->Fill(els_layer1_charge().at(els));
	  } else {

	    els_nValidPixelHits_incorCharge_forward->Fill(els_valid_pixelhits().at(els));
	    // veto forward entries with 0 valid pixelhits
	    if ( els_valid_pixelhits().at(els) == 0 ) continue;

	    // make following plots for 1,2 valid pixel hits
	    if ( els_valid_pixelhits().at(els) == 1 || els_valid_pixelhits().at(els) == 2 ) {
	      els_detIdFirstPixelHit_incorCharge_forward->Fill(els_layer1_det().at(els));
	      if ( els_layer1_det().at(els) == 2 ) {
		els_layerFirstPixelHit_incorCharge_forward->Fill(els_layer1_layer().at(els));
		els_chargeFirstPixelHit_incorCharge_forward->Fill(els_layer1_charge().at(els));
	      }
	    }
	  }

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
