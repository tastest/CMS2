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
//#include "../CORE/selections.cc"
#include "../CORE/electronSelections.cc"
//#include "../CORE/utilities.cc"
#include "tools.cc"

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

  float ptbins[10]= {10, 20, 30, 40, 50, 60, 70, 80, 90, 150};
  // float ptbins[5]= {10, 30, 40, 50, 200};
   float etabins[8] = {0, 0.5, 1.0, 1.28, 1.56, 1.84, 2.12, 2.5};
  
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
  TH2F *els_2Detavspt_recosim               = book2DVarHist ("els_eta_pt_recosim", 
						       "true Electron: #eta^{true} vs. p_{T}^{true}",
						       9, 
						       ptbins, 
						       7, 
						       etabins, 
						       "p_{T}",
						       "#eta",
						       "",
						       2);
  TH2 *els_2Detavspt_recosim_incorCharge =  book2DVarHist("els_eta_pt_recosim_incorCharge",
						     "true Electron with incorrect reconstructed Charge: #eta^{true} vs. p_{T}^{true}",
						     9, 
						     ptbins, 
						     7, 
						     etabins, 
						     "p_{T}",
						     "#eta",
						     "",
						     2);

//  els_2Detavspt_recosim               = book2DVarHist ("els_eta_pt_recosim", 
// 						       "true Electron: #eta^{true} vs. p_{T}^{true}",
// 						       4, 
// 						       ptbins, 
// 						       7, 
// 						       etabins, 
// 						       "p_{T}",
// 						       "#eta",
// 						       "",
// 						       2);
//   els_2Detavspt_recosim_incorCharge =  book2DVarHist("els_eta_pt_recosim_incorCharge",
// 						     "true Electron with incorrect reconstructed Charge: #eta^{true} vs. p_{T}^{true}",
// 						     4, 
// 						     ptbins, 
// 						     7, 
// 						     etabins, 
// 						     "p_{T}",
// 						     "#eta",
// 						     "",
// 						     2);
						    

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
	if(! (electronSelection_cand02(els) && electronId_extra(els))) continue;
	
	   // check how many electrons don't have an associated track
	els_trkId->Fill(els_trkidx().at(els));
	// tmp charge variable
	//double charge = els_charge().at(els); //majority logic
	//int charge = getChargeUsingMajorityLogic(els, 0.45);
	double charge = els_trk_charge().at(els);  //our method
	//double charge = els_sccharge().at(els);
	//if (charge !=els_charge().at(els)) std::cout<<"warning  " <<std::endl;
	
	if(isChargeFlip(els))continue ;
	

	//if ((cms2.els_trkidx().at(els) >= 0) && (cms2.els_trk_charge().at(els) != cms2.trks_charge().at(cms2.els_trkidx().at(els))) ) continue;

	
	//	else if ((els_trkidx().at(els) < 0) && (cms2.els_trk_charge().at(els) != els_sccharge().at(els)))continue;
	//if (isFromConversionHitPattern(els)) continue;
	

	  els_pt_reco->Fill(els_p4().at(els).Pt());
	  els_eta_reco->Fill(els_p4().at(els).eta());


	// fill reco_corCharge
	if ( (charge == -1 && els_mc_id().at(els) == 11) || (charge == 1 && els_mc_id().at(els) == -11) ) {
	  els_pt_reco_corCharge->Fill(els_p4().at(els).Pt());
	  els_eta_reco_corCharge->Fill(els_p4().at(els).eta());
	}

	// exclude reco which has no true electron match
	if ( abs(els_mc_id()[els]) != 11 ) continue;

	// fill recosim
	els_pt_recosim->Fill(els_mc_p4().at(els).Pt());
	els_eta_recosim->Fill(els_mc_p4().at(els).eta());
	
	els_2Detavspt_recosim->Fill(els_p4().at(els).Pt(), fabs(els_p4().at(els).eta())); 

	// correct charge identified
	if ( (charge == -1 && els_mc_id().at(els) == 11) || (charge == 1 && els_mc_id().at(els) == -11) ) {

	  els_pt_recosim_corCharge->Fill(els_mc_p4().at(els).Pt());
	  els_eta_recosim_corCharge->Fill(els_mc_p4().at(els).eta());

	  // incorrect charge identified
	} else {

	  els_pt_recosim_incorCharge->Fill(els_mc_p4().at(els).Pt());
	  els_eta_recosim_incorCharge->Fill(els_mc_p4().at(els).eta());

	  els_2Detavspt_recosim_incorCharge->Fill(els_p4().at(els).Pt(), fabs(els_p4().at(els).eta())); 
	}

      } //end of looping over reco eles
      
    }
    
  }

  if ( nEventsChain != nEventsTotal ) {
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  }

  return 0;
}
