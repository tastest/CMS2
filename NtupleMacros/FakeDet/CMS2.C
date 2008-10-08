#include <iostream>
#include <vector>
#include <set>

#include "TChain.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TROOT.h"

#include "TH1F.h"
#include "TH2F.h"

using namespace std;

#ifndef __CINT__
#include "CMS2.h"
CMS2 cms2;
#endif

#include "../Tools/selections.C"
#include "../Tools/utilities.C"

bool isDenominatorElectron(int index) {
  //
  // returns true if input fulfills certain cuts
  //

  // cut definition
  float pt_cut        		= 15.;
  float eta_cut       		= 2.5;
  float hOverE_cut    		= 0.2;
  unsigned int   njets_cut    	= 1;  // require at least N jets
  float HLT_jet_approx_uncorr 	= 30.0; // require leading jet to be larger than N GeV uncorrected v2_3
  float HLT_jet_approx 		= HLT_jet_approx_uncorr / cms2.jets_tq_noCorrF()[0]; // V00-05-01 Fake tuples have TQ corrected jets. Uncorreced cut is done here.

  bool result = true;

  if ( cms2.evt_njets() < njets_cut)                         result = false;
  if ( cms2.jets_p4()[0].Pt() < HLT_jet_approx )      	result = false;
  if ( cms2.els_p4()[index].Pt()  < pt_cut )                 result = false;
  if ( std::abs(cms2.els_p4()[index].Eta()) > eta_cut )      result = false;
  if ( !passElectronIsolation(index) )          	result = false;
  if ( cms2.els_hOverE()[index]   > hOverE_cut )             result = false;

  return result;

}

bool isNumeratorElectron(int index, int type=0) { 
  //
  // 1=loose, 2=tight
  //
  // returns true if input fulfills certain cuts
  //
  
  // cut definition
  float pt_cut        		= 15;
  float eta_cut       		= 2.5;
  unsigned int   njets_cut      = 1;  // require at least N jets
  float HLT_jet_approx_uncorr 	= 30.0; // require leading jet to be larger than N GeV uncorrected v2_3
  float HLT_jet_approx 		= HLT_jet_approx_uncorr / cms2.jets_tq_noCorrF()[0]; // V00-05-01 Fake tuples have TQ corrected jets. Uncorreced cut is done here.

  bool result = true;

  if ( cms2.evt_njets() < njets_cut)                         result = false;
  if ( cms2.jets_p4()[0].Pt() < HLT_jet_approx )      	result = false;
  if ( cms2.els_p4()[index].Pt()  < pt_cut )                 result = false;
  if ( std::abs(cms2.els_p4()[index].Eta()) > eta_cut )      result = false;
  if ( !passElectronIsolation(index) )          	result = false;

  if ( type == 1 ) {
    // loose
    if ( !goodLooseElectronWithoutIsolation(index) )   result = false;
  } else if ( type == 2 ) {
    // tight
    if ( !goodElectronWithoutIsolation(index) )   result = false;
  } else {
    cout << "WARNING: wrong electron type detected, please select loose (1) or tight (2)" << endl;
  }

  return result;
  
}

int ScanChain( TChain* chain, char * prefix = "", int nEvents = -1) {

  // event counting
  unsigned int nEventsChain=0;
  if(nEvents==-1) 
    nEventsChain = chain->GetEntries();
  else 
    nEventsChain = nEvents;
  unsigned int nEventsTotal = 0;
  unsigned int nEventsSelected = 0;

  // clear list of duplicates                                                                                                                         
  already_seen.clear();
  int duplicates_total_n = 0;
  int i_permille_old = 0;

  // book histograms

  vector<float> binsPt;
  binsPt.push_back(0.);
  //	binsPt.push_back(10.);
  binsPt.push_back(20.);
  //	binsPt.push_back(30.);
  binsPt.push_back(40.);
  //	binsPt.push_back(50.);
  binsPt.push_back(65.);
  //	binsPt.push_back(80.);
  binsPt.push_back(95.);
  //	binsPt.push_back(115.);
  binsPt.push_back(150.);
  vector<float> binsEta;
  binsEta.push_back(-2.5);
  //	binsEta.push_back(-2.1);
  binsEta.push_back(-1.9);
  //	binsEta.push_back(-1.7);
  binsEta.push_back(-1.5);
  //	binsEta.push_back(-1.3);
  binsEta.push_back(-1.1);
  //	binsEta.push_back(-0.9);
  binsEta.push_back(-0.7);
  //	binsEta.push_back(-0.5);
  binsEta.push_back(-0.3);
  // 	binsEta.push_back(-0.1);
  // 	binsEta.push_back(0.1);
  binsEta.push_back(0.3);
  //	binsEta.push_back(0.5);
  binsEta.push_back(0.7);
  //	binsEta.push_back(0.9);
  binsEta.push_back(1.1);
  //	binsEta.push_back(1.3);
  binsEta.push_back(1.5);
  //	binsEta.push_back(1.7);
  binsEta.push_back(1.9);
  //	binsEta.push_back(2.1);
  binsEta.push_back(2.5);

  TH1F *pt_num_ell  = book1DVarHist(Form("%s_hpt_num_ell",prefix),"p_{T} Numerator",binsPt,"p_{T} [GeV]","p_{T} Loose Electron Numerator");
  TH1F *pt_num_elt  = book1DVarHist(Form("%s_hpt_num_elt",prefix),"p_{T} Numerator",binsPt,"p_{T} [GeV]","p_{T} Tight Electron Numerator");
  TH1F *pt_den_ele  = book1DVarHist(Form("%s_hpt_den_ele",prefix),"p_{T} Denominator",binsPt,"p_{T} [GeV]","p_{T} Electron Denominator");

  TH1F *eta_num_ell = book1DVarHist(Form("%s_heta_num_ell",prefix),"#eta Numerator",binsEta,"#eta","#eta Loose Electron Numerator");
  TH1F *eta_num_elt = book1DVarHist(Form("%s_heta_num_elt",prefix),"#eta Numerator",binsEta,"#eta","#eta Tight Electron Numerator");
  TH1F *eta_den_ele = book1DVarHist(Form("%s_heta_den_ele",prefix),"#eta Denominator",binsEta,"#eta","#eta Electron Denominator");

  TH2F *num_ell = book2DVarHist(Form("%s_hnum_ell",prefix),"Numerator",binsEta,binsPt,"#eta","p_{T} [GeV]","Loose Electron Numerator");
  TH2F *num_elt = book2DVarHist(Form("%s_hnum_elt",prefix),"Numerator",binsEta,binsPt,"#eta","p_{T} [GeV]","Tight Electron Numerator");
  TH2F *den_ele = book2DVarHist(Form("%s_hden_ele",prefix),"Denominator",binsEta,binsPt,"#eta","p_{T} [GeV]","Electron Denominator");

  TH1F *pt_num_wo_leading_ell  = book1DVarHist(Form("%s_hpt_num_wo_leading_ell",prefix),"p_{T} Numerator without leading Jet",binsPt,"p_{T} [GeV]","p_{T} Loose Electron Numerator");
  TH1F *pt_num_wo_leading_elt  = book1DVarHist(Form("%s_hpt_num_wo_leading_elt",prefix),"p_{T} Numerator without leading Jet",binsPt,"p_{T} [GeV]","p_{T} Tight Electron Numerator");
  TH1F *pt_den_wo_leading_ele  = book1DVarHist(Form("%s_hpt_den_wo_leading_ele",prefix),"p_{T} Denominator without leading Jet",binsPt,"p_{T} [GeV]","p_{T} Electron Denominator");

  TH1F *eta_num_wo_leading_ell = book1DVarHist(Form("%s_heta_num_wo_leading_ell",prefix),"#eta Numerator without leading Jet",binsEta,"#eta","#eta Loose Electron Numerator");
  TH1F *eta_num_wo_leading_elt = book1DVarHist(Form("%s_heta_num_wo_leading_elt",prefix),"#eta Numerator without leading Jet",binsEta,"#eta","#eta Tight Electron Numerator");
  TH1F *eta_den_wo_leading_ele = book1DVarHist(Form("%s_heta_den_wo_leading_ele",prefix),"#eta Denominator without leading Jet",binsEta,"#eta","#eta Electron Denominator");

  TH2F *num_wo_leading_ell = book2DVarHist(Form("%s_hnum_wo_leading_ell",prefix),"Numerator without leading Jet",binsEta,binsPt,"#eta","p_{T} [GeV]","Loose Electron Numerator");
  TH2F *num_wo_leading_elt = book2DVarHist(Form("%s_hnum_wo_leading_elt",prefix),"Numerator without leading Jet",binsEta,binsPt,"#eta","p_{T} [GeV]","Tight Electron Numerator");
  TH2F *den_wo_leading_ele = book2DVarHist(Form("%s_hden_wo_leading_ele",prefix),"Denominator without leading Jet",binsEta,binsPt,"#eta","p_{T} [GeV]","Electron Denominator");

  TH1F *pt_num_wo_second_leading_ell  = book1DVarHist(Form("%s_hpt_num_wo_second_leading_ell",prefix),"p_{T} Numerator without second leading Jet",binsPt,"p_{T} [GeV]","p_{T} Loose Electron Numerator");
  TH1F *pt_num_wo_second_leading_elt  = book1DVarHist(Form("%s_hpt_num_wo_second_leading_elt",prefix),"p_{T} Numerator without second leading Jet",binsPt,"p_{T} [GeV]","p_{T} Tight Electron Numerator");
  TH1F *pt_den_wo_second_leading_ele  = book1DVarHist(Form("%s_hpt_den_wo_second_leading_ele",prefix),"p_{T} Denominator without second leading Jet",binsPt,"p_{T} [GeV]","p_{T} Electron Denominator");

  TH1F *eta_num_wo_second_leading_ell = book1DVarHist(Form("%s_heta_num_wo_second_leading_ell",prefix),"#eta Numerator without second leading Jet",binsEta,"#eta","#eta Loose Electron Numerator");
  TH1F *eta_num_wo_second_leading_elt = book1DVarHist(Form("%s_heta_num_wo_second_leading_elt",prefix),"#eta Numerator without second leading Jet",binsEta,"#eta","#eta Tight Electron Numerator");
  TH1F *eta_den_wo_second_leading_ele = book1DVarHist(Form("%s_heta_den_wo_second_leading_ele",prefix),"#eta Denominator without second leading Jet",binsEta,"#eta","#eta Electron Denominator");

  TH2F *num_wo_second_leading_ell = book2DVarHist(Form("%s_hnum_wo_second_leading_ell",prefix),"Numerator without second leading Jet",binsEta,binsPt,"#eta","p_{T} [GeV]","Loose Electron Numerator");
  TH2F *num_wo_second_leading_elt = book2DVarHist(Form("%s_hnum_wo_second_leading_elt",prefix),"Numerator without second leading Jet",binsEta,binsPt,"#eta","p_{T} [GeV]","Tight Electron Numerator");
  TH2F *den_wo_second_leading_ele = book2DVarHist(Form("%s_hden_wo_second_leading_ele",prefix),"Denominator without second leading Jet",binsEta,binsPt,"#eta","p_{T} [GeV]","Electron Denominator");
	
  TH1F *ele_n   = book1DHist(Form("%s_hele_n",prefix),"number of electrons per event",15,0.,15.,"n_{e}","Events");

  unsigned int nJetsBins = 50;
  float        nJetsLow  = 0.;
  float        nJetsHigh = 50.;

  TH1F *njets = book1DHist(Form("%s_hnjets",prefix),"Number of Jets",nJetsBins,nJetsLow,nJetsHigh,"n_{Jets}","Events");

  unsigned int jetPtBins = 4000;
  float        jetPtLow  = 0.;
  float        jetPtHigh = 4000.;

  TH1F *jetpt = book1DHist(Form("%s_hjetpt",prefix),"p_{T} of Jets",jetPtBins,jetPtLow,jetPtHigh,"p_{T} [GeV]","Events");

  unsigned int jetEtaBins = 30;
  float        jetEtaLow  = -3.;
  float        jetEtaHigh =  3.;

  TH1F *jeteta = book1DHist(Form("%s_hjeteta",prefix),"#eta of Jets",jetEtaBins,jetEtaLow,jetEtaHigh,"#eta","Events");

  unsigned int deltaRBins = 500;
  float        deltaRLow  = 0.;
  float        deltaRHigh = 10.;

  TH1F *deltaR = book1DHist(Form("%s_hdeltaR",prefix),"deltaR between jet and denominator",deltaRBins,deltaRLow,deltaRHigh,"#Delta R","Events");

  unsigned int tkIsoBins = 400;
  float        tkIsoLow  = 0.;
  float        tkIsoHigh = 200.;

  TH1F *tkIso = book1DHist(Form("%s_htkIso",prefix),"track isolation",tkIsoBins,tkIsoLow,tkIsoHigh,"tk_{iso}","Events");
  TH1F *tkIso_uncut = book1DHist(Form("%s_htkIso_uncut",prefix),"uncut track isolation",tkIsoBins,tkIsoLow,tkIsoHigh,"tk_{iso}","Events");

  unsigned int EOverpBins = 30;
  float        EOverpLow  = -1.5;
  float        EOverpHigh = 1.5;

  TH1F *EOverp = book1DHist(Form("%s_hEOverp",prefix),"E/p",EOverpBins,EOverpLow,EOverpHigh,"E/p","Events");
  TH1F *EOverp_uncut = book1DHist(Form("%s_hEOverp_uncut",prefix),"uncut E/p",EOverpBins,EOverpLow,EOverpHigh,"E/p","Events");

  unsigned int HOverEBins = 300;
  float        HOverELow  = -5.;
  float        HOverEHigh = 5.;

  TH1F *HOverE = book1DHist(Form("%s_hHOverE",prefix),"H/E",HOverEBins,HOverELow,HOverEHigh,"H/E","Events");
  TH1F *HOverE_uncut = book1DHist(Form("%s_hHOverE_uncut",prefix),"uncut H/E",HOverEBins,HOverELow,HOverEHigh,"H/E","Events");

  // cuts
  float deltaRCut = 0.2;

  // file loop
  TObjArray *listOfFiles = chain->GetListOfFiles();
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

      DorkyEventIdentifier id;
      if (cms2.els_d0().size() != 0) {
	id.run 		= cms2.evt_run();
	id.event 	= cms2.evt_event();
	id.trks_d0 	= cms2.els_d0()[0];
	id.hyp_lt_pt 	= cms2.els_p4()[0].pt();
	id.hyp_lt_eta 	= cms2.els_p4()[0].eta();
	id.hyp_lt_phi 	= cms2.els_p4()[0].phi();
      } else if (cms2.mus_d0().size() != 0) {
	id.run 		= cms2.evt_run();
	id.event 	= cms2.evt_event();
	id.trks_d0 	= cms2.mus_d0()[0];
	id.hyp_lt_pt 	= cms2.mus_p4()[0].pt();
	id.hyp_lt_eta 	= cms2.mus_p4()[0].eta();
	id.hyp_lt_phi 	= cms2.mus_p4()[0].phi();
      } else {
	continue;
      }

      if (is_duplicate(id)) {
	duplicates_total_n++;
	continue;
      }
      ++nEventsSelected;

      // Progress feedback to the user
      progressBar(i_permille_old, nEventsTotal, nEventsChain);

      // fill general histograms
      ele_n->Fill(cms2.evt_nels());
      njets->Fill(cms2.evt_njets());
      for ( unsigned int jet_counter = 0;
	    jet_counter < (unsigned int)cms2.evt_njets();
	    ++jet_counter ) {
	jetpt->Fill(cms2.jets_p4()[jet_counter].Pt());
	jeteta->Fill(cms2.jets_p4()[jet_counter].Eta());
      }
      for ( unsigned int electron_counter = 0;
	    electron_counter < (unsigned int)cms2.evt_nels();
	    ++electron_counter ) {
	tkIso_uncut->Fill(cms2.els_tkIso()[electron_counter]);
	EOverp_uncut->Fill(cms2.els_eOverPIn()[electron_counter]);
	HOverE_uncut->Fill(cms2.els_hOverE()[electron_counter]);
      }

      // electron loop
      for ( unsigned int electron_counter = 0;
	    electron_counter < (unsigned int)cms2.evt_nels();
	    ++electron_counter ) {

	// general electron histograms
	tkIso_uncut->Fill(cms2.els_tkIso()[electron_counter]);
	EOverp_uncut->Fill(cms2.els_eOverPIn()[electron_counter]);
	HOverE_uncut->Fill(cms2.els_hOverE()[electron_counter]);

	// fill every electron with pt > 150 GeV into last bin
	float pt = cms2.els_p4()[electron_counter].Pt();
	if (pt >= 150.) pt = 149.;

	if ( isDenominatorElectron(electron_counter)){
	  
	  pt_den_ele->Fill(pt);
	  eta_den_ele->Fill(cms2.els_p4()[electron_counter].Eta());
	  den_ele->Fill(cms2.els_p4()[electron_counter].Eta(),pt);

	  for ( unsigned int jet_counter = 0;
		jet_counter < (unsigned int)cms2.evt_njets();
		++jet_counter ) {
	    deltaR->Fill(CalculateDeltaR(cms2.els_p4()[electron_counter],cms2.jets_p4()[jet_counter]));
	  }
	  tkIso->Fill(cms2.els_tkIso()[electron_counter]);
	  EOverp->Fill(cms2.els_eOverPIn()[electron_counter]);
	  HOverE->Fill(cms2.els_hOverE()[electron_counter]);
	}

	// loose electrons
	if (isNumeratorElectron(electron_counter,1)){ 
	  // 1=loose, 2=tight
	  pt_num_ell->Fill(pt);
	  eta_num_ell->Fill(cms2.els_p4()[electron_counter].Eta());
	  num_ell->Fill(cms2.els_p4()[electron_counter].Eta(),pt);
	}

	// tight electrons
	if ( isNumeratorElectron(electron_counter,2)){ 
	  // 1=loose, 2=tight
	  pt_num_elt->Fill(pt);
	  eta_num_elt->Fill(cms2.els_p4()[electron_counter].Eta());
	  num_elt->Fill(cms2.els_p4()[electron_counter].Eta(),pt);
	}

	// exclude leading jet
	if ( CalculateDeltaR(cms2.els_p4()[electron_counter],cms2.jets_p4()[0]) >= deltaRCut ) {

	  if ( isDenominatorElectron(electron_counter)){
	    pt_den_wo_leading_ele->Fill(pt);
	    eta_den_wo_leading_ele->Fill(cms2.els_p4()[electron_counter].Eta());
	    den_wo_leading_ele->Fill(cms2.els_p4()[electron_counter].Eta(),pt);
	  }

	  // loose electrons
	  if (isNumeratorElectron(electron_counter,1)){ 
	    // 1=loose, 2=tight
	    pt_num_wo_leading_ell->Fill(pt);
	    eta_num_wo_leading_ell->Fill(cms2.els_p4()[electron_counter].Eta());
	    num_wo_leading_ell->Fill(cms2.els_p4()[electron_counter].Eta(),pt);
	  }

	  // tight electrons
	  if ( isNumeratorElectron(electron_counter,2)) { 
	    // 1=loose, 2=tight
	    pt_num_wo_leading_elt->Fill(pt);
	    eta_num_wo_leading_elt->Fill(cms2.els_p4()[electron_counter].Eta());
	    num_wo_leading_elt->Fill(cms2.els_p4()[electron_counter].Eta(),pt);
	  }

	}

	// exclude leading and second leading jet
	if ( cms2.evt_njets() > 1 ) {

	  if ( (CalculateDeltaR(cms2.els_p4()[electron_counter],cms2.jets_p4()[0]) >= deltaRCut) &&
	       (CalculateDeltaR(cms2.els_p4()[electron_counter],cms2.jets_p4()[1]) >= deltaRCut) ) {

	    if ( isDenominatorElectron(electron_counter)){
	      pt_den_wo_second_leading_ele->Fill(pt);
	      eta_den_wo_second_leading_ele->Fill(cms2.els_p4()[electron_counter].Eta());
	      den_wo_second_leading_ele->Fill(cms2.els_p4()[electron_counter].Eta(),pt);
	    }

	    // loose electrons
	    if (isNumeratorElectron(electron_counter,1)){ 
	      // 1=loose, 2=tight
	      pt_num_wo_second_leading_ell->Fill(pt);
	      eta_num_wo_second_leading_ell->Fill(cms2.els_p4()[electron_counter].Eta());
	      num_wo_second_leading_ell->Fill(cms2.els_p4()[electron_counter].Eta(),pt);
	    }

	    // tight electrons
	    if ( isNumeratorElectron(electron_counter,2)){ 
	      // 1=loose, 2=tight
	      pt_num_wo_second_leading_elt->Fill(pt);
	      eta_num_wo_second_leading_elt->Fill(cms2.els_p4()[electron_counter].Eta());
	      num_wo_second_leading_elt->Fill(cms2.els_p4()[electron_counter].Eta(),pt);
	    }
	  }
	}
      }
    }
  }

  std::cout << std::endl;
  if ( nEventsChain != nEventsTotal ) {
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  }
  
  std::cout << "Prefix: " << prefix << " processed: " << nEventsTotal
	    << " Events, found: " << duplicates_total_n 
	    << " Duplicates and selected: " << nEventsSelected
	    << " Events."
	    << endl << endl;

  cout<<" # Denom Electrons: " << eta_den_ele->GetEntries() << endl;
  cout<<" # Loose Electrons: " << eta_num_ell->GetEntries() << endl;
  cout<<" # Tight Electrons: " << eta_num_elt->GetEntries() << endl;
  cout<<" Loose Fraction: " << eta_num_ell->GetEntries()/eta_den_ele->GetEntries() << endl;
  cout<<" Tight Fraction: " << eta_num_elt->GetEntries()/eta_den_ele->GetEntries() << endl;

  cout<<" # Denom Electrons 20-40 GeV: " << pt_den_ele->GetBinContent(2) << endl;
  cout<<" # Loose Electrons 20-40 GeV: " << pt_num_ell->GetBinContent(2) << endl;
  cout<<" # Tight Electrons 20-40 GeV: " << pt_num_elt->GetBinContent(2) << endl;
  cout<<" Loose Fraction 20-40 GeV: " << pt_num_ell->GetBinContent(2)/pt_den_ele->GetBinContent(2) << endl;
  cout<<" Tight Fraction 20-40 GeV: " << pt_num_elt->GetBinContent(2)/pt_den_ele->GetBinContent(2) << endl;

  cout<<" # Denom Electrons without leading jet: " << eta_den_wo_leading_ele->GetEntries() << endl;
  cout<<" # Loose Electrons without leading jet: " << eta_num_wo_leading_ell->GetEntries() << endl;
  cout<<" # Tight Electrons without leading jet: " << eta_num_wo_leading_elt->GetEntries() << endl;
  cout<<" Loose Fraction without leading jet: " << eta_num_wo_leading_ell->GetEntries()/eta_den_wo_leading_ele->GetEntries() << endl;
  cout<<" Tight Fraction without leading jet: " << eta_num_wo_leading_elt->GetEntries()/eta_den_wo_leading_ele->GetEntries() << endl;

  cout<<" # Denom Electrons without second leading jet: " << eta_den_wo_second_leading_ele->GetEntries() << endl;
  cout<<" # Loose Electrons without second leading jet: " << eta_num_wo_second_leading_ell->GetEntries() << endl;
  cout<<" # Tight Electrons without second leading jet: " << eta_num_wo_second_leading_elt->GetEntries() << endl;
  cout<<" Loose Fraction without second leading jet: " << eta_num_wo_second_leading_ell->GetEntries()/eta_den_wo_second_leading_ele->GetEntries() << endl;
  cout<<" Tight Fraction without second leading jet: " << eta_num_wo_second_leading_elt->GetEntries()/eta_den_wo_second_leading_ele->GetEntries() << endl;

  return 0;
}
