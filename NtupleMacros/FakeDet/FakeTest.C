//now make the source file
#include <iostream>
#include <vector>
#include <set>

#include "TChain.h"
#include "TFile.h"

#include "TH1F.h"
#include "TH2F.h"

using namespace std;

#ifndef __CINT__
#include "CMS2.h"
CMS2 cms2;
#endif

#include "../Tools/selections.C"
#include "../Tools/utilities.C"

int ScanChain( TChain* chain, char * prefix="", int specDY=-1, float kFactor=1.0, int WjBG=0, TH2F* theFakeRate=0) {

  // Make sure the specDY flags is kosher                                                                                                             
  if (specDY < -1 || specDY > 2) {
    std::cout << "specDY flag is not allowed...quit" << std::endl;
    return 1;
  }

  // clear list of duplicates                                                                                                                         
  already_seen.clear();
  int duplicates_total_n = 0;
  double duplicates_total_weight = 0;
  int i_permille_old = 0;

  TObjArray *listOfFiles = chain->GetListOfFiles();

  unsigned int nEventsChain = chain->GetEntries();
  unsigned int nEventsTotal = 0;
  unsigned int nCandidatesSelected = 0;

  const unsigned int allBuckets = 3;
  char *suffix[allBuckets+1];
  suffix[0] = "ee";
  suffix[1] = "mm";
  suffix[2] = "em";
  suffix[3] = "all";

  // declare histograms

  TH1F* hnEleMCId[allBuckets+1];                // mc id of electron
  TH1F* hnEleMotherMCId[allBuckets+1];          // mc id of mother of electron

  TH1F* hnJet[allBuckets+1];       		// Njet distributions
  TH1F* helePt[allBuckets+1];      		// electron Pt
  TH1F* helePtTrue[allBuckets+1];      		// electron Pt, MC tagged electron (not matched to signal)
  TH1F* helePtFake[allBuckets+1];      		// electron Pt, MC anti-tagged electron (not matched to signal)
  TH1F* helePtTrueEW[allBuckets+1];      	// electron Pt, MC tagged electron (matched to W)
  						//	TH1F* helePtFakeEW[allBuckets+1];      // electron Pt, MC anti-tagged electron (matched to W)
  TH1F* helePtTrueFSR[allBuckets+1];      	// electron Pt, MC tagged electron (matched to Final state radiation from the lepton)
  						//plErr
  TH1F* hnJetplErr[allBuckets+1];       	// Njet distributions
  TH1F* helePtplErr[allBuckets+1];      	// electron Pt
  TH1F* helePtTrueplErr[allBuckets+1];      	// electron Pt, MC tagged electron (not matched to signal)
  TH1F* helePtFakeplErr[allBuckets+1];      	// electron Pt, MC anti-tagged electron (not matched to signal)
  TH1F* helePtTrueEWplErr[allBuckets+1];      	// electron Pt, MC tagged electron (matched to W)
  						//	TH1F* helePtFakeEWplErr[allBuckets+1];      // electron Pt, MC anti-tagged electron (matched to W)
  TH1F* helePtTrueFSRplErr[allBuckets+1];      	// electron Pt, MC tagged electron (matched to Final state radiation from the lepton)
  						//miErr
  TH1F* hnJetmiErr[allBuckets+1];       	// Njet distributions
  TH1F* helePtmiErr[allBuckets+1];      	// electron Pt
  TH1F* helePtTruemiErr[allBuckets+1];      	// electron Pt, MC tagged electron (not matched to signal)
  TH1F* helePtFakemiErr[allBuckets+1];      	// electron Pt, MC anti-tagged electron (not matched to signal)
  TH1F* helePtTrueEWmiErr[allBuckets+1];      	// electron Pt, MC tagged electron (matched to W)
  						//	TH1F* helePtFakeEWmiErr[allBuckets+1];      // electron Pt, MC anti-tagged electron (matched to W)
  TH1F* helePtTrueFSRmiErr[allBuckets+1];      	// electron Pt, MC tagged electron (matched to Final state radiation from the lepton)
  						//cont...
  TH1F* hmuPt[allBuckets+1];       		// muon Pt
  TH1F* hminLepPt[allBuckets+1];   		// minimum lepton Pt
  TH1F* hmaxLepPt[allBuckets+1];   		// maximum lepton Pt
  TH1F* helePhi[allBuckets+1];     		// electron phi
  TH1F* hmuPhi[allBuckets+1];      		// muon phi
  TH1F* hdphiLep[allBuckets+1];    		// delta phi between leptons
  TH1F* heleEta[allBuckets+1];     		// electron eta
  TH1F* heleEtaTrue[allBuckets+1];     		// electron eta
  TH1F* heleEtaTrueEW[allBuckets+1];     	// electron eta
  TH1F* heleEtaTrueFSR[allBuckets+1];     	// electron eta
  						//plErr
  TH1F* heleEtaplErr[allBuckets+1];     	// electron eta
  TH1F* heleEtaTrueplErr[allBuckets+1];     	// electron eta
  TH1F* heleEtaTrueEWplErr[allBuckets+1];     	// electron eta
  TH1F* heleEtaTrueFSRplErr[allBuckets+1];     	// electron eta
  						//miErr
  TH1F* heleEtamiErr[allBuckets+1];     	// electron eta
  TH1F* heleEtaTruemiErr[allBuckets+1];     	// electron eta
  TH1F* heleEtaTrueEWmiErr[allBuckets+1];     	// electron eta
  TH1F* heleEtaTrueFSRmiErr[allBuckets+1];     	// electron eta
  						//cont...
  TH1F* hmuEta[allBuckets+1];      		// muon eta
  TH1F* hdilMass[allBuckets+1];    		// dilepton mass
  TH1F* hdilMassTightWindow[allBuckets+1]; 	// dilepton mass, but zooming around Z
  TH1F* hdilPt[allBuckets+1];       		// dilepton Pt
  TH1F* hmet[allBuckets+1];       		// MET
  TH1F* hmetPhi[allBuckets+1];       		// MET phi
  TH2F* hmetVsDilepPt[allBuckets+1];  		// MET vs dilepton Pt

  TH2F* hmetOverPtVsDphi[allBuckets+1]; 	// MET/Lepton Pt vs DeltaPhi between MET and Lepton Pt
  TH1F* hptJet1[allBuckets+1];   		// Pt of 1st jet
  TH1F* hptJet2[allBuckets+1];   		// Pt of 2nd jet
  TH1F* hptJet3[allBuckets+1];   		// Pt of 3rd jet
  TH1F* hptJet4[allBuckets+1];   		// Pt of 4th jet
  TH1F* hetaJet1[allBuckets+1];   		// eta of 1st jet
  TH1F* hetaJet2[allBuckets+1];   		// eta of 2nd jet
  TH1F* hetaJet3[allBuckets+1];   		// eta of 3rd jet
  TH1F* hetaJet4[allBuckets+1];   		// eta of 4th jet

  for (unsigned int i=0; i<=allBuckets; ++i) {

    hnEleMCId[i] = new TH1F(Form("%s_hnEleMCId_%s",prefix,suffix[i]),Form("%s_hnEleMCId_%s",prefix,suffix[i]),
			    1000,0.,1000.);	
    hnEleMotherMCId[i] = new TH1F(Form("%s_hnEleMotherMCId_%s",prefix,suffix[i]),Form("%s_hnEleMotherMCId_%s",prefix,suffix[i]),
			    1000,0.,1000.);	

    hnJet[i] = new TH1F(Form("%s_hnJet_%s",prefix,suffix[i]),Form("%s_nJet_%s",prefix,suffix[i]),
			5,0.,5.);	
    helePt[i] = new TH1F(Form("%s_helePt_%s",prefix,suffix[i]),Form("%s_elePt_%s",prefix,suffix[i]),
			 150,0.,150.);
    helePtTrue[i] = new TH1F(Form("%s_helePtTrue_%s",prefix,suffix[i]),Form("%s_elePtTrue_%s",prefix,suffix[i]),
			     150,0.,150.);
    helePtFake[i] = new TH1F(Form("%s_helePtFake_%s",prefix,suffix[i]),Form("%s_elePtFake_%s",prefix,suffix[i]),
			     150,0.,150.);
    helePtTrueEW[i] = new TH1F(Form("%s_helePtTrueEW_%s",prefix,suffix[i]),Form("%s_elePtTrueEW_%s",prefix,suffix[i]),
			       150,0.,150.);
    //	  helePtFakeEW[i] = new TH1F(Form("%s_helePtFakeEW_%s",prefix,suffix[i]),Form("%s_elePtFakeEW_%s",prefix,suffix[i]),
    //			       150,0.,150.);
    helePtTrueFSR[i] = new TH1F(Form("%s_helePtTrueFSR_%s",prefix,suffix[i]),Form("%s_elePtTrueFSR_%s",prefix,suffix[i]),
				150,0.,150.);
    //plErr
    hnJetplErr[i] = new TH1F(Form("%s_hnJetplErr_%s",prefix,suffix[i]),Form("%s_nJetplErr_%s",prefix,suffix[i]),
			     5,0.,5.);	
    helePtplErr[i] = new TH1F(Form("%s_helePtplErr_%s",prefix,suffix[i]),Form("%s_elePtplErr_%s",prefix,suffix[i]),
			      150,0.,150.);
    helePtTrueplErr[i] = new TH1F(Form("%s_helePtTrueplErr_%s",prefix,suffix[i]),Form("%s_elePtTrueplErr_%s",prefix,suffix[i]),
				  150,0.,150.);
    helePtFakeplErr[i] = new TH1F(Form("%s_helePtFakeplErr_%s",prefix,suffix[i]),Form("%s_elePtFakeplErr_%s",prefix,suffix[i]),
				  150,0.,150.);
    helePtTrueEWplErr[i] = new TH1F(Form("%s_helePtTrueEWplErr_%s",prefix,suffix[i]),Form("%s_elePtTrueEWplErr_%s",prefix,suffix[i]),
				    150,0.,150.);
    //	  helePtFakeEWplErr[i] = new TH1F(Form("%s_helePtFakeEWplErr_%s",prefix,suffix[i]),Form("%s_elePtFakeEWplErr_%s",prefix,suffix[i]),
    //			       150,0.,150.);
    helePtTrueFSRplErr[i] = new TH1F(Form("%s_helePtTrueFSRplErr_%s",prefix,suffix[i]),Form("%s_elePtTrueFSRplErr_%s",prefix,suffix[i]),
				     150,0.,150.);
    //miErr
    hnJetmiErr[i] = new TH1F(Form("%s_hnJetmiErr_%s",prefix,suffix[i]),Form("%s_nJetmiErr_%s",prefix,suffix[i]),
			     5,0.,5.);	
    helePtmiErr[i] = new TH1F(Form("%s_helePtmiErr_%s",prefix,suffix[i]),Form("%s_elePtmiErr_%s",prefix,suffix[i]),
			      150,0.,150.);
    helePtTruemiErr[i] = new TH1F(Form("%s_helePtTruemiErr_%s",prefix,suffix[i]),Form("%s_elePtTruemiErr_%s",prefix,suffix[i]),
				  150,0.,150.);
    helePtFakemiErr[i] = new TH1F(Form("%s_helePtFakemiErr_%s",prefix,suffix[i]),Form("%s_elePtFakemiErr_%s",prefix,suffix[i]),
				  150,0.,150.);
    helePtTrueEWmiErr[i] = new TH1F(Form("%s_helePtTrueEWmiErr_%s",prefix,suffix[i]),Form("%s_elePtTrueEWmiErr_%s",prefix,suffix[i]),
				    150,0.,150.);
    //	  helePtFakeEWmiErr[i] = new TH1F(Form("%s_helePtFakeEWmiErr_%s",prefix,suffix[i]),Form("%s_elePtFakeEWmiErr_%s",prefix,suffix[i]),
    //			       150,0.,150.);
    helePtTrueFSRmiErr[i] = new TH1F(Form("%s_helePtTrueFSRmiErr_%s",prefix,suffix[i]),Form("%s_elePtTrueFSRmiErr_%s",prefix,suffix[i]),
				     150,0.,150.);
    //cont...
    hmuPt[i]  = new TH1F(Form("%s_hmuPt_%s",prefix,suffix[i]),Form("%s_muPt_%s",prefix,suffix[i]),
			 150,0.,150.);
    hminLepPt[i]  = new TH1F(Form("%s_hminLepPt_%s",prefix,suffix[i]),
			     Form("%s_minLepPt_%s",prefix,suffix[i]),150,0.,150.);
    hmaxLepPt[i]  = new TH1F(Form("%s_hmaxLepPt_%s",prefix,suffix[i]),
			     Form("%s_maxLepPt_%s",prefix,suffix[i]),150,0.,150.);
    helePhi[i] = new TH1F(Form("%s_helePhi_%s",prefix,suffix[i]),Form("%s_elePhi_%s",prefix,suffix[i]),
			  50,-1*TMath::Pi(), TMath::Pi());
    hmuPhi[i]  = new TH1F(Form("%s_hmuPhi_%s",prefix,suffix[i]),Form("%s_muPhi_%s",prefix,suffix[i]),
			  50,-1*TMath::Pi(), TMath::Pi());
    hdphiLep[i]  = new TH1F(Form("%s_hdphiLep_%s",prefix,suffix[i]),Form("%s_dphiLep_%s",prefix,suffix[i]),
			    50,0., TMath::Pi());
    heleEta[i] = new TH1F(Form("%s_heleEta_%s",prefix,suffix[i]),Form("%s_eleEta_%s",prefix,suffix[i]),
			  60, -3., 3.);
    heleEtaTrue[i] = new TH1F(Form("%s_heleEtaTrue_%s",prefix,suffix[i]),Form("%s_eleEtaTrue_%s",prefix,suffix[i]),
			      60, -3., 3.);
    heleEtaTrueEW[i] = new TH1F(Form("%s_heleEtaTrueEW_%s",prefix,suffix[i]),Form("%s_eleEtaTrueEW_%s",prefix,suffix[i]),
				60, -3., 3.);
    heleEtaTrueFSR[i] = new TH1F(Form("%s_heleEtaTrueFSR_%s",prefix,suffix[i]),Form("%s_eleEtaTrueFSR_%s",prefix,suffix[i]),
				 60, -3., 3.);
    //plErr
    heleEtaplErr[i] = new TH1F(Form("%s_heleEtaplErr_%s",prefix,suffix[i]),Form("%s_eleEtaplErr_%s",prefix,suffix[i]),
			       60, -3., 3.);
    heleEtaTrueplErr[i] = new TH1F(Form("%s_heleEtaTrueplErr_%s",prefix,suffix[i]),Form("%s_eleEtaTrueplErr_%s",prefix,suffix[i]),
				   60, -3., 3.);
    heleEtaTrueEWplErr[i] = new TH1F(Form("%s_heleEtaTrueEWplErr_%s",prefix,suffix[i]),Form("%s_eleEtaTrueEWplErr_%s",prefix,suffix[i]),
				     60, -3., 3.);
    heleEtaTrueFSRplErr[i] = new TH1F(Form("%s_heleEtaTrueFSRplErr_%s",prefix,suffix[i]),Form("%s_eleEtaTrueFSRplErr_%s",prefix,suffix[i]),
				      60, -3., 3.);
    //miErr
    heleEtamiErr[i] = new TH1F(Form("%s_heleEtamiErr_%s",prefix,suffix[i]),Form("%s_eleEtamiErr_%s",prefix,suffix[i]),
			       60, -3., 3.);
    heleEtaTruemiErr[i] = new TH1F(Form("%s_heleEtaTruemiErr_%s",prefix,suffix[i]),Form("%s_eleEtaTruemiErr_%s",prefix,suffix[i]),
				   60, -3., 3.);
    heleEtaTrueEWmiErr[i] = new TH1F(Form("%s_heleEtaTrueEWmiErr_%s",prefix,suffix[i]),Form("%s_eleEtaTrueEWmiErr_%s",prefix,suffix[i]),
				     60, -3., 3.);
    heleEtaTrueFSRmiErr[i] = new TH1F(Form("%s_heleEtaTrueFSRmiErr_%s",prefix,suffix[i]),Form("%s_eleEtaTrueFSRmiErr_%s",prefix,suffix[i]),
				      60, -3., 3.);
    //cont...
    hmuEta[i]  = new TH1F(Form("%s_hmuEta_%s",prefix,suffix[i]),Form("%s_muEta_%s",prefix,suffix[i]),
			  60, -3., 3.);
    hdilMass[i] = new TH1F(Form("%s_hdilMass_%s",prefix,suffix[i]),Form("%s_dilMass_%s",prefix,suffix[i]),
			   100, 0., 300.);
    hdilMassTightWindow[i] = new TH1F(Form("%s_hdilMassTightWindow_%s",prefix,suffix[i]),
				      Form("%s_dilMassTightWindow_%s",prefix,suffix[i]),
				      120, 60., 120.);
    hdilPt[i] = new TH1F(Form("%s_hdilPt_%s",prefix,suffix[i]),Form("%s_dilPt_%s",prefix,suffix[i]),
			 100, 0., 300.);
    hmet[i] = new TH1F(Form("%s_hmet_%s",prefix,suffix[i]),Form("%s_met_%s",prefix,suffix[i]),100,0.,200.);
    hmetPhi[i] = new TH1F(Form("%s_hmetPhi_%s",prefix,suffix[i]),Form("%s_metPhi_%s",prefix,suffix[i]),
			  50,-1*TMath::Pi(), TMath::Pi());
    hmetVsDilepPt[i] = new TH2F(Form("%s_hmetVsDilepPt_%s",prefix,suffix[i]),
				Form("%s_metVsDilepPt_%s",prefix,suffix[i]),
				100,0.,200.,100,0.,200.);
    hmetOverPtVsDphi[i] = new TH2F(Form("%s_hmetOverPtVsDphi_%s",prefix,suffix[i]),
				   Form("%s_metOverPtVsDphi_%s",prefix,suffix[i]),
				   100,0.,3.,50,0., TMath::Pi());
    hptJet1[i] = new TH1F(Form("%s_hptJet1_%s",prefix,suffix[i]),Form("%s_ptJet1_%s",prefix,suffix[i]),
			  100, 0., 300.);
    hptJet2[i] = new TH1F(Form("%s_hptJet2_%s",prefix,suffix[i]),Form("%s_ptJet2_%s",prefix,suffix[i]),
			  100, 0., 300.);
    hptJet3[i] = new TH1F(Form("%s_hptJet3_%s",prefix,suffix[i]),Form("%s_ptJet3_%s",prefix,suffix[i]),
			  100, 0., 300.);
    hptJet4[i] = new TH1F(Form("%s_hptJet4_%s",prefix,suffix[i]),Form("%s_ptJet4_%s",prefix,suffix[i]),
			  100, 0., 300.);

    hetaJet1[i] = new TH1F(Form("%s_hetaJet1_%s",prefix,suffix[i]),Form("%s_etaJet1_%s",prefix,suffix[i]),
			   50, -4., 4.);
    hetaJet2[i] = new TH1F(Form("%s_hetaJet2_%s",prefix,suffix[i]),Form("%s_etaJet2_%s",prefix,suffix[i]),
			   50, -4., 4.);
    hetaJet3[i] = new TH1F(Form("%s_hetaJet3_%s",prefix,suffix[i]),Form("%s_etaJet3_%s",prefix,suffix[i]),
			   50, -4., 4.);
    hetaJet4[i] = new TH1F(Form("%s_hetaJet4_%s",prefix,suffix[i]),Form("%s_etaJet4_%s",prefix,suffix[i]),
			   50, -4., 4.);


    hnJet[i]->Sumw2();
    helePt[i]->Sumw2();
    helePtTrue[i]->Sumw2();
    helePtFake[i]->Sumw2();
    helePtTrueEW[i]->Sumw2();
    //	  helePtFakeEW[i]->Sumw2();
    helePtTrueFSR[i]->Sumw2();
    //plErr
    hnJetplErr[i]->Sumw2();
    helePtplErr[i]->Sumw2();
    helePtTrueplErr[i]->Sumw2();
    helePtFakeplErr[i]->Sumw2();
    helePtTrueEWplErr[i]->Sumw2();
    //	  helePtFakeEWplErr[i]->Sumw2();
    helePtTrueFSRplErr[i]->Sumw2();
    //miErr
    hnJetmiErr[i]->Sumw2();
    helePtmiErr[i]->Sumw2();
    helePtTruemiErr[i]->Sumw2();
    helePtFakemiErr[i]->Sumw2();
    helePtTrueEWmiErr[i]->Sumw2();
    //	  helePtFakeEWmiErr[i]->Sumw2();
    helePtTrueFSRmiErr[i]->Sumw2();
    //cont...
    hmuPt[i]->Sumw2();
    hminLepPt[i]->Sumw2();
    hmaxLepPt[i]->Sumw2();
    helePhi[i]->Sumw2();
    hmuPhi[i]->Sumw2();
    hdphiLep[i]->Sumw2();
    heleEta[i]->Sumw2();
    heleEtaTrue[i]->Sumw2();
    heleEtaTrueEW[i]->Sumw2();
    heleEtaTrueFSR[i]->Sumw2();
    //plErr
    heleEtaplErr[i]->Sumw2();
    heleEtaTrueplErr[i]->Sumw2();
    heleEtaTrueEWplErr[i]->Sumw2();
    heleEtaTrueFSRplErr[i]->Sumw2();
    //miErr
    heleEtamiErr[i]->Sumw2();
    heleEtaTruemiErr[i]->Sumw2();
    heleEtaTrueEWmiErr[i]->Sumw2();
    heleEtaTrueFSRmiErr[i]->Sumw2();
    //cont...
    hmuEta[i]->Sumw2();
    hdilMass[i]->Sumw2();
    hdilMassTightWindow[i]->Sumw2();
    hdilPt[i]->Sumw2();
    hmet[i]->Sumw2();
    hmetPhi[i]->Sumw2();
    hptJet1[i]->Sumw2();
    hptJet2[i]->Sumw2();
    hptJet3[i]->Sumw2();
    hptJet4[i]->Sumw2();
    hetaJet1[i]->Sumw2();
    hetaJet2[i]->Sumw2();
    hetaJet3[i]->Sumw2();
    hetaJet4[i]->Sumw2();

  }

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

      if (cms2.trks_d0().size() == 0) continue;
      DorkyEventIdentifier id = { cms2.evt_run(), cms2.evt_event(), cms2.trks_d0()[0],
				  cms2.hyp_lt_p4()[0].pt(), cms2.hyp_lt_p4()[0].eta(), cms2.hyp_lt_p4()[0].phi() };
      if (is_duplicate(id)) {
	duplicates_total_n++;
	duplicates_total_weight += cms2.evt_scale1fb();
	continue;
      }

      // Progress feedback to the user
      progressBar(i_permille_old, nEventsTotal, nEventsChain);

      // The event weight including the kFactor (scaled to 1 fb-1)
      float weight = cms2.evt_scale1fb() * kFactor;

      // special handling for DY
      bool processEvent=true;
      if (specDY == 0) {
	if ( !isDYee() ) processEvent = false;
      } else if (specDY == 1) {
	if ( !isDYmm() ) processEvent = false;
      } else if (specDY == 2) {
	if ( !isDYtt() ) processEvent = false;
      }
      if (!processEvent) continue;
      
      // loop over trilepton candidates
      for ( unsigned int cand = 0; 
	    cand < cms2.hyp_ll_p4().size();
	    ++cand ) {

	unsigned int bucket = 99;
	if (cms2.hyp_type()[cand] == 3) bucket = 0;  // ee
	if (cms2.hyp_type()[cand] == 0) bucket = 1;  // mm
	if (cms2.hyp_type()[cand] == 1 || cms2.hyp_type()[cand] == 2) bucket=2; // em
	if (bucket == 99) {
	  std::cout << "WARNING: unknown dilepton type = " << cms2.hyp_type()[cand] << std::endl;
	  continue;
	}

	// Fake Probability....this will multiply the weight of the event
	// in filling the histograms.  It will be = 1.0 (ie no effect) 
	// unless we are soing something special with W+jets
	float fProb = 1.0;
	float	fProbplErr = 1.0;
	float	fProbmiErr = 1.0;

	// on top of this cut away all electrons below 20 GeV and outside +-2.4 in eta 080317:
	if( abs(cms2.hyp_ll_id()[cand]) == 11 ) {
	  if( cms2.hyp_ll_p4()[cand].Pt()  < 20. || abs(cms2.hyp_ll_p4()[cand].eta()) > 2.4 ) continue;
	}
	if(abs(cms2.hyp_lt_id()[cand]) == 11 ) {
	  if( cms2.hyp_lt_p4()[cand].Pt()  < 20. || abs(cms2.hyp_lt_p4()[cand].eta()) > 2.4 ) continue;
	}
	  

	// on top of this cut away all muons below 20 GeV and outside +-2.4 in eta 080409:
	if( abs(cms2.hyp_ll_id()[cand]) == 13 ) {
	  if( cms2.hyp_ll_p4()[cand].Pt()  < 20. || abs(cms2.hyp_ll_p4()[cand].eta()) > 2.4 ) continue;
	}
	if(abs(cms2.hyp_lt_id()[cand]) == 13 ) {
	  if( cms2.hyp_lt_p4()[cand].Pt()  < 20. || abs(cms2.hyp_lt_p4()[cand].eta()) > 2.4 ) continue;
	}

	// Apply all the selection cuts and the lepton quality cuts...unless it is special W+jet handling for ee and emu
	// Note that if it is a mu-mu event we dont do anything special for now
	if (WjBG == 0 || cms2.hyp_type()[cand]==0) {
	  // MET cut
	  if ( !pass4Met(cand) ) continue;

	  // opposite sign
	  if ( cms2.hyp_lt_id()[cand] * cms2.hyp_ll_id()[cand] > 0 ) continue;

	  // pt leptons cut
	  if ( cms2.hyp_ll_p4()[cand].Pt() < 20.0 ) return false;
	  if ( cms2.hyp_lt_p4()[cand].Pt() < 20.0 ) return false;

	  // Muon quality cuts
	  if (abs(cms2.hyp_lt_id()[cand]) == 13) {
	    if ( !goodMuonIsolated(cms2.hyp_lt_index()[cand]) ) continue;
	    // CHEAT - remove all but true W muons:
	    if ( !trueMuonFromW(cms2.hyp_lt_index()[cand]) ) continue;
	  }
	  if (abs(cms2.hyp_ll_id()[cand]) == 13) {
	    //	      if ( !goodGlobalMuon(cms2.hyp_ll_index()[cand]) ) continue;
	    if ( !goodMuonIsolated(cms2.hyp_ll_index()[cand]) ) continue;
	    // CHEAT - remove all but true W muons:
	    if ( !trueMuonFromW(cms2.hyp_ll_index()[cand]) ) continue;
	  }
	    
	  // Electron quality cuts
	  if (abs(cms2.hyp_lt_id()[cand]) == 11) {
	    if ( !goodElectronIsolated(cms2.hyp_lt_index()[cand],true) ) continue;
	  }
	  if (abs(cms2.hyp_ll_id()[cand]) == 11) {
	    if ( !goodElectronIsolated(cms2.hyp_ll_index()[cand],true) ) continue;
	  }

	} else {  // Special handling for W+jets case
	    
	  // MET cut
	  if ( !pass4Met(cand) ) continue;
	  
	  // opposite sign
	  if ( cms2.hyp_lt_id()[cand] * cms2.hyp_ll_id()[cand] > 0) continue;

	  // If it is an emu, apply muon ID and isolation on the muon
	  int ele_index=-1;
	  if (abs(cms2.hyp_lt_id()[cand]) == 13) {
	    if ( !goodMuonIsolated(cms2.hyp_lt_index()[cand]) ) continue;
	    // CHEAT - remove all but true W muons:
	    if ( !trueMuonFromW(cms2.hyp_lt_index()[cand]) ) continue;

	    ele_index = cms2.hyp_ll_index()[cand];   // the electron is the loose lepton
	  }
	  if (abs(cms2.hyp_ll_id()[cand]) == 13) {
	    if ( !goodMuonIsolated(cms2.hyp_ll_index()[cand]) ) continue;
	    // CHEAT - remove all but true W muons:
	    if ( !trueMuonFromW(cms2.hyp_ll_index()[cand]) ) continue;

	    ele_index = cms2.hyp_lt_index()[cand];   // the electron is the tight lepton
	  }

	  // If it is an emu, make sure that the electron is a fakeable object
	  // Use the value of ele_index to tag emu (sloppy, but works)
	  if (ele_index > -1) {
	    if ( !isFakeDenominatorElectron(ele_index) ) continue;
	    if (WjBG == 2) {
	      fProb      = fakeProb(ele_index,theFakeRate,0);
	      fProbplErr = fakeProb(ele_index,theFakeRate,1);
	      fProbmiErr = fakeProb(ele_index,theFakeRate,-1);
	    }
	  } 

	  // If it is an ee, see if one of them passes the cuts and take the other one
	  // as a fakeable object.  What happens is both pass cuts?  For now we will
	  // just take the 1st combination.  THIS IS VERY VERY VERY SLOPPY.....
	  if (ele_index == -1) {	      
	    if ( goodElectronIsolated(cms2.hyp_lt_index()[cand],true) ) {   // "tight" lepton is a good lepton
	      if (isFakeDenominatorElectron(cms2.hyp_ll_index()[cand])) {         
		ele_index = cms2.hyp_ll_index()[cand];
		if (WjBG == 2) {
		  fProb      = fakeProb(ele_index,theFakeRate,0);
		  fProbplErr = fakeProb(ele_index,theFakeRate,1);
		  fProbmiErr = fakeProb(ele_index,theFakeRate,-1);
		}
	      }
	    }
	    // looks like the previous combination did not work...try the other one
	    if (ele_index == -1) {
	      if ( goodElectronIsolated(cms2.hyp_ll_index()[cand],true) ) {   // "loose" lepton is a good lepton
		if (isFakeDenominatorElectron(cms2.hyp_lt_index()[cand])) {         
		  ele_index = cms2.hyp_lt_index()[cand];
		  if (WjBG == 2) {
		    fProb      = fakeProb(ele_index,theFakeRate,0);
		    fProbplErr = fakeProb(ele_index,theFakeRate,1);
		    fProbmiErr = fakeProb(ele_index,theFakeRate,-1);
		  }
		}
	      }
	    }
	    // If we did not find a valid combination, we are forced to quit
	    if (ele_index == -1) continue;
	  }
	}

	// passed all cuts, fill histograms
	weight = weight * fProb;    // histogram weight

	// jet count
	hnJet[bucket]->Fill(cms2.hyp_njets()[cand], weight);
	hnJet[3]->Fill(cms2.hyp_njets()[cand], weight);
	hnJetplErr[bucket]->Fill(cms2.hyp_njets()[cand], weight*fProbplErr/fProb);
	hnJetplErr[3]->Fill(cms2.hyp_njets()[cand], weight*fProbplErr/fProb);
	hnJetmiErr[bucket]->Fill(cms2.hyp_njets()[cand], weight*fProbmiErr/fProb);
	hnJetmiErr[3]->Fill(cms2.hyp_njets()[cand], weight*fProbmiErr/fProb);
	  
	// lepton Pt
	if (abs(cms2.hyp_lt_id()[cand]) == 11) helePt[bucket]->Fill(cms2.hyp_lt_p4()[cand].Pt(), weight);
	if (abs(cms2.hyp_ll_id()[cand]) == 11) helePt[bucket]->Fill(cms2.hyp_ll_p4()[cand].Pt(), weight);

	if (abs(cms2.hyp_lt_id()[cand]) == 11 && abs(cms2.hyp_lt_mc_id()[cand])  == 11 ) helePtTrue[bucket]->Fill(cms2.hyp_lt_p4()[cand].Pt(), weight);
	if (abs(cms2.hyp_ll_id()[cand]) == 11 && abs(cms2.hyp_ll_mc_id()[cand])  == 11 ) helePtTrue[bucket]->Fill(cms2.hyp_ll_p4()[cand].Pt(), weight);
	if (abs(cms2.hyp_lt_id()[cand]) == 11 && abs(cms2.hyp_lt_mc_id()[cand])  != 11 ) helePtFake[bucket]->Fill(cms2.hyp_lt_p4()[cand].Pt(), weight);
	if (abs(cms2.hyp_ll_id()[cand]) == 11 && abs(cms2.hyp_ll_mc_id()[cand])  != 11 ) helePtFake[bucket]->Fill(cms2.hyp_ll_p4()[cand].Pt(), weight);

	if (abs(cms2.hyp_lt_id()[cand]) == 11 && abs(cms2.hyp_lt_mc_id()[cand])  == 11 && abs(cms2.hyp_lt_mc_motherid()[cand]) == 24 ) helePtTrueEW[bucket]->Fill(cms2.hyp_lt_p4()[cand].Pt(), weight);
	if (abs(cms2.hyp_ll_id()[cand]) == 11 && abs(cms2.hyp_ll_mc_id()[cand])  == 11 && abs(cms2.hyp_ll_mc_motherid()[cand]) == 24 ) helePtTrueEW[bucket]->Fill(cms2.hyp_ll_p4()[cand].Pt(), weight);
	if (abs(cms2.hyp_lt_id()[cand]) == 11 && (abs(cms2.hyp_lt_mc_id()[cand])  == 22 ) && (abs(cms2.hyp_lt_mc_motherid()[cand]) == 11 || abs(cms2.hyp_lt_mc_motherid()[cand]) == 13 ) ) helePtTrueFSR[bucket]->Fill(cms2.hyp_lt_p4()[cand].Pt(), weight);
	if (abs(cms2.hyp_ll_id()[cand]) == 11 && (abs(cms2.hyp_ll_mc_id()[cand])  == 22 ) && (abs(cms2.hyp_ll_mc_motherid()[cand]) == 11 || abs(cms2.hyp_ll_mc_motherid()[cand]) == 13 ) ) helePtTrueFSR[bucket]->Fill(cms2.hyp_ll_p4()[cand].Pt(), weight);

	//ErrorBands
	// lepton Pt plErr
	if (abs(cms2.hyp_lt_id()[cand]) == 11) helePtplErr[bucket]->Fill(cms2.hyp_lt_p4()[cand].Pt(), weight*fProbplErr/fProb);
	if (abs(cms2.hyp_ll_id()[cand]) == 11) helePtplErr[bucket]->Fill(cms2.hyp_ll_p4()[cand].Pt(), weight*fProbplErr/fProb);

	if (abs(cms2.hyp_lt_id()[cand]) == 11 && abs(cms2.hyp_lt_mc_id()[cand])  == 11 ) helePtTrueplErr[bucket]->Fill(cms2.hyp_lt_p4()[cand].Pt(), weight*fProbplErr/fProb);
	if (abs(cms2.hyp_ll_id()[cand]) == 11 && abs(cms2.hyp_ll_mc_id()[cand])  == 11 ) helePtTrueplErr[bucket]->Fill(cms2.hyp_ll_p4()[cand].Pt(), weight*fProbplErr/fProb);
	if (abs(cms2.hyp_lt_id()[cand]) == 11 && abs(cms2.hyp_lt_mc_id()[cand])  != 11 ) helePtFakeplErr[bucket]->Fill(cms2.hyp_lt_p4()[cand].Pt(), weight*fProbplErr/fProb);
	if (abs(cms2.hyp_ll_id()[cand]) == 11 && abs(cms2.hyp_ll_mc_id()[cand])  != 11 ) helePtFakeplErr[bucket]->Fill(cms2.hyp_ll_p4()[cand].Pt(), weight*fProbplErr/fProb);

	if (abs(cms2.hyp_lt_id()[cand]) == 11 && abs(cms2.hyp_lt_mc_id()[cand])  == 11 && abs(cms2.hyp_lt_mc_motherid()[cand]) == 24 ) helePtTrueEWplErr[bucket]->Fill(cms2.hyp_lt_p4()[cand].Pt(), weight*fProbplErr/fProb);
	if (abs(cms2.hyp_ll_id()[cand]) == 11 && abs(cms2.hyp_ll_mc_id()[cand])  == 11 && abs(cms2.hyp_ll_mc_motherid()[cand]) == 24 ) helePtTrueEWplErr[bucket]->Fill(cms2.hyp_ll_p4()[cand].Pt(), weight*fProbplErr/fProb);
	// 	  if (abs(cms2.hyp_lt_id()[cand]) == 11 && abs(cms2.hyp_lt_mc_id()[cand])  != 11 && abs(cms2.hyp_lt_mc_motherid()[cand]) != 24 ) helePtFakeEWplErr[bucket]->Fill(cms2.hyp_lt_p4()[cand].Pt(), weight*fProbplErr/fProb);
	// 	  if (abs(cms2.hyp_ll_id()[cand]) == 11 && abs(cms2.hyp_ll_mc_id()[cand])  != 11 && abs(cms2.hyp_ll_mc_motherid()[cand]) != 24 ) helePtFakeEWplErr[bucket]->Fill(cms2.hyp_ll_p4()[cand].Pt(), weight*fProbplErr/fProb);
	if (abs(cms2.hyp_lt_id()[cand]) == 11 && (abs(cms2.hyp_lt_mc_id()[cand])  == 22 ) && (abs(cms2.hyp_lt_mc_motherid()[cand]) == 11 || abs(cms2.hyp_lt_mc_motherid()[cand]) == 13 ) ) helePtTrueFSRplErr[bucket]->Fill(cms2.hyp_lt_p4()[cand].Pt(), weight*fProbplErr/fProb);
	if (abs(cms2.hyp_ll_id()[cand]) == 11 && (abs(cms2.hyp_ll_mc_id()[cand])  == 22 ) && (abs(cms2.hyp_ll_mc_motherid()[cand]) == 11 || abs(cms2.hyp_ll_mc_motherid()[cand]) == 13 ) ) helePtTrueFSRplErr[bucket]->Fill(cms2.hyp_ll_p4()[cand].Pt(), weight*fProbplErr/fProb);
	// lepton Pt miErr
	if (abs(cms2.hyp_lt_id()[cand]) == 11) helePtmiErr[bucket]->Fill(cms2.hyp_lt_p4()[cand].Pt(), weight*fProbmiErr/fProb);
	if (abs(cms2.hyp_ll_id()[cand]) == 11) helePtmiErr[bucket]->Fill(cms2.hyp_ll_p4()[cand].Pt(), weight*fProbmiErr/fProb);

	if (abs(cms2.hyp_lt_id()[cand]) == 11 && abs(cms2.hyp_lt_mc_id()[cand])  == 11 ) helePtTruemiErr[bucket]->Fill(cms2.hyp_lt_p4()[cand].Pt(), weight*fProbmiErr/fProb);
	if (abs(cms2.hyp_ll_id()[cand]) == 11 && abs(cms2.hyp_ll_mc_id()[cand])  == 11 ) helePtTruemiErr[bucket]->Fill(cms2.hyp_ll_p4()[cand].Pt(), weight*fProbmiErr/fProb);
	if (abs(cms2.hyp_lt_id()[cand]) == 11 && abs(cms2.hyp_lt_mc_id()[cand])  != 11 ) helePtFakemiErr[bucket]->Fill(cms2.hyp_lt_p4()[cand].Pt(), weight*fProbmiErr/fProb);
	if (abs(cms2.hyp_ll_id()[cand]) == 11 && abs(cms2.hyp_ll_mc_id()[cand])  != 11 ) helePtFakemiErr[bucket]->Fill(cms2.hyp_ll_p4()[cand].Pt(), weight*fProbmiErr/fProb);

	if (abs(cms2.hyp_lt_id()[cand]) == 11 && abs(cms2.hyp_lt_mc_id()[cand])  == 11 && abs(cms2.hyp_lt_mc_motherid()[cand]) == 24 ) helePtTrueEWmiErr[bucket]->Fill(cms2.hyp_lt_p4()[cand].Pt(), weight*fProbmiErr/fProb);
	if (abs(cms2.hyp_ll_id()[cand]) == 11 && abs(cms2.hyp_ll_mc_id()[cand])  == 11 && abs(cms2.hyp_ll_mc_motherid()[cand]) == 24 ) helePtTrueEWmiErr[bucket]->Fill(cms2.hyp_ll_p4()[cand].Pt(), weight*fProbmiErr/fProb);
	// 	  if (abs(cms2.hyp_lt_id()[cand]) == 11 && abs(cms2.hyp_lt_mc_id()[cand])  != 11 && abs(cms2.hyp_lt_mc_motherid()[cand]) != 24 ) helePtFakeEWmiErr[bucket]->Fill(cms2.hyp_lt_p4()[cand].Pt(), weight*fProbmiErr/fProb);
	// 	  if (abs(cms2.hyp_ll_id()[cand]) == 11 && abs(cms2.hyp_ll_mc_id()[cand])  != 11 && abs(cms2.hyp_ll_mc_motherid()[cand]) != 24 ) helePtFakeEWmiErr[bucket]->Fill(cms2.hyp_ll_p4()[cand].Pt(), weight*fProbmiErr/fProb);
	if (abs(cms2.hyp_lt_id()[cand]) == 11 && (abs(cms2.hyp_lt_mc_id()[cand])  == 22 ) && (abs(cms2.hyp_lt_mc_motherid()[cand]) == 11 || abs(cms2.hyp_lt_mc_motherid()[cand]) == 13 ) ) helePtTrueFSRmiErr[bucket]->Fill(cms2.hyp_lt_p4()[cand].Pt(), weight*fProbmiErr/fProb);
	if (abs(cms2.hyp_ll_id()[cand]) == 11 && (abs(cms2.hyp_ll_mc_id()[cand])  == 22 ) && (abs(cms2.hyp_ll_mc_motherid()[cand]) == 11 || abs(cms2.hyp_ll_mc_motherid()[cand]) == 13 ) ) helePtTrueFSRmiErr[bucket]->Fill(cms2.hyp_ll_p4()[cand].Pt(), weight*fProbmiErr/fProb);


	//cont...
	if (abs(cms2.hyp_lt_id()[cand]) == 13) hmuPt[bucket]->Fill(cms2.hyp_lt_p4()[cand].Pt(), weight);
	if (abs(cms2.hyp_ll_id()[cand]) == 13) hmuPt[bucket]->Fill(cms2.hyp_ll_p4()[cand].Pt(), weight);
	hminLepPt[bucket]->Fill(min(cms2.hyp_ll_p4()[cand].Pt(), cms2.hyp_lt_p4()[cand].Pt()), weight);
	hmaxLepPt[bucket]->Fill(max(cms2.hyp_ll_p4()[cand].Pt(), cms2.hyp_lt_p4()[cand].Pt()), weight );

	if (abs(cms2.hyp_lt_id()[cand]) == 11) helePt[3]->Fill(cms2.hyp_lt_p4()[cand].Pt(), weight);
	if (abs(cms2.hyp_ll_id()[cand]) == 11) helePt[3]->Fill(cms2.hyp_ll_p4()[cand].Pt(), weight);

	if (abs(cms2.hyp_lt_id()[cand]) == 11 &&  abs(cms2.hyp_lt_mc_id()[cand])  == 11 ) helePtTrue[3]->Fill(cms2.hyp_lt_p4()[cand].Pt(), weight);
	if (abs(cms2.hyp_ll_id()[cand]) == 11 &&  abs(cms2.hyp_ll_mc_id()[cand])  == 11 ) helePtTrue[3]->Fill(cms2.hyp_ll_p4()[cand].Pt(), weight);
	if (abs(cms2.hyp_lt_id()[cand]) == 11 &&  abs(cms2.hyp_lt_mc_id()[cand])  != 11 ) helePtFake[3]->Fill(cms2.hyp_lt_p4()[cand].Pt(), weight);
	if (abs(cms2.hyp_ll_id()[cand]) == 11 &&  abs(cms2.hyp_ll_mc_id()[cand])  != 11 ) helePtFake[3]->Fill(cms2.hyp_ll_p4()[cand].Pt(), weight);

	if (abs(cms2.hyp_lt_id()[cand]) == 11 &&  abs(cms2.hyp_lt_mc_id()[cand])  == 11 && abs(cms2.hyp_lt_mc_motherid()[cand]) == 24 ) helePtTrueEW[3]->Fill(cms2.hyp_lt_p4()[cand].Pt(), weight);
	if (abs(cms2.hyp_ll_id()[cand]) == 11 &&  abs(cms2.hyp_ll_mc_id()[cand])  == 11 && abs(cms2.hyp_ll_mc_motherid()[cand]) == 24 ) helePtTrueEW[3]->Fill(cms2.hyp_ll_p4()[cand].Pt(), weight);
	// 	  if (abs(cms2.hyp_lt_id()[cand]) == 11 &&  abs(cms2.hyp_lt_mc_id()[cand])  != 11 && abs(cms2.hyp_lt_mc_motherid()[cand]) != 24 ) helePtFakeEW[3]->Fill(cms2.hyp_lt_p4()[cand].Pt(), weight);
	// 	  if (abs(cms2.hyp_ll_id()[cand]) == 11 &&  abs(cms2.hyp_ll_mc_id()[cand])  != 11 && abs(cms2.hyp_ll_mc_motherid()[cand]) != 24 ) helePtFakeEW[3]->Fill(cms2.hyp_ll_p4()[cand].Pt(), weight);
	if (abs(cms2.hyp_lt_id()[cand]) == 11 && (abs(cms2.hyp_lt_mc_id()[cand])  == 22 ) && (abs(cms2.hyp_lt_mc_motherid()[cand]) == 11 || abs(cms2.hyp_lt_mc_motherid()[cand]) == 13 ) ) helePtTrueFSR[3]->Fill(cms2.hyp_lt_p4()[cand].Pt(), weight);
	if (abs(cms2.hyp_ll_id()[cand]) == 11 && (abs(cms2.hyp_ll_mc_id()[cand])  == 22 ) && (abs(cms2.hyp_ll_mc_motherid()[cand]) == 11 || abs(cms2.hyp_ll_mc_motherid()[cand]) == 13 ) ) helePtTrueFSR[3]->Fill(cms2.hyp_ll_p4()[cand].Pt(), weight);
	//plErr
	if (abs(cms2.hyp_lt_id()[cand]) == 11) helePtplErr[3]->Fill(cms2.hyp_lt_p4()[cand].Pt(), weight*fProbplErr/fProb);
	if (abs(cms2.hyp_ll_id()[cand]) == 11) helePtplErr[3]->Fill(cms2.hyp_ll_p4()[cand].Pt(), weight*fProbplErr/fProb);

	if (abs(cms2.hyp_lt_id()[cand]) == 11 &&  abs(cms2.hyp_lt_mc_id()[cand])  == 11 ) helePtTrueplErr[3]->Fill(cms2.hyp_lt_p4()[cand].Pt(), weight*fProbplErr/fProb);
	if (abs(cms2.hyp_ll_id()[cand]) == 11 &&  abs(cms2.hyp_ll_mc_id()[cand])  == 11 ) helePtTrueplErr[3]->Fill(cms2.hyp_ll_p4()[cand].Pt(), weight*fProbplErr/fProb);
	if (abs(cms2.hyp_lt_id()[cand]) == 11 &&  abs(cms2.hyp_lt_mc_id()[cand])  != 11 ) helePtFakeplErr[3]->Fill(cms2.hyp_lt_p4()[cand].Pt(), weight*fProbplErr/fProb);
	if (abs(cms2.hyp_ll_id()[cand]) == 11 &&  abs(cms2.hyp_ll_mc_id()[cand])  != 11 ) helePtFakeplErr[3]->Fill(cms2.hyp_ll_p4()[cand].Pt(), weight*fProbplErr/fProb);

	if (abs(cms2.hyp_lt_id()[cand]) == 11 &&  abs(cms2.hyp_lt_mc_id()[cand])  == 11 && abs(cms2.hyp_lt_mc_motherid()[cand]) == 24 ) helePtTrueEWplErr[3]->Fill(cms2.hyp_lt_p4()[cand].Pt(), weight*fProbplErr/fProb);
	if (abs(cms2.hyp_ll_id()[cand]) == 11 &&  abs(cms2.hyp_ll_mc_id()[cand])  == 11 && abs(cms2.hyp_ll_mc_motherid()[cand]) == 24 ) helePtTrueEWplErr[3]->Fill(cms2.hyp_ll_p4()[cand].Pt(), weight*fProbplErr/fProb);
	// 	  if (abs(cms2.hyp_lt_id()[cand]) == 11 &&  abs(cms2.hyp_lt_mc_id()[cand])  != 11 && abs(cms2.hyp_lt_mc_motherid()[cand]) != 24 ) helePtFakeEWplErr[3]->Fill(cms2.hyp_lt_p4()[cand].Pt(), weight*fProbplErr/fProb);
	// 	  if (abs(cms2.hyp_ll_id()[cand]) == 11 &&  abs(cms2.hyp_ll_mc_id()[cand])  != 11 && abs(cms2.hyp_ll_mc_motherid()[cand]) != 24 ) helePtFakeEWplErr[3]->Fill(cms2.hyp_ll_p4()[cand].Pt(), weight*fProbplErr/fProb);
	if (abs(cms2.hyp_lt_id()[cand]) == 11 && (abs(cms2.hyp_lt_mc_id()[cand])  == 22 ) && (abs(cms2.hyp_lt_mc_motherid()[cand]) == 11 || abs(cms2.hyp_lt_mc_motherid()[cand]) == 13 ) ) helePtTrueFSRplErr[3]->Fill(cms2.hyp_lt_p4()[cand].Pt(), weight*fProbplErr/fProb);
	if (abs(cms2.hyp_ll_id()[cand]) == 11 && (abs(cms2.hyp_ll_mc_id()[cand])  == 22 ) && (abs(cms2.hyp_ll_mc_motherid()[cand]) == 11 || abs(cms2.hyp_ll_mc_motherid()[cand]) == 13 ) ) helePtTrueFSRplErr[3]->Fill(cms2.hyp_ll_p4()[cand].Pt(), weight*fProbplErr/fProb);
	//miErr
	if (abs(cms2.hyp_lt_id()[cand]) == 11) helePtmiErr[3]->Fill(cms2.hyp_lt_p4()[cand].Pt(), weight*fProbmiErr/fProb);
	if (abs(cms2.hyp_ll_id()[cand]) == 11) helePtmiErr[3]->Fill(cms2.hyp_ll_p4()[cand].Pt(), weight*fProbmiErr/fProb);

	if (abs(cms2.hyp_lt_id()[cand]) == 11 &&  abs(cms2.hyp_lt_mc_id()[cand])  == 11 ) helePtTruemiErr[3]->Fill(cms2.hyp_lt_p4()[cand].Pt(), weight*fProbmiErr/fProb);
	if (abs(cms2.hyp_ll_id()[cand]) == 11 &&  abs(cms2.hyp_ll_mc_id()[cand])  == 11 ) helePtTruemiErr[3]->Fill(cms2.hyp_ll_p4()[cand].Pt(), weight*fProbmiErr/fProb);
	if (abs(cms2.hyp_lt_id()[cand]) == 11 &&  abs(cms2.hyp_lt_mc_id()[cand])  != 11 ) helePtFakemiErr[3]->Fill(cms2.hyp_lt_p4()[cand].Pt(), weight*fProbmiErr/fProb);
	if (abs(cms2.hyp_ll_id()[cand]) == 11 &&  abs(cms2.hyp_ll_mc_id()[cand])  != 11 ) helePtFakemiErr[3]->Fill(cms2.hyp_ll_p4()[cand].Pt(), weight*fProbmiErr/fProb);

	if (abs(cms2.hyp_lt_id()[cand]) == 11 &&  abs(cms2.hyp_lt_mc_id()[cand])  == 11 && abs(cms2.hyp_lt_mc_motherid()[cand]) == 24 ) helePtTrueEWmiErr[3]->Fill(cms2.hyp_lt_p4()[cand].Pt(), weight*fProbmiErr/fProb);
	if (abs(cms2.hyp_ll_id()[cand]) == 11 &&  abs(cms2.hyp_ll_mc_id()[cand])  == 11 && abs(cms2.hyp_ll_mc_motherid()[cand]) == 24 ) helePtTrueEWmiErr[3]->Fill(cms2.hyp_ll_p4()[cand].Pt(), weight*fProbmiErr/fProb);
	// 	  if (abs(cms2.hyp_lt_id()[cand]) == 11 &&  abs(cms2.hyp_lt_mc_id()[cand])  != 11 && abs(cms2.hyp_lt_mc_motherid()[cand]) != 24 ) helePtFakeEWmiErr[3]->Fill(cms2.hyp_lt_p4()[cand].Pt(), weight*fProbmiErr/fProb);
	// 	  if (abs(cms2.hyp_ll_id()[cand]) == 11 &&  abs(cms2.hyp_ll_mc_id()[cand])  != 11 && abs(cms2.hyp_ll_mc_motherid()[cand]) != 24 ) helePtFakeEWmiErr[3]->Fill(cms2.hyp_ll_p4()[cand].Pt(), weight*fProbmiErr/fProb);
	if (abs(cms2.hyp_lt_id()[cand]) == 11 && (abs(cms2.hyp_lt_mc_id()[cand])  == 22 ) && (abs(cms2.hyp_lt_mc_motherid()[cand]) == 11 || abs(cms2.hyp_lt_mc_motherid()[cand]) == 13 ) ) helePtTrueFSRmiErr[3]->Fill(cms2.hyp_lt_p4()[cand].Pt(), weight*fProbmiErr/fProb);
	if (abs(cms2.hyp_ll_id()[cand]) == 11 && (abs(cms2.hyp_ll_mc_id()[cand])  == 22 ) && (abs(cms2.hyp_ll_mc_motherid()[cand]) == 11 || abs(cms2.hyp_ll_mc_motherid()[cand]) == 13 ) ) helePtTrueFSRmiErr[3]->Fill(cms2.hyp_ll_p4()[cand].Pt(), weight*fProbmiErr/fProb);


	//cont...
	if (abs(cms2.hyp_lt_id()[cand]) == 13) hmuPt[3]->Fill(cms2.hyp_lt_p4()[cand].Pt(), weight);
	if (abs(cms2.hyp_ll_id()[cand]) == 13) hmuPt[3]->Fill(cms2.hyp_ll_p4()[cand].Pt(), weight);
	hminLepPt[3]->Fill(min(cms2.hyp_ll_p4()[cand].Pt(), cms2.hyp_lt_p4()[cand].Pt()), weight);
	hmaxLepPt[3]->Fill(max(cms2.hyp_ll_p4()[cand].Pt(), cms2.hyp_lt_p4()[cand].Pt()), weight );

	// lepton Phi
	if (abs(cms2.hyp_lt_id()[cand]) == 11) helePhi[bucket]->Fill(cms2.hyp_lt_p4()[cand].phi(), weight);
	if (abs(cms2.hyp_ll_id()[cand]) == 11) helePhi[bucket]->Fill(cms2.hyp_ll_p4()[cand].phi(), weight);
	if (abs(cms2.hyp_lt_id()[cand]) == 13) hmuPhi[bucket]->Fill(cms2.hyp_lt_p4()[cand].phi(), weight);
	if (abs(cms2.hyp_ll_id()[cand]) == 13) hmuPhi[bucket]->Fill(cms2.hyp_ll_p4()[cand].phi(), weight);

	if (abs(cms2.hyp_lt_id()[cand]) == 11) helePhi[3]->Fill(cms2.hyp_lt_p4()[cand].phi(), weight);
	if (abs(cms2.hyp_ll_id()[cand]) == 11) helePhi[3]->Fill(cms2.hyp_ll_p4()[cand].phi(), weight);
	if (abs(cms2.hyp_lt_id()[cand]) == 13) hmuPhi[3]->Fill(cms2.hyp_lt_p4()[cand].phi(), weight);
	if (abs(cms2.hyp_ll_id()[cand]) == 13) hmuPhi[3]->Fill(cms2.hyp_ll_p4()[cand].phi(), weight);

	// delta phi btw leptons
	double dphi = fabs(cms2.hyp_lt_p4()[cand].phi() - cms2.hyp_ll_p4()[cand].phi());
	if (dphi > TMath::Pi()) dphi = TMath::TwoPi() - dphi;
	hdphiLep[bucket]->Fill(dphi, weight);
	hdphiLep[3]->Fill(dphi, weight);

	// electron mc information
	if (abs(cms2.hyp_lt_id()[cand]) == 11) {
	  hnEleMCId[bucket]->Fill(cms2.hyp_lt_mc_id()[cand], weight);
	  hnEleMotherMCId[bucket]->Fill(cms2.hyp_lt_mc_motherid()[cand], weight);
	}
	if (abs(cms2.hyp_ll_id()[cand]) == 11) {
	  hnEleMCId[bucket]->Fill(cms2.hyp_ll_mc_id()[cand], weight);
	  hnEleMotherMCId[bucket]->Fill(cms2.hyp_ll_mc_motherid()[cand], weight);
	}

	// lepton Eta
	if (abs(cms2.hyp_lt_id()[cand]) == 11) heleEta[bucket]->Fill(cms2.hyp_lt_p4()[cand].eta(), weight);
	if (abs(cms2.hyp_ll_id()[cand]) == 11) heleEta[bucket]->Fill(cms2.hyp_ll_p4()[cand].eta(), weight);

	if (abs(cms2.hyp_lt_id()[cand]) == 11) heleEtaTrue[bucket]->Fill(cms2.hyp_lt_p4()[cand].eta(), weight);
	if (abs(cms2.hyp_ll_id()[cand]) == 11) heleEtaTrue[bucket]->Fill(cms2.hyp_ll_p4()[cand].eta(), weight);
	if (abs(cms2.hyp_lt_id()[cand]) == 11 &&  abs(cms2.hyp_lt_mc_id()[cand])  == 11 && abs(cms2.hyp_lt_mc_motherid()[cand]) == 24 ) heleEtaTrueEW[bucket]->Fill(cms2.hyp_lt_p4()[cand].eta(), weight);
	if (abs(cms2.hyp_ll_id()[cand]) == 11 &&  abs(cms2.hyp_ll_mc_id()[cand])  == 11 && abs(cms2.hyp_ll_mc_motherid()[cand]) == 24 ) heleEtaTrueEW[bucket]->Fill(cms2.hyp_ll_p4()[cand].eta(), weight);
	if (abs(cms2.hyp_lt_id()[cand]) == 11 && (abs(cms2.hyp_lt_mc_id()[cand])  == 22 ) && (abs(cms2.hyp_lt_mc_motherid()[cand]) == 11 || abs(cms2.hyp_lt_mc_motherid()[cand]) == 13 ) ) heleEtaTrueFSR[bucket]->Fill(cms2.hyp_lt_p4()[cand].eta(), weight);
	if (abs(cms2.hyp_ll_id()[cand]) == 11 && (abs(cms2.hyp_ll_mc_id()[cand])  == 22 ) && (abs(cms2.hyp_ll_mc_motherid()[cand]) == 11 || abs(cms2.hyp_ll_mc_motherid()[cand]) == 13 ) ) heleEtaTrueFSR[bucket]->Fill(cms2.hyp_ll_p4()[cand].eta(), weight);

	if (abs(cms2.hyp_lt_id()[cand]) == 13) hmuEta[bucket]->Fill(cms2.hyp_lt_p4()[cand].eta(), weight);
	if (abs(cms2.hyp_ll_id()[cand]) == 13) hmuEta[bucket]->Fill(cms2.hyp_ll_p4()[cand].eta(), weight);

	if (abs(cms2.hyp_lt_id()[cand]) == 11) heleEta[3]->Fill(cms2.hyp_lt_p4()[cand].eta(), weight);
	if (abs(cms2.hyp_ll_id()[cand]) == 11) heleEta[3]->Fill(cms2.hyp_ll_p4()[cand].eta(), weight);

	if (abs(cms2.hyp_lt_id()[cand]) == 11) heleEtaTrue[3]->Fill(cms2.hyp_lt_p4()[cand].eta(), weight);
	if (abs(cms2.hyp_ll_id()[cand]) == 11) heleEtaTrue[3]->Fill(cms2.hyp_ll_p4()[cand].eta(), weight);
	if (abs(cms2.hyp_lt_id()[cand]) == 11 &&  abs(cms2.hyp_lt_mc_id()[cand])  == 11 && abs(cms2.hyp_lt_mc_motherid()[cand]) == 24 ) heleEtaTrueEW[3]->Fill(cms2.hyp_lt_p4()[cand].eta(), weight);
	if (abs(cms2.hyp_ll_id()[cand]) == 11 &&  abs(cms2.hyp_ll_mc_id()[cand])  == 11 && abs(cms2.hyp_ll_mc_motherid()[cand]) == 24 ) heleEtaTrueEW[3]->Fill(cms2.hyp_ll_p4()[cand].eta(), weight);
	if (abs(cms2.hyp_lt_id()[cand]) == 11 && (abs(cms2.hyp_lt_mc_id()[cand])  == 22 ) && (abs(cms2.hyp_lt_mc_motherid()[cand]) == 11 || abs(cms2.hyp_lt_mc_motherid()[cand]) == 13 ) ) heleEtaTrueFSR[3]->Fill(cms2.hyp_lt_p4()[cand].eta(), weight);
	if (abs(cms2.hyp_ll_id()[cand]) == 11 && (abs(cms2.hyp_ll_mc_id()[cand])  == 22 ) && (abs(cms2.hyp_ll_mc_motherid()[cand]) == 11 || abs(cms2.hyp_ll_mc_motherid()[cand]) == 13 ) ) heleEtaTrueFSR[3]->Fill(cms2.hyp_ll_p4()[cand].eta(), weight);
	//plErr
	// lepton Eta
	if (abs(cms2.hyp_lt_id()[cand]) == 11) heleEtaplErr[bucket]->Fill(cms2.hyp_lt_p4()[cand].eta(), weight*fProbplErr/fProb);
	if (abs(cms2.hyp_ll_id()[cand]) == 11) heleEtaplErr[bucket]->Fill(cms2.hyp_ll_p4()[cand].eta(), weight*fProbplErr/fProb);

	if (abs(cms2.hyp_lt_id()[cand]) == 11) heleEtaTrueplErr[bucket]->Fill(cms2.hyp_lt_p4()[cand].eta(), weight*fProbplErr/fProb);
	if (abs(cms2.hyp_ll_id()[cand]) == 11) heleEtaTrueplErr[bucket]->Fill(cms2.hyp_ll_p4()[cand].eta(), weight*fProbplErr/fProb);
	if (abs(cms2.hyp_lt_id()[cand]) == 11 &&  abs(cms2.hyp_lt_mc_id()[cand])  == 11 && abs(cms2.hyp_lt_mc_motherid()[cand]) == 24 ) heleEtaTrueEWplErr[bucket]->Fill(cms2.hyp_lt_p4()[cand].eta(), weight*fProbplErr/fProb);
	if (abs(cms2.hyp_ll_id()[cand]) == 11 &&  abs(cms2.hyp_ll_mc_id()[cand])  == 11 && abs(cms2.hyp_ll_mc_motherid()[cand]) == 24 ) heleEtaTrueEWplErr[bucket]->Fill(cms2.hyp_ll_p4()[cand].eta(), weight*fProbplErr/fProb);
	if (abs(cms2.hyp_lt_id()[cand]) == 11 && (abs(cms2.hyp_lt_mc_id()[cand])  == 22 ) && (abs(cms2.hyp_lt_mc_motherid()[cand]) == 11 || abs(cms2.hyp_lt_mc_motherid()[cand]) == 13 ) ) heleEtaTrueFSRplErr[bucket]->Fill(cms2.hyp_lt_p4()[cand].eta(), weight*fProbplErr/fProb);
	if (abs(cms2.hyp_ll_id()[cand]) == 11 && (abs(cms2.hyp_ll_mc_id()[cand])  == 22 ) && (abs(cms2.hyp_ll_mc_motherid()[cand]) == 11 || abs(cms2.hyp_ll_mc_motherid()[cand]) == 13 ) ) heleEtaTrueFSRplErr[bucket]->Fill(cms2.hyp_ll_p4()[cand].eta(), weight*fProbplErr/fProb);

	if (abs(cms2.hyp_lt_id()[cand]) == 11) heleEtaplErr[3]->Fill(cms2.hyp_lt_p4()[cand].eta(), weight*fProbplErr/fProb);
	if (abs(cms2.hyp_ll_id()[cand]) == 11) heleEtaplErr[3]->Fill(cms2.hyp_ll_p4()[cand].eta(), weight*fProbplErr/fProb);

	if (abs(cms2.hyp_lt_id()[cand]) == 11) heleEtaTrueplErr[3]->Fill(cms2.hyp_lt_p4()[cand].eta(), weight*fProbplErr/fProb);
	if (abs(cms2.hyp_ll_id()[cand]) == 11) heleEtaTrueplErr[3]->Fill(cms2.hyp_ll_p4()[cand].eta(), weight*fProbplErr/fProb);
	if (abs(cms2.hyp_lt_id()[cand]) == 11 &&  abs(cms2.hyp_lt_mc_id()[cand])  == 11 && abs(cms2.hyp_lt_mc_motherid()[cand]) == 24 ) heleEtaTrueEWplErr[3]->Fill(cms2.hyp_lt_p4()[cand].eta(), weight*fProbplErr/fProb);
	if (abs(cms2.hyp_ll_id()[cand]) == 11 &&  abs(cms2.hyp_ll_mc_id()[cand])  == 11 && abs(cms2.hyp_ll_mc_motherid()[cand]) == 24 ) heleEtaTrueEWplErr[3]->Fill(cms2.hyp_ll_p4()[cand].eta(), weight*fProbplErr/fProb);
	if (abs(cms2.hyp_lt_id()[cand]) == 11 && (abs(cms2.hyp_lt_mc_id()[cand])  == 22 ) && (abs(cms2.hyp_lt_mc_motherid()[cand]) == 11 || abs(cms2.hyp_lt_mc_motherid()[cand]) == 13 ) ) heleEtaTrueFSRplErr[3]->Fill(cms2.hyp_lt_p4()[cand].eta(), weight*fProbplErr/fProb);
	if (abs(cms2.hyp_ll_id()[cand]) == 11 && (abs(cms2.hyp_ll_mc_id()[cand])  == 22 ) && (abs(cms2.hyp_ll_mc_motherid()[cand]) == 11 || abs(cms2.hyp_ll_mc_motherid()[cand]) == 13 ) ) heleEtaTrueFSRplErr[3]->Fill(cms2.hyp_ll_p4()[cand].eta(), weight*fProbplErr/fProb);
	//miErr
	// lepton Eta
	if (abs(cms2.hyp_lt_id()[cand]) == 11) heleEtamiErr[bucket]->Fill(cms2.hyp_lt_p4()[cand].eta(), weight*fProbmiErr/fProb);
	if (abs(cms2.hyp_ll_id()[cand]) == 11) heleEtamiErr[bucket]->Fill(cms2.hyp_ll_p4()[cand].eta(), weight*fProbmiErr/fProb);

	if (abs(cms2.hyp_lt_id()[cand]) == 11) heleEtaTruemiErr[bucket]->Fill(cms2.hyp_lt_p4()[cand].eta(), weight*fProbmiErr/fProb);
	if (abs(cms2.hyp_ll_id()[cand]) == 11) heleEtaTruemiErr[bucket]->Fill(cms2.hyp_ll_p4()[cand].eta(), weight*fProbmiErr/fProb);
	if (abs(cms2.hyp_lt_id()[cand]) == 11 &&  abs(cms2.hyp_lt_mc_id()[cand])  == 11 && abs(cms2.hyp_lt_mc_motherid()[cand]) == 24 ) heleEtaTrueEWmiErr[bucket]->Fill(cms2.hyp_lt_p4()[cand].eta(), weight*fProbmiErr/fProb);
	if (abs(cms2.hyp_ll_id()[cand]) == 11 &&  abs(cms2.hyp_ll_mc_id()[cand])  == 11 && abs(cms2.hyp_ll_mc_motherid()[cand]) == 24 ) heleEtaTrueEWmiErr[bucket]->Fill(cms2.hyp_ll_p4()[cand].eta(), weight*fProbmiErr/fProb);
	if (abs(cms2.hyp_lt_id()[cand]) == 11 && (abs(cms2.hyp_lt_mc_id()[cand])  == 22 ) && (abs(cms2.hyp_lt_mc_motherid()[cand]) == 11 || abs(cms2.hyp_lt_mc_motherid()[cand]) == 13 ) ) heleEtaTrueFSRmiErr[bucket]->Fill(cms2.hyp_lt_p4()[cand].eta(), weight*fProbmiErr/fProb);
	if (abs(cms2.hyp_ll_id()[cand]) == 11 && (abs(cms2.hyp_ll_mc_id()[cand])  == 22 ) && (abs(cms2.hyp_ll_mc_motherid()[cand]) == 11 || abs(cms2.hyp_ll_mc_motherid()[cand]) == 13 ) ) heleEtaTrueFSRmiErr[bucket]->Fill(cms2.hyp_ll_p4()[cand].eta(), weight*fProbmiErr/fProb);

	if (abs(cms2.hyp_lt_id()[cand]) == 11) heleEtamiErr[3]->Fill(cms2.hyp_lt_p4()[cand].eta(), weight*fProbmiErr/fProb);
	if (abs(cms2.hyp_ll_id()[cand]) == 11) heleEtamiErr[3]->Fill(cms2.hyp_ll_p4()[cand].eta(), weight*fProbmiErr/fProb);

	if (abs(cms2.hyp_lt_id()[cand]) == 11) heleEtaTruemiErr[3]->Fill(cms2.hyp_lt_p4()[cand].eta(), weight*fProbmiErr/fProb);
	if (abs(cms2.hyp_ll_id()[cand]) == 11) heleEtaTruemiErr[3]->Fill(cms2.hyp_ll_p4()[cand].eta(), weight*fProbmiErr/fProb);
	if (abs(cms2.hyp_lt_id()[cand]) == 11 &&  abs(cms2.hyp_lt_mc_id()[cand])  == 11 && abs(cms2.hyp_lt_mc_motherid()[cand]) == 24 ) heleEtaTrueEWmiErr[3]->Fill(cms2.hyp_lt_p4()[cand].eta(), weight*fProbmiErr/fProb);
	if (abs(cms2.hyp_ll_id()[cand]) == 11 &&  abs(cms2.hyp_ll_mc_id()[cand])  == 11 && abs(cms2.hyp_ll_mc_motherid()[cand]) == 24 ) heleEtaTrueEWmiErr[3]->Fill(cms2.hyp_ll_p4()[cand].eta(), weight*fProbmiErr/fProb);
	if (abs(cms2.hyp_lt_id()[cand]) == 11 && (abs(cms2.hyp_lt_mc_id()[cand])  == 22 ) && (abs(cms2.hyp_lt_mc_motherid()[cand]) == 11 || abs(cms2.hyp_lt_mc_motherid()[cand]) == 13 ) ) heleEtaTrueFSRmiErr[3]->Fill(cms2.hyp_lt_p4()[cand].eta(), weight*fProbmiErr/fProb);
	if (abs(cms2.hyp_ll_id()[cand]) == 11 && (abs(cms2.hyp_ll_mc_id()[cand])  == 22 ) && (abs(cms2.hyp_ll_mc_motherid()[cand]) == 11 || abs(cms2.hyp_ll_mc_motherid()[cand]) == 13 ) ) heleEtaTrueFSRmiErr[3]->Fill(cms2.hyp_ll_p4()[cand].eta(), weight*fProbmiErr/fProb);


	//cont...
	if (abs(cms2.hyp_lt_id()[cand]) == 13) hmuEta[3]->Fill(cms2.hyp_lt_p4()[cand].eta(), weight);
	if (abs(cms2.hyp_ll_id()[cand]) == 13) hmuEta[3]->Fill(cms2.hyp_ll_p4()[cand].eta(), weight);

	// dilepton mass
	hdilMass[bucket]->Fill(cms2.hyp_p4()[cand].mass(), weight);
	hdilMassTightWindow[bucket]->Fill(cms2.hyp_p4()[cand].mass(), weight);
	hdilMass[3]->Fill(cms2.hyp_p4()[cand].mass(), weight);
	hdilMassTightWindow[3]->Fill(cms2.hyp_p4()[cand].mass(), weight);

	// dilepton pt
	hdilPt[bucket]->Fill(cms2.hyp_p4()[cand].Pt(), weight);
	hdilPt[3]->Fill(cms2.hyp_p4()[cand].Pt(), weight);

	// Met and Met phi
	hmet[bucket]->Fill(cms2.hyp_met()[cand], weight);      
	hmetPhi[bucket]->Fill(cms2.hyp_metPhi()[cand], weight);      
	hmet[3]->Fill(cms2.hyp_met()[cand], weight);      
	hmetPhi[3]->Fill(cms2.hyp_metPhi()[cand], weight);      

	// Met vs dilepton Pt
	hmetVsDilepPt[bucket]->Fill(cms2.hyp_met()[cand], cms2.hyp_p4()[cand].Pt(), weight);
	hmetVsDilepPt[3]->Fill(cms2.hyp_met()[cand], cms2.hyp_p4()[cand].Pt(), weight);
    
	// Met over dilepton Pt vs deltaphi btw the two
	double dphi2 = fabs(cms2.hyp_p4()[cand].phi() - cms2.hyp_metPhi()[cand]);
	if (dphi2 > TMath::Pi()) dphi2 = TMath::TwoPi() - dphi2;
	hmetOverPtVsDphi[bucket]->Fill(cms2.hyp_met()[cand]/cms2.hyp_p4()[cand].Pt(), dphi2, weight);
	hmetOverPtVsDphi[3]->Fill(cms2.hyp_met()[cand]/cms2.hyp_p4()[cand].Pt(), dphi2, weight);

      }
    }
  }

  std::cout << std::endl;
  if ( nEventsChain != nEventsTotal ) {
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  }
  
  std::cout << "Prefix: " << prefix << " processed: " << nEventsTotal
	    << " Events, found: " << duplicates_total_n 
	    << " Duplicates and selected: " << nCandidatesSelected
	    << " Candidates."
	    << endl << endl;

  return 0;
}
