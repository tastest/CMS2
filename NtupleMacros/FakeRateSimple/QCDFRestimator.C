/* Usage:
   root[0] .L ttDilCounts_looper.C++
   root [1] TFile *_file0 = TFile::Open("ntuple_file.root")
   root [2] TChain *chain = new TChain("Events")
   root [3] chain->Add("ntuple_file.root")
   root [4] ttDilCounts_looper a 
   root [5] a.ScanChain(chain) // will give the same results
*/
#include <iostream>
#include <vector>
#include <map>

#include "TChain.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TROOT.h"

#include "TChain.h"
#include "TFile.h"
#include "TChainElement.h"
#include "TDirectory.h"
#include "TROOT.h"

//from slava
#include "TH1F.h"
#include "TH2F.h"
#include "Math/LorentzVector.h"
#include "TMath.h"
#include <algorithm>
#include "TRandom2.h"
#include <fstream>
#include "TChain.h"

#include "CMS2.cc"
#include "../CORE/selections.cc"
#include "../CORE/utilities.cc"
#include "../CORE/electronSelections.cc"
#include "../CORE/muonSelections.cc"

#include "../Tools/fakerates.cc"
#include "QCDFRestimator.h"

// dbarge
ofstream logfile("category4_printout.txt");


//------------------------------------------------------------

typedef vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >  VofP4;

// this is Jake's magic to sort jets by Pt
//Bool_t comparePt(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > lv1,
//                 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > lv2) {
//  return lv1.pt() > lv2.pt();
//}


//------------------------------------------------------------
Float_t QCDFRestimator::GetValueTH2F(Float_t x, Float_t y, TH2F* h) {
 
  Int_t binx = h->GetXaxis()->FindBin(x);
  Int_t biny = h->GetYaxis()->FindBin(y);
  return h->GetBinContent(binx, biny);
 
} //------------------------------------------------------------


// dbarge
bool isNumEl(int iEl){
  Double_t pt = cms2.els_p4()[iEl].Pt();
  Double_t eta = cms2.els_p4()[iEl].Eta();

  if(
     (pt >= 10.) &&
     (fabs(eta)<=2.4) &&
     // dEtaIn < 0.007, 0.010 (EB, EE)
     // dPhiIn < 0.020, 0.025 (EB, EE)
     // hoe < 0.01, 0.01 (EB, EE)
     // sigmaIEtaIEta < N/A , 0.03 (EB, EE)
     // E2x5Max /E5x5 > 0.90, N/A (EB, EE) 
     electronId_cand01(iEl) &&
     electronImpact_cand01(iEl) &&             	    // d0corr < .02
     electronIsolation_relsusy_cand1(iEl,true) < 0.1 &&   // relative isolation < .1
     !isFromConversionPartnerTrack(iEl) &&        	// dist < .02 dcot < .02
		 electronId_noMuon(iEl)
     ){
    return true;
  } else {
    return false;
  }
}
bool isDenomEl(int iEl){
  Double_t pt = cms2.els_p4()[iEl].Pt();
  Double_t eta = cms2.els_p4()[iEl].Eta();
  if(
     (pt >= 10.) &&
     (fabs(eta)<=2.4) &&
     // electronId_cand01(iEl) &&
     electronImpact_cand01(iEl) &&             	    // d0corr < .02
     electronIsolation_relsusy_cand1(iEl,true) < 0.1 &&   // relative isolation < .1
     !isFromConversionPartnerTrack(iEl) &&       	  // dist < .02 dcot < .02
		 electronId_noMuon(iEl)
     ){
    return true;
  } else{
    return false;
  }
}

// dbarge
/* muon numerator aka selection definition */
bool isNumMu(int iMu){
  // Pt & Eta cuts
  Double_t pt   = cms2.mus_p4().at(iMu).Pt();
  Double_t eta  = cms2.mus_p4().at(iMu).Eta();
  if( pt < 10.0 || fabs(eta) > 2.4 ) return false;
  /* muonId cuts - numerator */
  if(
      // abs(eta) <= 2.4
      // global chisq per degree of freedom < 10
      // global muon
      // tracker muon
      // # of valid hits in silicon fit >= 11
      // ECalE < 4
      // HcalE < 6
      // glb fit has hits
      // d0 from beamspot <= .02
      // ISO < 0.2
      muonId(iMu, Nominal)
     ){
    return true;
  } else {
    return false;
  }
}
  
/* muon denominator aka fakeable object definition */
bool isDenomMu(int iMu){
  // Pt & Eta cuts
  Double_t pt   = cms2.mus_p4().at(iMu).Pt();
  Double_t eta  = cms2.mus_p4().at(iMu).Eta();
  if( pt < 10.0 || fabs(eta) > 2.4 ) return false;
  /* muonId cuts - denomintor */
  if(

      /* dummy muon denominator definition */
      // what was done in AN 2009/041:
      // relax chisq/N to 20
      // remove requirement on # silicon hits
      // relax isolation

      TMath::Abs(cms2.mus_p4()[iMu].eta()) <= 2.4 &&								        //eta
      //cms2.mus_gfit_chi2().at(iMu)/cms2.mus_gfit_ndof().at(iMu) < 10 &&		//glb fit chisq
      cms2.mus_gfit_chi2().at(iMu)/cms2.mus_gfit_ndof().at(iMu) < 20 &&		  //glb fit chisq
      (	( (cms2.mus_type().at(iMu) ) & (1<<1) ) != 0	) &&						      // global muon
      (	( (cms2.mus_type().at(iMu) ) & (1<<2)) != 0 ) &&						        // tracker muon
      //cms2.mus_validHits().at(iMu) >= 11 &&											          // # of tracker hits
      cms2.mus_iso_ecalvetoDep().at(iMu) <= 4 &&									          // ECalE < 4 
      cms2.mus_iso_hcalvetoDep().at(iMu) <= 6 &&									          // HCalE < 6 
      cms2.mus_gfit_validSTAHits().at(iMu) > 0 &&									          // Glb fit must have hits in mu chambers
      TMath::Abs(cms2.mus_d0corr().at(iMu)) <= 0.02 &&							        // d0 from beamspot
      muonIsoValue(iMu) <= 0.1															                // Isolation cut
     ){
    return true;
  } else{
    return false;
  }
}



bool QCDFRestimator::testJetsForElectrons( vector<LorentzVector>& jetP4, 
                                           const LorentzVector& elP4) {

  for(unsigned int i = 0; i < jetP4.size() ; i++) {
    
    float elphi  = elP4.Phi();
    float jetphi = jetP4[i].Phi();
    
    float eleta  = elP4.Eta();
    float jeteta = jetP4[i].Eta();
    
    float dphi = elphi - jetphi;
    float deta = eleta - jeteta;
    if(fabs(dphi) > TMath::Pi() ) dphi = 2*TMath::Pi() - fabs(dphi);
   
    double dR = sqrt(dphi*dphi + deta*deta);
    if (dR < 0.4) 
      return true;
  }
    
  return false;
}

//------------------------------------------------------------
			    
int QCDFRestimator::ScanChainAppTest ( TChain* chain, TString prefix, float kFactor, int prescale){

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();
  TH2F *h_FR[2];
  TH2F *h_FRErr[2];
  h_FR[0]     = (TH2F*)rootdir->Get("QCD_FRptvseta_el");
  if(h_FR[0] == NULL)
    cout << "h_FR0 not found " << endl;
  h_FRErr[0]  = (TH2F*)rootdir->Get("QCD_FRErrptvseta_el");
  h_FR[1]     = (TH2F*)rootdir->Get("QCD_FRptvseta_mu");
  h_FRErr[1]  = (TH2F*)rootdir->Get("QCD_FRErrptvseta_mu");

  using namespace std;
  
  // bins
  Float_t pt[17] = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160};
  Float_t eta[13] = {-3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3};
  Float_t nbins[7] = {0,1,2,3,4,5,6};

  //book Histograms
  char *suffix[2] =  {"el", "mu"};
  for ( unsigned int suf = 0; suf < 2; ++suf ) {
    // nJet plot
    h_predictednJets[suf]  = new TH1F( Form("%s_predictednJets_%s",prefix.Data(), suffix[suf]),
                                       "predicted NJet distribution, FO object", 6, nbins);
    h_predictednJets[suf]->Sumw2();
    h_actualnJets[suf]     = new TH1F(	Form("%s_actualnJets_%s",prefix.Data(), suffix[suf]), 
                                        "actual NJet distribution", 6, nbins);
    h_actualnJets[suf]->Sumw2();
    h_nJets3D[suf]          = new TH3F(Form("%s_nJets3D_%s",prefix.Data(), suffix[suf]),
                                       "3D histo to store error info", 6, nbins, 12, eta, 16, pt);
    h_nJets3D[suf]->Sumw2();

    // MC category plot
    h_predictedTrueCat[suf]  = new TH1F( Form("%s_predictedTrueCat_%s",prefix.Data(), suffix[suf]),
                                         "predicted MC category distribution, FO object", 6, nbins);
    h_predictedTrueCat[suf]->Sumw2();
    h_actualTrueCat[suf]     = new TH1F(	Form("%s_actualTrueCat_%s",prefix.Data(), suffix[suf]), 
                                          "actual MC category distribution", 6, nbins);
    h_actualTrueCat[suf]->Sumw2();
    h_TrueCat3D[suf]          = new TH3F(Form("%s_TrueCat3D_%s",prefix.Data(), suffix[suf]),
                                         "3D histo to store error info", 6, nbins, 12, eta, 16, pt);
    h_TrueCat3D[suf]->Sumw2();

    h_truecomposition_num[suf] = new TH1F(  Form("%s_truecomposition_num_%s",prefix.Data(), suffix[suf]),
                                            Form("%s_truecomposition_num_%s",prefix.Data(), suffix[suf]), 5, -0.5, 4.5);
    h_truecomposition_num[suf]->Sumw2();
    h_truecomposition_denom[suf] = new TH1F(  Form("%s_truecomposition_denom_%s",prefix.Data(), suffix[suf]),
                                              Form("%s_truecomposition_denom_%s",prefix.Data(), suffix[suf]), 5, -0.5, 4.5);
    h_truecomposition_denom[suf]->Sumw2();
    h_truecomposition_ratio[suf] = new TH1F(  Form("%s_truecomposition_ratio_%s",prefix.Data(), suffix[suf]),
                                              Form("%s_truecomposition_ratio_%s",prefix.Data(), suffix[suf]), 5, -0.5, 4.5);
    h_truecomposition_ratio[suf]->Sumw2();
  }

  //Event Loop
  int nAfterPrescale = 0;
  TObjArray *listOfFiles = chain->GetListOfFiles();
  unsigned int nEventsChain=chain->GetEntries();
  unsigned int nEventsTotal = 0;
  // file loop
  TIter fileIter(listOfFiles);
  int nAllEvents = 0;
  map<int,int> m_events;
	
  while(TChainElement *currentFile = (TChainElement*)fileIter.Next() ) {
    TFile f(currentFile->GetTitle());
    TTree *tree = (TTree*)f.Get("Events");
    cms2.Init(tree);
    unsigned int nEntries = tree->GetEntries();
    unsigned int nLoop = nEntries;
    unsigned int z;
    for( z = 0; z < nLoop; z++) {
      nAllEvents++;
      // Progress feedback to the user
      int mod = 100000;
      int iz = nAllEvents/mod;
      if (nAllEvents-mod*iz == 0) cout << "Processing event " << nAllEvents+1 << " of sample " << prefix << endl;
      cms2.GetEntry(z);
      ++nEventsTotal;

      //weight
      float weight = kFactor*cms2.evt_scale1fb();
      weight = 1.0;
		
      for(unsigned int iHyp = 0; iHyp < cms2.hyp_p4().size(); iHyp++) {
	
        //opposite sign - Oli + Ingo
        if( cms2.hyp_lt_id()[iHyp]*cms2.hyp_ll_id()[iHyp] < 0) continue;
		        
        int nJets = cms2.hyp_jets_p4()[iHyp].size();
        nJets = min(nJets, 5);
        int hyp_type = cms2.hyp_type()[iHyp];
		
        //emu case
        if(hyp_type == 0 || hyp_type == 3) continue;
        int iEl = 0;
        int iMu = 0;
        if(hyp_type == 2) {
          iEl = cms2.hyp_lt_index()[iHyp];
          iMu = cms2.hyp_ll_index()[iHyp];
        } 
        if (hyp_type == 1) {
          iEl = cms2.hyp_ll_index()[iHyp];
          iMu = cms2.hyp_lt_index()[iHyp];
        } 
		
        /* electrons */

        //look only in true W->mu events. ONLY FOR W+JETS SAMPLE!!!
        //Do this by looking at the gen-level info only
        //if(isTrueLeptonfromW(13)) {
        //look only at global mus!
        //don't look at cases where the electron comes from a mu (mu->mu+gamma->electron)
        //	if(trueMuonFromW_WJets(iMu) && (2 & cms2.mus_type()[iMu]) && cms2.mus_p4()[iMu].Pt() > 10.) {
        if( prefix.Contains("WJets") && trueMuonFromW(iMu) || 
            prefix.Contains("TTbar") && ttbarconstituents( iHyp) == 2 && trueMuonFromW(iMu) ) {
          if( TMath::Max(cms2.hyp_lt_p4()[iHyp].pt(),cms2.hyp_ll_p4()[iHyp].pt()) > 20. ) {
            if( TMath::Min(cms2.hyp_lt_p4()[iHyp].pt(),cms2.hyp_ll_p4()[iHyp].pt()) > 10.) {
              Double_t pt = cms2.els_p4()[iEl].Pt();
              Double_t eta = cms2.els_p4()[iEl].Eta();
              if(trueGammaFromMuon(iEl))	continue;
		
              //dbarge
              if( isNumEl(iEl) ){	// if a numerator electron
                h_actualnJets[0]->Fill(nJets, weight);
                int cat = elFakeMCCategory(iEl);
                h_actualTrueCat[0]->Fill(cat, weight);
                h_truecomposition_num[0]->Fill( cat, weight);
                if( cat == 4 ){
                  logfile << "WJnum:\t" << cms2.els_mc_id().at(iEl) << "\t" << cms2.els_mc_motherid().at(iEl) << endl;
                }
              } else { // not a numerator electron
                if( isDenomEl(iEl) ){ // is a fakeable object
                  Float_t FR    = GetValueTH2F(fabs(eta), min(pt,149.0), h_FR[0]);
                  Float_t FRErr = GetValueTH2F(fabs(eta), pt, h_FRErr[0]);
                  h_predictednJets[0]->Fill(nJets, weight*FR/(1-FR));
                  h_nJets3D[0]        ->Fill(nJets, fabs(eta), pt, weight*FRErr);
                  int cat = elFakeMCCategory(iEl);
                  h_predictedTrueCat[0]->Fill(cat, weight*FR/(1-FR));
                  h_TrueCat3D[0]        ->Fill(cat, fabs(eta), pt, weight*FRErr);
                  if( cat == 4 ){
                    logfile << "WJden:\t" << cms2.els_mc_id().at(iEl) << "\t" << cms2.els_mc_motherid().at(iEl) << endl;
                  }
                }
              }
							// Derek's plots should not filled using the FR/(1-FR) logic, confusing)
							if( isDenomEl(iEl) ){ // is a fakeable object
								h_truecomposition_denom[0]->Fill( cat, weight);
							}
		
            } // ll pt > 10
          } // lt pt > 20 
        } // is true muon from W
        /* end electrons */

        /* muons */
        if(trueElectronFromW(iEl) && cms2.els_p4().at(iEl).Pt() > 10.) {  // electron from W with Pt >= 10
          //if(!isNumElSUSY09(iEl))  continue; // Require the electron to be isolated and good
          if( !isNumEl(iEl) )  continue;  // electron passes selection 

          //if(!isFakeableMuSUSY09(iMu)) continue;
          if( !isDenomMu(iMu) ) continue; // denomintaor aka fakeable object muons
          if(trueMuonFromW(iMu) )	cout << "SHOULD NEVER GET HERE!!!!!!" << endl;

          // fill denominators aka fakeable objects
          Double_t pt = cms2.mus_p4().at(iMu).Pt();
          Double_t eta = cms2.mus_p4().at(iMu).Eta();
          Float_t FR    = GetValueTH2F(fabs(eta), pt, h_FR[1]);
          Float_t FRErr = GetValueTH2F(fabs(eta), pt, h_FRErr[1]);
          h_predictednJets[1]->Fill(nJets, weight*FR);
          h_nJets3D[1]->Fill(nJets, fabs(eta), pt, weight*FRErr);
          h_predictedTrueCat[1]->Fill(elFakeMCCategory(iEl), weight*FR);
          h_TrueCat3D[1]->Fill(elFakeMCCategory(iEl), fabs(eta), pt, weight*FRErr);
          //h_truecomposition_denom[1]->Fill( muFakeMCCategory(iMu), weight); // true composition

          //if( !isNumMuSUSY09(iMu) ) continue;
          if( !isNumMu(iMu) ) continue; // numerator aka selected muons

          // fill
          h_actualnJets[1]->Fill(nJets, weight);
          h_actualTrueCat[1]->Fill(elFakeMCCategory(iEl), weight);
          //h_truecomposition_num[1]->Fill( muFakeMCCategory(iMu), weight); // true composition
        } // is true electron from W
        /* end muons */


      } // end hyp loop
    } // end event loop
  } // end file loop

  if ( nEventsChain != nEventsTotal ) {
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  }

  for(unsigned int i = 0; i < 2; i++) {
    //do the errors for the WJets
    Float_t totalErr = 0.;
    for(unsigned int ieta = 1; ieta < h_nJets3D[i]->GetNbinsY() + 1; ieta++) {
      for(unsigned int ipt = 1; ipt < h_nJets3D[i]->GetNbinsZ() + 1; ipt++) {
        Float_t temp23 = 0.;  
        for(unsigned int iJet = 1; iJet < h_nJets3D[i]->GetNbinsX() + 1; iJet++) {
          temp23 = temp23 + h_nJets3D[i]->GetBinContent(iJet, ieta, ipt);
        }
        totalErr = pow(temp23,2) + totalErr;
      }
    }
    cout << "****** Error for  " << suffix[i] << sqrt(totalErr) << endl;
    for(unsigned int iJet = 1; iJet < h_nJets3D[i]->GetNbinsX() + 1; iJet++) {
      Float_t err2 = 0.;
      for(unsigned int ieta = 1; ieta < h_nJets3D[i]->GetNbinsY() + 1; ieta++) {
        for(unsigned int ipt = 1; ipt < h_nJets3D[i]->GetNbinsZ() + 1; ipt++) {
          Float_t temp = h_nJets3D[i]->GetBinContent(iJet, ieta, ipt);
          err2 = err2 + pow(temp,2);
        }
      }
      Float_t err = sqrt(err2);
      h_predictednJets[i]->SetBinError(iJet, err);
    }
	
    for(unsigned int ieta = 1; ieta < h_TrueCat3D[i]->GetNbinsY() + 1; ieta++) {
      for(unsigned int ipt = 1; ipt < h_TrueCat3D[i]->GetNbinsZ() + 1; ipt++) {
        Float_t temp23 = 0.;  
        for(unsigned int iMainBin = 1; iMainBin < h_TrueCat3D[i]->GetNbinsX() + 1; iMainBin++) {
          temp23 = temp23 + h_TrueCat3D[i]->GetBinContent(iMainBin, ieta, ipt);
        }
        totalErr = pow(temp23,2) + totalErr;
      }
    }
    cout << "****** Error for  " << suffix[i] << sqrt(totalErr) << endl;
    for(unsigned int iMainBin = 1; iMainBin < h_TrueCat3D[i]->GetNbinsX() + 1; iMainBin++) {
      Float_t err2 = 0.;
      for(unsigned int ieta = 1; ieta < h_TrueCat3D[i]->GetNbinsY() + 1; ieta++) {
        for(unsigned int ipt = 1; ipt < h_TrueCat3D[i]->GetNbinsZ() + 1; ipt++) {
          Float_t temp = h_TrueCat3D[i]->GetBinContent(iMainBin, ieta, ipt);
          err2 = err2 + pow(temp,2);
        }
      }
      Float_t err = sqrt(err2);
      h_predictedTrueCat[i]->SetBinError(iMainBin, err);
    }
  }//lepton flavor loop
      
  
  std::cout<<"Done with "<<prefix<<std::endl;
  rootdir = gDirectory->GetDirectory("Rint:"); 
  rootdir->cd(); 
  

  //dbarge
  h_truecomposition_ratio[0]->Divide( 	h_truecomposition_num[0], h_truecomposition_denom[0], 1, 1, "B");
  h_truecomposition_ratio[1]->Divide( 	h_truecomposition_num[1], h_truecomposition_denom[1], 1, 1, "B");
  
  return 0;
}


//--------------------------------------------------------------------------------

int QCDFRestimator::ScanChainQCD ( TChain* chain, TString prefix, float kFactor, int prescale,
                                   float pthatmin, float pthatmax){
  
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  using namespace std;
  
  //book Histograms
  bookHistos(prefix.Data());

  //Event Loop
  int nAfterPrescale = 0;
  TObjArray *listOfFiles = chain->GetListOfFiles();

  unsigned int nEventsChain=chain->GetEntries();
  unsigned int nEventsTotal = 0;
  
  
  // file loop
  TIter fileIter(listOfFiles);
  int nAllEvents = 0;
  map<int,int> m_events;
  while(TChainElement *currentFile = (TChainElement*)fileIter.Next() ) {
		//
		// hardcoded QCD bin pthat treatment
		// based on filename
		// kind of a hack
		// if more qcd bins are added
		// enlarge treatment here
		TString filename = currentFile->GetTitle();
		if ( filename.Contains("_Pt30_") ) {
			pthatmin = 0.0;
			pthatmax = 80.0;
		} else if ( filename.Contains("_Pt80_") ) {
			pthatmin = 30.0;
			pthatmax = 99999999999999999999999999999.;
		}

    TFile f(filename.Data());
    TTree *tree = (TTree*)f.Get("Events");
    cms2.Init(tree);
    unsigned int nEntries = tree->GetEntries();
    unsigned int nLoop = nEntries;
    
    unsigned int z;
    float maxaa = 0.;

    for( z = 0; z < nLoop; z++) {
      nAllEvents++;
      // Progress feedback to the user
      int mod = 100000;
      int iz = nAllEvents/mod;
      if (nAllEvents-mod*iz == 0) cout << "Processing event " << nAllEvents+1 << " of sample " << prefix << endl;
      cms2.GetEntry(z);
      ++nEventsTotal;


      //cout << pthatmin << " " << cms2.genps_pthat() << " " << pthatmax << endl;

      //if the pthat of the event is not within acceptable limits, quit
      //carefull.....WJets doesn't have the genps_pthat variable
      //if(cms2.genps_pthat() < pthatmin || cms2.genps_pthat() > pthatmax){
      if(cms2.genps_qScale() < pthatmin || cms2.genps_qScale() > pthatmax){
        continue;
      }
      
      float weight = kFactor*cms2.evt_scale1fb()*0.1;
			// deactivated 2/5/10
      // weight = 1.0; //weights are the same, screws the Binomial errors up
      
			for(int iEl = 0 ; iEl < cms2.els_p4().size(); iEl++) {

        if(!isDenomEl(iEl)) continue;
        Double_t pt = cms2.els_p4()[iEl].Pt();
        Double_t eta = cms2.els_p4()[iEl].Eta();
        Double_t phi = cms2.els_p4()[iEl].Phi();

        int cat = elFakeMCCategory(iEl);
        h_truecomposition_denom[0]->Fill( cat, weight);
        if( cat == 4 ){
          logfile 	<< "QCDden:\t" << cms2.els_mc_id().at(iEl) << "\t" 
                    << cms2.els_mc_motherid().at(iEl) << endl;
        }
        h_FOptvseta[0]->Fill(fabs(eta), min(pt,149.0), weight);
        h_FOpt[0]     ->Fill(min(pt,149.0), weight);
        h_FOeta[0]    ->Fill(fabs(eta), weight);


        //figure out the mcid of the closest status==3 particle
        Double_t dR = 5.0;
        unsigned int matchedId = 99999;
        for(int genIdx = 0; genIdx < cms2.genps_id().size(); genIdx++) {
          int mcid = cms2.genps_id()[genIdx];
          if( abs(mcid)==2212) continue;
          if( abs(mcid) == 6 || abs(mcid)==24) continue;
          Double_t deta = eta - cms2.genps_p4()[genIdx].Eta();
          Double_t dphi = fabs(phi - cms2.genps_p4()[genIdx].Phi());
          if( dphi > TMath::Pi()) 
            dphi = 2*TMath::Pi() - dphi;
          if(sqrt(deta*deta + dphi*dphi) < dR) {
            dR = sqrt(deta*deta + dphi*dphi);
            matchedId = abs(mcid);
          }
        }
        if(matchedId == 5 || matchedId == 4 ) {
          h_FOhfpt[0] ->Fill(min(pt,149.0), weight);
          h_FOhfeta[0]->Fill(fabs(eta), weight);
        }
        if(matchedId < 4) {
          h_FOlqpt[0]->Fill(min(pt,149.0), weight);
          h_FOlqeta[0]->Fill(fabs(eta), weight);
        }
        if(matchedId == 21) {
          h_FOgpt[0]->Fill(min(pt,149.0), weight);
          h_FOgeta[0]->Fill(fabs(eta), weight);
        }
        if(dR < 1.0) 
          h_FOmc3Id[0]->Fill(matchedId, weight);
        h_FOmc3dR[0]->Fill(dR, weight);
	  

        if(!isNumEl(iEl)) continue;

        // int cat = elFakeMCCategory(iEl);
        // cat already defined above
        h_truecomposition_num[0]->Fill( cat, weight);
        if( cat == 4 ){
          logfile << "QCDnum:\t" << cms2.els_mc_id().at(iEl) << "\t" << cms2.els_mc_motherid().at(iEl) << endl;
        }

        h_numptvseta[0]->Fill(fabs(eta), min(pt,149.0), weight);
        h_numpt[0]     ->Fill(min(pt,149.0), weight);
        h_numeta[0]    ->Fill(fabs(eta), weight);
        if(matchedId == 5 || matchedId == 4 ) {
          h_numhfpt[0] ->Fill(min(pt,149.0), weight);
          h_numhfeta[0]->Fill(fabs(eta), weight);
        }
        if(matchedId < 4) {
          h_numlqpt[0]->Fill(min(pt,149.0), weight);
          h_numlqeta[0]->Fill(fabs(eta), weight);
        }
        if(matchedId == 21) {
          h_numgpt[0]->Fill(min(pt,149.0), weight);
          h_numgeta[0]->Fill(fabs(eta), weight);
        }
        if(dR < 1.0)
          h_nummc3Id[0]->Fill(matchedId, weight);
        h_nummc3dR[0]->Fill(dR, weight);
	
      } //electron loop

      /* muons */
      for(int iMu = 0 ; iMu < cms2.mus_p4().size(); iMu++) {  // loop on muons
        //if(!isFakeableMuSUSY09(iMu)) continue;
        if( !isDenomMu(iMu) ) continue; // denominator muons aka fakeable objects
        // fill denominators aka fakeable objects
        Double_t pt = cms2.mus_p4()[iMu].Pt();
        Double_t eta = cms2.mus_p4()[iMu].Eta();
        Double_t phi = cms2.mus_p4()[iMu].Phi();
        h_FOptvseta[1]->Fill(fabs(eta), min(pt,149.0), weight);
        h_FOpt[1]     ->Fill(min(pt,149.0), weight);
        h_FOeta[1]    ->Fill(fabs(eta), weight);
        int cat = muFakeMCCategory(iMu);
        h_truecomposition_denom[1]->Fill( cat, weight); // true composition of muon denominator aka fakeable object

        // figure out the mcid of the closest status==3 particle
        Double_t dR = 5.0;
        unsigned int matchedId = 99999;
        for(int genIdx = 0; genIdx < cms2.genps_id().size(); genIdx++) {
          int mcid = cms2.genps_id()[genIdx];
          if( abs(mcid)==2212)
            continue;
          if( abs(mcid) == 6 || abs(mcid) == 24)
            continue;
          Double_t deta = eta - cms2.genps_p4()[genIdx].Eta();
          Double_t dphi = fabs(phi - cms2.genps_p4()[genIdx].Phi());
          if( dphi > TMath::Pi()) 
            dphi = 2*TMath::Pi() - dphi;
          if(sqrt(deta*deta + dphi*dphi) < dR) {
            dR = sqrt(deta*deta + dphi*dphi);
            matchedId = abs(mcid);
          }
        }
        if(matchedId == 5 || matchedId == 4 ) {
          h_FOhfpt[1] ->Fill(min(pt,149.0), weight);
          h_FOhfeta[1]->Fill(fabs(eta), weight);
        }
        if(matchedId < 4) {
          h_FOlqpt[1]->Fill(min(pt,149.0), weight);
          h_FOlqeta[1]->Fill(fabs(eta), weight);
        }
        if(matchedId == 21) {
          h_FOgpt[1]->Fill(min(pt,149.0), weight);
          h_FOgeta[1]->Fill(fabs(eta), weight);
        }
        if(dR < 1.0) 
          h_FOmc3Id[1]->Fill(matchedId, weight);
        h_FOmc3dR[1]->Fill(dR, weight);
	
        //if( !isNumMuSUSY09(iMu) ) continue;
        if( !isNumMu(iMu) ) continue; // numerator muons aka selection muons

        // fill numerators aka selections
        h_numptvseta[1]->Fill(fabs(eta), min(pt,149.0), weight);
        h_numpt[1]     ->Fill(min(pt,149.0), weight);
        h_numeta[1]    ->Fill(fabs(eta), weight);
        if(matchedId == 5 || matchedId == 4 ) {
          h_numhfpt[1] ->Fill(min(pt,149.0), weight);
          h_numhfeta[1]->Fill(fabs(eta), weight);
        }
        if(matchedId < 4) {
          h_numlqpt[1]->Fill(min(pt,149.0), weight);
          h_numlqeta[1]->Fill(fabs(eta), weight);
        }
        if(matchedId == 21) {
          h_numgpt[1]->Fill(min(pt,149.0), weight);
          h_numgeta[1]->Fill(fabs(eta), weight);
        }
        if(dR < 1.0) h_nummc3Id[1]->Fill(matchedId, weight);
        h_nummc3dR[1]->Fill(dR, weight);
        h_truecomposition_num[1]->Fill( cat, weight); // true composition of muon numerator aka selection muons
      } // end loop on muons

    }//event loop
  }//file loop
  
  if ( nEventsChain != nEventsTotal ) {
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  }
  

  for(unsigned int i = 0; i < 2; i++) {
    
    h_FOptvseta[i]->Sumw2();
    h_FOeta[i]->Sumw2();
    h_FOpt[i]->Sumw2();
    h_FOhfpt[i]->Sumw2();
    h_FOhfeta[i]->Sumw2();
    h_FOlqpt[i]->Sumw2();
    h_FOlqeta[i]->Sumw2();
    h_FOgpt[i]->Sumw2();
    h_FOgeta[i]->Sumw2();
    
    h_FOmc3Id[i]->Sumw2();
    
    h_numptvseta[i]->Sumw2();
    h_numeta[i]->Sumw2();
    h_numpt[i]->Sumw2();
    h_numhfpt[i]->Sumw2();
    h_numhfeta[i]->Sumw2();
    h_numlqpt[i]->Sumw2();
    h_numlqeta[i]->Sumw2();
    h_numgpt[i]->Sumw2();
    h_numgeta[i]->Sumw2();
    h_nummc3Id[i]->Sumw2();
    
    
    
    h_FRptvseta[i]->Divide(h_numptvseta[i], h_FOptvseta[i], 1.,1.,"B");
    h_FRpt[i]->Divide(h_numpt[i], h_FOpt[i], 1., 1., "B");
    h_FReta[i]->Divide(h_numeta[i], h_FOeta[i], 1., 1., "B");
    h_FRhfpt[i]->Divide(h_numhfpt[i], h_FOhfpt[i], 1., 1., "B");
    h_FRhfeta[i]->Divide(h_numhfeta[i], h_FOhfeta[i], 1., 1., "B");
    h_FRlqpt[i]->Divide(h_numlqpt[i], h_FOlqpt[i], 1., 1., "B");
    h_FRlqeta[i]->Divide(h_numlqeta[i], h_FOlqeta[i], 1., 1., "B");
    h_FRgpt[i]->Divide(h_numgpt[i], h_FOgpt[i], 1., 1., "B");
    h_FRgeta[i]->Divide(h_numgeta[i], h_FOgeta[i], 1., 1., "B");


    h_FRmc3Id[i]->Divide(h_nummc3Id[i], h_FOmc3Id[i], 1., 1., "B");
    
    
    //fill the FR errors 
    for(unsigned int ptbin = 1; ptbin < h_FRpt[i]->GetNbinsX()+1; ptbin++) {
      for(unsigned int etabin = 1; etabin < h_FReta[i]->GetNbinsX() + 1; etabin++) {
        //Float_t err2 = pow(h_FRpt[i]->GetBinError(ptbin),2) + pow(h_FReta[i]->GetBinError(etabin),2);
        Float_t err = h_FRptvseta[i]->GetBinError(etabin, ptbin);
        h_FRErrptvseta[i]->SetBinContent(etabin, ptbin, err);
        if(etabin == 1 && ptbin == 1)
          cout << err << endl;
      }//eta loop
    }//pt loop
  }//lepton flavor loop
      
  
  std::cout<<"Done with "<<prefix<<std::endl;
  rootdir = gDirectory->GetDirectory("Rint:"); 
  rootdir->cd(); 
  

  //dbarge
  h_truecomposition_ratio[0]->Divide(	h_truecomposition_num[0], h_truecomposition_denom[0], 1, 1, "B");
  h_truecomposition_ratio[1]->Divide(	h_truecomposition_num[1], h_truecomposition_denom[1], 1, 1, "B");

  return 0;
}

//------------------------------------------------------------

void QCDFRestimator::bookHistos(const char *sample) {

  char *flavor[2] = {"el", "mu"};
  //FO --> electron

  int suf = 0;
  h_truecomposition_num[suf] = new TH1F(  Form("%s_truecomposition_num_%s",sample, flavor[suf]),
                                          Form("%s_truecomposition_num_%s",sample, flavor[suf]), 5, -0.5, 4.5);
  h_truecomposition_num[suf]->Sumw2();
  h_truecomposition_denom[suf] = new TH1F(  Form("%s_truecomposition_denom_%s",sample, flavor[suf]),
                                            Form("%s_truecomposition_denom_%s",sample, flavor[suf]), 5, -0.5, 4.5);
  h_truecomposition_denom[suf]->Sumw2();
  h_truecomposition_ratio[suf] = new TH1F(  Form("%s_truecomposition_ratio_%s",sample, flavor[suf]),
                                            Form("%s_truecomposition_ratio_%s",sample, flavor[suf]), 5, -0.5, 4.5);
  h_truecomposition_ratio[suf]->Sumw2();

  suf = 1;
  h_truecomposition_num[suf] = new TH1F(  Form("%s_truecomposition_num_%s",sample, flavor[suf]),
                                          Form("%s_truecomposition_num_%s",sample, flavor[suf]), 5, -0.5, 4.5);
  h_truecomposition_num[suf]->Sumw2();
  h_truecomposition_denom[suf] = new TH1F(  Form("%s_truecomposition_denom_%s",sample, flavor[suf]),
                                            Form("%s_truecomposition_denom_%s",sample, flavor[suf]), 5, -0.5, 4.5);
  h_truecomposition_denom[suf]->Sumw2();
  h_truecomposition_ratio[suf] = new TH1F(  Form("%s_truecomposition_ratio_%s",sample, flavor[suf]),
                                            Form("%s_truecomposition_ratio_%s",sample, flavor[suf]), 5, -0.5, 4.5);
  h_truecomposition_ratio[suf]->Sumw2();

  Float_t ptel[4] = {10,20,60,150};
  Float_t etael[3] = {0, 1.479, 2.4};
  h_FOptvseta[0] = new TH2F(Form("%s_FOptvseta_%s", sample, flavor[0]),
                            Form("%s pt vs eta of FO, %s", flavor[0], sample),
                            2, etael, 3, ptel);
  h_FOpt[0] = new TH1F(Form("%s_FOpt_%s", sample, flavor[0]), 
                       Form("%s pt of FO, %s", flavor[0], sample), 
                       3, ptel);
  h_FOeta[0] = new TH1F(Form("%s_FOeta_%s", sample, flavor[0]), 
                        Form("%s eta of FO, %s", flavor[0], sample), 
                        2, etael);
  h_FOhfpt[0] = new TH1F(Form("%s_FOhfpt_%s", sample, flavor[0]), 
                         Form("%s pt of FO matched to HF, %s", flavor[0], sample), 
                         3, ptel);
  h_FOhfeta[0] = new TH1F(Form("%s_FOhfeta_%s", sample, flavor[0]), 
                          Form("%s eta of FO matched to HF, %s", flavor[0], sample), 
                          2, etael);
  h_FOlqpt[0] = new TH1F(Form("%s_FOlqpt_%s", sample, flavor[0]), 
                         Form("%s pt of FO matched to Light quarks, %s", flavor[0], sample), 
                         3, ptel);
  h_FOlqeta[0] = new TH1F(Form("%s_FOlqeta_%s", sample, flavor[0]), 
                          Form("%s eta of FO matched to Light quarks, %s", flavor[0], sample), 
                          2, etael);
  h_FOgpt[0] = new TH1F(Form("%s_FOgpt_%s", sample, flavor[0]), 
                        Form("%s pt of FO matched to gluons, %s", flavor[0], sample), 
                        3, ptel);
  h_FOgeta[0] = new TH1F(Form("%s_FOgeta_%s", sample, flavor[0]), 
                         Form("%s eta of FO matched to gluons, %s", flavor[0], sample), 
                         2, etael);


  
  h_FOmc3Id[0] = new TH1F(Form("%s_FOmc3Id_%s", sample, flavor[0]),
                          Form("Id of closest mc3 particle and %s in %s", flavor[0], sample),
                          30, 0,30);
  h_FOmc3dR[0] = new TH1F(Form("%s_FOmc3dR_%s", sample, flavor[0]),
                          Form("dR between closest mc3 particle and %s in %s", flavor[0], sample),
                          100, 0, 2.5);
   
  //tight selection info
  h_numptvseta[0] = new TH2F(Form("%s_numptvseta_%s", sample, flavor[0]), 
                             Form("%s pt vs eta of num, %s",flavor[0], sample),
                             2, etael, 3, ptel);
  h_numpt[0] = new TH1F(Form("%s_numpt_%s", sample, flavor[0]), 
                        Form("%s pt of num, %s",flavor[0], sample),
                        3, ptel);
  h_numeta[0] = new TH1F(Form("%s_numeta_%s", sample, flavor[0]), 
                         Form("%s eta of num, %s", flavor[0], sample),
                         2, etael);
  h_numhfpt[0] = new TH1F(Form("%s_numhfpt_%s", sample, flavor[0]), 
                          Form("%s pt of num matched to HF, %s", flavor[0], sample), 
                          3, ptel);
  h_numhfeta[0] = new TH1F(Form("%s_numhfeta_%s", sample, flavor[0]), 
                           Form("%s eta of num matched to HF, %s", flavor[0], sample), 
                           2, etael);
  h_numlqpt[0] = new TH1F(Form("%s_numlqpt_%s", sample, flavor[0]), 
                          Form("%s pt of num matched to Light quarks, %s", flavor[0], sample), 
                          3, ptel);
  h_numlqeta[0] = new TH1F(Form("%s_numlqeta_%s", sample, flavor[0]), 
                           Form("%s eta of num matched to Light quarks, %s", flavor[0], sample), 
                           2, etael);
  h_numgpt[0] = new TH1F(Form("%s_numgpt_%s", sample, flavor[0]), 
                         Form("%s pt of num matched to gluons, %s", flavor[0], sample), 
                         3, ptel);
  h_numgeta[0] = new TH1F(Form("%s_numgeta_%s", sample, flavor[0]), 
                          Form("%s eta of num matched to gluons, %s", flavor[0], sample), 
                          2, etael);
  h_nummc3Id[0] = new TH1F(Form("%s_nummc3Id_%s", sample, flavor[0]),
                           Form("Id of closest mc3 particle and %s in %s", flavor[0], sample),
                           30, 0,30);
  h_nummc3dR[0] = new TH1F(Form("%s_nummc3dR_%s", sample, flavor[0]),
                           Form("dR between closest mc3 particle and %s in %s", flavor[0], sample),
                           100, 0, 2.5);
    

    
  //FR --> electron
  h_FRptvseta[0] = new TH2F(Form("%s_FRptvseta_%s", sample, flavor[0]), 
                            Form("%s pt vs eta of num, %s", flavor[0], sample),
                            2, etael, 3, ptel);
  h_FRErrptvseta[0] = new TH2F(Form("%s_FRErrptvseta_%s", sample, flavor[0]), 
                               Form("%s, FRErr(pt, eta) of num, %s", flavor[0], sample),
                               2, etael, 3, ptel);
  h_FRpt[0] = new TH1F(Form("%s_FRpt_%s", sample, flavor[0]), 
                       Form("%s pt of FR, %s", flavor[0], sample),
                       3, ptel);
  h_FReta[0] = new TH1F(Form("%s_FReta_%s", sample, flavor[0]), 
                        Form("%s eta of FR, %s", flavor[0], sample),
                        2, etael);
  h_FRhfpt[0] = new TH1F(Form("%s_FRhfpt_%s", sample, flavor[0]), 
                         Form("%s pt of FR matched to HF, %s", flavor[0], sample), 
                         3, ptel);
  h_FRhfeta[0] = new TH1F(Form("%s_FRhfeta_%s", sample, flavor[0]), 
                          Form("%s eta of FR matched to HF, %s", flavor[0], sample), 
                          2, etael);
  h_FRlqpt[0] = new TH1F(Form("%s_FRlqpt_%s", sample, flavor[0]), 
                         Form("%s pt of FR matched to Light quarks, %s", flavor[0], sample), 
                         3, ptel);
  h_FRlqeta[0] = new TH1F(Form("%s_FRlqeta_%s", sample, flavor[0]), 
                          Form("%s eta of FR matched to Light quarks, %s", flavor[0], sample), 
                          2, etael);
  h_FRgpt[0] = new TH1F(Form("%s_FRgpt_%s", sample, flavor[0]), 
                        Form("%s pt of FR matched to gluons, %s", flavor[0], sample), 
                        3, ptel);
  h_FRgeta[0] = new TH1F(Form("%s_FRgeta_%s", sample, flavor[0]), 
                         Form("%s eta of FR matched to gluons, %s", flavor[0], sample), 
                         2, etael);
  h_FRmc3Id[0] = new TH1F(Form("%s_FRmc3Id_%s", sample, flavor[0]),
                          Form("%s FR as a function of mcId (status==3), %s", flavor[0], sample),
                          30, 0, 30);
    

  

  //FO --> muons
  Float_t pt[4] = {10,20,60,150};
  Float_t eta[3] = {0, 1.5, 2.4};
  h_FOptvseta[1] = new TH2F(Form("%s_FOptvseta_%s", sample, flavor[1]),
                            Form("%s pt vs eta of FO, %s", flavor[1], sample),
                            2, eta, 3, pt);
  h_FOpt[1] = new TH1F(Form("%s_FOpt_%s", sample, flavor[1]), 
                       Form("%s pt of FO, %s", flavor[1], sample), 
                       3, pt);
  h_FOeta[1] = new TH1F(Form("%s_FOeta_%s", sample, flavor[1]), 
                        Form("%s eta of FO, %s", flavor[1], sample), 
                        2, eta);
  h_FOhfpt[1] = new TH1F(Form("%s_FOhfpt_%s", sample, flavor[1]), 
                         Form("%s pt of FO matched to HF, %s", flavor[1], sample), 
                         3, pt);
  h_FOhfeta[1] = new TH1F(Form("%s_FOhfeta_%s", sample, flavor[1]), 
                          Form("%s eta of FO matched to HF, %s", flavor[1], sample), 
                          2, eta);
  h_FOlqpt[1] = new TH1F(Form("%s_FOlqpt_%s", sample, flavor[1]), 
                         Form("%s pt of FO matched to Light quarks, %s", flavor[1], sample), 
                         3, pt);
  h_FOlqeta[1] = new TH1F(Form("%s_FOlqeta_%s", sample, flavor[1]), 
                          Form("%s eta of FO matched to Light quarks, %s", flavor[1], sample), 
                          2, eta);
  h_FOgpt[1] = new TH1F(Form("%s_FOgpt_%s", sample, flavor[1]), 
                        Form("%s pt of FO matched to gluons, %s", flavor[1], sample), 
                        3, pt);
  h_FOgeta[1] = new TH1F(Form("%s_FOgeta_%s", sample, flavor[1]), 
                         Form("%s eta of FO matched to gluons, %s", flavor[1], sample), 
                         2, eta);
  h_FOmc3Id[1] = new TH1F(Form("%s_FOmc3Id_%s", sample, flavor[1]),
                          Form("Id of closest mc3 particle and %s in %s", flavor[1], sample),
                          30, 0,30);
  h_FOmc3dR[1] = new TH1F(Form("%s_FOmc3dR_%s", sample, flavor[1]),
                          Form("dR between closest mc3 particle and %s in %s", flavor[1], sample),
                          100, 0, 2.5);
  

  //tight selection info
  h_numptvseta[1] = new TH2F(Form("%s_numptvseta_%s", sample, flavor[1]), 
                             Form("%s pt vs eta of num, %s",flavor[1], sample),
                             2, eta, 3, pt);
  h_numpt[1] = new TH1F(Form("%s_numpt_%s", sample, flavor[1]), 
                        Form("%s pt of num, %s",flavor[1], sample),
                        3, pt);
  h_numeta[1] = new TH1F(Form("%s_numeta_%s", sample, flavor[1]), 
                         Form("%s eta of num, %s", flavor[1], sample),
                         2, eta);
  h_numhfpt[1] = new TH1F(Form("%s_numhfpt_%s", sample, flavor[1]), 
                          Form("%s pt of num matched to HF, %s", flavor[1], sample), 
                          3, pt);
  h_numhfeta[1] = new TH1F(Form("%s_numhfeta_%s", sample, flavor[1]), 
                           Form("%s eta of num matched to HF, %s", flavor[1], sample), 
                           2, eta);
  h_numlqpt[1] = new TH1F(Form("%s_numlqpt_%s", sample, flavor[1]), 
                          Form("%s pt of num matched to Light quarks, %s", flavor[1], sample), 
                          3, pt);
  h_numlqeta[1] = new TH1F(Form("%s_numlqeta_%s", sample, flavor[1]), 
                           Form("%s eta of num matched to Light quarks, %s", flavor[1], sample), 
                           2, eta);
  h_numgpt[1] = new TH1F(Form("%s_numgpt_%s", sample, flavor[1]), 
                         Form("%s pt of num matched to gluons, %s", flavor[1], sample), 
                         3, pt);
  h_numgeta[1] = new TH1F(Form("%s_numgeta_%s", sample, flavor[1]), 
                          Form("%s eta of num matched to gluons, %s", flavor[1], sample), 
                          2, eta);
  h_nummc3Id[1] = new TH1F(Form("%s_nummc3Id_%s", sample, flavor[1]),
                           Form("Id of closest mc3 particle and %s in %s", flavor[1], sample),
                           30, 0,30);
  h_nummc3dR[1] = new TH1F(Form("%s_nummc3dR_%s", sample, flavor[1]),
                           Form("dR between closest mc3 particle and %s in %s", flavor[1], sample),
                           100, 0, 2.5);
    

  
  //FR --> muon
  h_FRptvseta[1] = new TH2F(Form("%s_FRptvseta_%s", sample, flavor[1]), 
                            Form("%s pt vs eta of num, %s", flavor[1], sample),
                            2, eta, 3, pt);
  h_FRErrptvseta[1] = new TH2F(Form("%s_FRErrptvseta_%s", sample, flavor[1]), 
                               Form("%s, FRErr(pt, eta) of num, %s", flavor[1], sample),
                               2, eta, 3, pt);
  h_FRpt[1] = new TH1F(Form("%s_FRpt_%s", sample, flavor[1]), 
                       Form("%s pt of FR, %s", flavor[1], sample),
                       3, pt);
  h_FReta[1] = new TH1F(Form("%s_FReta_%s", sample, flavor[1]), 
                        Form("%s eta of FR, %s", flavor[1], sample),
                        2, eta);
  h_FRhfpt[1] = new TH1F(Form("%s_FRhfpt_%s", sample, flavor[1]), 
                         Form("%s pt of FR matched to HF, %s", flavor[1], sample), 
                         3, pt);
  h_FRhfeta[1] = new TH1F(Form("%s_FRhfeta_%s", sample, flavor[1]), 
                          Form("%s eta of FR matched to HF, %s", flavor[1], sample), 
                          2, eta);
  h_FRlqpt[1] = new TH1F(Form("%s_FRlqpt_%s", sample, flavor[1]), 
                         Form("%s pt of FR matched to Light quarks, %s", flavor[1], sample), 
                         3, pt);
  h_FRlqeta[1] = new TH1F(Form("%s_FRlqeta_%s", sample, flavor[1]), 
                          Form("%s eta of FR matched to Light quarks, %s", flavor[1], sample), 
                          2, eta);
  h_FRgpt[1] = new TH1F(Form("%s_FRgpt_%s", sample, flavor[1]), 
                        Form("%s pt of FR matched to gluons, %s", flavor[1], sample), 
                        3, pt);
  h_FRgeta[1] = new TH1F(Form("%s_FRgeta_%s", sample, flavor[1]), 
                         Form("%s eta of FR matched to gluons, %s", flavor[1], sample), 
                         2, eta);
  h_FRmc3Id[1] = new TH1F(Form("%s_FRmc3Id_%s", sample, flavor[1]),
                          Form("%s FR as a function of mcId (status==3), %s", flavor[1], sample),
                          30, 0, 30);
  
}


 
//------------------------------------------------------------  

bool QCDFRestimator::passFakeJetTrigger(float unCorrJetPtCut) {
  
  for(int iJet = 0 ; iJet < cms2.jets_p4().size(); iJet++) {
    if(cms2.jets_p4()[iJet].Pt()*cms2.jets_pat_noCorrF()[iJet] > unCorrJetPtCut ) 
      return true;
  }
  return false;
}
