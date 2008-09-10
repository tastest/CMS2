//now make the source file
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include "TChain.h"
#include "TChainElement.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "Math/LorentzVector.h"
#include "TMath.h"
#include "TStopwatch.h"
#include <algorithm>
#include <set>
#include "TCanvas.h"
#include "TRegexp.h"
#include "TLorentzVector.h"

using namespace std;

#ifndef __CINT__
//#include "../Tools/CMS2_soup_messedwith_Class.h"
#include "CMS2_Class.h"
CMS2 cms2;
#include "../Tools/selections.C"
#include "../Tools/matchTools.C"
#endif

static int hypos_total_n[4];
static double hypos_total_weight[4];
static unsigned int bGen;
static unsigned int bbarGen;
static unsigned int lep;
static unsigned int alep;

// static double evt_scale1fb;

enum Sample {WW, WZ, ZZ, Wjets, DYee, DYmm, DYtt, ttbar, tW}; // signal samples
enum Hypothesis {MM, EM, EE, ALL}; // hypothesis types (em and me counted as same) and all

// this is Jake's magic to sort jets by Pt
Bool_t comparePt(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > lv1, 
                 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > lv2) {
   return lv1.pt() > lv2.pt();
}

struct DorkyEventIdentifier {
     // this is a workaround for not having unique event id's in MC
     unsigned long int run, event;
     float trks_d0;
     float hyp_lt_pt, hyp_lt_eta, hyp_lt_phi;
     bool operator < (const DorkyEventIdentifier &) const;
     bool operator == (const DorkyEventIdentifier &) const;
};

bool DorkyEventIdentifier::operator < (const DorkyEventIdentifier &other) const
{
     if (run != other.run)
	  return run < other.run;
     if (event != other.event)
	  return event < other.event;
     // the floating point numbers are not easy, because we're
     // comapring ones that are truncated (because they were written
     // to file and read back in) with ones that are not truncated.
     if (fabs(trks_d0 - other.trks_d0) > 1e-6 * trks_d0)
	  return trks_d0 < other.trks_d0;
     if (fabs(hyp_lt_pt - other.hyp_lt_pt) > 1e-6 * hyp_lt_pt)
	  return hyp_lt_pt < other.hyp_lt_pt;
     if (fabs(hyp_lt_eta - other.hyp_lt_eta) > 1e-6 * hyp_lt_eta)
	  return hyp_lt_eta < other.hyp_lt_eta;
     if (fabs(hyp_lt_phi - other.hyp_lt_phi) > 1e-6 * hyp_lt_phi)
	  return hyp_lt_phi < other.hyp_lt_phi;
     // if the records are exactly the same, then r1 is not less than
     // r2.  Duh!
     return false;
}

bool DorkyEventIdentifier::operator == (const DorkyEventIdentifier &other) const
{
     if (run != other.run)
	  return false;
     if (event != other.event)
	  return false;
     // the floating point numbers are not easy, because we're
     // comapring ones that are truncated (because they were written
     // to file and read back in) with ones that are not truncated.
     if (fabs(trks_d0 - other.trks_d0) > 1e-6 * trks_d0)
	  return false;
     if (fabs(hyp_lt_pt - other.hyp_lt_pt) > 1e-6 * hyp_lt_pt)
	  return false;
     if (fabs(hyp_lt_eta - other.hyp_lt_eta) > 1e-6 * hyp_lt_eta)
	  return false;
     if (fabs(hyp_lt_phi - other.hyp_lt_phi) > 1e-6 * hyp_lt_phi)
	  return false;
     return true;
}

static std::set<DorkyEventIdentifier> already_seen;
bool is_duplicate (const DorkyEventIdentifier &id)
{
     std::pair<std::set<DorkyEventIdentifier>::const_iterator, bool> ret = 
	  already_seen.insert(id);
     return !ret.second;
}

// filter events by process
bool filterByProcess( enum Sample sample ) {
  switch (sample) {
  case WW: case WZ: case ZZ: case tW:
    return true;
  case Wjets:
    return cms2.evt_CSA07Process() < 11;
  case DYee: 
    return (cms2.evt_CSA07Process() > 10 && cms2.evt_CSA07Process() < 22 && isDYee() );
  case DYmm:
    return (cms2.evt_CSA07Process() > 10 && cms2.evt_CSA07Process() < 22 && isDYmm() );
  case DYtt:
    return (cms2.evt_CSA07Process() > 10 && cms2.evt_CSA07Process() < 22 && isDYtt() );
  case ttbar:
    return (cms2.evt_CSA07Process() > 21 && cms2.evt_CSA07Process() < 27);
  }
  return false;
}

// filter candidates by hypothesis
Hypothesis filterByHypothesis( int candidate ) {
  switch (candidate) {
  case 0:
    return MM;
  case 1: case 2:
    return EM;
  case 3:
    return EE;
  }
}

//  Book histograms...
//  Naming Convention:
//  Prefix comes from the sample and it is passed to the scanning function
//  Suffix is "ee" "em" "em" "all" which depends on the final state
//  For example: histogram named tt_hnJet_ee would be the Njet distribution
//  for the ee final state in the ttbar sample.

// MAKE SURE TO CAL SUMW2 FOR EACH 1D HISTOGRAM BEFORE FILLING!!!!!!

TH1F* hnJet[4];       // Njet distributions
TH1F* helePt[4];      // electron Pt
TH1F* hmuPt[4];       // muon Pt
TH1F* hmuPtFromSilicon[4];    // muon Pt (from tracker)
TH1F* hminLepPt[4];   // minimum lepton Pt
TH1F* hmaxLepPt[4];   // maximum lepton Pt
TH1F* helePhi[4];     // electron phi
TH1F* hmuPhi[4];      // muon phi
TH1F* hdphiLep[4];    // delta phi between leptons
TH1F* heleEta[4];     // electron eta
TH1F* hmuEta[4];      // muon eta
TH1F* hdilMass[4];    // dilepton mass
TH1F* hdilMassTightWindow[4]; // dilepton mass, but zooming around Z
TH1F* hdilPt[4];       // dilepton Pt
TH1F* hmet[4];       // MET
TH1F* hmetPhi[4];       // MET phi
TH2F* hmetVsDilepPt[4];  // MET vs dilepton Pt

TH2F* hmetOverPtVsDphi[4]; // MET/Lepton Pt vs DeltaPhi between MET and Lepton Pt
TH2F* hdphillvsmll[4]; // delta phi between leptons vs dilepton mass
TH1F* hptJet1[4];   // Pt of 1st jet
TH1F* hptJet2[4];   // Pt of 2nd jet
TH1F* hptJet3[4];   // Pt of 3rd jet
TH1F* hptJet4[4];   // Pt of 4th jet
TH1F* hetaJet1[4];   // eta of 1st jet
TH1F* hetaJet2[4];   // eta of 2nd jet
TH1F* hetaJet3[4];   // eta of 3rd jet
TH1F* hetaJet4[4];   // eta of 4th jet
TH1F* numTightLep[4]; // number of tight leptons per event.
TH1F* heleSumPt[4];   // sumPt for electron isolation
TH1F* hmuSumPt[4];   // sumPt for muon isolation
TH1F* hmuSumIso[4];  // sum of trk pt, em et, had et in cone of 0.3  
TH1F* heleRelIso[4]; //  Iso variable defined as pt/(pt+sum) for electron
TH1F* hmuRelIso[4]; //  Iso variable defined as pt/(pt+sum) for muons
TH1F* hnJetLepVeto[4]; //njet distribution after requiring numTightLep < 3.

// fkw September 2008 hists:
TH1F* hbfoundEta[4]; //eta distribution of gen b or bbar that is within dR<0.2 of found jet.
TH1F* hbnotfoundEta[4]; //eta distribution of gen b or bbar with no jet found within dR<0.2.
TH1F* hbfoundPt[4]; //pT distribution of gen b or bbar that is within dR<0.2 of found jet.
TH1F* hbnotfoundPt[4]; //pT distribution of gen b or bbar with no jet found within dR<0.2.
TH1F* hbfoundmuondR[4]; // dR between gen b/bbar and reco muon that is not candidate.
TH1F* hbnotfoundmuondR[4]; // dR between gen b/bbar and reco muon that is not candidate.
TH1F* hbfoundmuonEta[4]; // eta of closest reco muon that is not candidate and gen b/bbar.
TH1F* hbnotfoundmuonEta[4]; // eta of closest reco muon that is not candidate and gen b/bbar.
TH1F* hbfoundmuonPt[4]; // pT of closest reco muon that is not candidate and gen b/bbar.
TH1F* hbnotfoundmuonPt[4]; // pT of closest reco muon that is not candidate and gen b/bbar.
TH1F* hbfoundmuonD0sig[4]; // d0 significance of closest reco muon that is not candidate and gen b/bbar.
TH1F* hbnotfoundmuonD0sig[4]; // d0 significance of closest reco muon that is not candidate and gen b/bbar.

// fkw September 2nd round hists:
TH1F* hbmuonLowdRPt[4];
TH1F* hbmuonHidRPt[4];
TH1F* hbmuonLowdRd0sig[4];
TH1F* hbmuonHidRd0sig[4];
TH1F* hbfoundmuonancestor[4];
TH1F* hbnotfoundmuonancestor[4];

// fkw September 2008 final hist used for muon tags estimate of top bkg
TH2F* hextramuonsvsnjet[4];

// fkw September 2008 hists for jetVeto and impact parameter study
TH2F* hntrkjetvsncalojet[4];
TH2F* hnalltrkjetvsncalojet[4];
TH1F* hntrkpertrkjet[4];
TH1F* hmaxd0sigoftrksintrkjet[4];
TH1F* hsumd0sigoftrksintrkjet[4];
TH1F* hprobd0sigoftrksintrkjet[4];
TH1F* htrkjetpt[4];
TH1F* halltrkjetpt[4];

// Sanjay histogram start

TH1F* hnTrkJet[4];
TH1F* hnCaloJet[4];
TH1F* hnTotalTrkJet[4];
TH1F* hnTotalCaloJet[4];
TH1F* h_itemTrkJet[4];
TH1F* h_itemCaloJet[4];
TH1F* h_itemCaloandTrkJet[27][4];

std::vector<TLorentzVector>* trk_jets = new std::vector<TLorentzVector>();
std::vector<TLorentzVector>* calo_jets = new std::vector<TLorentzVector>();

  static const Double_t etbin[27] =
    { 0., 5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100., 105., 110., 115., 120., 125., 130.};

// Sanjay histogram end

void hypo (int i_hyp, double kFactor) 
{
     // Cut on lepton Pt
     if (cms2.hyp_lt_p4()[i_hyp].pt() < 20.0) return;
     if (cms2.hyp_ll_p4()[i_hyp].pt() < 20.0) return;
     
     // Require opposite sign
     if ( cms2.hyp_lt_id()[i_hyp] * cms2.hyp_ll_id()[i_hyp] > 0 ) return;
     
     // Z mass veto using hyp_leptons for ee and mumu final states
     if (cms2.hyp_type()[i_hyp] == 0 || cms2.hyp_type()[i_hyp] == 3) {
	  if (inZmassWindow(cms2.hyp_p4()[i_hyp].mass())) return;
     }
     
     // special handling for DY --- moved to event loop
#if 0
     bool processEvent=true;
     if (specDY == 0) {
	  if ( !isDYee() ) processEvent = false;
     } else if (specDY == 1) {
	  if ( !isDYmm() ) processEvent = false;
     } else if (specDY == 2) {
	  if ( !isDYtt() ) processEvent = false;
     }
     if (!processEvent) return;
#endif

     // fkw: There are 2 cuts that were not applied in this:
     //     if ( cms2.hyp_njets()[i_hyp] != 0) return;  //only 0-jet bin
     if (!pass4Met(i_hyp)) return; //metspecial cut, and extra tight MET if ee or mumu.
     // fkw: now we go back to the standard ttbar selection.

     // Dima's MET requirement
     if (!pass2Met(i_hyp)) return;
     
     // Muon quality cuts, including isolation
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 13 && !goodMuonIsolated(cms2.hyp_lt_index()[i_hyp]) ) return;
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 13 && !goodMuonIsolated(cms2.hyp_ll_index()[i_hyp]) ) return;
     
     // Electron quality cuts, including isolation
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11 && !goodElectronIsolated(cms2.hyp_lt_index()[i_hyp]) ) return;
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 11 && !goodElectronIsolated(cms2.hyp_ll_index()[i_hyp]) ) return;
     
     // Z veto using additional leptons in the event
     if (additionalZveto()) return;

     //count the number of tight leptons:
     // This is an FKW variable which I turned off since I do not understand it 
     // and also uses simpleIdPlus, so it should be fixed up before being turned back on.
     //    int inumTightLep = numTightLeptons();
     int inumTightLep = 0;    

     // The event weight including the kFactor (scaled to 1 fb-1)
     float weight = cms2.evt_scale1fb() * kFactor;

     // For top group political reasons, we rescale to 10 pb-1
     //  weight = weight/100.

     // If we made it to here, we passed all cuts and we are ready to fill
     int myType = 99;
     if (cms2.hyp_type()[i_hyp] == 3) myType = 0;  // ee
     if (cms2.hyp_type()[i_hyp] == 0) myType = 1;  // mm
     if (cms2.hyp_type()[i_hyp] == 1 || cms2.hyp_type()[i_hyp] == 2) myType=2; // em
     if (myType == 99) {
	  std::cout << "YUK:  unknown dilepton type = " << cms2.hyp_type()[i_hyp] << std::endl;
	  return;
     }

     //cout << " passed all cuts for WW " << endl;

     trk_jets->clear();
     calo_jets->clear();

     double jetet = 0;
     double jeteta = 3.0;

     // TrkJets & CaloJet save it after the lepton subtraction

     for ( unsigned int itrkjet=0; itrkjet<cms2.trkjets_p4().size(); ++itrkjet) {
       if ((abs(cms2.hyp_lt_id()[i_hyp]) == 11 && dRbetweenVectors(cms2.hyp_lt_p4()[i_hyp],cms2.trkjets_p4()[itrkjet]) < 0.4)||
           (abs(cms2.hyp_ll_id()[i_hyp]) == 11 && dRbetweenVectors(cms2.hyp_ll_p4()[i_hyp],cms2.trkjets_p4()[itrkjet]) < 0.4)
           ) continue;
       TLorentzVector p(cms2.trkjets_p4()[itrkjet].Px(), cms2.trkjets_p4()[itrkjet].Py(), cms2.trkjets_p4()[itrkjet].Pz(), cms2.trkjets_p4()[itrkjet].E());
       if (p.Perp() < jetet) continue;
       if (fabs(p.Eta()) > jeteta) continue;
       trk_jets->push_back(p);
     }

     vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > calo_jets_p4(cms2.jets_p4());

     for (int icalojet=0; icalojet<calo_jets_p4.size(); ++icalojet) {
       if ((abs(cms2.hyp_lt_id()[i_hyp]) == 11 && dRbetweenVectors(cms2.hyp_lt_p4()[i_hyp],calo_jets_p4[icalojet]) < 0.4)||
	   (abs(cms2.hyp_ll_id()[i_hyp]) == 11 && dRbetweenVectors(cms2.hyp_ll_p4()[i_hyp],calo_jets_p4[icalojet]) < 0.4)
	   ) continue;
       TLorentzVector p(calo_jets_p4[icalojet].Px()*cms2.jets_tq_noCorrF()[icalojet], 
			calo_jets_p4[icalojet].Py()*cms2.jets_tq_noCorrF()[icalojet], 
			calo_jets_p4[icalojet].Pz()*cms2.jets_tq_noCorrF()[icalojet],  
			calo_jets_p4[icalojet].E()*cms2.jets_tq_noCorrF()[icalojet]);//p is now uncorrected jet energy
       if (p.Perp() < jetet) continue;
       if (fabs(p.Eta()) > jeteta) continue;
       //	 if (p.Perp() < 28 ) cout << p.Perp() << endl;
       calo_jets->push_back(p);
     }

     //     get the data
     std::vector<TLorentzVector> trkjets = *trk_jets;
     std::vector<TLorentzVector> cjets = *calo_jets;

     hnTrkJet[myType]->Fill(NjetVeto(trkjets, 15),weight);
     hnCaloJet[myType]->Fill(NjetVeto(cjets, 15),weight);
     hnTrkJet[3]->Fill(NjetVeto(trkjets, 15),weight);
     hnCaloJet[3]->Fill(NjetVeto(cjets, 15),weight);

     //     cout << cjets.size() << endl;

     for(unsigned int k=0; k<26; ++k) {
       if (NjetVeto(trkjets, etbin[k]) == 0 ) {
	 h_itemTrkJet[myType]->Fill(etbin[k], weight);
	 h_itemTrkJet[3]->Fill(etbin[k], weight);
       }

       if (NjetVeto(cjets, etbin[k]) == 0 ) {
	 h_itemCaloJet[myType]->Fill(etbin[k], weight);
	 h_itemCaloJet[3]->Fill(etbin[k], weight);
       }
     }


     for(unsigned int k=0; k<26; ++k) {
       if (NjetVeto(trkjets, 10000.) == 0 ) {
	 hnTotalTrkJet[myType]->Fill(etbin[k], weight);
	 hnTotalTrkJet[3]->Fill(etbin[k], weight);
       }
       if (NjetVeto(cjets, 10000.) == 0 ) {
	 hnTotalCaloJet[myType]->Fill(etbin[k], weight);
	 hnTotalCaloJet[3]->Fill(etbin[k], weight);
       }
     }

     for(unsigned int m=0; m<26; ++m) {
       for(unsigned int k=0; k<26; ++k) {
	 if ((NjetVeto(trkjets, etbin[m]) == 0 ) && (NjetVeto(cjets, etbin[k]) == 0 )) {
	   h_itemCaloandTrkJet[m][myType]->Fill(etbin[k], weight);
	   h_itemCaloandTrkJet[m][3]->Fill(etbin[k], weight);
	 }
       }
     }

     hypos_total_n[myType]++;
     hypos_total_n[3]++;
     hypos_total_weight[myType] += weight;
     hypos_total_weight[3] += weight;

     // jet count
     hnJet[myType]->Fill(cms2.hyp_njets()[i_hyp], weight);
     hnJet[3]->Fill(cms2.hyp_njets()[i_hyp], weight);
     if( inumTightLep < 3) {
	  hnJetLepVeto[myType]->Fill(cms2.hyp_njets()[i_hyp], weight);
	  hnJetLepVeto[3]->Fill(cms2.hyp_njets()[i_hyp], weight);
     }
     numTightLep[myType]->Fill(inumTightLep,weight);
     numTightLep[3]->Fill(inumTightLep,weight);
    
     // lepton Pt
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11) helePt[myType]->Fill(cms2.hyp_lt_p4()[i_hyp].pt(), weight);
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 11) helePt[myType]->Fill(cms2.hyp_ll_p4()[i_hyp].pt(), weight);
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 13) hmuPt[myType]->Fill(cms2.hyp_lt_p4()[i_hyp].pt(), weight);
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 13) hmuPt[myType]->Fill(cms2.hyp_ll_p4()[i_hyp].pt(), weight);
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 13) hmuPtFromSilicon[myType]->Fill(cms2.mus_trk_p4().at(cms2.hyp_lt_index()[i_hyp]).pt(), weight);
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 13) hmuPtFromSilicon[myType]->Fill(cms2.mus_trk_p4().at(cms2.hyp_ll_index()[i_hyp]).pt(), weight);
     hminLepPt[myType]->Fill(min(cms2.hyp_ll_p4()[i_hyp].pt(), cms2.hyp_lt_p4()[i_hyp].pt()), weight);
     hmaxLepPt[myType]->Fill(max(cms2.hyp_ll_p4()[i_hyp].pt(), cms2.hyp_lt_p4()[i_hyp].pt()), weight );
    
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11) helePt[3]->Fill(cms2.hyp_lt_p4()[i_hyp].pt(), weight);
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 11) helePt[3]->Fill(cms2.hyp_ll_p4()[i_hyp].pt(), weight);
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 13) hmuPt[3]->Fill(cms2.hyp_lt_p4()[i_hyp].pt(), weight);
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 13) hmuPt[3]->Fill(cms2.hyp_ll_p4()[i_hyp].pt(), weight);
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 13) hmuPtFromSilicon[3]->Fill(cms2.mus_trk_p4().at(cms2.hyp_lt_index()[i_hyp]).pt(), weight);
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 13) hmuPtFromSilicon[3]->Fill(cms2.mus_trk_p4().at(cms2.hyp_ll_index()[i_hyp]).pt(), weight);
     hminLepPt[3]->Fill(min(cms2.hyp_ll_p4()[i_hyp].pt(), cms2.hyp_lt_p4()[i_hyp].pt()), weight);
     hmaxLepPt[3]->Fill(max(cms2.hyp_ll_p4()[i_hyp].pt(), cms2.hyp_lt_p4()[i_hyp].pt()), weight );
    
     // lepton Phi
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11) helePhi[myType]->Fill(cms2.hyp_lt_p4()[i_hyp].phi(), weight);
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 11) helePhi[myType]->Fill(cms2.hyp_ll_p4()[i_hyp].phi(), weight);
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 13) hmuPhi[myType]->Fill(cms2.hyp_lt_p4()[i_hyp].phi(), weight);
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 13) hmuPhi[myType]->Fill(cms2.hyp_ll_p4()[i_hyp].phi(), weight);
    
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11) helePhi[3]->Fill(cms2.hyp_lt_p4()[i_hyp].phi(), weight);
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 11) helePhi[3]->Fill(cms2.hyp_ll_p4()[i_hyp].phi(), weight);
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 13) hmuPhi[3]->Fill(cms2.hyp_lt_p4()[i_hyp].phi(), weight);
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 13) hmuPhi[3]->Fill(cms2.hyp_ll_p4()[i_hyp].phi(), weight);
    
     // dilepton mass
     hdilMass[myType]->Fill(cms2.hyp_p4()[i_hyp].mass(), weight);
     hdilMassTightWindow[myType]->Fill(cms2.hyp_p4()[i_hyp].mass(), weight);
     hdilMass[3]->Fill(cms2.hyp_p4()[i_hyp].mass(), weight);
     hdilMassTightWindow[3]->Fill(cms2.hyp_p4()[i_hyp].mass(), weight);
    
     // delta phi btw leptons
     double dphi = fabs(cms2.hyp_lt_p4()[i_hyp].phi() - cms2.hyp_ll_p4()[i_hyp].phi());
     if (dphi > TMath::Pi()) dphi = TMath::TwoPi() - dphi;
     hdphiLep[myType]->Fill(dphi, weight);
     hdphiLep[3]->Fill(dphi, weight);
    
     // dphill vs mll, i.e. the 2d correlation between the previous two variables
     hdphillvsmll[myType]->Fill(cms2.hyp_p4()[i_hyp].mass(), dphi, weight);
     hdphillvsmll[3]->Fill(cms2.hyp_p4()[i_hyp].mass(), dphi, weight);
    
     // lepton Eta
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11) heleEta[myType]->Fill(cms2.hyp_lt_p4()[i_hyp].eta(), weight);
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 11) heleEta[myType]->Fill(cms2.hyp_ll_p4()[i_hyp].eta(), weight);
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 13) hmuEta[myType]->Fill(cms2.hyp_lt_p4()[i_hyp].eta(), weight);
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 13) hmuEta[myType]->Fill(cms2.hyp_ll_p4()[i_hyp].eta(), weight);
    
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11) heleEta[3]->Fill(cms2.hyp_lt_p4()[i_hyp].eta(), weight);
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 11) heleEta[3]->Fill(cms2.hyp_ll_p4()[i_hyp].eta(), weight);
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 13) hmuEta[3]->Fill(cms2.hyp_lt_p4()[i_hyp].eta(), weight);
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 13) hmuEta[3]->Fill(cms2.hyp_ll_p4()[i_hyp].eta(), weight);


     // electron trk isolation 
     double temp_lt_iso = cms2.hyp_lt_iso()[i_hyp];  // so that min works
     double temp_ll_iso = cms2.hyp_ll_iso()[i_hyp];  // so that min works
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11) heleSumPt[myType]->Fill(min(temp_lt_iso,24.99),weight);
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11) heleSumPt[3]->Fill(min(temp_lt_iso,24.99),weight);
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 11) heleSumPt[myType]->Fill(min(temp_ll_iso,24.99),weight);
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 11) heleSumPt[3]->Fill(min(temp_ll_iso,24.99),weight);

     // muon trk isolation
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 13) hmuSumPt[myType]->Fill(min(temp_lt_iso,24.99),weight);
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 13) hmuSumPt[3]->Fill(min(temp_lt_iso,24.99),weight);
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 13) hmuSumPt[myType]->Fill(min(temp_ll_iso,24.99),weight);
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 13) hmuSumPt[3]->Fill(min(temp_ll_iso,24.99),weight);

     // muon trk+calo isolation
     double combIso_lt = -1.;
     double combIso_ll = -1.;
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 13)
	  combIso_lt = cms2.mus_iso03_sumPt().at(cms2.hyp_lt_index()[i_hyp])+cms2.mus_iso03_emEt().at(cms2.hyp_lt_index()[i_hyp])+cms2.mus_iso03_hadEt().at(cms2.hyp_lt_index()[i_hyp]);
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 13)
	  combIso_ll = cms2.mus_iso03_sumPt().at(cms2.hyp_ll_index()[i_hyp])+cms2.mus_iso03_emEt().at(cms2.hyp_ll_index()[i_hyp])+cms2.mus_iso03_hadEt().at(cms2.hyp_ll_index()[i_hyp]);
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 13) hmuSumIso[myType]->Fill(min(combIso_lt,24.99),weight);
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 13) hmuSumIso[3]->Fill(min(combIso_lt,24.99),weight);
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 13) hmuSumIso[myType]->Fill(min(combIso_ll,24.99),weight);
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 13) hmuSumIso[3]->Fill(min(combIso_ll,24.99),weight);
    

     // Relative isolation... muons
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 13) {
	  double thisSum =  cms2.mus_iso03_sumPt().at(cms2.hyp_lt_index()[i_hyp]) +  
	       cms2.mus_iso03_emEt().at(cms2.hyp_lt_index()[i_hyp])  +
	       cms2.mus_iso03_hadEt().at(cms2.hyp_lt_index()[i_hyp]);
	  double thisPt  = cms2.mus_p4().at(cms2.hyp_lt_index()[i_hyp]).pt();
	  double temp    = thisPt / (thisPt+thisSum);
	  hmuRelIso[myType]->Fill(temp, weight);
	  hmuRelIso[3]->Fill(temp, weight);
     }
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 13) {
	  double thisSum =  cms2.mus_iso03_sumPt().at(cms2.hyp_ll_index()[i_hyp]) +  
	       cms2.mus_iso03_emEt().at(cms2.hyp_ll_index()[i_hyp])  +
	       cms2.mus_iso03_hadEt().at(cms2.hyp_ll_index()[i_hyp]);
	  double thisPt  = cms2.mus_p4().at(cms2.hyp_ll_index()[i_hyp]).pt();
	  double temp    = thisPt / (thisPt+thisSum);
	  hmuRelIso[myType]->Fill(temp, weight);
	  hmuRelIso[3]->Fill(temp, weight);
     }


     // Relative isolation... electrons
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11) {
	  double thisSum =  cms2.hyp_lt_iso()[i_hyp];
	  double thisPt  = cms2.hyp_lt_p4()[i_hyp].pt();
	  double temp    = thisPt / (thisPt+thisSum);
	  heleRelIso[myType]->Fill(temp, weight);
	  heleRelIso[3]->Fill(temp, weight);
     }
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 11) {
	  double thisSum =  cms2.hyp_ll_iso()[i_hyp];
	  double thisPt  = cms2.hyp_ll_p4()[i_hyp].pt();
	  double temp    = thisPt / (thisPt+thisSum);
	  heleRelIso[myType]->Fill(temp, weight);
	  heleRelIso[3]->Fill(temp, weight);
     }

     // dilepton pt
     hdilPt[myType]->Fill(cms2.hyp_p4()[i_hyp].pt(), weight);
     hdilPt[3]->Fill(cms2.hyp_p4()[i_hyp].pt(), weight);
    
     // Met and Met phi
     hmet[myType]->Fill(cms2.hyp_met()[i_hyp], weight);      
     hmetPhi[myType]->Fill(cms2.hyp_metPhi()[i_hyp], weight);      
     hmet[3]->Fill(cms2.hyp_met()[i_hyp], weight);      
     hmetPhi[3]->Fill(cms2.hyp_metPhi()[i_hyp], weight);      
    
     // Met vs dilepton Pt
     hmetVsDilepPt[myType]->Fill(cms2.hyp_met()[i_hyp], cms2.hyp_p4()[i_hyp].pt(), weight);
     hmetVsDilepPt[3]->Fill(cms2.hyp_met()[i_hyp], cms2.hyp_p4()[i_hyp].pt(), weight);
    
     // Met over dilepton Pt vs deltaphi btw the two
     double dphi2 = fabs(cms2.hyp_p4()[i_hyp].phi() - cms2.hyp_metPhi()[i_hyp]);
     if (dphi2 > TMath::Pi()) dphi2 = TMath::TwoPi() - dphi2;
     hmetOverPtVsDphi[myType]->Fill(cms2.hyp_met()[i_hyp]/cms2.hyp_p4()[i_hyp].pt(), dphi2, weight);
     hmetOverPtVsDphi[3]->Fill(cms2.hyp_met()[i_hyp]/cms2.hyp_p4()[i_hyp].pt(), dphi2, weight);
    
     // Make a vector of sorted jets, fill jet histograms
     if (cms2.hyp_njets()[i_hyp] > 0) {
	  vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > my_hyp_jets_p4(cms2.hyp_jets_p4()[i_hyp]);
	  sort(my_hyp_jets_p4.begin(), my_hyp_jets_p4.end(), comparePt);   // sort them by Pt
	  hptJet1[myType]->Fill(my_hyp_jets_p4[0].Pt(), weight);
	  hptJet1[3]->Fill(my_hyp_jets_p4[0].Pt(), weight);
	  hetaJet1[myType]->Fill(my_hyp_jets_p4[0].Eta(), weight);
	  hetaJet1[3]->Fill(my_hyp_jets_p4[0].Eta(), weight);
	  if (cms2.hyp_njets()[i_hyp] > 1) {
	       hptJet2[myType]->Fill(my_hyp_jets_p4[1].Pt(), weight);
	       hptJet2[3]->Fill(my_hyp_jets_p4[1].Pt(), weight);
	       hetaJet2[myType]->Fill(my_hyp_jets_p4[1].Eta(), weight);
	       hetaJet2[3]->Fill(my_hyp_jets_p4[1].Eta(), weight);
	  }
	  if (cms2.hyp_njets()[i_hyp] > 2) {
	       hptJet3[myType]->Fill(my_hyp_jets_p4[2].Pt(), weight);
	       hptJet3[3]->Fill(my_hyp_jets_p4[2].Pt(), weight);
	       hetaJet3[myType]->Fill(my_hyp_jets_p4[2].Eta(), weight);
	       hetaJet3[3]->Fill(my_hyp_jets_p4[2].Eta(), weight);
	  }
	  if (cms2.hyp_njets()[i_hyp] > 3) {
	       hptJet4[myType]->Fill(my_hyp_jets_p4[3].Pt(), weight);
	       hptJet4[3]->Fill(my_hyp_jets_p4[3].Pt(), weight);
	       hetaJet4[myType]->Fill(my_hyp_jets_p4[3].Eta(), weight);
	       hetaJet4[3]->Fill(my_hyp_jets_p4[3].Eta(), weight);
	  }
     }//end of if-clause requiring at least one jet

     //fkw September 2008 histograms
     if ( bGen != 0) {//this event has at least one b-quark in genps block
       int iJet = match4vector(cms2.genps_p4()[bGen],cms2.hyp_jets_p4()[i_hyp],0.2);
       int iMuon = match4vector(cms2.genps_p4()[bGen],cms2.mus_p4());
       //Note: The algorithm here is somewhat flawed. I will miss the extra muon if the candidate muon is
       //      the closest muon to the b-quark instead of the extra muon.
       //      For an accurate count of the muon tag, I thus have a 2d hist of # of muons vs # of jets/ 
       if ( iMuon != -1 && ( 
	    myType == 0 ||
	    ( myType == 1 && cms2.hyp_lt_index()[i_hyp] != iMuon && cms2.hyp_ll_index()[i_hyp] != iMuon ) ||
	    ( myType == 2 && abs(cms2.hyp_lt_id()[i_hyp]) == 11 && cms2.hyp_ll_index()[i_hyp] != iMuon ) ||
	    ( myType == 2 && abs(cms2.hyp_ll_id()[i_hyp]) == 11 && cms2.hyp_lt_index()[i_hyp] != iMuon )
	    ) ) {//found a muon other than the candidate leptons
	 double dR = dRbetweenVectors(cms2.genps_p4()[bGen],cms2.mus_p4()[iMuon]);
	 double eta = cms2.mus_p4()[iMuon].Eta();
	 double pt = cms2.mus_p4()[iMuon].Pt();
	 double d0sig = cms2.mus_d0()[iMuon]/cms2.mus_d0Err()[iMuon];
	 //fill hists irrespective of b-quark found
	 if( dR < 0.4) { 
	   hbmuonLowdRPt[myType]->Fill(pt,weight);
	   hbmuonLowdRPt[3]->Fill(pt,weight);
	   hbmuonLowdRd0sig[myType]->Fill(d0sig,weight);
	   hbmuonLowdRd0sig[3]->Fill(d0sig,weight);
	 } else {
	   hbmuonHidRPt[myType]->Fill(pt,weight);
	   hbmuonHidRPt[3]->Fill(pt,weight);
	   hbmuonHidRd0sig[myType]->Fill(d0sig,weight);
	   hbmuonHidRd0sig[3]->Fill(d0sig,weight);
	 }
	 // now seperate into found and notfound
	 if( iJet != -1) {//found a matching jet
	   hbfoundEta[myType]->Fill(cms2.genps_p4()[bGen].Eta(),weight);
	   hbfoundPt[myType]->Fill(cms2.genps_p4()[bGen].Pt(),weight);
	   hbfoundEta[3]->Fill(cms2.genps_p4()[bGen].Eta(),weight);
	   hbfoundPt[3]->Fill(cms2.genps_p4()[bGen].Pt(),weight);
	   hbfoundmuondR[myType]->Fill(dR,weight);
	   hbfoundmuonEta[myType]->Fill(eta,weight);
	   hbfoundmuonPt[myType]->Fill(pt,weight);
	   hbfoundmuonD0sig[myType]->Fill(d0sig,weight);
	   hbfoundmuondR[3]->Fill(dR,weight);
	   hbfoundmuonEta[3]->Fill(eta,weight);
	   hbfoundmuonPt[3]->Fill(pt,weight);
	   hbfoundmuonD0sig[3]->Fill(d0sig,weight);
      	   //hbfoundmuonancestor[myType]->Fill(cms2.mus_mcidx()[iMuon],weight);
	   //hbfoundmuonancestor[3]->Fill(cms2.mus_mcidx()[iMuon],weight);
	 }else {//no matching jet found
	   hbnotfoundEta[myType]->Fill(cms2.genps_p4()[bGen].Eta(),weight);
	   hbnotfoundPt[myType]->Fill(cms2.genps_p4()[bGen].Pt(),weight);
	   hbnotfoundEta[3]->Fill(cms2.genps_p4()[bGen].Eta(),weight);
	   hbnotfoundPt[3]->Fill(cms2.genps_p4()[bGen].Pt(),weight);
	   hbnotfoundmuondR[myType]->Fill(dR,weight);
	   hbnotfoundmuonEta[myType]->Fill(eta,weight);
	   hbnotfoundmuonPt[myType]->Fill(pt,weight);
	   hbnotfoundmuonD0sig[myType]->Fill(d0sig,weight);
	   hbnotfoundmuondR[3]->Fill(dR,weight);
	   hbnotfoundmuonEta[3]->Fill(eta,weight);
	   hbnotfoundmuonPt[3]->Fill(pt,weight);
	   hbnotfoundmuonD0sig[3]->Fill(d0sig,weight);
	   //hbnotfoundmuonancestor[myType]->Fill(cms2.mus_mcidx()[iMuon],weight);
	   //hbnotfoundmuonancestor[3]->Fill(cms2.mus_mcidx()[iMuon],weight);
	 }
       } else {//no muon found
	 if( iJet != -1) {//found a matching jet
	   hbfoundEta[myType]->Fill(cms2.genps_p4()[bGen].Eta(),weight);
	   hbfoundPt[myType]->Fill(cms2.genps_p4()[bGen].Pt(),weight);
	   hbfoundEta[3]->Fill(cms2.genps_p4()[bGen].Eta(),weight);
	   hbfoundPt[3]->Fill(cms2.genps_p4()[bGen].Pt(),weight);
	 }else {//no matching jet found
	   hbnotfoundEta[myType]->Fill(cms2.genps_p4()[bGen].Eta(),weight);
	   hbnotfoundPt[myType]->Fill(cms2.genps_p4()[bGen].Pt(),weight);
	   hbnotfoundEta[3]->Fill(cms2.genps_p4()[bGen].Eta(),weight);
	   hbnotfoundPt[3]->Fill(cms2.genps_p4()[bGen].Pt(),weight);
	 }
       }//end of if muon found
     }//end of if bgen
     if ( bbarGen != 0) {//this event has at least one bbar-quark in genps block
       int iJet = match4vector(cms2.genps_p4()[bbarGen],cms2.hyp_jets_p4()[i_hyp],0.2);
       int iMuon = match4vector(cms2.genps_p4()[bbarGen],cms2.mus_p4());
       //Note: The algorithm here is somewhat flawed. I will miss the extra muon if the candidate muon is
       //      the closest muon to the b-quark instead of the extra muon.
       //      For an accurate count of the muon tag, I thus have a 2d hist of # of muons vs # of jets/ 
       if ( iMuon != -1 && ( 
	    myType == 0 ||
	    ( myType == 1 && cms2.hyp_lt_index()[i_hyp] != iMuon && cms2.hyp_ll_index()[i_hyp] != iMuon ) ||
	    ( myType == 2 && abs(cms2.hyp_lt_id()[i_hyp]) == 11 && cms2.hyp_ll_index()[i_hyp] != iMuon ) ||
	    ( myType == 2 && abs(cms2.hyp_ll_id()[i_hyp]) == 11 && cms2.hyp_lt_index()[i_hyp] != iMuon )
	    )) {//found a muon
	 double dR = dRbetweenVectors(cms2.genps_p4()[bbarGen],cms2.mus_p4()[iMuon]);
	 double eta = cms2.mus_p4()[iMuon].Eta();
	 double pt = cms2.mus_p4()[iMuon].Pt();
	 double d0sig = cms2.mus_d0()[iMuon]/cms2.mus_d0Err()[iMuon];
	 //fill hists irrespective of b-quark found
	 if( dR < 0.4) { 
	   hbmuonLowdRPt[myType]->Fill(pt,weight);
	   hbmuonLowdRPt[3]->Fill(pt,weight);
	   hbmuonLowdRd0sig[myType]->Fill(d0sig,weight);
	   hbmuonLowdRd0sig[3]->Fill(d0sig,weight);
	 } else {
	   hbmuonHidRPt[myType]->Fill(pt,weight);
	   hbmuonHidRPt[3]->Fill(pt,weight);
	   hbmuonHidRd0sig[myType]->Fill(d0sig,weight);
	   hbmuonHidRd0sig[3]->Fill(d0sig,weight);
	 }
	 // now seperate into found and notfound
	 if( iJet != -1) {//found a matching jet
	   hbfoundEta[myType]->Fill(cms2.genps_p4()[bbarGen].Eta(),weight);
	   hbfoundPt[myType]->Fill(cms2.genps_p4()[bbarGen].Pt(),weight);
	   hbfoundEta[3]->Fill(cms2.genps_p4()[bbarGen].Eta(),weight);
	   hbfoundPt[3]->Fill(cms2.genps_p4()[bbarGen].Pt(),weight);
	   hbfoundmuondR[myType]->Fill(dR,weight);
	   hbfoundmuonEta[myType]->Fill(eta,weight);
	   hbfoundmuonPt[myType]->Fill(pt,weight);
	   hbfoundmuonD0sig[myType]->Fill(d0sig,weight);
	   hbfoundmuondR[3]->Fill(dR,weight);
	   hbfoundmuonEta[3]->Fill(eta,weight);
	   hbfoundmuonPt[3]->Fill(pt,weight);
	   hbfoundmuonD0sig[3]->Fill(d0sig,weight);
	   //hbfoundmuonancestor[myType]->Fill(cms2.mus_mcidx()[iMuon],weight);
	   //hbfoundmuonancestor[3]->Fill(cms2.mus_mcidx()[iMuon],weight);
	 }else {//no matching jet found
	   hbnotfoundEta[myType]->Fill(cms2.genps_p4()[bbarGen].Eta(),weight);
	   hbnotfoundPt[myType]->Fill(cms2.genps_p4()[bbarGen].Pt(),weight);
	   hbnotfoundEta[3]->Fill(cms2.genps_p4()[bbarGen].Eta(),weight);
	   hbnotfoundPt[3]->Fill(cms2.genps_p4()[bbarGen].Pt(),weight);
	   hbnotfoundmuondR[myType]->Fill(dR,weight);
	   hbnotfoundmuonEta[myType]->Fill(eta,weight);
	   hbnotfoundmuonPt[myType]->Fill(pt,weight);
	   hbnotfoundmuonD0sig[myType]->Fill(d0sig,weight);
	   hbnotfoundmuondR[3]->Fill(dR,weight);
	   hbnotfoundmuonEta[3]->Fill(eta,weight);
	   hbnotfoundmuonPt[3]->Fill(pt,weight);
	   hbnotfoundmuonD0sig[3]->Fill(d0sig,weight);
	   //hbnotfoundmuonancestor[myType]->Fill(cms2.mus_mcidx()[iMuon],weight);
	   //hbnotfoundmuonancestor[3]->Fill(cms2.mus_mcidx()[iMuon],weight);
	 }
       } else {//no muon found
	 if( iJet != -1) {//found a matching jet
	   hbfoundEta[myType]->Fill(cms2.genps_p4()[bbarGen].Eta(),weight);
	   hbfoundPt[myType]->Fill(cms2.genps_p4()[bbarGen].Pt(),weight);
	   hbfoundEta[3]->Fill(cms2.genps_p4()[bbarGen].Eta(),weight);
	   hbfoundPt[3]->Fill(cms2.genps_p4()[bbarGen].Pt(),weight);
	 }else {//no matching jet found
	   hbnotfoundEta[myType]->Fill(cms2.genps_p4()[bbarGen].Eta(),weight);
	   hbnotfoundPt[myType]->Fill(cms2.genps_p4()[bbarGen].Pt(),weight);
	   hbnotfoundEta[3]->Fill(cms2.genps_p4()[bbarGen].Eta(),weight);
	   hbnotfoundPt[3]->Fill(cms2.genps_p4()[bbarGen].Pt(),weight);
	 }
       }//end of if muon found
     }//end of if bbargen

     //2d hist for muon tag counting
     float countmus = 0;
     for (int imu=0; imu < cms2.mus_charge().size(); ++imu) {
       if ( myType == 0 ||
	    ( myType == 1 && cms2.hyp_lt_index()[i_hyp] != imu && cms2.hyp_ll_index()[i_hyp] != imu ) ||
	    ( myType == 2 && abs(cms2.hyp_lt_id()[i_hyp]) == 11 && cms2.hyp_ll_index()[i_hyp] != imu ) ||
	    ( myType == 2 && abs(cms2.hyp_ll_id()[i_hyp]) == 11 && cms2.hyp_lt_index()[i_hyp] != imu )
	    ) ++countmus;       
     }
     hextramuonsvsnjet[myType]->Fill(countmus, cms2.hyp_njets()[i_hyp],weight);
     hextramuonsvsnjet[3]->Fill(countmus, cms2.hyp_njets()[i_hyp],weight);

     //2d hist for trkjet vs njet counting
     //I probably got this backwards which version of trkjets needs to be muon vetoed and which doesn't!
     //fkw checked this. I think the alltrkjets includes all trks 
     //    while the trkjets excludes isolated high pt muons that pass our nominal cuts.
     float ptcut = 10.0;
     int count = 0;
     float maxpt = 0.0;

     //deal with trkjets first
     for ( unsigned int itrkjet=0; itrkjet<cms2.trkjets_p4().size(); ++itrkjet){//loop over trkjets
       if ((abs(cms2.hyp_lt_id()[i_hyp]) == 11 && dRbetweenVectors(cms2.hyp_lt_p4()[i_hyp],cms2.trkjets_p4()[itrkjet]) < 0.4)||
	   (abs(cms2.hyp_ll_id()[i_hyp]) == 11 && dRbetweenVectors(cms2.hyp_ll_p4()[i_hyp],cms2.trkjets_p4()[itrkjet]) < 0.4)
	   ) continue;//veto trkjets within dR<0.4 of electrons in hypothesis.
       //muons in hypothesis are allowed to be as close as they want to a trkjet because they are excluded from trkjet anyway.
       if ( cms2.trkjets_p4()[itrkjet].pt() > ptcut ) {
	 ++count;
	 vector<unsigned int> idxIntrkJet = idxInCone( cms2.trkjets_p4()[itrkjet],cms2.trks_trk_p4());
	 hntrkpertrkjet[myType]->Fill(idxIntrkJet.size(),weight);
	 hntrkpertrkjet[3]->Fill(idxIntrkJet.size(),weight);
	 double maxd0sig = 0.0;
	 double sumd0sig = 0.0;
	 double probd0sig = 0.0;
	 for (unsigned int itrk=0; itrk<idxIntrkJet.size(); ++itrk){//loop over trks in trkjet cone
	   double d0sig = abs(cms2.trks_d0()[idxIntrkJet[itrk]]/cms2.trks_d0Err()[idxIntrkJet[itrk]]);
	   if ( d0sig > maxd0sig ) maxd0sig = d0sig;
	   sumd0sig = sumd0sig+d0sig;
	   probd0sig = probd0sig + d0sig*d0sig;
	 }
	 hmaxd0sigoftrksintrkjet[myType]->Fill(maxd0sig,weight);
	 hsumd0sigoftrksintrkjet[myType]->Fill(sumd0sig,weight);
	 hprobd0sigoftrksintrkjet[myType]->Fill(probd0sig,weight);
	 hmaxd0sigoftrksintrkjet[3]->Fill(maxd0sig,weight);
	 hsumd0sigoftrksintrkjet[3]->Fill(sumd0sig,weight);
	 hprobd0sigoftrksintrkjet[3]->Fill(probd0sig,weight);
       }//end of ptcut on trkjets
       if ( cms2.trkjets_p4()[itrkjet].pt() > maxpt ) maxpt = cms2.trkjets_p4()[itrkjet].pt();
     }//end of loop over trkjets
     htrkjetpt[myType]->Fill(maxpt,weight);
     htrkjetpt[3]->Fill(maxpt,weight);
     hntrkjetvsncalojet[myType]->Fill(count, cms2.hyp_njets()[i_hyp],weight);
     hntrkjetvsncalojet[3]->Fill(count, cms2.hyp_njets()[i_hyp],weight);

     //deal with alltrkjets second
     //don't bother with lifetime hists for alltrkjets for now.
     count = 0;
     maxpt = 0.0;
     for ( unsigned int itrkjet=0; itrkjet<cms2.alltrkjets_p4().size(); ++itrkjet){//loop over alltrkjets
       if ( dRbetweenVectors(cms2.hyp_lt_p4()[i_hyp],cms2.alltrkjets_p4()[itrkjet]) < 0.4 ||
	    dRbetweenVectors(cms2.hyp_ll_p4()[i_hyp],cms2.alltrkjets_p4()[itrkjet]) < 0.4
	    ) continue; //veto alltrkjets for both electrons and muons!
       if ( cms2.alltrkjets_p4()[itrkjet].pt() > ptcut ) ++count;
       if ( cms2.alltrkjets_p4()[itrkjet].pt() > maxpt ) maxpt = cms2.alltrkjets_p4()[itrkjet].pt();
     }//end of loop over alltrkjets
     halltrkjetpt[myType]->Fill(maxpt,weight);
     halltrkjetpt[3]->Fill(maxpt,weight);
     hnalltrkjetvsncalojet[myType]->Fill(count, cms2.hyp_njets()[i_hyp],weight);
     hnalltrkjetvsncalojet[3]->Fill(count, cms2.hyp_njets()[i_hyp],weight);


}//end of void hypo

int ScanChain( TChain* chain, enum Sample sample ) {

  unsigned int nEventsChain = chain->GetEntries();  // number of entries in chain --> number of events from all files
  unsigned int nEventsTotal = 0;

  const unsigned int numHypTypes = 4;  // number of hypotheses: MM, EM, EE, ALL

 // declare and create array of histograms
  const char sample_names[][1024] = { "ww", "wz", "zz", "wjets", "dyee", "dymm", "dytt", "ttbar", "tw" };
  const char *prefix = sample_names[sample];

  // DY samples are supposed to get an additional k-factor of 1.2
  double kFactor = 1;
  switch (sample) {
  case DYee: case DYmm: case DYtt:
       kFactor = 1.12;
       break;
  case ttbar:
       kFactor = 1.85;
       break;
  default:
       break;
  }

//   switch (sample) {
//   case WW:
//        evt_scale1fb = 0.1538;
//        break;
//   default:
//        break;
//   }
  char *suffix[3];
  suffix[0] = "ee";
  suffix[1] = "mm";
  suffix[2] = "em";
  suffix[3] = "all";
  
  // The statement below should work but does not work due to bug in root when TH2 are also used
  // Rene Brun promised a fix.
  //TH1::SetDefaultSumw2(kTRUE); // do errors properly based on weights
  
  for (int i=0; i<4; i++) {

    hnTrkJet[i] = new TH1F(Form("%s_hnTrkJet_%s",prefix,suffix[i]),Form("%s_hnTrkJet_%s",prefix,suffix[i]),
			   10,0.,10.);	
    hnCaloJet[i] = new TH1F(Form("%s_hnCaloJet_%s",prefix,suffix[i]),Form("%s_hnCaloJet_%s",prefix,suffix[i]),
			    10,0.,10.);	
    
    hnTotalTrkJet[i] = new TH1F(Form("%s_hnTotalTrkJet_%s",prefix,suffix[i]),Form("%s_hnTotalTrkJet_%s",prefix,suffix[i]), 26,0.,130.);
    hnTotalCaloJet[i] = new TH1F(Form("%s_hnTotalCaloJet_%s",prefix,suffix[i]),Form("%s_hnTotalCaloJet_%s",prefix,suffix[i]), 26,0.,130.);

    h_itemTrkJet[i] = new TH1F(Form("%s_h_itemTrkJet_%s",prefix,suffix[i]),Form("%s_h_itemTrkJet_%s",prefix,suffix[i]), 26,0.,130.);
    h_itemCaloJet[i] = new TH1F(Form("%s_h_itemCaloJet_%s",prefix,suffix[i]),Form("%s_h_itemCaloJet_%s",prefix,suffix[i]), 26,0.,130.);

    h_itemCaloandTrkJet[0][i] = new TH1F(Form("%s_h_itemCaloandTrkJet_0_%s",prefix,suffix[i]),Form("%s_h_itemCaloandTrkJet_0_%s",prefix,suffix[i]), 26,0.,130.);
    h_itemCaloandTrkJet[1][i] = new TH1F(Form("%s_h_itemCaloandTrkJet_1_%s",prefix,suffix[i]),Form("%s_h_itemCaloandTrkJet_1_%s",prefix,suffix[i]), 26,0.,130.);
    h_itemCaloandTrkJet[2][i] = new TH1F(Form("%s_h_itemCaloandTrkJet_2_%s",prefix,suffix[i]),Form("%s_h_itemCaloandTrkJet_2_%s",prefix,suffix[i]), 26,0.,130.);
    h_itemCaloandTrkJet[3][i] = new TH1F(Form("%s_h_itemCaloandTrkJet_3_%s",prefix,suffix[i]),Form("%s_h_itemCaloandTrkJet_3_%s",prefix,suffix[i]), 26,0.,130.);
    h_itemCaloandTrkJet[4][i] = new TH1F(Form("%s_h_itemCaloandTrkJet_4_%s",prefix,suffix[i]),Form("%s_h_itemCaloandTrkJet_4_%s",prefix,suffix[i]), 26,0.,130.);
    h_itemCaloandTrkJet[5][i] = new TH1F(Form("%s_h_itemCaloandTrkJet_5_%s",prefix,suffix[i]),Form("%s_h_itemCaloandTrkJet_5_%s",prefix,suffix[i]), 26,0.,130.);
    h_itemCaloandTrkJet[6][i] = new TH1F(Form("%s_h_itemCaloandTrkJet_6_%s",prefix,suffix[i]),Form("%s_h_itemCaloandTrkJet_6_%s",prefix,suffix[i]), 26,0.,130.);
    h_itemCaloandTrkJet[7][i] = new TH1F(Form("%s_h_itemCaloandTrkJet_7_%s",prefix,suffix[i]),Form("%s_h_itemCaloandTrkJet_7_%s",prefix,suffix[i]), 26,0.,130.);
    h_itemCaloandTrkJet[8][i] = new TH1F(Form("%s_h_itemCaloandTrkJet_8_%s",prefix,suffix[i]),Form("%s_h_itemCaloandTrkJet_8_%s",prefix,suffix[i]), 26,0.,130.);
    h_itemCaloandTrkJet[9][i] = new TH1F(Form("%s_h_itemCaloandTrkJet_9_%s",prefix,suffix[i]),Form("%s_h_itemCaloandTrkJet_9_%s",prefix,suffix[i]), 26,0.,130.);
    h_itemCaloandTrkJet[10][i] = new TH1F(Form("%s_h_itemCaloandTrkJet_10_%s",prefix,suffix[i]),Form("%s_h_itemCaloandTrkJet_10_%s",prefix,suffix[i]), 26,0.,130.);
    h_itemCaloandTrkJet[11][i] = new TH1F(Form("%s_h_itemCaloandTrkJet_11_%s",prefix,suffix[i]),Form("%s_h_itemCaloandTrkJet_11_%s",prefix,suffix[i]), 26,0.,130.);
    h_itemCaloandTrkJet[12][i] = new TH1F(Form("%s_h_itemCaloandTrkJet_12_%s",prefix,suffix[i]),Form("%s_h_itemCaloandTrkJet_12_%s",prefix,suffix[i]), 26,0.,130.);
    h_itemCaloandTrkJet[13][i] = new TH1F(Form("%s_h_itemCaloandTrkJet_13_%s",prefix,suffix[i]),Form("%s_h_itemCaloandTrkJet_13_%s",prefix,suffix[i]), 26,0.,130.);
    h_itemCaloandTrkJet[14][i] = new TH1F(Form("%s_h_itemCaloandTrkJet_14_%s",prefix,suffix[i]),Form("%s_h_itemCaloandTrkJet_14_%s",prefix,suffix[i]), 26,0.,130.);
    h_itemCaloandTrkJet[15][i] = new TH1F(Form("%s_h_itemCaloandTrkJet_15_%s",prefix,suffix[i]),Form("%s_h_itemCaloandTrkJet_15_%s",prefix,suffix[i]), 26,0.,130.);
    h_itemCaloandTrkJet[16][i] = new TH1F(Form("%s_h_itemCaloandTrkJet_16_%s",prefix,suffix[i]),Form("%s_h_itemCaloandTrkJet_16_%s",prefix,suffix[i]), 26,0.,130.);
    h_itemCaloandTrkJet[17][i] = new TH1F(Form("%s_h_itemCaloandTrkJet_17_%s",prefix,suffix[i]),Form("%s_h_itemCaloandTrkJet_17_%s",prefix,suffix[i]), 26,0.,130.);
    h_itemCaloandTrkJet[18][i] = new TH1F(Form("%s_h_itemCaloandTrkJet_18_%s",prefix,suffix[i]),Form("%s_h_itemCaloandTrkJet_18_%s",prefix,suffix[i]), 26,0.,130.);
    h_itemCaloandTrkJet[19][i] = new TH1F(Form("%s_h_itemCaloandTrkJet_19_%s",prefix,suffix[i]),Form("%s_h_itemCaloandTrkJet_19_%s",prefix,suffix[i]), 26,0.,130.);
    h_itemCaloandTrkJet[20][i] = new TH1F(Form("%s_h_itemCaloandTrkJet_20_%s",prefix,suffix[i]),Form("%s_h_itemCaloandTrkJet_20_%s",prefix,suffix[i]), 26,0.,130.);
    h_itemCaloandTrkJet[21][i] = new TH1F(Form("%s_h_itemCaloandTrkJet_21_%s",prefix,suffix[i]),Form("%s_h_itemCaloandTrkJet_21_%s",prefix,suffix[i]), 26,0.,130.);
    h_itemCaloandTrkJet[22][i] = new TH1F(Form("%s_h_itemCaloandTrkJet_22_%s",prefix,suffix[i]),Form("%s_h_itemCaloandTrkJet_22_%s",prefix,suffix[i]), 26,0.,130.);
    h_itemCaloandTrkJet[23][i] = new TH1F(Form("%s_h_itemCaloandTrkJet_23_%s",prefix,suffix[i]),Form("%s_h_itemCaloandTrkJet_23_%s",prefix,suffix[i]), 26,0.,130.);
    h_itemCaloandTrkJet[24][i] = new TH1F(Form("%s_h_itemCaloandTrkJet_24_%s",prefix,suffix[i]),Form("%s_h_itemCaloandTrkJet_24_%s",prefix,suffix[i]), 26,0.,130.);
    h_itemCaloandTrkJet[25][i] = new TH1F(Form("%s_h_itemCaloandTrkJet_25_%s",prefix,suffix[i]),Form("%s_h_itemCaloandTrkJet_25_%s",prefix,suffix[i]), 26,0.,130.);
    h_itemCaloandTrkJet[26][i] = new TH1F(Form("%s_h_itemCaloandTrkJet_26_%s",prefix,suffix[i]),Form("%s_h_itemCaloandTrkJet_26_%s",prefix,suffix[i]), 26,0.,130.);
    h_itemCaloandTrkJet[27][i] = new TH1F(Form("%s_h_itemCaloandTrkJet_27_%s",prefix,suffix[i]),Form("%s_h_itemCaloandTrkJet_27_%s",prefix,suffix[i]), 26,0.,130.);

    hnJet[i] = new TH1F(Form("%s_hnJet_%s",prefix,suffix[i]),Form("%s_nJet_%s",prefix,suffix[i]),
			5,0.,5.);	
    hnJetLepVeto[i] = new TH1F(Form("%s_hnJetLepVeto_%s",prefix,suffix[i]),Form("%s_nJetLepVeto_%s",prefix,suffix[i]),
			       5,0.,5.);	
    numTightLep[i] = new TH1F(Form("%s_numTightLep_%s",prefix,suffix[i]),Form("%s_numTightLep_%s",prefix,suffix[i]),
			      10,0.,10.);	
    helePt[i] = new TH1F(Form("%s_helePt_%s",prefix,suffix[i]),Form("%s_elePt_%s",prefix,suffix[i]),
			       150,0.,150.);
    hmuPt[i]  = new TH1F(Form("%s_hmuPt_%s",prefix,suffix[i]),Form("%s_muPt_%s",prefix,suffix[i]),
			 150,0.,150.);
    hmuPtFromSilicon[i]  = new TH1F(Form("%s_hmuPtFromSilicon_%s",prefix,suffix[i]),
                                    Form("%s_muPtFromSilicon_%s",prefix,suffix[i]),150,0.,150.);
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
    hdphillvsmll[i] = new TH2F(Form("%s_dphillvsmll_%s",prefix,suffix[i]),
			       Form("%s_dphillvsmll_%s",prefix,suffix[i]),
			       100,10.,210.,50,0., TMath::Pi());
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
    
    heleSumPt[i] = new TH1F(Form("%s_heleSumPt_%s",prefix,suffix[i]),Form("%s_heleSumPt_%s",prefix,suffix[i]),
			    100, 0., 25.);
    hmuSumPt[i] = new TH1F(Form("%s_hmuSumPt_%s",prefix,suffix[i]),Form("%s_hmuSumPt_%s",prefix,suffix[i]),
			    100, 0., 25.);
    hmuSumIso[i] = new TH1F(Form("%s_hmuIsoSum_%s",prefix,suffix[i]),Form("%s_hmuIsoSum_%s",prefix,suffix[i]),
			    100, 0., 25.);
    heleRelIso[i] = new TH1F(Form("%s_heleRelIso_%s",prefix,suffix[i]),Form("%s_heleRelIso_%s",prefix,suffix[i]),
			     100, 0., 1.);
    hmuRelIso[i] = new TH1F(Form("%s_hmuRelIso_%s",prefix,suffix[i]),Form("%s_hmuRelIso_%s",prefix,suffix[i]),
			     100, 0., 1.);

    //fkw new hists September 2008
    hbfoundEta[i] = new TH1F(Form("%s_hbfoundEta_%s",prefix,suffix[i]),Form("%s_hbfoundEta_%s",prefix,suffix[i]),100, -5., 5.);			 
    hbfoundPt[i] = new TH1F(Form("%s_hbfoundPt_%s",prefix,suffix[i]),Form("%s_hbfoundPt_%s",prefix,suffix[i]),100, 0., 100.);			 
    hbfoundmuondR[i] = new TH1F(Form("%s_hbfoundmuondR_%s",prefix,suffix[i]),Form("%s_hbfoundmuondR_%s",prefix,suffix[i]),100, 0., 2.);			 
    hbfoundmuonPt[i] = new TH1F(Form("%s_hbfoundmuonPt_%s",prefix,suffix[i]),Form("%s_hbfoundmuonPt_%s",prefix,suffix[i]),100, 0., 100.);			 
    hbfoundmuonEta[i] = new TH1F(Form("%s_hbfoundmuonEta_%s",prefix,suffix[i]),Form("%s_hbfoundmuonEta_%s",prefix,suffix[i]),100, -5., 5.);			 
    hbfoundmuonD0sig[i] = new TH1F(Form("%s_hbfoundmuonD0sig_%s",prefix,suffix[i]),Form("%s_hbfoundmuonD0sig_%s",prefix,suffix[i]),100, -5., 5.);			 

    hbnotfoundEta[i] = new TH1F(Form("%s_hbnotfoundEta_%s",prefix,suffix[i]),Form("%s_hbnotfoundEta_%s",prefix,suffix[i]),100, -5., 5.);			 
    hbnotfoundPt[i] = new TH1F(Form("%s_hbnotfoundPt_%s",prefix,suffix[i]),Form("%s_hbnotfoundPt_%s",prefix,suffix[i]),100, 0., 100.);			 
    hbnotfoundmuondR[i] = new TH1F(Form("%s_hbnotfoundmuondR_%s",prefix,suffix[i]),Form("%s_hbnotfoundmuondR_%s",prefix,suffix[i]),100, 0., 2.);			 
    hbnotfoundmuonPt[i] = new TH1F(Form("%s_hbnotfoundmuonPt_%s",prefix,suffix[i]),Form("%s_hbnotfoundmuonPt_%s",prefix,suffix[i]),100, 0., 100.);			 
    hbnotfoundmuonEta[i] = new TH1F(Form("%s_hbnotfoundmuonEta_%s",prefix,suffix[i]),Form("%s_hbnotfoundmuonEta_%s",prefix,suffix[i]),100, -5., 5.);			 
    hbnotfoundmuonD0sig[i] = new TH1F(Form("%s_hbnotfoundmuonD0sig_%s",prefix,suffix[i]),Form("%s_hbnotfoundmuonD0sig_%s",prefix,suffix[i]),100, -5., 5.);			 

    // fkw September 2008 2nd round hists:
    hbmuonLowdRPt[i] = new TH1F(Form("%s_hbmuonLowdRPt_%s",prefix,suffix[i]),Form("%s_hbmuonLowdRPt_%s",prefix,suffix[i]),100, 0., 100.);		 
    hbmuonHidRPt[i] = new TH1F(Form("%s_hbmuonHidRPt_%s",prefix,suffix[i]),Form("%s_hbmuonHidRPt_%s",prefix,suffix[i]),100, 0., 100.);	 
    hbmuonLowdRd0sig[i] = new TH1F(Form("%s_hbmuonLowdRd0sig_%s",prefix,suffix[i]),Form("%s_hbmuonLowdRd0sig_%s",prefix,suffix[i]),100, -5., 5.);			 
    hbmuonHidRd0sig[i] = new TH1F(Form("%s_hbmuonHidRd0sig_%s",prefix,suffix[i]),Form("%s_hbmuonHidRd0sig_%s",prefix,suffix[i]),100, -5., 5.);			 
    hbfoundmuonancestor[i] = new TH1F(Form("%s_hbfoundmuonancestor_%s",prefix,suffix[i]),Form("%s_hbfoundmuonancestor_%s",prefix,suffix[i]),100, 0, 100); 
    hbnotfoundmuonancestor[i] = new TH1F(Form("%s_hbnotfoundmuonancestor_%s",prefix,suffix[i]),Form("%s_hbnotfoundmuonancestor_%s",prefix,suffix[i]),100, 0, 100); 

    // fkw September 2008 final hist used for muon tag estimate of top bkg
    hextramuonsvsnjet[i] = new TH2F(Form("%s_extramuonsvsnjet_%s",
					 prefix,suffix[i]),
			       Form("%s_extramuonsvsnjet_%s",prefix,suffix[i]),
			       10,0.0,10.0,10,0.0,10.0);

    // fkw September 2008 hists for jetveto and impact parameter study
    hntrkjetvsncalojet[i] = new TH2F(Form("%s_trkjetvsncalojet_%s",
					 prefix,suffix[i]),
			       Form("%s_trkjetvsncalojet_%s",prefix,suffix[i]),
			       10,0.0,10.0,10,0.0,10.0);
    hnalltrkjetvsncalojet[i] = new TH2F(Form("%s_alltrkjetvsncalojet_%s",
					 prefix,suffix[i]),
			       Form("%s_alltrkjetvsncalojet_%s",prefix,suffix[i]),
			       10,0.0,10.0,10,0.0,10.0);
    hntrkpertrkjet[i] = new TH1F(Form("%s_ntrkpertrkjet_%s",prefix,suffix[i]),
				 Form("%s_ntrkpertrkjet_%s",prefix,suffix[i]),
				 50,0.0,50.0);
    hmaxd0sigoftrksintrkjet[i] = new TH1F(Form("%s_maxd0sigoftrksintrkjet_%s",prefix,suffix[i]),
				 Form("%s_maxd0sigoftrksintrkjet_%s",prefix,suffix[i]),
					  100,0.0,10.0);
    hsumd0sigoftrksintrkjet[i] = new TH1F(Form("%s_sumd0sigoftrksintrkjet_%s",prefix,suffix[i]),
				 Form("%s_sumd0sigoftrksintrkjet_%s",prefix,suffix[i]),
					  100,0.0,50.0);
    hprobd0sigoftrksintrkjet[i] = new TH1F(Form("%s_probd0sigoftrksintrkjet_%s",prefix,suffix[i]),
				 Form("%s_probd0sigoftrksintrkjet_%s",prefix,suffix[i]),
					   100,0.0,100.0);
    htrkjetpt[i] = new TH1F(Form("%s_trkjetpt_%s",prefix,suffix[i]),
				 Form("%s_trkjetpt_%s",prefix,suffix[i]),
					   100,0.0,100.0);
    halltrkjetpt[i] = new TH1F(Form("%s_alltrkjetpt_%s",prefix,suffix[i]),
				 Form("%s_alltrkjetpt_%s",prefix,suffix[i]),
					   100,0.0,100.0);


    hnJet[i]->Sumw2();
    hnJetLepVeto[i]->Sumw2();
    numTightLep[i]->Sumw2();
    helePt[i]->Sumw2();
    hmuPt[i]->Sumw2();
    hmuPtFromSilicon[i]->Sumw2();
    hminLepPt[i]->Sumw2();
    hmaxLepPt[i]->Sumw2();
    helePhi[i]->Sumw2();
    hmuPhi[i]->Sumw2();
    hdphiLep[i]->Sumw2();
    heleEta[i]->Sumw2();
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
    heleSumPt[i]->Sumw2();
    hmuSumPt[i]->Sumw2();
    hmuSumIso[i]->Sumw2();
    heleRelIso[i]->Sumw2(); 
    hmuRelIso[i]->Sumw2(); 
    hbfoundEta[i]->Sumw2();
    hbfoundPt[i]->Sumw2();
    hbfoundmuondR[i]->Sumw2();
    hbfoundmuonEta[i]->Sumw2();
    hbfoundmuonPt[i]->Sumw2();
    hbfoundmuonD0sig[i]->Sumw2();
    hbnotfoundEta[i]->Sumw2();
    hbnotfoundPt[i]->Sumw2();
    hbnotfoundmuondR[i]->Sumw2();
    hbnotfoundmuonEta[i]->Sumw2();
    hbnotfoundmuonPt[i]->Sumw2();
    hbnotfoundmuonD0sig[i]->Sumw2();
    hbmuonLowdRPt[i]->Sumw2();
    hbmuonHidRPt[i]->Sumw2();
    hbmuonLowdRd0sig[i]->Sumw2();
    hbmuonHidRd0sig[i]->Sumw2();
    hbfoundmuonancestor[i]->Sumw2();
    hbnotfoundmuonancestor[i]->Sumw2();
    hextramuonsvsnjet[i]->Sumw2();
    hntrkjetvsncalojet[i]->Sumw2();
    hnalltrkjetvsncalojet[i]->Sumw2();
    hntrkpertrkjet[i]->Sumw2();
    hmaxd0sigoftrksintrkjet[i]->Sumw2();
    hsumd0sigoftrksintrkjet[i]->Sumw2();
    hprobd0sigoftrksintrkjet[i]->Sumw2();
    htrkjetpt[i]->Sumw2();
    halltrkjetpt[i]->Sumw2();
  }

  memset(hypos_total_n, 0, sizeof(hypos_total_n));
  memset(hypos_total_weight, 0, sizeof(hypos_total_weight));

  // clear list of duplicates
  already_seen.clear();
  int duplicates_total_n = 0;
  double duplicates_total_weight = 0;

  int i_permille_old = 0;
  // file loop
  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  while (TChainElement *currentFile = (TChainElement*)fileIter.Next()) {
       // need to call TFile::Open(), since the file is not
       // necessarily a plain TFile (TNetFile, TDcacheFile, etc)
//        printf("current file: %s (%s), %s\n", currentFile->GetName(), 
// 	      currentFile->GetTitle(), currentFile->IsA()->GetName());
       TFile *f = TFile::Open(currentFile->GetTitle()); 
       TTree *tree = (TTree*)f->Get("Events");
       
       cms2.Init(tree);  // set branch addresses for TTree tree

       TStopwatch t;
       //Event Loop
       unsigned int nEvents = tree->GetEntries();
       for( unsigned int event = 0; event < nEvents; ++event) {
	    cms2.GetEntry(event);  // get entries for Event number event from branches of TTree tree
	    ++nEventsTotal;
	    if (cms2.trks_d0().size() == 0)
		 continue;
	    DorkyEventIdentifier id = { cms2.evt_run(), cms2.evt_event(), cms2.trks_d0()[0], 
					cms2.hyp_lt_p4()[0].pt(), cms2.hyp_lt_p4()[0].eta(), cms2.hyp_lt_p4()[0].phi() };
	    if (is_duplicate(id)) {
		 duplicates_total_n++;
		 duplicates_total_weight += cms2.evt_scale1fb();
		 continue;
	    }

	    int i_permille = (int)floor(1000 * nEventsTotal / float(nEventsChain));
	    if (i_permille != i_permille_old) {
		 // xterm magic from L. Vacavant and A. Cerri
		 printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
			"\033[0m\033[32m <---\033[0m\015", i_permille/10.);
		 fflush(stdout);
		 i_permille_old = i_permille;
	    }
	    
	    // filter by process
	    if ( !filterByProcess(sample) ) continue;
	    
	    // fkw, here go all histos that should be filled per event instead of per hyp:
	    // loop over generator particles:
	    //Note: top = +-6, W = +-24, b = +-5
	    //cout << " Event = " << event << endl;
	    bGen = 0;
	    bbarGen = 0;
	    lep = 0;
	    alep = 0;
	    unsigned int nGen = cms2.genps_id().size();
	    for( unsigned int iGen = 0; iGen < nGen; ++iGen){
	      if ( cms2.genps_id()[iGen] == 5 && cms2.genps_id_mother()[iGen] == 6 ) bGen = iGen;
	      if ( cms2.genps_id()[iGen] == -5 && cms2.genps_id_mother()[iGen] == -6 ) bbarGen = iGen;
	      if ( cms2.genps_id()[iGen] == 11 || cms2.genps_id()[iGen] == 13 ) lep = iGen;
	      if ( cms2.genps_id()[iGen] == -11 || cms2.genps_id()[iGen] == -13 ) alep = iGen;
	    } 
	    //cout << " bgen= " << bGen << " bbarGen= " << bbarGen << " lep= " << lep << " alep= " << alep << endl;
	    // fkw, end of per event filling of histos.
	    
	    // loop over hypothesis candidates
	    unsigned int nHyps = cms2.hyp_type().size();
	    for( unsigned int i_hyp = 0; i_hyp < nHyps; ++i_hyp ) {
	      hypo(i_hyp, kFactor);
	    }
       }
       t.Stop();
       printf("Real time: %u events / %f s = %e event/s\n", nEvents, 
	      t.RealTime(), nEvents / t.RealTime());
       printf("CPU time: %u events / %f s = %e event/s\n", nEvents, 
	      t.CpuTime(), nEvents / t.CpuTime());
       delete f;
  }
  if ( nEventsChain != nEventsTotal ) {
       printf("ERROR: number of events from files (%d) is not equal to total number"
	      " of events (%d)\n", nEventsChain, nEventsTotal);
  }

  printf("Total candidate count (ee em mm all): %d %d %d %d.  Total weight %f %f %f %f\n",   
	 hypos_total_n[0], hypos_total_n[1], hypos_total_n[2], hypos_total_n[3], 
	 hypos_total_weight[0], hypos_total_weight[1], hypos_total_weight[2], hypos_total_weight[3]);
  printf("Total duplicate count: %d.  Total weight %f\n",   
	 duplicates_total_n, duplicates_total_weight);
  
  return 0;
}
