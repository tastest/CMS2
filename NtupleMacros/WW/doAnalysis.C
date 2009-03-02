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
#include "TDatabasePDG.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooCategory.h"

using namespace std;

#ifndef __CINT__
#include "../CORE/CMS2.cc"
#include "../CORE/utilities.cc"
#include "../CORE/selections.cc"
#endif

TH1F* hypos_total;
TH1F* hypos_total_weighted;

enum Sample {WW, WZ, ZZ, Wjets, DYee, DYmm, DYtt, ttbar, tW}; // signal samples
enum Hypothesis {MM, EM, EE, ALL}; // hypothesis types (em and me counted as same) and all

// this is Jake's magic to sort jets by Pt
Bool_t comparePt(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > lv1, 
                 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > lv2) {
   return lv1.pt() > lv2.pt();
}

struct DorkyEventIdentifier {
     // this is a workaround for not having unique event id's in MC
     unsigned long int run, event, lumi;
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
  case DYee: 
    return isDYee();
  case DYmm:
    return isDYmm();
  case DYtt:
    return isDYtt();
  default:
    return true;
  }
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
  cout << "Unknown type: " << candidate << "Abort" << endl;
  assert(0);
  return MM;
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
TH1F* heleSumPt[4];   // sumPt for electron isolation
TH1F* hmuSumPt[4];   // sumPt for muon isolation
TH1F* hmuSumIso[4];  // sum of trk pt, em et, had et in cone of 0.3  
TH1F* heleRelIso[4]; //  Iso variable defined as pt/(pt+sum) for electron
TH1F* hmuRelIso[4]; //  Iso variable defined as pt/(pt+sum) for muons

// Histograms with only basic cuts
TH1F* helTrkIsoPassId;      // electron trk isolation passed el ID
TH1F* helTrkPatIsoPassId;   // electron trk PAT isolation passed  el ID
TH1F* helTrkIsoNoId;        // electron trk isolation no el ID
TH1F* helTrkPatIsoNoId;     // electron trk PAT isolation no el ID
TH1F* helTrkIsoFailId;      // electron trk isolation failed el ID
TH1F* helTrkPatIsoFailId;   // electron trk PAT failed el ID

TH1F* helEcalJuraIsoPassId;  // electron ECAL Jurassic isolation based on basic clusters with el ID
TH1F* helEcalPatIsoPassId;   // electron ECAL PAT Jurassic isolation based on rec hits with el ID
TH1F* helEcalJuraIsoNoId;    // electron ECAL Jurassic isolation based on basic clusters no el ID
TH1F* helEcalPatIsoNoId;     // electron ECAL PAT Jurassic isolation based on rec hits no el ID
TH1F* helEcalJuraIsoFailId;    // electron ECAL Jurassic isolation based on basic clusters failed el ID
TH1F* helEcalPatIsoFailId;     // electron ECAL PAT Jurassic isolation based on rec hits failed el ID

TH1F* helHcalConeIsoPassId;  // electron HCAL tower isolation based on basic clusters with el ID
TH1F* helHcalPatIsoPassId;   // electron HCAL PAT isolation based on rec hits with el ID
TH1F* helHcalConeIsoNoId;    // electron HCAL tower isolation based on basic clusters no el ID
TH1F* helHcalPatIsoNoId;     // electron HCAL PAT isolation based on rec hits no el ID
TH1F* helHcalConeIsoFailId;    // electron HCAL tower isolation based on basic clusters failed el ID
TH1F* helHcalPatIsoFailId;     // electron HCAL PAT isolation based on rec hits failed el ID

TH1F* helRelIsoPassId;  // electron relative iso (trk+ecal+hcal) CMS2 with weights 1,1,1 passed el ID
TH1F* helRelIsoNoId;    // electron relative iso (trk+ecal+hcal) CMS2 with weights 1,1,1 no el ID
TH1F* helRelIsoFailId;  // electron relative iso (trk+ecal+hcal) CMS2 with weights 1,1,1 failed el ID

TH1F* helRelPatIsoPassId;  // electron relative iso (trk+ecal+hcal) PAT with weights 1,1,1 passed el ID
TH1F* helRelPatIsoNoId;    // electron relative iso (trk+ecal+hcal) PAT with weights 1,1,1 no el ID
TH1F* helRelPatIsoFailId;  // electron relative iso (trk+ecal+hcal) PAT with weights 1,1,1 failed el ID

TH1F* hemElRelIso;  // electron relative iso for emu final selection
TH1F* hemMuRelIso;  // muon relative iso for emu final selection

TH1F* hmaxJPTEt;         // energy distribution for the most energetic jet
TH1F* hmaxCaloJetEt;     // energy distribution for the most energetic jet
TH1F* hmaxTrkJetEt;      // energy distribution for the most energetic jet
TH1F* hmaxCaloTrkJetEt;  // energy distribution for the most energetic jet
TH1F* hmaxCaloTrkJet2Et; // energy distribution for the most energetic jet
TH1F* hmaxGenJetEt;      // energy distribution for the most energetic jet


// fkw September 2008 final hist used for muon tags estimate of top bkg
TH2F* hextramuonsvsnjet[4];

struct hypo_monitor{
  std::vector<std::pair<std::string,unsigned int> > counters;
  void count(unsigned int index, const char* name){
    unsigned int current_size = counters.size();
    for ( unsigned int i=current_size; i<=index; ++i ) 
      counters.push_back( std::pair<std::string,unsigned int>("",0) );
    counters[index].first = name;
    counters[index].second++;
  }
  void print(){
    for ( unsigned int i=0; i<counters.size(); ++i ) 
      std::cout << counters[i].first << "\t" << counters[i].second << std::endl;
  }
};
    
hypo_monitor monitor;

void checkIsolation(int i_hyp, double weight){
  // LT
  if ( abs(cms2.hyp_lt_id()[i_hyp]) == 11 ) { 
    helTrkIsoNoId->Fill(      cms2.els_tkIso()[cms2.hyp_lt_index()[i_hyp]], weight );
    helTrkPatIsoNoId->Fill(   cms2.els_pat_trackIso()[cms2.hyp_lt_index()[i_hyp]], weight );
    helEcalJuraIsoNoId->Fill( cms2.els_ecalJuraIso()[cms2.hyp_lt_index()[i_hyp]], weight );
    helEcalPatIsoNoId->Fill(  cms2.els_pat_ecalIso()[cms2.hyp_lt_index()[i_hyp]], weight );
    helHcalConeIsoNoId->Fill( cms2.els_hcalConeIso()[cms2.hyp_lt_index()[i_hyp]], weight );
    helHcalPatIsoNoId->Fill(  cms2.els_pat_hcalIso()[cms2.hyp_lt_index()[i_hyp]], weight );
    double sum = cms2.els_tkIso()[cms2.hyp_lt_index()[i_hyp]] +
      cms2.els_ecalJuraIso()[cms2.hyp_lt_index()[i_hyp]] +
      cms2.els_hcalConeIso()[cms2.hyp_lt_index()[i_hyp]];
    helRelIsoNoId->Fill(cms2.hyp_lt_p4()[i_hyp].pt()/(cms2.hyp_lt_p4()[i_hyp].pt()+sum) , weight);
    double sum2 = cms2.els_pat_trackIso()[cms2.hyp_lt_index()[i_hyp]] +
      cms2.els_pat_ecalIso()[cms2.hyp_lt_index()[i_hyp]] +
      cms2.els_pat_hcalIso()[cms2.hyp_lt_index()[i_hyp]];
    helRelPatIsoNoId->Fill(cms2.hyp_lt_p4()[i_hyp].pt()/(cms2.hyp_lt_p4()[i_hyp].pt()+sum2) , weight);
    if ( cms2.els_robustId()[cms2.hyp_lt_index()[i_hyp]] ) {
      helTrkIsoPassId->Fill(      cms2.els_tkIso()[cms2.hyp_lt_index()[i_hyp]], weight );
      helTrkPatIsoPassId->Fill(   cms2.els_pat_trackIso()[cms2.hyp_lt_index()[i_hyp]], weight );
      helEcalJuraIsoPassId->Fill( cms2.els_ecalJuraIso()[cms2.hyp_lt_index()[i_hyp]], weight );
      helEcalPatIsoPassId->Fill(  cms2.els_pat_ecalIso()[cms2.hyp_lt_index()[i_hyp]], weight );
      helHcalConeIsoPassId->Fill( cms2.els_hcalConeIso()[cms2.hyp_lt_index()[i_hyp]], weight );
      helHcalPatIsoPassId->Fill(  cms2.els_pat_hcalIso()[cms2.hyp_lt_index()[i_hyp]], weight );
      double sum = cms2.els_tkIso()[cms2.hyp_lt_index()[i_hyp]] +
	cms2.els_ecalJuraIso()[cms2.hyp_lt_index()[i_hyp]] +
	cms2.els_hcalConeIso()[cms2.hyp_lt_index()[i_hyp]];
      helRelIsoPassId->Fill(cms2.hyp_lt_p4()[i_hyp].pt()/(cms2.hyp_lt_p4()[i_hyp].pt()+sum) , weight);
      double sum2 = cms2.els_pat_trackIso()[cms2.hyp_lt_index()[i_hyp]] +
	cms2.els_pat_ecalIso()[cms2.hyp_lt_index()[i_hyp]] +
	cms2.els_pat_hcalIso()[cms2.hyp_lt_index()[i_hyp]];
      helRelPatIsoPassId->Fill(cms2.hyp_lt_p4()[i_hyp].pt()/(cms2.hyp_lt_p4()[i_hyp].pt()+sum2) , weight);
    } else {
      helTrkIsoFailId->Fill(      cms2.els_tkIso()[cms2.hyp_lt_index()[i_hyp]], weight );
      helTrkPatIsoFailId->Fill(   cms2.els_pat_trackIso()[cms2.hyp_lt_index()[i_hyp]], weight );
      helEcalJuraIsoFailId->Fill( cms2.els_ecalJuraIso()[cms2.hyp_lt_index()[i_hyp]], weight );
      helEcalPatIsoFailId->Fill(  cms2.els_pat_ecalIso()[cms2.hyp_lt_index()[i_hyp]], weight );
      helHcalConeIsoFailId->Fill( cms2.els_hcalConeIso()[cms2.hyp_lt_index()[i_hyp]], weight );
      helHcalPatIsoFailId->Fill(  cms2.els_pat_hcalIso()[cms2.hyp_lt_index()[i_hyp]], weight );
      double sum = cms2.els_tkIso()[cms2.hyp_lt_index()[i_hyp]] +
	cms2.els_ecalJuraIso()[cms2.hyp_lt_index()[i_hyp]] +
	cms2.els_hcalConeIso()[cms2.hyp_lt_index()[i_hyp]];
      helRelIsoFailId->Fill(cms2.hyp_lt_p4()[i_hyp].pt()/(cms2.hyp_lt_p4()[i_hyp].pt()+sum) , weight);
      double sum2 = cms2.els_pat_trackIso()[cms2.hyp_lt_index()[i_hyp]] +
	cms2.els_pat_ecalIso()[cms2.hyp_lt_index()[i_hyp]] +
	cms2.els_pat_hcalIso()[cms2.hyp_lt_index()[i_hyp]];
      helRelPatIsoFailId->Fill(cms2.hyp_lt_p4()[i_hyp].pt()/(cms2.hyp_lt_p4()[i_hyp].pt()+sum2) , weight);
    }
  }
  // LL
  if ( abs(cms2.hyp_ll_id()[i_hyp]) == 11 ) { 
    helTrkIsoNoId->Fill(      cms2.els_tkIso()[cms2.hyp_ll_index()[i_hyp]], weight );
    helTrkPatIsoNoId->Fill(   cms2.els_pat_trackIso()[cms2.hyp_ll_index()[i_hyp]], weight );
    helEcalJuraIsoNoId->Fill( cms2.els_ecalJuraIso()[cms2.hyp_ll_index()[i_hyp]], weight );
    helEcalPatIsoNoId->Fill(  cms2.els_pat_ecalIso()[cms2.hyp_ll_index()[i_hyp]], weight );
    helHcalConeIsoNoId->Fill( cms2.els_hcalConeIso()[cms2.hyp_ll_index()[i_hyp]], weight );
    helHcalPatIsoNoId->Fill(  cms2.els_pat_hcalIso()[cms2.hyp_ll_index()[i_hyp]], weight );
    double sum = cms2.els_tkIso()[cms2.hyp_ll_index()[i_hyp]] +
      cms2.els_ecalJuraIso()[cms2.hyp_ll_index()[i_hyp]] +
      cms2.els_hcalConeIso()[cms2.hyp_ll_index()[i_hyp]];
    helRelIsoNoId->Fill(cms2.hyp_ll_p4()[i_hyp].pt()/(cms2.hyp_ll_p4()[i_hyp].pt()+sum) , weight);
    double sum2 = cms2.els_pat_trackIso()[cms2.hyp_ll_index()[i_hyp]] +
      cms2.els_pat_ecalIso()[cms2.hyp_ll_index()[i_hyp]] +
      cms2.els_pat_hcalIso()[cms2.hyp_ll_index()[i_hyp]];
    helRelPatIsoNoId->Fill(cms2.hyp_ll_p4()[i_hyp].pt()/(cms2.hyp_ll_p4()[i_hyp].pt()+sum2) , weight);
    if ( cms2.els_robustId()[cms2.hyp_ll_index()[i_hyp]] ) {
      helTrkIsoPassId->Fill(      cms2.els_tkIso()[cms2.hyp_ll_index()[i_hyp]], weight );
      helTrkPatIsoPassId->Fill(   cms2.els_pat_trackIso()[cms2.hyp_ll_index()[i_hyp]], weight );
      helEcalJuraIsoPassId->Fill( cms2.els_ecalJuraIso()[cms2.hyp_ll_index()[i_hyp]], weight );
      helEcalPatIsoPassId->Fill(  cms2.els_pat_ecalIso()[cms2.hyp_ll_index()[i_hyp]], weight );
      helHcalConeIsoPassId->Fill( cms2.els_hcalConeIso()[cms2.hyp_ll_index()[i_hyp]], weight );
      helHcalPatIsoPassId->Fill(  cms2.els_pat_hcalIso()[cms2.hyp_ll_index()[i_hyp]], weight );
      double sum = cms2.els_tkIso()[cms2.hyp_ll_index()[i_hyp]] +
	cms2.els_ecalJuraIso()[cms2.hyp_ll_index()[i_hyp]] +
	cms2.els_hcalConeIso()[cms2.hyp_ll_index()[i_hyp]];
      helRelIsoPassId->Fill(cms2.hyp_ll_p4()[i_hyp].pt()/(cms2.hyp_ll_p4()[i_hyp].pt()+sum) , weight);
      double sum2 = cms2.els_pat_trackIso()[cms2.hyp_ll_index()[i_hyp]] +
	cms2.els_pat_ecalIso()[cms2.hyp_ll_index()[i_hyp]] +
	cms2.els_pat_hcalIso()[cms2.hyp_ll_index()[i_hyp]];
      helRelPatIsoPassId->Fill(cms2.hyp_ll_p4()[i_hyp].pt()/(cms2.hyp_ll_p4()[i_hyp].pt()+sum2) , weight);
    } else {
      helTrkIsoFailId->Fill(      cms2.els_tkIso()[cms2.hyp_ll_index()[i_hyp]], weight );
      helTrkPatIsoFailId->Fill(   cms2.els_pat_trackIso()[cms2.hyp_ll_index()[i_hyp]], weight );
      helEcalJuraIsoFailId->Fill( cms2.els_ecalJuraIso()[cms2.hyp_ll_index()[i_hyp]], weight );
      helEcalPatIsoFailId->Fill(  cms2.els_pat_ecalIso()[cms2.hyp_ll_index()[i_hyp]], weight );
      helHcalConeIsoFailId->Fill( cms2.els_hcalConeIso()[cms2.hyp_ll_index()[i_hyp]], weight );
      helHcalPatIsoFailId->Fill(  cms2.els_pat_hcalIso()[cms2.hyp_ll_index()[i_hyp]], weight );
      double sum = cms2.els_tkIso()[cms2.hyp_ll_index()[i_hyp]] +
	cms2.els_ecalJuraIso()[cms2.hyp_ll_index()[i_hyp]] +
	cms2.els_hcalConeIso()[cms2.hyp_ll_index()[i_hyp]];
      helRelIsoFailId->Fill(cms2.hyp_ll_p4()[i_hyp].pt()/(cms2.hyp_ll_p4()[i_hyp].pt()+sum) , weight);
      double sum2 = cms2.els_pat_trackIso()[cms2.hyp_ll_index()[i_hyp]] +
	cms2.els_pat_ecalIso()[cms2.hyp_ll_index()[i_hyp]] +
	cms2.els_pat_hcalIso()[cms2.hyp_ll_index()[i_hyp]];
      helRelPatIsoFailId->Fill(cms2.hyp_ll_p4()[i_hyp].pt()/(cms2.hyp_ll_p4()[i_hyp].pt()+sum2) , weight);
    }
  }
}

void getIsolationSidebandsAfterSelections(int i_hyp, double weight, RooDataSet* dataset, bool passedAllLeptonRequirements){
  // em case
  // both leptons are relaxed 
  RooArgSet set( *(dataset->get()) );
  set.setCatIndex("selected",passedAllLeptonRequirements?1:0);
  set.setRealValue("event",cms2.evt_event());
  set.setRealValue("run",cms2.evt_run());
  set.setRealValue("lumi",cms2.evt_lumiBlock());
  set.setCatLabel("sample_type","data_relaxed_iso");
  
  if ( cms2.hyp_type()[i_hyp] == 1 || cms2.hyp_type()[i_hyp] == 2 ){
    unsigned int imu = cms2.hyp_lt_index()[i_hyp];
    unsigned int iel = cms2.hyp_ll_index()[i_hyp];
    if ( TMath::Abs(cms2.hyp_lt_id()[i_hyp])==11 ){
      imu = cms2.hyp_ll_index()[i_hyp];
      iel = cms2.hyp_lt_index()[i_hyp];
    }
    if ( goodElectronWithoutIsolation(iel) && goodMuonIsolated(imu) ) {
      hemElRelIso->Fill( el_rel_iso(iel,true), weight );
      set.setCatLabel("hyp_type","em");
      set.setCatLabel("fake_type","electron");
      set.setRealValue("iso",el_rel_iso(iel,true));
      dataset->add(set,weight);
    
    }
    if ( goodElectronIsolated(iel,true) && goodMuonWithoutIsolation(imu) ) {
      hemMuRelIso->Fill( mu_rel_iso(imu), weight );
      set.setCatLabel("hyp_type","em");
      set.setCatLabel("fake_type","muon");
      set.setRealValue("iso",mu_rel_iso(imu));
      dataset->add(set,weight);
    }
  }
  // mumu case - only one muon is relaxed
  // taking the first pair that works, since looking for the least
  // isolated lepton would increase the non-isolated contribution
  // from signal sources by double counting the tail
  if ( cms2.hyp_type()[i_hyp] == 0){
    unsigned int imu1 = cms2.hyp_lt_index()[i_hyp];
    unsigned int imu2 = cms2.hyp_ll_index()[i_hyp];
    bool relaxed_muon_found(false);
    if ( goodMuonWithoutIsolation(imu1) && goodMuonIsolated(imu2) ) {
      set.setCatLabel("hyp_type","mm");
      set.setCatLabel("fake_type","muon");
      set.setRealValue("iso",mu_rel_iso(imu1));
      dataset->add(set,weight);
      relaxed_muon_found = true;
    }
    if ( !relaxed_muon_found &&
	 goodMuonIsolated(imu1) && goodMuonWithoutIsolation(imu2) ) {
      set.setCatLabel("hyp_type","mm");
      set.setCatLabel("fake_type","muon");
      set.setRealValue("iso",mu_rel_iso(imu2));
      dataset->add(set,weight);
    }
  }
  // ee case - only one muon is relaxed
  // taking the first pair that works, since looking for the least
  // isolated lepton would increase the non-isolated contribution
  // from signal sources by double counting the tail
  if ( cms2.hyp_type()[i_hyp] == 3){
    unsigned int iel1 = cms2.hyp_lt_index()[i_hyp];
    unsigned int iel2 = cms2.hyp_ll_index()[i_hyp];
    bool relaxed_electron_found(false);
    if ( goodElectronWithoutIsolation(iel1) && goodElectronIsolated(iel2,true) ) {
      set.setCatLabel("hyp_type","ee");
      set.setCatLabel("fake_type","electron");
      set.setRealValue("iso",el_rel_iso(iel1,true));
      dataset->add(set,weight);
      relaxed_electron_found = true;
    }
    if ( !relaxed_electron_found &&
	 goodElectronIsolated(iel1) && goodElectronWithoutIsolation(iel2) ) {
      set.setCatLabel("hyp_type","ee");
      set.setCatLabel("fake_type","electron");
      set.setRealValue("iso",el_rel_iso(iel2,true));
      dataset->add(set,weight);
    }
  }
}

void find_most_energetic_jets(int i_hyp, double weight)
{
  {
    double jptMax(0.);
    for ( unsigned int i=0; i < cms2.jpts_p4().size(); ++i) {
      if ( cms2.jpts_p4()[i].Et() < jptMax ) continue;
      if ( TMath::Abs(cms2.jpts_p4()[i].eta()) > 3.0 ) continue;
      if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.jpts_p4()[i])) < 0.4 ||
	   TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.jpts_p4()[i])) < 0.4 ) continue;
      jptMax = cms2.jpts_p4()[i].Et();
    }
    hmaxJPTEt->Fill(jptMax, weight);
  }
  {
    double genJetMax(0.);
    for ( unsigned int i=0; i < cms2.genjets_p4().size(); ++i) {
      if ( cms2.genjets_p4()[i].Et() < genJetMax ) continue;
      if ( TMath::Abs(cms2.genjets_p4()[i].eta()) > 3.0 ) continue;
      if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.genjets_p4()[i])) < 0.4 ||
	   TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.genjets_p4()[i])) < 0.4 ) continue;
      genJetMax = cms2.genjets_p4()[i].Et();
    }
    hmaxGenJetEt->Fill(genJetMax, weight);
  }
  {
    double caloJetMax(0.);
    for ( unsigned int i=0; i < cms2.jets_pat_jet_p4().size(); ++i) {
      if ( cms2.jets_pat_jet_p4()[i].Et()*cms2.jets_pat_noCorrF()[i] < caloJetMax ) continue;
      if ( TMath::Abs(cms2.jets_pat_jet_p4()[i].eta()) > 3.0 ) continue;
      if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.jets_pat_jet_p4()[i])) < 0.4 ||
	   TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.jets_pat_jet_p4()[i])) < 0.4 ) continue;
      caloJetMax = cms2.jets_pat_jet_p4()[i].Et()*cms2.jets_pat_noCorrF()[i];
    }
    hmaxCaloJetEt->Fill(caloJetMax, weight);
    double trkJetMax(0.);
    for ( unsigned int i=0; i < cms2.trkjets_p4().size(); ++i) {
	 if ( cms2.trkjets_p4()[i].Et() < trkJetMax ) continue;
	 if ( TMath::Abs(cms2.trkjets_p4()[i].eta()) > 3.0 ) continue;
	 if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.trkjets_p4()[i])) < 0.4 ||
	      TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.trkjets_p4()[i])) < 0.4 ) continue;
	 trkJetMax = cms2.trkjets_p4()[i].Et();
    }
    hmaxTrkJetEt->Fill(trkJetMax, weight);
    hmaxCaloTrkJetEt->Fill((caloJetMax+trkJetMax)/2, weight);
    hmaxCaloTrkJet2Et->Fill(std::max(caloJetMax,trkJetMax), weight);
  }
}  

void hypo (int i_hyp, double kFactor, RooDataSet* dataset = 0) 
{
     int myType = 99;
     if (cms2.hyp_type()[i_hyp] == 3) myType = 0;  // ee
     if (cms2.hyp_type()[i_hyp] == 0) myType = 1;  // mm
     if (cms2.hyp_type()[i_hyp] == 1 || cms2.hyp_type()[i_hyp] == 2) myType=2; // em
     if (myType == 99) {
	  std::cout << "YUK:  unknown dilepton type = " << cms2.hyp_type()[i_hyp] << std::endl;
	  return;
     }
     // The event weight including the kFactor (scaled to 1 fb-1)
     float weight = cms2.evt_scale1fb() * kFactor;

     unsigned int icounter(0);
     monitor.count(icounter++,"Total number of hypothesis: ");
     
     // Cut on lepton Pt
     if (cms2.hyp_lt_p4()[i_hyp].pt() < 20.0) return;
     if (cms2.hyp_ll_p4()[i_hyp].pt() < 20.0) return;
     monitor.count(icounter++,"Total number of hypothesis after lepton pt cut: ");
     
     // Require opposite sign
     if ( cms2.hyp_lt_id()[i_hyp] * cms2.hyp_ll_id()[i_hyp] > 0 ) return;
     
     // check electron isolation and id (no selection at this point)
     checkIsolation(i_hyp, weight);

     // Z mass veto using hyp_leptons for ee and mumu final states
     if (cms2.hyp_type()[i_hyp] == 0 || cms2.hyp_type()[i_hyp] == 3) {
       if (inZmassWindow(cms2.hyp_p4()[i_hyp].mass())) return;
     }

     // Z veto using additional leptons in the event
     // if (additionalZveto()) return;
     monitor.count(icounter++,"Total number of hypothesis after lepton pt + z vetos: ");
     
     // track corrected MET
     // const TVector3 trkCorr = correctMETforTracks();
     const TVector3 trkCorr; // no tcMET correction
     
     if (!pass2Met(i_hyp, trkCorr)) return;
     if (!pass4Met(i_hyp, trkCorr)) return;
     monitor.count(icounter++,"Total number of hypothesis after lepton pt + z vetos + MET cuts: ");
     
     bool goodEvent = true;

     unsigned int nJPT = nJPTs(i_hyp);
     if (nJPT>0) goodEvent = false;
     int countmus = numberOfExtraMuons(i_hyp,true);
     // int countmus = numberOfExtraMuons(i_hyp);
     int nExtraVetoMuons = numberOfExtraMuons(i_hyp,false);;
     if (nExtraVetoMuons) goodEvent = false;

     bool passedAllLeptonRequirements = true;
     // Muon quality cuts, including isolation
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 13 && !goodMuonIsolated(cms2.hyp_lt_index()[i_hyp]) ) passedAllLeptonRequirements = false;
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 13 && !goodMuonIsolated(cms2.hyp_ll_index()[i_hyp]) ) passedAllLeptonRequirements = false;
     
     // Electron quality cuts, including isolation
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11 && !goodElectronIsolated(cms2.hyp_lt_index()[i_hyp],true) ) passedAllLeptonRequirements = false;
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 11 && !goodElectronIsolated(cms2.hyp_ll_index()[i_hyp],true) ) passedAllLeptonRequirements = false;
     
     if (goodEvent && dataset )
       getIsolationSidebandsAfterSelections(i_hyp, weight, dataset, passedAllLeptonRequirements);
     
     if ( !passedAllLeptonRequirements ) return;
     monitor.count(icounter++,"Total number of hypothesis after full lepton selection + z vetos + MET cuts: ");
     
     // trkjet veto
     // if ( !passTrkJetVeto(i_hyp) ) return;
     
     // find most energetic jets
     find_most_energetic_jets(i_hyp,weight);
     // 2D hist for muon tag counting
     hextramuonsvsnjet[myType]->Fill(countmus, nJPT, weight);
     hextramuonsvsnjet[3]->Fill(countmus, nJPT, weight);
     
     if ( ! goodEvent ) return;

     // -------------------------------------------------------------------//
     // If we made it to here, we passed all cuts and we are ready to fill //
     // -------------------------------------------------------------------//
     
     hypos_total->Fill(myType);
     hypos_total->Fill(3);
     hypos_total_weighted->Fill(myType,weight);
     hypos_total_weighted->Fill(3,weight);

     // jet count
     hnJet[myType]->Fill(cms2.hyp_njets()[i_hyp], weight);
     hnJet[3]->Fill(cms2.hyp_njets()[i_hyp], weight);
     
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


}//end of void hypo

RooDataSet* MakeNewDataset(const char* name)
{
  RooRealVar set_iso("iso","iso",0.,1.);
  RooRealVar set_event("event","event",0);
  RooRealVar set_run("run","run",0);
  RooRealVar set_lumi("lumi","lumi",0);
  RooRealVar set_weight("weight","weight",0);
  RooCategory set_selected("selected","Passed final WW selection requirements");
  set_selected.defineType("true",1);
  set_selected.defineType("false",0);

  RooCategory set_hyp_type("hyp_type","Hypothesis type");
  set_hyp_type.defineType("ee",0);
  set_hyp_type.defineType("mm",1);
  set_hyp_type.defineType("em",2);
  
  RooCategory set_fake_type("fake_type","Define type of lepton for which isolation is extracted");
  set_fake_type.defineType("electron",0);
  set_fake_type.defineType("muon",1);

  RooCategory set_sample_type("sample_type","Sample type");
  set_sample_type.defineType("data_relaxed_iso",0);  // full sample with final selection 
  set_sample_type.defineType("control_sample_signal_iso",1);

  RooDataSet* dataset = new RooDataSet(name, name,
				       RooArgSet(set_event,set_run,set_lumi,
						 set_iso,set_selected,set_weight,
						 set_hyp_type,set_fake_type,set_sample_type) );
  dataset->setWeightVar(set_weight);
  return dataset;
}

void AddIsoSignalControlSample( int i_hyp, double kFactor, RooDataSet* dataset = 0 ) {
  if ( !dataset ) return;
  // The event weight including the kFactor (scaled to 1 fb-1)
  float weight = cms2.evt_scale1fb() * kFactor;
  // Cut on lepton Pt
  if (cms2.hyp_lt_p4()[i_hyp].pt() < 20.0) return;
  if (cms2.hyp_ll_p4()[i_hyp].pt() < 20.0) return;
  // Require opposite sign
  if ( cms2.hyp_lt_id()[i_hyp] * cms2.hyp_ll_id()[i_hyp] > 0 ) return;
  // Z mass veto using hyp_leptons for ee and mumu final states
  if ( cms2.hyp_type()[i_hyp] == 1 || cms2.hyp_type()[i_hyp] == 2 ) return;
  if (! inZmassWindow(cms2.hyp_p4()[i_hyp].mass())) return;
  RooArgSet set( *(dataset->get()) );
  set.setCatIndex("selected",0);
  set.setRealValue("event",cms2.evt_event());
  set.setRealValue("run",cms2.evt_run());
  set.setRealValue("lumi",cms2.evt_lumiBlock());
  set.setCatLabel("sample_type","control_sample_signal_iso");
	
  if ( cms2.hyp_type()[i_hyp] == 3 ){
    set.setCatLabel("hyp_type","ee");
    set.setCatLabel("fake_type","electron");
    if ( goodElectronIsolated(cms2.hyp_lt_index()[i_hyp],true) &&
	 goodElectronWithoutIsolation(cms2.hyp_ll_index()[i_hyp]) ){
      set.setRealValue("iso",el_rel_iso(cms2.hyp_ll_index()[i_hyp],true));
      dataset->add(set,weight);
    }
    if ( goodElectronIsolated(cms2.hyp_ll_index()[i_hyp],true) &&
	 goodElectronWithoutIsolation(cms2.hyp_lt_index()[i_hyp]) ){
      set.setRealValue("iso",el_rel_iso(cms2.hyp_lt_index()[i_hyp],true));
      dataset->add(set,weight);
    }
  } else {
    set.setCatLabel("hyp_type","mm");
    set.setCatLabel("fake_type","muon");
    if ( goodMuonIsolated(cms2.hyp_lt_index()[i_hyp]) &&
	 goodMuonWithoutIsolation(cms2.hyp_ll_index()[i_hyp]) ){
      set.setRealValue("iso",mu_rel_iso(cms2.hyp_ll_index()[i_hyp]));
      dataset->add(set,weight);
    }
    if ( goodMuonIsolated(cms2.hyp_ll_index()[i_hyp]) &&
	 goodMuonWithoutIsolation(cms2.hyp_lt_index()[i_hyp]) ){
      set.setRealValue("iso",mu_rel_iso(cms2.hyp_lt_index()[i_hyp]));
      dataset->add(set,weight);
    }
  }
}

RooDataSet* ScanChain( TChain* chain, enum Sample sample ) {
  
  unsigned int nEventsChain = chain->GetEntries();  // number of entries in chain --> number of events from all files
  unsigned int nEventsTotal = 0;

  // const unsigned int numHypTypes = 4;  // number of hypotheses: MM, EM, EE, ALL

 // declare and create array of histograms
  const char sample_names[][1024] = { "ww", "wz", "zz", "wjets", "dyee", "dymm", "dytt", "ttbar", "tw" };
  const char *prefix = sample_names[sample];
  RooDataSet* dataset = MakeNewDataset(sample_names[sample]);
  // DY samples are supposed to get an additional k-factor of 1.2
  double kFactor = 1;
  /* 
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
  */

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
  
  hypos_total          = new TH1F(Form("%s_hypos_total",prefix),"Total number of hypothesis counts",4,0,4);
  hypos_total_weighted = new TH1F(Form("%s_hypos_total_weighted",prefix),"Total number of hypotheses (weighted)",4,0,4);
  hypos_total_weighted->Sumw2();

  for (unsigned int i=0; i<4; ++i){
    hypos_total->GetXaxis()->SetBinLabel(i+1,suffix[i]);
    hypos_total_weighted->GetXaxis()->SetBinLabel(i+1,suffix[i]);
  }
  
  // The statement below should work but does not work due to bug in root when TH2 are also used
  // Rene Brun promised a fix.
  //TH1::SetDefaultSumw2(kTRUE); // do errors properly based on weights
  
  for (int i=0; i<4; i++) {

    hnJet[i] = new TH1F(Form("%s_hnJet_%s",prefix,suffix[i]),Form("%s_nJet_%s",prefix,suffix[i]),
			5,0.,5.);	
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

    // fkw September 2008 final hist used for muon tag estimate of top bkg
    hextramuonsvsnjet[i] = new TH2F(Form("%s_extramuonsvsnjet_%s",
					 prefix,suffix[i]),
			       Form("%s_extramuonsvsnjet_%s",prefix,suffix[i]),
			       10,0.0,10.0,10,0.0,10.0);


    hnJet[i]->Sumw2();
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
    hextramuonsvsnjet[i]->Sumw2();

  }
  
  helTrkIsoPassId = new TH1F(Form("%s_helTrkIsoPassId",prefix),
				  Form("%s - electron trk isolation passed robust el id",prefix),
				  100, 0., 20.);
  helTrkIsoPassId->Sumw2();
  helTrkIsoFailId = new TH1F(Form("%s_helTrkIsoFailId",prefix),
				  Form("%s - electron trk isolation failed robust el id",prefix),
				  100, 0., 20.);
  helTrkIsoFailId->Sumw2();
  helTrkIsoNoId   = new TH1F(Form("%s_helTrkIsoNoId",prefix),
				  Form("%s - electron trk isolation without el id",prefix),
				  100, 0., 20.);
  helTrkIsoNoId->Sumw2();
  helTrkPatIsoPassId  = new TH1F(Form("%s_helTrkPatIsoPassId",prefix),
				  Form("%s - electron trk PAT isolation passed robust el id",prefix),
				  100, 0., 20.);
  helTrkPatIsoPassId->Sumw2();
  helTrkPatIsoFailId  = new TH1F(Form("%s_helTrkPatIsoFailId",prefix),
				  Form("%s - electron trk PAT isolation failed robust el id",prefix),
				  100, 0., 20.);
  helTrkPatIsoFailId->Sumw2();
  helTrkPatIsoNoId    = new TH1F(Form("%s_helTrkPatIsoNoId",prefix),
				  Form("%s - electron trk PAT isolation without el id",prefix),
				  100, 0., 20.);
  helTrkPatIsoNoId->Sumw2();
  
  helEcalJuraIsoPassId = new TH1F(Form("%s_helEcalJuraIsoPassId",prefix),
				  Form("%s - electron ecal jurassic isolation based on basic clusters passed robust el id",prefix),
				  100, 0., 20.);
  helEcalJuraIsoPassId->Sumw2();
  helEcalJuraIsoFailId = new TH1F(Form("%s_helEcalJuraIsoFailId",prefix),
				  Form("%s - electron ecal jurassic isolation based on basic clusters failed robust el id",prefix),
				  100, 0., 20.);
  helEcalJuraIsoFailId->Sumw2();
  helEcalJuraIsoNoId   = new TH1F(Form("%s_helEcalJuraIsoNoId",prefix),
				  Form("%s - electron ecal jurassic isolation based on basic clusters without el id",prefix),
				  100, 0., 20.);
  helEcalJuraIsoNoId->Sumw2();
  helEcalPatIsoPassId  = new TH1F(Form("%s_helEcalPatIsoPassId",prefix),
				  Form("%s - electron ecal jurassic isolation based on rechits pass robust el id",prefix),
				  100, 0., 20.);
  helEcalPatIsoPassId->Sumw2();
  helEcalPatIsoFailId  = new TH1F(Form("%s_helEcalPatIsoFailId",prefix),
				  Form("%s - electron ecal jurassic isolation based on rechits failed robust el id",prefix),
				  100, 0., 20.);
  helEcalPatIsoFailId->Sumw2();
  helEcalPatIsoNoId    = new TH1F(Form("%s_helEcalPatIsoNoId",prefix),
				  Form("%s - electron ecal jurassic isolation based on rechits without el id",prefix),
				  100, 0., 20.);
  helEcalPatIsoNoId->Sumw2();

  helHcalConeIsoPassId = new TH1F(Form("%s_helHcalConeIsoPassId",prefix),
				  Form("%s - electron hcal cone isolation based on basic clusters passed robust el id",prefix),
				  100, 0., 20.);
  helHcalConeIsoPassId->Sumw2();
  helHcalConeIsoFailId = new TH1F(Form("%s_helHcalConeIsoFailId",prefix),
				  Form("%s - electron hcal cone isolation based on basic clusters failed robust el id",prefix),
				  100, 0., 20.);
  helHcalConeIsoFailId->Sumw2();
  helHcalConeIsoNoId   = new TH1F(Form("%s_helHcalConeIsoNoId",prefix),
				  Form("%s - electron hcal cone isolation based on basic clusters without el id",prefix),
				  100, 0., 20.);
  helHcalConeIsoNoId->Sumw2();
  helHcalPatIsoPassId  = new TH1F(Form("%s_helHcalPatIsoPassId",prefix),
				  Form("%s - electron hcal cone isolation based on rechits pass robust el id",prefix),
				  100, 0., 20.);
  helHcalPatIsoPassId->Sumw2();
  helHcalPatIsoFailId  = new TH1F(Form("%s_helHcalPatIsoFailId",prefix),
				  Form("%s - electron hcal cone isolation based on rechits failed robust el id",prefix),
				  100, 0., 20.);
  helHcalPatIsoFailId->Sumw2();
  helHcalPatIsoNoId    = new TH1F(Form("%s_helHcalPatIsoNoId",prefix),
				  Form("%s - electron hcal cone isolation based on rechits without el id",prefix),
				  100, 0., 20.);
  helHcalPatIsoNoId->Sumw2();
  
  helRelIsoPassId = new TH1F(Form("%s_helRelIsoPassId",prefix),
				  Form("%s - electron relative iso (trk+ecal+hcal) CMS2 with weights 1,1,1 passed robust el id",prefix),
				  120, 0., 1.2);
  helRelIsoPassId->Sumw2();
  helRelIsoFailId = new TH1F(Form("%s_helRelIsoFailId",prefix),
				  Form("%s - electron relative iso (trk+ecal+hcal) CMS2 with weights 1,1,1 failed robust el id",prefix),
				  120, 0., 1.2);
  helRelIsoFailId->Sumw2();
  helRelIsoNoId   = new TH1F(Form("%s_helRelIsoNoId",prefix),
				  Form("%s - electron relative iso (trk+ecal+hcal) CMS2 with weights 1,1,1 without el id",prefix),
				  120, 0., 1.2);
  helRelIsoNoId->Sumw2();
  helRelPatIsoPassId  = new TH1F(Form("%s_helRelPatIsoPassId",prefix),
				  Form("%s - electron relative iso (trk+ecal+hcal) PAT with weights 1,1,1 passed robust el id",prefix),
				  120, 0., 1.2);
  helRelPatIsoPassId->Sumw2();
  helRelPatIsoFailId  = new TH1F(Form("%s_helRelPatIsoFailId",prefix),
				  Form("%s - electron relative iso (trk+ecal+hcal) PAT with weights 1,1,1 failed robust el id",prefix),
				  120, 0., 1.2);
  helRelPatIsoFailId->Sumw2();
  helRelPatIsoNoId    = new TH1F(Form("%s_helRelPatIsoNoId",prefix),
				  Form("%s - electron relative iso (trk+ecal+hcal) PAT with weights 1,1,1 without el id",prefix),
				  120, 0., 1.2);
  helRelPatIsoNoId->Sumw2();

  hemElRelIso = new TH1F(Form("%s_hemElRelIso",prefix),
			 Form("%s - electron relative iso for emu final selection",prefix),
			 120, 0., 1.2);
  hemElRelIso->Sumw2();
  hemMuRelIso = new TH1F(Form("%s_hemMuRelIso",prefix),
			 Form("%s - muon relative iso for emu final selection",prefix),
			 120, 0., 1.2);
  hemMuRelIso->Sumw2();

  hmaxJPTEt = new TH1F(Form("%s_hmaxJPTEt",prefix), Form("%s - most energetic jet Et (JPT)",prefix), 100, 0., 100);
  hmaxJPTEt->Sumw2();
  hmaxCaloJetEt = new TH1F(Form("%s_hmaxCaloJetEt",prefix), Form("%s - most energetic jet Et (CaloJet)",prefix), 100, 0., 100);
  hmaxCaloJetEt->Sumw2();
  hmaxTrkJetEt = new TH1F(Form("%s_hmaxTrkJetEt",prefix), Form("%s - most energetic jet Et (TrkJet)",prefix), 100, 0., 100);
  hmaxTrkJetEt->Sumw2();
  hmaxCaloTrkJetEt = new TH1F(Form("%s_hmaxCaloTrkJetEt",prefix), Form("%s - most energetic jet Et (average of Calo + Trk Jets)",prefix), 100, 0., 100);
  hmaxCaloTrkJetEt->Sumw2();
  hmaxCaloTrkJet2Et = new TH1F(Form("%s_hmaxCaloTrkJet2Et",prefix), Form("%s - most energetic jet Et (Max Calo and Trk Jets)",prefix), 100, 0., 100);
  hmaxCaloTrkJet2Et->Sumw2();
  hmaxGenJetEt = new TH1F(Form("%s_hmaxGenJetEt",prefix), Form("%s - most energetic jet Et (GenJet)",prefix), 100, 0., 100);
  hmaxGenJetEt->Sumw2();

  // clear list of duplicates
  already_seen.clear();
  int duplicates_total_n = 0;
  double duplicates_total_weight = 0;

  int i_permille_old = 0;
  // file loop
  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  monitor.counters.clear();
  while (TChainElement *currentFile = (TChainElement*)fileIter.Next()) {
       // need to call TFile::Open(), since the file is not
       // necessarily a plain TFile (TNetFile, TDcacheFile, etc)
//        printf("current file: %s (%s), %s\n", currentFile->GetName(), 
// 	      currentFile->GetTitle(), currentFile->IsA()->GetName());
       TFile *f = TFile::Open(currentFile->GetTitle()); 
       assert(f);
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
	    DorkyEventIdentifier id = { cms2.evt_run(), cms2.evt_event(), cms2.evt_lumiBlock(), cms2.trks_d0()[0], 
					cms2.hyp_lt_p4()[0].pt(), cms2.hyp_lt_p4()[0].eta(), cms2.hyp_lt_p4()[0].phi() };
	    if (is_duplicate(id)) {
		 duplicates_total_n++;
		 duplicates_total_weight += cms2.evt_scale1fb();
		 // cout << "Duplicate event found. Run: " << cms2.evt_run() << ", Event:" << cms2.evt_event() << ", Lumi: " << cms2.evt_lumiBlock() << endl;
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
	    // fkw, end of per event filling of histos.
	    
	    // loop over hypothesis candidates
	    unsigned int nHyps = cms2.hyp_type().size();
	    for( unsigned int i_hyp = 0; i_hyp < nHyps; ++i_hyp ) {
	      hypo(i_hyp, kFactor, dataset);
	      AddIsoSignalControlSample(i_hyp, kFactor, dataset);
	    }
       }
       t.Stop();
       printf("Finished processing file: %s\n",currentFile->GetTitle());
       printf("Real time: %u events / %f s = %e event/s\n", nEvents, 
	      t.RealTime(), nEvents / t.RealTime());
       printf("CPU time: %u events / %f s = %e event/s\n", nEvents, 
	      t.CpuTime(), nEvents / t.CpuTime());
       printf("Total duplicate count: %d.  Total weight %f\n",   
	      duplicates_total_n, duplicates_total_weight);
       delete f;
  }
  monitor.print();
  if ( nEventsChain != nEventsTotal ) {
       printf("ERROR: number of events from files (%d) is not equal to total number"
	      " of events (%d)\n", nEventsChain, nEventsTotal);
  }

  printf("Total candidate count (ee mm em all): %.0f %.0f %.0f %0.f.  Total weight %f %f %f %f\n",   
	 hypos_total->GetBinContent(1), hypos_total->GetBinContent(2), 
	 hypos_total->GetBinContent(3), hypos_total->GetBinContent(4),
	 hypos_total_weighted->GetBinContent(1), hypos_total_weighted->GetBinContent(2), 
	 hypos_total_weighted->GetBinContent(3), hypos_total_weighted->GetBinContent(4));
  
  return dataset;
}

