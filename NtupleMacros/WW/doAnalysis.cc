//now make the source file
#include "doAnalysis.h"
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
#include "Math/VectorUtil.h"
#include "TSystem.h"
#include "TPRegexp.h"
#include "monitor.h"

using namespace std;

#ifndef __CINT__
#include "CORE/CMS2.h"
#include "CORE/electronSelections.h"
#include "CORE/muonSelections.h"
#endif

//
// Key analysis method implementation
//

bool goodElectronWithoutIsolation(unsigned int i){
  return ww_elBase(i) && ww_elId(i) && ww_eld0PV(i);
}

bool goodElectronIsolated(unsigned int i){
  return ww_elBase(i) && ww_elId(i) && ww_eld0PV(i) && ww_elIso(i);
  // return ww_eld0(i) && ww_elIso(i)<0.1;
}

bool fakableElectron(unsigned int i){
  // extrapolate in id
  return ww_elBase(i) && ww_eld0(i) && ww_elIso(i);
}

bool goodMuonWithoutIsolation(unsigned int i){
  return ww_muBase(i) && ww_mud0PV(i) && ww_muId(i);
}

bool goodMuonIsolated(unsigned int i){
  return ww_muBase(i) && ww_mud0PV(i) && ww_muId(i) && ww_muIso(i); 
}

// double metValue(){    return cms2.evt_tcmet(); }
// double metPhiValue(){ return cms2.evt_tcmetPhi(); }
double metValue(){    return cms2.evt35X_tcmet(); }
double metPhiValue(){ return cms2.evt35X_tcmetPhi(); }

bool passedMetRequirements(unsigned int i_hyp){
  // if ( cms2.hyp_p4().at(i_hyp).mass()>130 ) return true;
  HypothesisType type = getHypothesisType(cms2.hyp_type()[i_hyp]);
  double pMet = projectedMet(i_hyp);
  // if ( type == EM && cms2.hyp_p4().at(i_hyp).mass()>90 ) return true;
  if ( pMet < 20 ) return false;
  if (type == EE || type == MM) {
    // double dmass = fabs(cms2.hyp_p4()[i_hyp].mass()-91);
    if ( metValue() < 45 ) return false;
    // if ( !metBalance(i_hyp) ) return false;
  }
  return true;
}

JetType jetType(){
  return pfJet;
  // return jptJet;
}

unsigned int numberOfJets(unsigned int i_hyp){
  return getJets(jetType(), i_hyp, 20, 3.0).size();
}


//
// Electron Id
//

bool ww_elBase(unsigned int index){
  if (cms2.els_p4().at(index).pt() < 20.0) return false;
  if (fabs(cms2.els_p4().at(index).eta()) > 2.5) return false;
  return true;
}
bool ww_elId(unsigned int index){
  if( fabs(cms2.els_conv_dist().at(index)) < 0.02 &&
      fabs(cms2.els_conv_dcot().at(index)) < 0.02) return false;
  if (! (electronId_VBTF(index, VBTF_35X_80) & (1<<ELEID_ID)) ) return false;
  // if (! (electronId_VBTF(index, VBTF_35X_70) & (1<<ELEID_ID)) ) return false;
  // if (! (electronId_CIC(index, 4, CIC_SUPERTIGHT) & (1<<ELEID_ID)) ) return false;
  
  // conversion rejection - hit based
  //if ( cms2.els_exp_innerlayers().at(index) > 0 ) return false;
  //  int ctfIndex = cms2.els_trkidx().at(index);
  // if ( ctfIndex >=0 && 
  //     cms2.els_charge().at(index)!=cms2.trks_charge().at(ctfIndex) ) return false;
  return true;
}

bool ww_eld0(unsigned int index){
  return fabs(cms2.els_d0corr()[index]) < 0.02;
}

bool ww_eld0PV(unsigned int index){
  if ( cms2.vtxs_sumpt().empty() ) return false;
  unsigned int iMax = 0;
  double sumPtMax = cms2.vtxs_sumpt().at(0);
  for ( unsigned int i = iMax+1; i < cms2.vtxs_sumpt().size(); ++i )
    if ( cms2.vtxs_sumpt().at(i) > sumPtMax ){
      iMax = i;
      sumPtMax = cms2.vtxs_sumpt().at(i);
    }
  double dxyPV = cms2.els_d0()[index]-
    cms2.vtxs_position()[iMax].x()*sin(cms2.els_trk_p4()[index].phi())+
    cms2.vtxs_position()[iMax].y()*cos(cms2.els_trk_p4()[index].phi());
  return fabs(dxyPV) < 0.02;
}

double ww_elIsoVal(unsigned int index){
  float sum = cms2.els_tkIso().at(index);
  sum += max(0., (cms2.els_ecalIso().at(index) -1.));
  sum += cms2.els_hcalIso().at(index);
  double pt = cms2.els_p4().at(index).pt();
  return sum/pt;
}

bool ww_elIso(unsigned int index){
  return ww_elIsoVal(index)<0.1;
}

//
// Muon Id
//

bool ww_muBase(unsigned int index){
  if (cms2.mus_p4().at(index).pt() < 20.0) return false;
  if (fabs(cms2.mus_p4().at(index).eta()) > 2.4) return false;
  return true;
}
bool ww_mud0(unsigned int index){
  return fabs(cms2.mus_d0corr()[index]) < 0.02;
}
bool ww_mud0PV(unsigned int index){
  if ( cms2.vtxs_sumpt().empty() ) return false;
  unsigned int iMax = 0;
  double sumPtMax = cms2.vtxs_sumpt().at(0);
  for ( unsigned int i = iMax+1; i < cms2.vtxs_sumpt().size(); ++i )
    if ( cms2.vtxs_sumpt().at(i) > sumPtMax ){
      iMax = i;
      sumPtMax = cms2.vtxs_sumpt().at(i);
    }
  double dxyPV = cms2.mus_d0()[index]-
    cms2.vtxs_position()[iMax].x()*sin(cms2.mus_trk_p4()[index].phi())+
    cms2.vtxs_position()[iMax].y()*cos(cms2.mus_trk_p4()[index].phi());
  return fabs(dxyPV) < 0.02;
}
bool ww_muId(unsigned int index){
  if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) return false; //glb fit chisq
  if (((cms2.mus_type().at(index)) & (1<<1)) == 0)    return false; // global muon
  if (((cms2.mus_type().at(index)) & (1<<2)) == 0)    return false; // tracker muon
  if (cms2.mus_validHits().at(index) < 11)            return false; // # of tracker hits
  if (cms2.mus_gfit_validSTAHits().at(index)==0 ) return false;
  return true;
}

double ww_muIsoVal(unsigned int index){
  double sum =  cms2.mus_iso03_sumPt().at(index) +
    cms2.mus_iso03_emEt().at(index)  +
    cms2.mus_iso03_hadEt().at(index);
  double pt  = cms2.mus_p4().at(index).pt();
  return sum/pt;
}
bool ww_muIso(unsigned int index){
  return ww_muIsoVal(index)<0.15;
}
unsigned int numberOfSoftMuons(int i_hyp, bool nonisolated){
  unsigned int nMuons = 0;
  for (int imu=0; imu < int(cms2.mus_charge().size()); ++imu) {
    // quality cuts
    // if (  ((cms2.mus_goodmask()[imu]) & (1<<14)) == 0 ) continue; // TMLastStationOptimizedLowPtTight
    if (  ((cms2.mus_goodmask()[imu]) & (1<<19)) == 0 ) continue; // TMLastStationAngTight
    if ( cms2.mus_p4()[imu].pt() < 3 ) continue;
    if ( TMath::Abs(cms2.mus_d0corr()[imu]) > 0.2) continue;
    if ( cms2.mus_validHits()[imu] < 11) continue;
    if ( TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && cms2.hyp_lt_index()[i_hyp] == imu ) continue;
    if ( TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && cms2.hyp_ll_index()[i_hyp] == imu ) continue;
    if ( nonisolated && ww_muIsoVal(imu)>0.1 && cms2.mus_p4()[imu].pt()>20 ) continue;
    ++nMuons;
  }
  return nMuons;
}

unsigned int numberOfExtraLeptons(int i_hyp, double minPt){
  unsigned int nMuons = 0;
  for (int i=0; i < int(cms2.mus_charge().size()); ++i) {
    // printf("Muon: %u, pt: %0.2f\n",i,cms2.mus_p4().at(i).pt());
    if ( cms2.mus_p4()[i].pt() < minPt ) continue;
    // printf("\tpassed minPt\n");
    if ( TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && cms2.hyp_lt_index()[i_hyp] == i ) continue;
    if ( TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && cms2.hyp_ll_index()[i_hyp] == i ) continue;
    // printf("\tpassed hyp letpons\n");
    if ( ! (ww_mud0PV(i) && ww_muId(i) && ww_muIso(i)&&
	    fabs(cms2.mus_p4().at(i).eta()) <2.4) ) continue;
    // printf("\tpassed all\n");
    ++nMuons;
  }
  unsigned int nElectrons = 0;
  for (int i=0; i < int(cms2.els_charge().size()); ++i) {
    // printf("Electron: %u, pt: %0.2f\n",i,cms2.els_p4().at(i).pt());
    if ( cms2.els_p4()[i].pt() < minPt ) continue;
    // printf("\tpassed minPt\n");
    if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.els_p4().at(i)) <0.1) ) continue;
    if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.els_p4().at(i)) <0.1) ) continue;
    // printf("\tpassed hyp letpons\n");
    if ( !(ww_elId(i) && ww_eld0PV(i) && ww_elIso(i) && 
	   fabs(cms2.els_p4().at(i).eta()) < 2.5) ) continue;
    // printf("\tpassed all\n");
    ++nElectrons;
  }
  return nMuons+nElectrons;
}


//
// Triger
//

bool passedTriggerRequirements(HypTypeInNtuples type) {
  bool hlt_ele15_lw_l1r = cms2.passHLTTrigger("HLT_Ele15_LW_L1R");
  bool hltMu9           = cms2.passHLTTrigger("HLT_Mu9");
  return hlt_ele15_lw_l1r || hltMu9;
  /*
  if (type == MuMu && ! (hltMu9) ) return false;
  if ((type == ElMu || type == MuEl) && ! (hltMu9 || hlt_ele15_lw_l1r)) return false;
  if (type == ElEl && ! hlt_ele15_lw_l1r) return false;     

  return true;
  */
}

//
// MET
//
double nearestDeltaPhi(double Phi, int i_hyp)
{
  double tightDPhi = fabs(cms2.hyp_lt_p4()[i_hyp].Phi() - Phi);
  tightDPhi = std::min(2*TMath::Pi() - tightDPhi, tightDPhi);
  double looseDPhi = fabs(cms2.hyp_ll_p4()[i_hyp].Phi() - Phi);
  looseDPhi = std::min(2*TMath::Pi() - looseDPhi, looseDPhi);
  return TMath::Min(tightDPhi, looseDPhi);
}

double projectedMet(unsigned int i_hyp)
{
  double DeltaPhi = nearestDeltaPhi(metPhiValue(),i_hyp);
  if (DeltaPhi < TMath::Pi()/2) return metValue()*TMath::Sin(DeltaPhi);
  return metValue();
}

bool metBalance (unsigned int i_hyp) {
  if( metValue()/cms2.hyp_p4()[i_hyp].pt() < 0.9 ) return false;
  // if( metValue()/cms2.hyp_p4()[i_hyp].pt() < 0.6 &&
  //acos(cos(metPhiValue()-cms2.hyp_p4()[i_hyp].phi() - 3.1416)) < 0.25 ) return false;
  return true;
}

HypTypeInNtuples hypType(unsigned int i_hyp){
  HypTypeInNtuples type = HypTypeInNtuples(cms2.hyp_type().at(i_hyp));
  return type;
}

bool ww2009_met(unsigned int i_hyp){
  HypTypeInNtuples type = hypType(i_hyp);
  double pMet = projectedMet(i_hyp);

  if ( pMet < 20 ) return false;

  if (type == ElEl || type == MuMu) {
    if ( metValue() < 45 ) return false;
    if ( !metBalance(i_hyp) ) return false;
  }
  return true;
}

//
// Jets
//

Bool_t comparePt(LorentzVector lv1, LorentzVector lv2) {
   return lv1.pt() > lv2.pt();
}

std::vector<LorentzVector> 
getJets(JetType type, int i_hyp, double etThreshold, double etaMax, bool sortJets)
{
     std::vector<LorentzVector> jets;
     const double vetoCone    = 0.4;
     
     switch ( type ){
     case jptJet:
       for ( unsigned int i=0; i < cms2.jpts_p4().size(); ++i) {
	 if ( cms2.jpts_p4()[i].Et() < etThreshold ) continue;
	 if ( TMath::Abs(cms2.jpts_p4()[i].eta()) > etaMax ) continue;
	 if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.jpts_p4()[i])) < vetoCone ||
	      TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.jpts_p4()[i])) < vetoCone ) continue;
	 jets.push_back(cms2.jpts_p4()[i]);
       }
       break;
     case pfJet:
       for ( unsigned int i=0; i < cms2.pfjets_p4().size(); ++i) {
	 if ( cms2.pfjets_p4()[i].pt() < etThreshold ) continue;
	 if ( TMath::Abs(cms2.pfjets_p4()[i].eta()) > etaMax ) continue;
	 if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.pfjets_p4()[i])) < vetoCone ||
	      TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.pfjets_p4()[i])) < vetoCone ) continue;
	 jets.push_back(cms2.pfjets_p4()[i]);
       }
       break;
     case GenJet:
       for ( unsigned int i=0; i < cms2.genjets_p4().size(); ++i) {
	 if ( cms2.genjets_p4()[i].Et() < etThreshold ) continue;
	 if ( TMath::Abs(cms2.genjets_p4()[i].eta()) > etaMax ) continue;
	 if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.genjets_p4()[i])) < vetoCone ||
	      TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.genjets_p4()[i])) < vetoCone ) continue;
	 jets.push_back(cms2.genjets_p4()[i]);
       }
       break;
     case CaloJet:
       for ( unsigned int i=0; i < cms2.jets_pat_jet_p4().size(); ++i) {
	 if ( cms2.jets_pat_jet_p4()[i].Et() < etThreshold ) continue;
	 if ( TMath::Abs(cms2.jets_pat_jet_p4()[i].eta()) > etaMax ) continue;
	 if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.jets_pat_jet_p4()[i])) < vetoCone ||
	      TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.jets_pat_jet_p4()[i])) < vetoCone ) continue;
	 jets.push_back(cms2.jets_pat_jet_p4()[i]);
       }
       break;
     case TrkJet:
       for ( unsigned int i=0; i < cms2.trkjets_p4().size(); ++i) {
	 if ( cms2.trkjets_p4()[i].Et() < etThreshold ) continue;
	 if ( TMath::Abs(cms2.trkjets_p4()[i].eta()) > etaMax ) continue;
	 if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.trkjets_p4()[i])) < vetoCone ||
	      TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.trkjets_p4()[i])) < vetoCone ) continue;
	 jets.push_back(cms2.trkjets_p4()[i]);
       }
       break;
     default:
       std::cout << "ERROR: not supported jet type is requested: " << type << " FixIt!" << std::endl;
     }
     if ( sortJets ) std::sort(jets.begin(), jets.end(), comparePt);
     return jets;
}

//
// Various other cuts
//

bool inZmassWindow(float mass){
  // return ( mass > 76. && mass < 106. );
  return fabs(mass - 91.1876) < 15;
}

double BTag(JetType type, unsigned int iJet){
  // find a jet in the jets list
  // that matches the current type jet
  LorentzVector jetP4;
  switch ( type ) {
  case jptJet:
    jetP4 = cms2.jpts_p4().at(iJet);
    break;
  case CaloJet:
    // return cms2.jets_jetProbabilityBJetTag().at(iJet);
    return cms2.jets_combinedSecondaryVertexBJetTag().at(iJet);
    break;
  default:
    std::cout << "ERROR: not supported jet type is requested: " << type << " FixIt!" << std::endl;
    assert(0);
  }
  int refJet = -1;
  for ( unsigned int i=0; i < cms2.jets_p4().size(); ++i) {
    if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(jetP4,cms2.jets_p4()[i])) > 0.3 ) continue;
    refJet = i;
  }
  if (refJet == -1){
    // std::cout << "Warning: failed to find a matching jet for b-tagging." << std::endl; 
    return 0.0;
  }
  // return cms2.jets_jetProbabilityBJetTag().at(refJet);
  return cms2.jets_combinedSecondaryVertexBJetTag().at(refJet);
}

//
// Tools
//

static std::set<EventIdentifier> already_seen;
bool is_duplicate (const EventIdentifier &id)
{
     std::pair<std::set<EventIdentifier>::const_iterator, bool> ret = 
	  already_seen.insert(id);
     return !ret.second;
}

bool isDYee() {
  if (getDrellYanType() == 0) return true;
  return false;
}
bool isDYmm() {
  if (getDrellYanType() == 1) return true;
  return false;
}
bool isDYtt() {
  if (getDrellYanType() == 2) return true;
  return false;
}

bool isWW() {
  if (getVVType() == 0) return true;
  return false;
}

bool isWZ() {
  if (getVVType() == 1) return true;
  return false;
}

bool isZZ() {
  if (getVVType() == 2) return true;
  return false;
}

//-------------------------------------------------
// Auxiliary function to scan the doc line and 
// identify DY-> ee vs mm vs tt
//-------------------------------------------------
unsigned int getDrellYanType() {
  bool foundEP = false;
  bool foundEM = false;
  bool foundMP = false;
  bool foundMM = false;
  bool foundTP = false;
  bool foundTM = false;
  for (unsigned int i = 0; i < cms2.genps_id().size(); ++i) {
    if ( cms2.genps_id_mother().at(i) == 23 ){
      switch ( TMath::Abs(cms2.genps_id().at(i)) ){
      case 11:
	return 0;
	break;
      case 13:
	return 1;
	break;
      case 15:
	return 2;
	break;
      default:
	break;
      }
    }
    switch ( cms2.genps_id().at(i) ){
    case 11:
      foundEM = true;
      break;
    case -11:
      foundEP = true;
      break;
    case 13:
      foundMM = true;
      break;
    case -13:
      foundMP = true;
      break;
    case 15:
      foundTM = true;
      break;
    case -15:
      foundTP = true;
      break;
    default:
      break;
    }
  }
  
  if ( foundEP && foundEM ) return 0;  //DY->ee
  if ( foundMP && foundMM ) return 1;  //DY->mm
  if ( foundTP && foundTM ) return 2;  //DY->tautau
  std::cout << "Does not look like a DY event" << std::endl;
  return 999;
}

unsigned int getVVType() {
  // types:
  //   0 - WW
  //   1 - WZ
  //   2 - ZZ
  unsigned int nZ(0);
  unsigned int nW(0);
  std::vector<std::vector<int> > leptons;
  std::vector<int> mothers;

  bool verbose = false;

  for (unsigned int i = 0; i < cms2.genps_id().size(); ++i) {
    int pid = cms2.genps_id().at(i);
    int mid = cms2.genps_id_mother().at(i);
    if ( verbose ) std::cout << "Gen particle id: " << pid << ",\t mother id: " << mid <<std::endl;
    if ( abs(pid)<11 || abs(pid)>16 ) continue;
    if ( mid == 23 ) ++nZ;
    if ( abs(mid) == 24 ) ++nW;
    // now we need to really understand the pattern.
    unsigned int mIndex = 0;
    while ( mIndex < mothers.size() && mid != mothers[mIndex] ) ++mIndex;
    if ( mIndex == mothers.size() ) {
      mothers.push_back(mid);
      leptons.push_back(std::vector<int>());
    }
    leptons[mIndex].push_back(pid);
    if (mothers.size()>3){
      if (verbose) std::cout << "WARNING: failed to identify event (too many mothers)" << std::endl;
      return 999;
    }
  }

  if ( nZ == 4 ) {
    if ( verbose ) std::cout << "Event type ZZ" << std::endl;
    return 2;
  }
  if ( nW == 4 ) {
    if ( verbose ) std::cout << "Event type WW" << std::endl;
    return 0;
  }
  if ( nW == 2 && nZ == 2 ) {
    if ( verbose ) std::cout << "Event type WZ" << std::endl;
    return 1;
  }
  unsigned int nNus(0);
  for ( unsigned int i=0; i<mothers.size(); ++i ){
      nNus += leptons[i].size();
  }
  if ( mothers.size() < 3 && nNus == 4){
    for ( unsigned int i=0; i<mothers.size(); ++i ){
      if ( mothers[i] != 23 && abs(mothers[i]) != 24 ){
	if( leptons[i].size() != 2 && leptons[i].size() != 4){
	  if (verbose) std::cout << "WARNING: failed to identify event (unexpected number of daughters)" << std::endl;
	  if (verbose) std::cout << "\tnumber of daughters for first mother: " <<  leptons[0].size() << std::endl;
	  if (verbose) std::cout << "\tnumber of daughters for second mother: " <<  leptons[1].size() << std::endl;
	  return 999;
	}
	if ( abs(leptons[i][0]) == abs(leptons[i][1]) )
	  nZ += 2;
	else
	  nW += 2;
	if ( leptons[i].size()==4 ){
	  // now it's a wild guess, it's fraction should be small
	  if ( abs(leptons[i][2]) == abs(leptons[i][3]) )
	    nZ += 2;
	  else
	    nW += 2;
	}
      }
    }
  } else {
    // here be dragons
    
    // if we have 2 leptons and 3 neutrinos and they all of the same
    // generation, we assume it's ZZ (can be WZ also), but if
    // neutrinos are from different generations, than we conclude it's
    // WZ. 
    
    std::set<int> nus;
    for ( unsigned int i=0; i<mothers.size(); ++i )
      for ( unsigned int j=0; j<leptons[i].size(); ++j ) 
	if ( abs(leptons[i][j]) == 12 ||
	     abs(leptons[i][j]) == 14 ||
	     abs(leptons[i][j]) == 16 )
	  nus.insert(abs(leptons[i][j]));
    
    if ( nNus == 5 ){
      if ( nus.size() == 1 ) return 2;
      if ( nus.size() == 2 ) return 1;
    }
    
    if ( verbose ) std::cout << "WARNING: failed to identify event" << std::endl;
    return 999;
  }

  if ( nZ+nW != 4 ){
    if (verbose) std::cout << "WARNING: failed to identify event (wrong number of bosons)" << std::endl;
    if (verbose) std::cout << "\tfirst mother id: " << mothers[0] << std::endl;
    if (verbose) std::cout << "\tsecond mother id: " << mothers[1] << std::endl;
    if (verbose) std::cout << "\tnumber of daughters for first mother: " << leptons[0].size() << std::endl;
    if (verbose) std::cout << "\tnumber of daughters for second mother: " << leptons[1].size() << std::endl;
    if (verbose) std::cout << "\tnumber of Zs: " << nZ << std::endl;
    if (verbose) std::cout << "\tnumber of Ws: " << nW << std::endl;
    return 999;
  }

  if ( nZ == 4 ) {
    if ( verbose ) std::cout << "Event type ZZ" << std::endl;
    return 2;
  }
  if ( nW == 4 ) {
    if ( verbose ) std::cout << "Event type WW" << std::endl;
    return 0;
  }
  // this covers screws in logic, i.e. most hard to identify events end up being WZ
  if ( verbose ) std::cout << "Event type WZ (can be wrong)" << std::endl;
  return 1;
}

//
// Histograms
// 
//  Book histograms...
//  Naming Convention:
//  Prefix comes from the sample and it is passed to the scanning function
//  Suffix is "ee" "em" "em" "all" which depends on the final state
//  For example: histogram named tt_hnJet_ee would be the Njet distribution
//  for the ee final state in the ttbar sample.

// MAKE SURE TO CAL SUMW2 FOR EACH 1D HISTOGRAM BEFORE FILLING!!!!!!

// Histogram with all cuts applied

TH1F* hnJet[4];      // Njet distributions
TH1F* helePt[4];     // electron Pt
TH1F* hmuPt[4];      // muon Pt
TH1F* hminLepPt[4];  // minimum lepton Pt
TH1F* hmaxLepPt[4];  // maximum lepton Pt
TH1F* helePhi[4];    // electron phi
TH1F* hmuPhi[4];     // muon phi
TH1F* hdphiLep[4];   // delta phi between leptons
TH1F* heleEta[4];    // electron eta
TH1F* hmuEta[4];     // muon eta
TH1F* hdilMass[4];   // dilepton mass
TH1F* hdilPt[4];     // dilepton Pt
TH1F* hmet[4];       // MET
TH1F* hmetPhi[4];    // MET phi
TH1F* hptJet1[4];    // Pt of 1st jet
TH1F* hptJet2[4];    // Pt of 2nd jet
TH1F* hptJet3[4];    // Pt of 3rd jet
TH1F* hetaJet1[4];   // eta of 1st jet
TH1F* hetaJet2[4];   // eta of 2nd jet
TH1F* hetaJet3[4];   // eta of 3rd jet
TH1F* hElRelIso[4];  // electron relative iso (final selection, but one electron iso is relaxed)
TH1F* hMuRelIso[4];  // muon relative iso (final selection, but one muon iso is relaxed

TH2F* hmetVsDilepPt[4];    // MET vs dilepton Pt
TH2F* hmetOverPtVsDphi[4]; // MET/Lepton Pt vs DeltaPhi between MET and Lepton Pt
TH2F* hdphillvsmll[4];     // delta phi between leptons vs dilepton mass
TH2F* helFRfakable[4];     // Fake rate study: rate of fakable objects
TH2F* hFakableRateSingleElectron; // Fake rate study: rate of fakable objects
TH2F* hFinalRateSingleElectron;   // Fake rate study: rate of final objects
TH1F* hIsoSingleMuon;             // isolation background 
TH1F* hIsoSingleElectron;         // isolation background 

TH1F* hmaxJPTEt;         // energy distribution for the most energetic jet
TH2F* hmaxBtagVsJPTEt;   // btag vs energy distribution for the most energetic jet

//
// Not cleaned area
//

TH2F* helFRfakable_fakerate[4];     // Fake rate study: rate of fakable objects using the fakerate.cc
TH1F* hypos_total;
TH1F* hypos_total_weighted;

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

TH1F* hmaxCaloJetEt;     // energy distribution for the most energetic jet
TH1F* hmaxTrkJetEt;      // energy distribution for the most energetic jet
TH1F* hmaxCaloTrkJetEt;  // energy distribution for the most energetic jet
TH1F* hmaxCaloTrkJet2Et; // energy distribution for the most energetic jet
TH1F* hmaxGenJetEt;      // energy distribution for the most energetic jet

TH1F* hCentralBquarkEtaAfterVeto;   
TH1F* hForwardBquarkEtaAfterVeto;   

// fkw September 2008 final hist used for muon tags estimate of top bkg
TH2F* hextramuonsvsnjet[4];

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
    if ( electronId_classBasedLoose(cms2.hyp_lt_index()[i_hyp]) ) {
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
    if ( electronId_classBasedLoose(cms2.hyp_ll_index()[i_hyp]) ) {
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
  HypothesisType type = getHypothesisType(cms2.hyp_type()[i_hyp]);
  RooArgSet set( *(dataset->get()) );
  set.setCatIndex("selected",passedAllLeptonRequirements?1:0);
  set.setRealValue("event",cms2.evt_event());
  set.setRealValue("run",cms2.evt_run());
  set.setRealValue("lumi",cms2.evt_lumiBlock());
  set.setCatLabel("sample_type","data_relaxed_iso");
  
  // em case
  if ( type == EM ){
    unsigned int imu = cms2.hyp_lt_index()[i_hyp];
    unsigned int iel = cms2.hyp_ll_index()[i_hyp];
    if ( TMath::Abs(cms2.hyp_lt_id()[i_hyp])==11 ){
      imu = cms2.hyp_ll_index()[i_hyp];
      iel = cms2.hyp_lt_index()[i_hyp];
    }
    if ( goodElectronWithoutIsolation(iel) && 
	 goodMuonIsolated(imu) ) {
      hElRelIso[type]->Fill( ww_elIsoVal(iel), weight );
      set.setCatLabel("hyp_type","em");
      set.setCatLabel("fake_type","electron");
      set.setRealValue("iso", ww_elIsoVal(iel) );
      dataset->add(set,weight);
    
    }
    if ( goodElectronIsolated(iel) && goodMuonWithoutIsolation(imu) ) {
      hMuRelIso[type]->Fill( ww_muIsoVal(imu), weight );
      set.setCatLabel("hyp_type","em");
      set.setCatLabel("fake_type","muon");
      set.setRealValue("iso",ww_muIsoVal(imu));
      dataset->add(set,weight);
    }
  }

  // mm case
  if ( type == MM ){
    unsigned int imu1 = cms2.hyp_lt_index()[i_hyp];
    unsigned int imu2 = cms2.hyp_ll_index()[i_hyp];
    if ( goodMuonWithoutIsolation(imu1) && goodMuonIsolated(imu2) ) {
      hMuRelIso[type]->Fill( ww_muIsoVal(imu1), weight );
      set.setCatLabel("hyp_type","mm");
      set.setCatLabel("fake_type","muon");
      set.setRealValue("iso",ww_muIsoVal(imu1));
      dataset->add(set,weight);
    }
    if ( goodMuonIsolated(imu1) && goodMuonWithoutIsolation(imu2) ) {
      hMuRelIso[type]->Fill( ww_muIsoVal(imu2), weight );
      set.setCatLabel("hyp_type","mm");
      set.setCatLabel("fake_type","muon");
      set.setRealValue("iso",ww_muIsoVal(imu2));
      dataset->add(set,weight);
    }
  }

  // ee case
  if ( type == EE){
    unsigned int iel1 = cms2.hyp_lt_index()[i_hyp];
    unsigned int iel2 = cms2.hyp_ll_index()[i_hyp];
    if ( goodElectronWithoutIsolation(iel1) && goodElectronIsolated(iel2) ) {
      hElRelIso[type]->Fill( ww_elIsoVal(iel1), weight );
      set.setCatLabel("hyp_type","ee");
      set.setCatLabel("fake_type","electron");
      set.setRealValue("iso",ww_elIsoVal(iel1));
      dataset->add(set,weight);
    }
    if ( goodElectronIsolated(iel1) && goodElectronWithoutIsolation(iel2) ) {
      hElRelIso[type]->Fill( ww_elIsoVal(iel2), weight );
      set.setCatLabel("hyp_type","ee");
      set.setCatLabel("fake_type","electron");
      set.setRealValue("iso",ww_elIsoVal(iel2));
      dataset->add(set,weight);
    }
  }
}

void find_most_energetic_jets(int i_hyp, double weight)
{
  {
    double jptMax(0.);
    int jptMaxIndex(-1);
    // double 
    for ( unsigned int i=0; i < cms2.jpts_p4().size(); ++i) {
      if ( cms2.jpts_p4()[i].Et() < jptMax ) continue;
      if ( TMath::Abs(cms2.jpts_p4()[i].eta()) > 3.0 ) continue;
      if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.jpts_p4()[i])) < 0.4 ||
	   TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.jpts_p4()[i])) < 0.4 ) continue;
      jptMax = cms2.jpts_p4()[i].Et();
      jptMaxIndex = i;
    }
    hmaxJPTEt->Fill(jptMax, weight);
    if (jptMaxIndex >= 0)
      hmaxBtagVsJPTEt->Fill(jptMax, BTag(jptJet,jptMaxIndex), weight);
    else
      hmaxBtagVsJPTEt->Fill(jptMax, 0.0, weight);
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

void extractFakeRateSingleLepton(){
  // weights are ignore. see no reason to complicate stuff
  for ( unsigned int i=0; i < cms2.els_p4().size(); ++i ){
    if ( cms2.els_p4().at(i).pt() < 20 ) continue;
    if ( fakableElectron(i) ) 
      hFakableRateSingleElectron->Fill(cms2.els_p4().at(i).pt(), fabs(cms2.els_p4().at(i).eta()) );
    if ( goodElectronIsolated(i) ) 
      hFinalRateSingleElectron->Fill(cms2.els_p4().at(i).pt(), fabs(cms2.els_p4().at(i).eta()) );  
  }
}

void extractIsoSingleLepton(){
  // weights are ignore. see no reason to complicate stuff
  for ( unsigned int i=0; i < cms2.mus_p4().size(); ++i ){
    if ( cms2.mus_p4().at(i).pt() < 20 ) continue;
    if ( goodMuonWithoutIsolation(i) ) hIsoSingleMuon->Fill(ww_muIsoVal(i));
  }
  for ( unsigned int i=0; i < cms2.els_p4().size(); ++i ){
    if ( cms2.els_p4().at(i).pt() < 20 ) continue;
    if ( goodElectronWithoutIsolation(i) ) hIsoSingleElectron->Fill(ww_elIsoVal(i));
  }
}

void countFakableObjectsAfterAllSelections(unsigned int i_hyp, 
					   double weight, 
					   bool passedLTElFakableRequirements, 
					   bool passedLLElFakableRequirements,
					   bool passedLTFinalRequirements,
					   bool passedLLFinalRequirements)
{
  HypothesisType type = getHypothesisType(cms2.hyp_type()[i_hyp]);
  if ( passedLTElFakableRequirements && !passedLTFinalRequirements && passedLLFinalRequirements ){
    helFRfakable[type]->Fill(cms2.hyp_lt_p4().at(i_hyp).pt(),
			     fabs(cms2.hyp_lt_p4().at(i_hyp).eta()),
			     weight);
    helFRfakable[3]->Fill(cms2.hyp_lt_p4().at(i_hyp).pt(),
			  fabs(cms2.hyp_lt_p4().at(i_hyp).eta()),
			  weight);
    
    helFRfakable_fakerate[type]->Fill(fabs(cms2.hyp_lt_p4().at(i_hyp).eta()), 
				      cms2.hyp_lt_p4().at(i_hyp).pt(),
				      weight);
    helFRfakable_fakerate[3]->Fill( fabs(cms2.hyp_lt_p4().at(i_hyp).eta()),
				    cms2.hyp_lt_p4().at(i_hyp).pt(),
				    weight);
    

  }
  if ( passedLLElFakableRequirements && passedLTFinalRequirements && !passedLLFinalRequirements ){
    helFRfakable[type]->Fill(cms2.hyp_ll_p4().at(i_hyp).pt(),
			     fabs(cms2.hyp_ll_p4().at(i_hyp).eta()),
			     weight);
    helFRfakable[3]->Fill(cms2.hyp_ll_p4().at(i_hyp).pt(),
			  fabs(cms2.hyp_ll_p4().at(i_hyp).eta()),
			  weight);
    
    helFRfakable_fakerate[type]->Fill( fabs(cms2.hyp_ll_p4().at(i_hyp).eta()),
				       cms2.hyp_ll_p4().at(i_hyp).pt(),
				       weight);
    helFRfakable_fakerate[3]->Fill( fabs(cms2.hyp_ll_p4().at(i_hyp).eta()),
				    cms2.hyp_ll_p4().at(i_hyp).pt(),
				    weight);

  }
}

void hypo (int i_hyp, double kFactor, RooDataSet* dataset) 
{
  /*
  unsigned int nGenLeptons = 0;
  for ( unsigned int i=0; i<cms2.genps_id().size(); ++i)
    if ( abs(cms2.genps_id().at(i)) == 11 || abs(cms2.genps_id().at(i)) == 13 )
      nGenLeptons++;
  if ( nGenLeptons < 2 ) return;
  */
  // if (cms2.evt_event()!=101838) return;
  HypothesisType type = getHypothesisType(cms2.hyp_type()[i_hyp]);
  
  // The event weight including the kFactor (scaled to 1 fb-1)
  float weight = cms2.evt_scale1fb() * kFactor;

  monitor.nEvtProcessed = cms2.evt_nEvts();
  monitor.count(cms2, type, "Total number before cuts");
     
  // if ( cms2.hyp_FVFit_prob()[i_hyp] < 0.005 ) return;
  // monitor.count(cms2, type, "after vertex cut");

  if ( ! passedTriggerRequirements( hypType(i_hyp) ) )return;
  monitor.count(cms2, type, "after trigger requirements");

  // Require same sign
  if ( cms2.hyp_lt_id()[i_hyp] * cms2.hyp_ll_id()[i_hyp] > 0 ) return;

  // Baseline cuts
  if (abs(cms2.hyp_lt_id()[i_hyp]) == 13 && !ww_muBase(cms2.hyp_lt_index()[i_hyp]) ) return;
  if (abs(cms2.hyp_ll_id()[i_hyp]) == 13 && !ww_muBase(cms2.hyp_ll_index()[i_hyp]) ) return;
  if (abs(cms2.hyp_lt_id()[i_hyp]) == 11 && !ww_elBase(cms2.hyp_lt_index()[i_hyp]) ) return;
  if (abs(cms2.hyp_ll_id()[i_hyp]) == 11 && !ww_elBase(cms2.hyp_ll_index()[i_hyp]) ) return;
  
  monitor.count(cms2,type,"after previous + baseline cuts");

  // TEMPORARY
  /*{
    if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && !ww_mud0PV(cms2.hyp_lt_index()[i_hyp]) ) return;
    if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && !ww_mud0PV(cms2.hyp_ll_index()[i_hyp]) ) return;
    if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11 && !ww_eld0PV(cms2.hyp_lt_index()[i_hyp]) ) return;
    if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11 && !ww_eld0PV(cms2.hyp_ll_index()[i_hyp]) ) return;

    monitor.count(cms2,type,"after previous + d0");

    if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && ww_muIsoVal(cms2.hyp_lt_index()[i_hyp])>0.15 ) return;
    if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && ww_muIsoVal(cms2.hyp_ll_index()[i_hyp])>0.15 ) return;
    if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11 && ww_elIsoVal(cms2.hyp_lt_index()[i_hyp])>0.1 ) return;
    if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11 && ww_elIsoVal(cms2.hyp_ll_index()[i_hyp])>0.1 ) return;

    monitor.count(cms2,type,"after previous + iso");
    
    if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && !ww_muId(cms2.hyp_lt_index()[i_hyp]) ) return;
    if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && !ww_muId(cms2.hyp_ll_index()[i_hyp]) ) return;
    if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11 && 
	! (electronId_VBTF(cms2.hyp_lt_index()[i_hyp], VBTF_35X_80) & (1<<ELEID_ID))  ) return;
    if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11 && 
	! (electronId_VBTF(cms2.hyp_ll_index()[i_hyp], VBTF_35X_80) & (1<<ELEID_ID))  ) return;

    monitor.count(cms2,type,"after previous + lepton id");

    if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11 && 
	(fabs(cms2.els_conv_dist().at(cms2.hyp_lt_index()[i_hyp])) < 0.02 &&
	 fabs(cms2.els_conv_dcot().at(cms2.hyp_lt_index()[i_hyp])) < 0.02 )) return;
    if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11 && 
	(fabs(cms2.els_conv_dist().at(cms2.hyp_ll_index()[i_hyp])) < 0.02 &&
	 fabs(cms2.els_conv_dcot().at(cms2.hyp_ll_index()[i_hyp])) < 0.02 )) return;

    monitor.count(cms2,type,"after previous + conv rejection");

    if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && !goodMuonIsolated(cms2.hyp_lt_index()[i_hyp]) ) return;
    if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && !goodMuonIsolated(cms2.hyp_ll_index()[i_hyp]) ) return;
    if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11 && !goodElectronIsolated(cms2.hyp_lt_index()[i_hyp]) ) return;
    if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11 && !goodElectronIsolated(cms2.hyp_ll_index()[i_hyp]) ) return;

    monitor.count(cms2,type,"after previous + lepton id/iso");
    if ( metValue()<20 ) return;
    monitor.count(cms2,type,"after previous + met>20");
  }*/

  if (cms2.hyp_p4()[i_hyp].mass() < 12) return;
  monitor.count(cms2,type,"after previous + lepton pt and mll cuts");
     
  // check electron isolation and id (no selection at this point)
  checkIsolation(i_hyp, weight);

  // Z mass veto using hyp_leptons for ee and mumu final states
  if ( type == EE || type == MM) {
    if (inZmassWindow(cms2.hyp_p4()[i_hyp].mass())) return;
  }

  // Z veto using additional leptons in the event
  // if (additionalZveto()) return;
  monitor.count(cms2,type,"after previous + z veto cuts");
     
  // MET
  // if (cms2.evt_tcmet()<20) return;
  // monitor.count(cms2,type,"after previous + MET>20 cuts: ");
  
  if (!passedMetRequirements(i_hyp)) return;
  monitor.count(cms2,type,"after previous + Full MET cuts: ");
     
  bool goodEvent = true;
  bool passedJetVeto = true;

  unsigned int nJets = numberOfJets(i_hyp);
  if (nJets>0) {
    goodEvent = false;
    passedJetVeto = false;
  }
  int countmus = numberOfSoftMuons(i_hyp,true);
  int nExtraVetoMuons = numberOfSoftMuons(i_hyp,false);
  if (nExtraVetoMuons) goodEvent = false;
  if (numberOfExtraLeptons(i_hyp,10)) goodEvent = false;

  bool passedLTFinalRequirements = true;
  bool passedLLFinalRequirements = true;
  bool passedLTElFakableRequirements = true;
  bool passedLLElFakableRequirements = true;

  // Muon quality cuts, including isolation
  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13){
    if ( !goodMuonIsolated(cms2.hyp_lt_index()[i_hyp]) ) passedLTFinalRequirements = false;
    passedLTElFakableRequirements = false;
  }
  if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13){
    if ( !goodMuonIsolated(cms2.hyp_ll_index()[i_hyp]) ) passedLLFinalRequirements = false;
    passedLLElFakableRequirements = false;
  }  
  // Electron quality cuts, including isolation
  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11){
    if ( !goodElectronIsolated(cms2.hyp_lt_index()[i_hyp]) ) passedLTFinalRequirements = false;
    if ( !fakableElectron(cms2.hyp_lt_index()[i_hyp]) ) passedLTElFakableRequirements = false;
  }
  if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11){
    if ( !goodElectronIsolated(cms2.hyp_ll_index()[i_hyp]) ) passedLLFinalRequirements = false;
    if ( !fakableElectron(cms2.hyp_ll_index()[i_hyp]) ) passedLLElFakableRequirements = false;
  }
  if (goodEvent && dataset )
    getIsolationSidebandsAfterSelections(i_hyp, weight, dataset, passedLTFinalRequirements && passedLLFinalRequirements);
  
  if (goodEvent) {
    countFakableObjectsAfterAllSelections(i_hyp, weight, 
					  passedLTElFakableRequirements, passedLLElFakableRequirements, 
					  passedLTFinalRequirements, passedLLFinalRequirements);
  }
     
  if ( !passedLTFinalRequirements || !passedLLFinalRequirements ) return;
  monitor.count(cms2,type,"after previous + lepton id/iso cuts");
     
  // trkjet veto
  // if ( !passTrkJetVeto(i_hyp) ) return;
     
  // find most energetic jets
  find_most_energetic_jets(i_hyp,weight);
  
  // 2D hist for muon tag counting
  hextramuonsvsnjet[type]->Fill(countmus, nJets, weight);
  hextramuonsvsnjet[3]->Fill(countmus, nJets, weight);
     
  if ( passedJetVeto ) {
    // loop over gen particles
    float centralBQuarkEta(100);
    float forwardBQuarkEta(0);
    unsigned int nBQuarks(0);
    for ( unsigned int i=0; i<cms2.genps_id().size(); ++i )
      if ( abs(cms2.genps_id()[i])==5 ) {
	++nBQuarks;
	float abs_quark_eta = fabs(cms2.genps_p4()[i].eta());
	if ( centralBQuarkEta > abs_quark_eta ) centralBQuarkEta = abs_quark_eta;
	if ( forwardBQuarkEta < abs_quark_eta ) forwardBQuarkEta = abs_quark_eta;
      }
    if ( nBQuarks>0 ) hCentralBquarkEtaAfterVeto->Fill(centralBQuarkEta);
    if ( nBQuarks>1 ) hForwardBquarkEtaAfterVeto->Fill(forwardBQuarkEta);
  }
  if ( ! goodEvent ) return;
  
  monitor.count(cms2,type,"after all cuts (including soft and extra lepton)");

  // -------------------------------------------------------------------//
  // If we made it to here, we passed all cuts and we are ready to fill //
  // -------------------------------------------------------------------//

  hypos_total->Fill(type);
  hypos_total->Fill(3);
  hypos_total_weighted->Fill(type,weight);
  hypos_total_weighted->Fill(3,weight);

  // jet count
  hnJet[type]->Fill(cms2.hyp_njets()[i_hyp], weight);
  hnJet[3]->Fill(cms2.hyp_njets()[i_hyp], weight);
     
  // lepton Pt
  if (abs(cms2.hyp_lt_id()[i_hyp]) == 11) helePt[type]->Fill(cms2.hyp_lt_p4()[i_hyp].pt(), weight);
  if (abs(cms2.hyp_ll_id()[i_hyp]) == 11) helePt[type]->Fill(cms2.hyp_ll_p4()[i_hyp].pt(), weight);
  if (abs(cms2.hyp_lt_id()[i_hyp]) == 13) hmuPt[type]->Fill(cms2.hyp_lt_p4()[i_hyp].pt(), weight);
  if (abs(cms2.hyp_ll_id()[i_hyp]) == 13) hmuPt[type]->Fill(cms2.hyp_ll_p4()[i_hyp].pt(), weight);
     hminLepPt[type]->Fill(min(cms2.hyp_ll_p4()[i_hyp].pt(), cms2.hyp_lt_p4()[i_hyp].pt()), weight);
     hmaxLepPt[type]->Fill(max(cms2.hyp_ll_p4()[i_hyp].pt(), cms2.hyp_lt_p4()[i_hyp].pt()), weight );
    
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11) helePt[3]->Fill(cms2.hyp_lt_p4()[i_hyp].pt(), weight);
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 11) helePt[3]->Fill(cms2.hyp_ll_p4()[i_hyp].pt(), weight);
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 13) hmuPt[3]->Fill(cms2.hyp_lt_p4()[i_hyp].pt(), weight);
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 13) hmuPt[3]->Fill(cms2.hyp_ll_p4()[i_hyp].pt(), weight);
     hminLepPt[3]->Fill(min(cms2.hyp_ll_p4()[i_hyp].pt(), cms2.hyp_lt_p4()[i_hyp].pt()), weight);
     hmaxLepPt[3]->Fill(max(cms2.hyp_ll_p4()[i_hyp].pt(), cms2.hyp_lt_p4()[i_hyp].pt()), weight );
    
     // lepton Phi
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11) helePhi[type]->Fill(cms2.hyp_lt_p4()[i_hyp].phi(), weight);
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 11) helePhi[type]->Fill(cms2.hyp_ll_p4()[i_hyp].phi(), weight);
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 13) hmuPhi[type]->Fill(cms2.hyp_lt_p4()[i_hyp].phi(), weight);
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 13) hmuPhi[type]->Fill(cms2.hyp_ll_p4()[i_hyp].phi(), weight);
    
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11) helePhi[3]->Fill(cms2.hyp_lt_p4()[i_hyp].phi(), weight);
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 11) helePhi[3]->Fill(cms2.hyp_ll_p4()[i_hyp].phi(), weight);
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 13) hmuPhi[3]->Fill(cms2.hyp_lt_p4()[i_hyp].phi(), weight);
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 13) hmuPhi[3]->Fill(cms2.hyp_ll_p4()[i_hyp].phi(), weight);
    
     // dilepton mass
     hdilMass[type]->Fill(cms2.hyp_p4()[i_hyp].mass(), weight);
     hdilMass[3]->Fill(cms2.hyp_p4()[i_hyp].mass(), weight);
    
     // delta phi btw leptons
     double dphi = fabs(cms2.hyp_lt_p4()[i_hyp].phi() - cms2.hyp_ll_p4()[i_hyp].phi());
     if (dphi > TMath::Pi()) dphi = TMath::TwoPi() - dphi;
     hdphiLep[type]->Fill(dphi, weight);
     hdphiLep[3]->Fill(dphi, weight);
    
     // dphill vs mll, i.e. the 2d correlation between the previous two variables
     hdphillvsmll[type]->Fill(cms2.hyp_p4()[i_hyp].mass(), dphi, weight);
     hdphillvsmll[3]->Fill(cms2.hyp_p4()[i_hyp].mass(), dphi, weight);
    
     // lepton Eta
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11) heleEta[type]->Fill(cms2.hyp_lt_p4()[i_hyp].eta(), weight);
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 11) heleEta[type]->Fill(cms2.hyp_ll_p4()[i_hyp].eta(), weight);
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 13) hmuEta[type]->Fill(cms2.hyp_lt_p4()[i_hyp].eta(), weight);
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 13) hmuEta[type]->Fill(cms2.hyp_ll_p4()[i_hyp].eta(), weight);
    
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11) heleEta[3]->Fill(cms2.hyp_lt_p4()[i_hyp].eta(), weight);
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 11) heleEta[3]->Fill(cms2.hyp_ll_p4()[i_hyp].eta(), weight);
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 13) hmuEta[3]->Fill(cms2.hyp_lt_p4()[i_hyp].eta(), weight);
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 13) hmuEta[3]->Fill(cms2.hyp_ll_p4()[i_hyp].eta(), weight);


     // dilepton pt
     hdilPt[type]->Fill(cms2.hyp_p4()[i_hyp].pt(), weight);
     hdilPt[3]->Fill(cms2.hyp_p4()[i_hyp].pt(), weight);
    
     // Met and Met phi
     hmet[type]->Fill(metValue(), weight);      
     hmetPhi[type]->Fill(metPhiValue(), weight);      
     hmet[3]->Fill(metValue(), weight);      
     hmetPhi[3]->Fill(metPhiValue(), weight);      
    
     // Met vs dilepton Pt
     hmetVsDilepPt[type]->Fill(metValue(), cms2.hyp_p4()[i_hyp].pt(), weight);
     hmetVsDilepPt[3]->Fill(metValue(), cms2.hyp_p4()[i_hyp].pt(), weight);
    
     // Met over dilepton Pt vs deltaphi btw the two
     double dphi2 = fabs(cms2.hyp_p4()[i_hyp].phi() - metPhiValue());
     if (dphi2 > TMath::Pi()) dphi2 = TMath::TwoPi() - dphi2;
     hmetOverPtVsDphi[type]->Fill(metValue()/cms2.hyp_p4()[i_hyp].pt(), dphi2, weight);
     hmetOverPtVsDphi[3]->Fill(metValue()/cms2.hyp_p4()[i_hyp].pt(), dphi2, weight);
    
     // get a vector of sorted jets, fill jet histograms
     std::vector<LorentzVector> sortedJets = getJets(jetType(), i_hyp, 0, 5.0, true);
     if ( !sortedJets.empty() ) {
	  hptJet1[type]->Fill(sortedJets[0].Pt(), weight);
	  hptJet1[3]->Fill(sortedJets[0].Pt(), weight);
	  hetaJet1[type]->Fill(sortedJets[0].Eta(), weight);
	  hetaJet1[3]->Fill(sortedJets[0].Eta(), weight);
	  if (sortedJets.size() > 1) {
	    hptJet2[type]->Fill(sortedJets[0].Pt(), weight);
	    hptJet2[3]->Fill(sortedJets[0].Pt(), weight);
	    hetaJet2[type]->Fill(sortedJets[0].Eta(), weight);
	    hetaJet2[3]->Fill(sortedJets[0].Eta(), weight);
	  }
	  if (sortedJets.size() > 2) {
	    hptJet3[type]->Fill(sortedJets[0].Pt(), weight);
	    hptJet3[3]->Fill(sortedJets[0].Pt(), weight);
	    hetaJet3[type]->Fill(sortedJets[0].Eta(), weight);
	    hetaJet3[3]->Fill(sortedJets[0].Eta(), weight);
	  }
     }

}//end of void hypo

RooDataSet* MakeNewDataset(const char* name)
{
  RooRealVar set_iso("iso","iso",0.,10.);
  RooRealVar set_event("event","event",0);
  RooRealVar set_run("run","run",0);
  RooRealVar set_lumi("lumi","lumi",0);
  RooRealVar set_weight("weight","weight",0);
  RooCategory set_selected("selected","Passed final WW selection requirements");
  set_selected.defineType("true",1);
  set_selected.defineType("false",0);

  RooCategory set_hyp_type("hyp_type","Hypothesis type");
  set_hyp_type.defineType(HypothesisTypeName(MM),MM);
  set_hyp_type.defineType(HypothesisTypeName(EM),EM);
  set_hyp_type.defineType(HypothesisTypeName(EE),EE);
  
  RooCategory set_fake_type("fake_type","Define type of lepton for which isolation is extracted");
  set_fake_type.defineType("electron",0);
  set_fake_type.defineType("muon",1);

  RooCategory set_sample_type("sample_type","Sample type");
  set_sample_type.defineType("data_relaxed_iso",0);  // full sample with final selection 
  set_sample_type.defineType("control_sample_signal_iso",1);

  RooDataSet* dataset = new RooDataSet(name, "N-1 dataset",
				       RooArgSet(set_event,set_run,set_lumi,
						 set_iso,set_selected,set_weight,
						 set_hyp_type,set_fake_type,set_sample_type),
				       RooFit::WeightVar(set_weight) );
  return dataset;
}

void AddIsoSignalControlSample( int i_hyp, double kFactor, RooDataSet* dataset) {
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
    if ( goodElectronIsolated(cms2.hyp_lt_index()[i_hyp]) &&
	 goodElectronWithoutIsolation(cms2.hyp_ll_index()[i_hyp]) ){
      set.setRealValue("iso",ww_elIsoVal(cms2.hyp_ll_index()[i_hyp]));
      dataset->add(set,weight);
    }
    if ( goodElectronIsolated(cms2.hyp_ll_index()[i_hyp]) &&
	 goodElectronWithoutIsolation(cms2.hyp_lt_index()[i_hyp]) ){
      set.setRealValue("iso",ww_elIsoVal(cms2.hyp_lt_index()[i_hyp]));
      dataset->add(set,weight);
    }
  } else {
    set.setCatLabel("hyp_type","mm");
    set.setCatLabel("fake_type","muon");
    if ( goodMuonIsolated(cms2.hyp_lt_index()[i_hyp]) &&
	 goodMuonWithoutIsolation(cms2.hyp_ll_index()[i_hyp]) ){
      set.setRealValue("iso",ww_muIsoVal(cms2.hyp_ll_index()[i_hyp]));
      dataset->add(set,weight);
    }
    if ( goodMuonIsolated(cms2.hyp_ll_index()[i_hyp]) &&
	 goodMuonWithoutIsolation(cms2.hyp_lt_index()[i_hyp]) ){
      set.setRealValue("iso",ww_muIsoVal(cms2.hyp_lt_index()[i_hyp]));
      dataset->add(set,weight);
    }
  }
}

RooDataSet* ScanChain( TChain* chain, 
		       enum Sample sample, 
		       double integratedLumi, // in unit of pb^-1
		       double xsec,
		       bool identifyEvents, 
		       bool qcdBackground) 
{
  // chain->SetParallelUnzip(kTRUE);
  // gErrorIgnoreLevel = 3000; // suppress warnings about missing dictionaries 
  unsigned int nEventsChain = chain->GetEntries();  // number of entries in chain --> number of events from all files
  gErrorIgnoreLevel = -1;
  unsigned int nEventsTotal = 0;

 // declare and create array of histograms
  const char *prefix = SampleName(sample);
  RooDataSet* dataset = MakeNewDataset(prefix);
  
  hypos_total          = new TH1F(Form("%s_hypos_total",prefix),"Total number of hypothesis counts",4,0,4);
  hypos_total_weighted = new TH1F(Form("%s_hypos_total_weighted",prefix),"Total number of hypotheses (weighted)",4,0,4);
  hypos_total_weighted->Sumw2();

  for (unsigned int i=0; i<4; ++i){
    hypos_total->GetXaxis()->SetBinLabel(i+1,HypothesisTypeName(i));
    hypos_total_weighted->GetXaxis()->SetBinLabel(i+1,HypothesisTypeName(i));
  }
  
  const Double_t ptbins[4] = {20,30,80,200};
  const Double_t ptbins_fakerate[4] = {10,20,60,150};
  const Double_t etabins_fakerate[4] = {0,1.479,2.5};
  
  for (unsigned int i=0; i<4; i++) {

    hnJet[i]      = new TH1F(Form("%s_hnJet_%s",     prefix,HypothesisTypeName(i)), "Number of jets after all cuts" , 5,0.,5.);	
    helePt[i]     = new TH1F(Form("%s_helePt_%s",    prefix,HypothesisTypeName(i)), "Electron Pt after all cuts", 150,0.,150.);
    hmuPt[i]      = new TH1F(Form("%s_hmuPt_%s",     prefix,HypothesisTypeName(i)), "Muon Pt after all cuts", 150,0.,150.);
    hminLepPt[i]  = new TH1F(Form("%s_hminLepPt_%s", prefix,HypothesisTypeName(i)), "Minimum lepton Pt after all cuts", 150,0.,150.);
    hmaxLepPt[i]  = new TH1F(Form("%s_hmaxLepPt_%s", prefix,HypothesisTypeName(i)), "Maximum lepton Pt after all cuts", 150,0.,150.);
    helePhi[i]    = new TH1F(Form("%s_helePhi_%s",   prefix,HypothesisTypeName(i)), "Electron phi after all cuts", 64, -3.2, 3.2);
    hmuPhi[i]     = new TH1F(Form("%s_hmuPhi_%s",    prefix,HypothesisTypeName(i)), "Muon phi after all cuts", 64, -3.2, 3.2);
    hdphiLep[i]   = new TH1F(Form("%s_hdphiLep_%s",  prefix,HypothesisTypeName(i)), "Delta Phi of the two leptons after all cuts", 64, 0, 3.2);
    heleEta[i]    = new TH1F(Form("%s_heleEta_%s",   prefix,HypothesisTypeName(i)), "Electron eta after all cuts", 60, -3., 3.);
    hmuEta[i]     = new TH1F(Form("%s_hmuEta_%s",    prefix,HypothesisTypeName(i)), "Muon eta after all cuts", 60, -3., 3.);
    hdilMass[i]   = new TH1F(Form("%s_hdilMass_%s",  prefix,HypothesisTypeName(i)), "Di-lepton mass after all cuts", 300, 0., 300.);
    hdilPt[i]     = new TH1F(Form("%s_hdilPt_%s",    prefix,HypothesisTypeName(i)), "Di-lepton pt", 300, 0., 300.);
    hmet[i]       = new TH1F(Form("%s_hmet_%s",      prefix,HypothesisTypeName(i)), "MET after all cuts", 100,0.,200.);
    hmetPhi[i]    = new TH1F(Form("%s_hmetPhi_%s",   prefix,HypothesisTypeName(i)), "MET phi after all cuts", 64, -3.2, 3.2);
    hptJet1[i]    = new TH1F(Form("%s_hptJet1_%s",   prefix,HypothesisTypeName(i)), "Leading jet Pt after all cuts", 300, 0., 300.);
    hptJet2[i]    = new TH1F(Form("%s_hptJet2_%s",   prefix,HypothesisTypeName(i)), "Second jet Pt after all cuts", 300, 0., 300.);
    hptJet3[i]    = new TH1F(Form("%s_hptJet3_%s",   prefix,HypothesisTypeName(i)), "Third jet Pt after all cuts", 300, 0., 300.);
    hetaJet1[i]   = new TH1F(Form("%s_hetaJet1_%s",  prefix,HypothesisTypeName(i)), "Leading jet Eta after all cuts", 50, -5., 5.);
    hetaJet2[i]   = new TH1F(Form("%s_hetaJet2_%s",  prefix,HypothesisTypeName(i)), "Second jet Eta after all cuts", 50, -5., 5.);
    hetaJet3[i]   = new TH1F(Form("%s_hetaJet3_%s",  prefix,HypothesisTypeName(i)), "Third jet Eta after all cuts", 50, -5., 5.);
    hElRelIso[i]  = new TH1F(Form("%s_hElRelIso_%s", prefix,HypothesisTypeName(i)), "Electron relative isolation after all cuts but electron iso",200, 0., 2.);
    hMuRelIso[i]  = new TH1F(Form("%s_hMuRelIso_%s", prefix,HypothesisTypeName(i)), "Muon relative isolation after all cuts but muon iso", 200, 0., 2.);
    hmetVsDilepPt[i]    = new TH2F(Form("%s_hmetVsDilepPt_%s",   prefix,HypothesisTypeName(i)), "MET vs di-lepton Pt after all cuts", 100,0.,200.,100,0.,200.);
    hmetVsDilepPt[i]->Sumw2();
    hmetOverPtVsDphi[i] = new TH2F(Form("%s_hmetOverPtVsDphi_%s",prefix,HypothesisTypeName(i)), "MET/Pt vs di-lepton phi after all cuts", 100,0.,3.,32,0., 3.2);
    hmetOverPtVsDphi[i]->Sumw2();
    hdphillvsmll[i]     = new TH2F(Form("%s_dphillvsmll_%s",     prefix,HypothesisTypeName(i)), "dPhi of the two leptons vs di-lepton mass", 100,10.,210.,32,0.,3.2);
    hdphillvsmll[i]->Sumw2();
    hextramuonsvsnjet[i]= new TH2F(Form("%s_extramuonsvsnjet_%s",prefix,HypothesisTypeName(i)), "Number of extra muon vs number of jets", 10,0.0,10.0,10,0.0,10.0);
    hextramuonsvsnjet[i]->Sumw2();
    helFRfakable[i]     = new TH2F(Form("%s_helFRfakable_%s",    prefix,HypothesisTypeName(i)), "FR study: rate of fakable objects", 3,ptbins,2,0,3.0);
    helFRfakable[i]->Sumw2();
    // fakable object by the fakerate.cc Revision 1.15
    helFRfakable_fakerate[i] = new TH2F(Form("%s_helFRfakable_fakerate_%s", prefix,HypothesisTypeName(i)), "FR study: rate of fakable objects",2,etabins_fakerate,3,ptbins_fakerate);
    helFRfakable_fakerate[i]->Sumw2();

    hnJet[i]->Sumw2();
    helePt[i]->Sumw2();
    hmuPt[i]->Sumw2();
    hminLepPt[i]->Sumw2();
    hmaxLepPt[i]->Sumw2();
    helePhi[i]->Sumw2();
    hmuPhi[i]->Sumw2();
    hdphiLep[i]->Sumw2();
    heleEta[i]->Sumw2();
    hmuEta[i]->Sumw2();
    hdilMass[i]->Sumw2();
    hdilPt[i]->Sumw2();
    hmet[i]->Sumw2();
    hmetPhi[i]->Sumw2();
    hptJet1[i]->Sumw2();
    hptJet2[i]->Sumw2();
    hptJet3[i]->Sumw2();
    hetaJet1[i]->Sumw2();
    hetaJet2[i]->Sumw2();
    hetaJet3[i]->Sumw2();
    hElRelIso[i]->Sumw2(); 
    hMuRelIso[i]->Sumw2(); 

  }
  
  helTrkIsoPassId = new TH1F(Form("%s_helTrkIsoPassId",prefix),        Form("%s - electron trk isolation passed robust el id",prefix),  100, 0., 20.);
  helTrkIsoPassId->Sumw2();
  helTrkIsoFailId = new TH1F(Form("%s_helTrkIsoFailId",prefix),        Form("%s - electron trk isolation failed robust el id",prefix),  100, 0., 20.);
  helTrkIsoFailId->Sumw2();
  helTrkIsoNoId   = new TH1F(Form("%s_helTrkIsoNoId",prefix),          Form("%s - electron trk isolation without el id",prefix), 100, 0., 20.);
  helTrkIsoNoId->Sumw2();
  helTrkPatIsoPassId  = new TH1F(Form("%s_helTrkPatIsoPassId",prefix), Form("%s - electron trk PAT isolation passed robust el id",prefix), 100, 0., 20.);
  helTrkPatIsoPassId->Sumw2();
  helTrkPatIsoFailId  = new TH1F(Form("%s_helTrkPatIsoFailId",prefix), Form("%s - electron trk PAT isolation failed robust el id",prefix), 100, 0., 20.);
  helTrkPatIsoFailId->Sumw2();
  helTrkPatIsoNoId    = new TH1F(Form("%s_helTrkPatIsoNoId",prefix),   Form("%s - electron trk PAT isolation without el id",prefix), 100, 0., 20.);
  helTrkPatIsoNoId->Sumw2();
  
  helEcalJuraIsoPassId = new TH1F(Form("%s_helEcalJuraIsoPassId",prefix), Form("%s - electron ecal jurassic isolation based on basic clusters passed robust el id",prefix), 100, 0., 20.);
  helEcalJuraIsoPassId->Sumw2();
  helEcalJuraIsoFailId = new TH1F(Form("%s_helEcalJuraIsoFailId",prefix), Form("%s - electron ecal jurassic isolation based on basic clusters failed robust el id",prefix), 100, 0., 20.);
  helEcalJuraIsoFailId->Sumw2();
  helEcalJuraIsoNoId   = new TH1F(Form("%s_helEcalJuraIsoNoId",prefix), Form("%s - electron ecal jurassic isolation based on basic clusters without el id",prefix), 100, 0., 20.);
  helEcalJuraIsoNoId->Sumw2();
  helEcalPatIsoPassId  = new TH1F(Form("%s_helEcalPatIsoPassId",prefix),Form("%s - electron ecal jurassic isolation based on rechits pass robust el id",prefix), 100, 0., 20.);
  helEcalPatIsoPassId->Sumw2();
  helEcalPatIsoFailId  = new TH1F(Form("%s_helEcalPatIsoFailId",prefix),Form("%s - electron ecal jurassic isolation based on rechits failed robust el id",prefix), 100, 0., 20.);
  helEcalPatIsoFailId->Sumw2();
  helEcalPatIsoNoId    = new TH1F(Form("%s_helEcalPatIsoNoId",prefix),  Form("%s - electron ecal jurassic isolation based on rechits without el id",prefix), 100, 0., 20.);
  helEcalPatIsoNoId->Sumw2();

  helHcalConeIsoPassId = new TH1F(Form("%s_helHcalConeIsoPassId",prefix), Form("%s - electron hcal cone isolation based on basic clusters passed robust el id",prefix), 100, 0., 20.);
  helHcalConeIsoPassId->Sumw2();
  helHcalConeIsoFailId = new TH1F(Form("%s_helHcalConeIsoFailId",prefix), Form("%s - electron hcal cone isolation based on basic clusters failed robust el id",prefix), 100, 0., 20.);
  helHcalConeIsoFailId->Sumw2();
  helHcalConeIsoNoId   = new TH1F(Form("%s_helHcalConeIsoNoId",prefix),   Form("%s - electron hcal cone isolation based on basic clusters without el id",prefix), 100, 0., 20.);
  helHcalConeIsoNoId->Sumw2();
  helHcalPatIsoPassId  = new TH1F(Form("%s_helHcalPatIsoPassId",prefix),  Form("%s - electron hcal cone isolation based on rechits pass robust el id",prefix), 100, 0., 20.);
  helHcalPatIsoPassId->Sumw2();
  helHcalPatIsoFailId  = new TH1F(Form("%s_helHcalPatIsoFailId",prefix),  Form("%s - electron hcal cone isolation based on rechits failed robust el id",prefix), 100, 0., 20.);
  helHcalPatIsoFailId->Sumw2();
  helHcalPatIsoNoId    = new TH1F(Form("%s_helHcalPatIsoNoId",prefix),	  Form("%s - electron hcal cone isolation based on rechits without el id",prefix), 100, 0., 20.);
  helHcalPatIsoNoId->Sumw2();
  
  helRelIsoPassId = new TH1F(Form("%s_helRelIsoPassId",prefix),  Form("%s - electron relative iso (trk+ecal+hcal) CMS2 with weights 1,1,1 passed robust el id",prefix), 120, 0., 1.2);
  helRelIsoPassId->Sumw2();
  helRelIsoFailId = new TH1F(Form("%s_helRelIsoFailId",prefix),  Form("%s - electron relative iso (trk+ecal+hcal) CMS2 with weights 1,1,1 failed robust el id",prefix), 120, 0., 1.2);
  helRelIsoFailId->Sumw2();
  helRelIsoNoId   = new TH1F(Form("%s_helRelIsoNoId",prefix),    Form("%s - electron relative iso (trk+ecal+hcal) CMS2 with weights 1,1,1 without el id",prefix), 120, 0., 1.2);
  helRelIsoNoId->Sumw2();
  helRelPatIsoPassId  = new TH1F(Form("%s_helRelPatIsoPassId",prefix), Form("%s - electron relative iso (trk+ecal+hcal) PAT with weights 1,1,1 passed robust el id",prefix), 120, 0., 1.2);
  helRelPatIsoPassId->Sumw2();
  helRelPatIsoFailId  = new TH1F(Form("%s_helRelPatIsoFailId",prefix), Form("%s - electron relative iso (trk+ecal+hcal) PAT with weights 1,1,1 failed robust el id",prefix), 120, 0., 1.2);
  helRelPatIsoFailId->Sumw2();
  helRelPatIsoNoId    = new TH1F(Form("%s_helRelPatIsoNoId",prefix),   Form("%s - electron relative iso (trk+ecal+hcal) PAT with weights 1,1,1 without el id",prefix), 120, 0., 1.2);
  helRelPatIsoNoId->Sumw2();

  hmaxJPTEt = new TH1F(Form("%s_hmaxJPTEt",prefix),               Form("%s - most energetic jet Et (JPT)",prefix), 200, 0., 200);
  hmaxJPTEt->Sumw2();
  hmaxBtagVsJPTEt = new TH2F(Form("%s_hmaxBtagVsJPTEt",prefix),   Form("%s - most energetic jet Et (JPT) vs b-tagger output",prefix), 20, 0., 100, 10, 0, 1);
  hmaxBtagVsJPTEt->Sumw2();
  hmaxCaloJetEt = new TH1F(Form("%s_hmaxCaloJetEt",prefix),       Form("%s - most energetic jet Et (CaloJet)",prefix), 200, 0., 200);
  hmaxCaloJetEt->Sumw2();
  hmaxTrkJetEt = new TH1F(Form("%s_hmaxTrkJetEt",prefix),         Form("%s - most energetic jet Et (TrkJet)",prefix), 200, 0., 200);
  hmaxTrkJetEt->Sumw2();
  hmaxCaloTrkJetEt = new TH1F(Form("%s_hmaxCaloTrkJetEt",prefix), Form("%s - most energetic jet Et (average of Calo + Trk Jets)",prefix), 200, 0., 200);
  hmaxCaloTrkJetEt->Sumw2();
  hmaxCaloTrkJet2Et = new TH1F(Form("%s_hmaxCaloTrkJet2Et",prefix), Form("%s - most energetic jet Et (Max Calo and Trk Jets)",prefix), 200, 0., 200);
  hmaxCaloTrkJet2Et->Sumw2();
  hmaxGenJetEt = new TH1F(Form("%s_hmaxGenJetEt",prefix),         Form("%s - most energetic jet Et (GenJet)",prefix), 200, 0., 200);
  hmaxGenJetEt->Sumw2();
  hCentralBquarkEtaAfterVeto = new TH1F(Form("%s_centralBQuarkEtaAfterVeto",prefix), Form("%s - central b quark eta distribution after jet veto",prefix), 20, 0, 10);
  hForwardBquarkEtaAfterVeto = new TH1F(Form("%s_forwardBQuarkEtaAfterVeto",prefix), Form("%s - forward b quark eta distribution after jet veto",prefix), 20, 0, 10);

  if (qcdBackground) {
    hFakableRateSingleElectron = new TH2F(Form("%s_hFakableRateSingleElectron",prefix), "FR study: rate of fakable objects", 3,ptbins,2,0,3.0);
    hFakableRateSingleElectron->Sumw2();
    hFinalRateSingleElectron   = new TH2F(Form("%s_hFinalRateSingleElectron",  prefix), "FR study: rate of final objects", 3,ptbins,2,0,3.0);
    hFinalRateSingleElectron->Sumw2();
    hIsoSingleMuon             = new TH1F(Form("%s_hIsoSingleMuon",prefix),             "Muon isolation distribution for goodMuonWithoutIsolation", 100,0,10);
    hIsoSingleMuon->Sumw2();
    hIsoSingleElectron         = new TH1F(Form("%s_hIsoSingleElectron",prefix),         "Electron isolation distribution for goodElectronWithoutIsolation", 100,0,10);
    hIsoSingleElectron->Sumw2();
  }

  // clear list of duplicates
  already_seen.clear();
  int duplicates_total_n = 0;
  double duplicates_total_weight = 0;
  int nFailedIdentification = 0;
  int nFilteredOut = 0;

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
       assert(tree);

       cms2.Init(tree);  // set branch addresses for TTree tree
       
       TStopwatch t;
       //Event Loop
       unsigned int nEvents = tree->GetEntries();
       for( unsigned int event = 0; event < nEvents; ++event) {
	 cms2.GetEntry(event);  // get entries for Event number event from branches of TTree tree
	 
	 // Reset the kFactor if the xsec argument is specified
	 double kFactor = integratedLumi/1000.0;
	 if( xsec > 0) kFactor = (integratedLumi/1000.0) * xsec / (cms2.evt_xsec_excl()*cms2.evt_kfactor());
       
	 ++nEventsTotal;
	 if (qcdBackground) {
	   // get fake rates
	   extractFakeRateSingleLepton();
	   // isolation
	   extractIsoSingleLepton();
	 }
	 
	 if (cms2.trks_d0().size() == 0) continue;  // needed to get rid of back Monte Carlo events in CMSSW_2_X analysis
	 if (cms2.hyp_type().size() == 0) continue; // skip events without hypothesis
	 EventIdentifier id = { cms2.evt_run(), cms2.evt_event(), cms2.evt_lumiBlock(), cms2.trks_d0()[0], 
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
	 
	 if ( identifyEvents ){
	   // check if we know what we are looking at
	   if ( ! isIdentified(sample) ) nFailedIdentification++;
	   
	   // filter by process
	   if ( ! filterByProcess(sample) ) {
	     nFilteredOut++;
	     continue;
	   }
	 }
	 
	 // loop over hypothesis candidates
	 unsigned int nHyps = cms2.hyp_type().size();
	 for( unsigned int i_hyp = 0; i_hyp < nHyps; ++i_hyp ) {
	   if(cms2.hyp_p4().at(i_hyp).mass2() < 0 ) break;
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
  // monitor.printEvents(3);
  if ( nEventsChain != nEventsTotal ) {
       printf("ERROR: number of events from files (%d) is not equal to total number"
	      " of events (%d)\n", nEventsChain, nEventsTotal);
  }
  printf("Total number of skipped events due to bad identification: %d (%0.0f %%)\n",   
	   nFailedIdentification, nFailedIdentification*100.0/(nEventsChain+1e-5));
  printf("Total number of filtered out events: %d (%0.0f %%)\n",   
	   nFilteredOut, nFilteredOut*100.0/(nEventsChain+1e-5));
  printf("Total candidate count (%s %s %s %s): %.0f %.0f %.0f %0.f.\n",
	 HypothesisTypeName(0), HypothesisTypeName(1), 
	 HypothesisTypeName(2), HypothesisTypeName(3),
	 hypos_total->GetBinContent(1), hypos_total->GetBinContent(2), 
	 hypos_total->GetBinContent(3), hypos_total->GetBinContent(4));
  printf("Total weighted candidate yeild (%s %s %s %s): %f %f %f %f\n",   
	 HypothesisTypeName(0), HypothesisTypeName(1), 
	 HypothesisTypeName(2), HypothesisTypeName(3),
	 hypos_total_weighted->GetBinContent(1), hypos_total_weighted->GetBinContent(2), 
	 hypos_total_weighted->GetBinContent(3), hypos_total_weighted->GetBinContent(4));
  
  return dataset;
}

bool EventIdentifier::operator < (const EventIdentifier &other) const
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

bool EventIdentifier::operator == (const EventIdentifier &other) const
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

// filter events by process
bool filterByProcess( enum Sample sample ) {
  switch (sample) {
  case DYee: 
    return isDYee();
  case DYmm:
    return isDYmm();
  case DYtt:
    return isDYtt();
  case WW:
    return isWW();
  case WZ:
    return isWZ();
  case ZZ:
    return isZZ();
  default:
    return true;
  }
}
  
bool isIdentified( enum Sample sample ) {
  switch (sample) {
  case DYee:
  case DYmm:
  case DYtt:
    return getDrellYanType()!=999;
  case WW:
  case WZ:
  case ZZ:
    return getVVType()!=999;
  default:
    return true;
  }
}

void ProcessSample( std::string file_pattern, 
		    Sample sample, 
		    double integratedLumi,
		    double xsec,
		    RooDataSet* output_dataset, 
		    Color_t color, 
		    bool identifyEvents,
		    bool qcdBackground)
{
  std::vector<string> vec;
  vec.push_back(file_pattern);
  ProcessSample(vec,sample,integratedLumi,xsec,output_dataset,color,identifyEvents,qcdBackground);
}

void ProcessSample( std::vector<std::string> file_patterns, 
		    Sample sample, 
		    double integratedLumi,
		    double xsec,
		    RooDataSet* output_dataset, 
		    Color_t color, 
		    bool identifyEvents,
		    bool qcdBackground)
{
  TChain *tchain = new TChain("Events");
  for ( std::vector<std::string>::const_iterator pattern = file_patterns.begin();
	pattern != file_patterns.end(); ++pattern )
    tchain->Add(pattern->c_str());
  if ( gSystem->Getenv("SkimSamples") ){
    std::cout << "Skimming " << SampleName(sample) << " with skim name " << 
      gSystem->Getenv("SkimName") << " .." << std::endl;
    SkimChain(tchain);
  } else {
    std::cout << "Processing " << SampleName(sample) << ".." << std::endl;
    RooDataSet* data = ScanChain(tchain,sample,integratedLumi,xsec,identifyEvents,qcdBackground);
    if( data ){
      if ( output_dataset )
	output_dataset->append(*data);
      else
	output_dataset=data;
    }
    
    const char* sampleName = SampleName(sample);
    TRegexp reg(sampleName, kFALSE);
    
    TList* list = gDirectory->GetList() ;
    if (!list) {
      cout << "Failed to set color for " << sampleName << endl;
      return;
    }
    TIterator* iter = list->MakeIterator();

    while (TObject* obj = iter->Next()) {
      if (! obj->InheritsFrom(TH1::Class())) continue;
      
      TString name = obj->GetName();

      if (TString(sampleName).MaybeRegexp()) {
	if (TString(obj->GetName()).Index(reg) < 0 ) continue;
      } else if (! name.BeginsWith(sampleName)) continue;
      
      ((TH1*)obj)->SetFillColor(color);
      ((TH1*)obj)->SetLineColor(color);
      ((TH1*)obj)->SetMarkerColor(color);
    }
  }
}

void SkimChain(TChain* chain){
  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  TPRegexp fileNameMatchPattern("(.*?)/([^/]*)$");
  TPRegexp dataDirMatchPattern("^data");
  if ( gSystem->Getenv("DataDir") ) dataDirMatchPattern = TPRegexp(gSystem->Getenv("DataDir"));
  unsigned int nEventsTotal = 0;
  unsigned int nEventsSelected = 0;
  while (TChainElement *currentFile = (TChainElement*)fileIter.Next()) {
    cout << "file: " << currentFile->GetTitle() << endl;
    TObjArray* matches = fileNameMatchPattern.MatchS(currentFile->GetTitle());
    assert(matches);
    cout << "matches->GetLast(): " << matches->GetLast() << endl;
    assert(matches && matches->GetLast()==2);
    TString inputFileName(currentFile->GetTitle());
    TString fileName(((TObjString*)matches->At(2))->GetString());
    TString outputDirectory(((TObjString*)matches->At(1))->GetString());
    if ( gSystem->Getenv("SkimDir") ) 
      dataDirMatchPattern.Substitute(outputDirectory,gSystem->Getenv("SkimDir"));
    assert(gSystem->Getenv("SkimName"));
    TString skimname = gSystem->Getenv("SkimName");
    outputDirectory += "/" + skimname;
    // make output directory if it doesn't exist yet
    if ( gSystem->AccessPathName(outputDirectory.Data()) ){
      gSystem->mkdir(outputDirectory.Data(),true);
      assert( !gSystem->AccessPathName(outputDirectory.Data()) );
    }
    TString outputFileName(outputDirectory+"/"+fileName);
    cout << "Skimming " << inputFileName << " -> " << outputFileName << endl;

    TFile *output = TFile::Open(outputFileName.Data(), "RECREATE");
    assert(output);
    TFile *input = TFile::Open(inputFileName.Data());
    assert(input);
    TTree *tree = (TTree*)input->Get("Events");
    tree->SetBranchStatus("EventAuxiliary",0);
    TTree *newtree = tree->CloneTree(0);
    newtree->SetDirectory(output);
    
    cms2.Init(newtree);
    cms2.Init(tree);
    
    // Event Loop
    const unsigned int nEvents = tree->GetEntries();
    int i_permille_old = 0;
    for (unsigned int event = 0; event < nEvents; ++event, ++nEventsTotal) {
      int i_permille = (int)floor(10000 * event / float(nEvents));
      if (i_permille != i_permille_old) {
	// xterm magic from L. Vacavant and A. Cerri
	if (isatty(1)) {
	  printf("\015\033[32m ---> \033[1m\033[31m%5.2f%%"
		 "\033[0m\033[32m <---\033[0m\015", i_permille/100.);
	  fflush(stdout);
	}
	i_permille_old = i_permille;
      }
      cms2.GetEntry(event);
      //set condition to skip event
      if ( not passedSkimSelection() ) continue;
      
      ++nEventsSelected;
      cms2.LoadAllBranches();
      // fill the new tree
      newtree->Fill();
    }
    output->cd();
    newtree->Write();
    output->Close();
  }
  cout << Form("Processed events: %u, \tselected: %u\n",nEventsTotal,nEventsSelected) << endl;
}

bool passedSkimSelection()
{
  unsigned int nHyps = cms2.hyp_type().size();
  for( unsigned int i_hyp = 0; i_hyp < nHyps; ++i_hyp ) {
    // min lepton Pt
    if ( cms2.hyp_lt_p4().at(i_hyp).pt() < 10 || cms2.hyp_ll_p4().at(i_hyp).pt() < 10 ) continue;
    // charge
    if ( cms2.hyp_lt_id()[i_hyp] * cms2.hyp_ll_id()[i_hyp] > 0 ) continue;
    // met
    if ( cms2.evt35X_tcmet() < 20 &&
	 cms2.evt_tcmet() < 20 &&
	 cms2.evt_pfmet() < 20 ) continue;
    // id & iso
    if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && !goodMuonIsolated(cms2.hyp_lt_index()[i_hyp]) ) continue;
    if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && !goodMuonIsolated(cms2.hyp_ll_index()[i_hyp]) ) continue;
    if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11 && !goodElectronIsolated(cms2.hyp_lt_index()[i_hyp]) ) continue;
    if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11 && !goodElectronIsolated(cms2.hyp_ll_index()[i_hyp]) ) continue;
    
    return true;
  }
  return false;
}
