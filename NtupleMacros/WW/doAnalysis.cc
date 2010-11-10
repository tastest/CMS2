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
#include "../Tools/goodrun.cc"
using namespace std;

#ifndef __CINT__
#include "CORE/CMS2.h"
#include "CORE/electronSelections.h"
#include "CORE/muonSelections.h"
#include "CORE/jetSelections.h"
#endif

enum jetregion { HCAL, HF, ALLJET};

enum hyp_selection {
  PASS_ZSEL,
  PASS_MET,
  PASS_JETVETO,  
  PASS_LT_FINAL,
  PASS_LL_FINAL,
  PASS_SOFTMUVETO,
  PASS_EXTRALEPTONVETO,
  PASS_TOPVETO,
  PASS_1BJET,
  PASS_1JET
};

const cuts_t pass_all = (1<< PASS_ZSEL) | (1<<PASS_MET) | (1<<PASS_JETVETO) | (1<<PASS_LT_FINAL) 
  | (1<<PASS_LL_FINAL) | (1<<PASS_SOFTMUVETO) | (1<<PASS_EXTRALEPTONVETO) | (1<<PASS_TOPVETO);

bool applyJEC = true;

std::vector<std::string> jetcorr_filenames_jpt;
FactorizedJetCorrector *jet_corrector_jpt;

std::vector<std::string> jetcorr_filenames_pf;
FactorizedJetCorrector *jet_corrector_pf;

std::vector<std::string> jetcorr_filenames_calo;
FactorizedJetCorrector *jet_corrector_calo;

std::vector<std::string> jetcorr_filenames_trk;
FactorizedJetCorrector *jet_corrector_trk;


//
// Key analysis method implementation
//

bool goodElectronWithoutIsolation(unsigned int i){
  return ww_elBase(i) && ww_elId(i) && ww_eld0PV(i);
}

bool goodElectronIsolated(unsigned int i){
  bool ptcut = cms2.els_p4().at(i).pt() >= 20.0;
  bool core = ptcut && pass_electronSelection( i, electronSelection_wwV1);
  bool internal = ww_elBase(i) && ww_elId(i) && ww_eld0PV(i) && ww_elIso(i);
  assert(core==internal);
  return core;
}

bool fakableElectron(unsigned int i){
  bool ptcut = cms2.els_p4().at(i).pt() >= 20.0;
  // extrapolate in partial id, iso and d0
  return ptcut && pass_electronSelection( i, electronSelectionFO_el_wwV1_v2);
  // extrapolate in id
  // return ww_elBase(i) && ww_eld0(i) && ww_elIso(i);
}

bool goodMuonWithoutIsolation(unsigned int i){
  return ww_muBase(i) && ww_mud0PV(i) && ww_muId(i);
}

bool goodMuonIsolated(unsigned int i){
  bool ptcut = cms2.mus_p4().at(i).pt() >= 20.0;
  bool core = ptcut && muonId(i, NominalWWV1);
  bool internal = ww_muBase(i) && ww_mud0PV(i) && ww_muId(i) && ww_muIso(i); 
  assert(core==internal);
  return core;
}

bool fakableMuon(unsigned int i){
  bool ptcut = cms2.mus_p4().at(i).pt() >= 20.0;
  // extrapolate in iso
  // return muonId(i, muonSelectionFO_mu_ww);
  return ptcut && muonId(i, muonSelectionFO_mu_wwV1_iso10);
  // return ww_muBase(i) && ww_muId(i) && ww_muIsoVal(i)<1.0 && fabs(cms2.mus_d0corr()[i]) < 2; 
}

double metValue(){    return cms2.evt_tcmet(); }
double metPhiValue(){ return cms2.evt_tcmetPhi(); }

bool passedMetRequirements(unsigned int i_hyp){
  // if ( cms2.hyp_p4().at(i_hyp).mass()>130 ) return true;
  HypothesisType type = getHypothesisType(cms2.hyp_type()[i_hyp]);
  double pMet = projectedMet(i_hyp);
  // if ( type == EM && cms2.hyp_p4().at(i_hyp).mass()>90 ) return true;
  if ( pMet < 20 ) return false;
  if (type == EE || type == MM) {
    // double dmass = fabs(cms2.hyp_p4()[i_hyp].mass()-91);
    // if ( metValue() < 45 ) return false;
    if ( pMet < 35 ) return false;
    // if ( !metBalance(i_hyp) ) return false;
  }
  return true;
}


WWJetType jetType(){
  return pfJet;
  // return jptJet;
}

std::vector<LorentzVector> getDefaultJets(unsigned int i_hyp, bool btagged=false){
  return getJets(jetType(), i_hyp, 25, 5.0, false, btagged); // V1
}

unsigned int numberOfJets(unsigned int i_hyp){
  return getDefaultJets(i_hyp, false).size(); 
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
  // if ( cms2.els_exp_innerlayers().at(index) > 0 ) return false;
  if ( cms2.els_exp_innerlayers39X().at(index) > 0 ) return false;
  //  int ctfIndex = cms2.els_trkidx().at(index);
  // if ( ctfIndex >=0 && 
  //     cms2.els_charge().at(index)!=cms2.trks_charge().at(ctfIndex) ) return false;
  return true;
}
 
bool ww_eld0(unsigned int index){
  return fabs(cms2.els_d0corr()[index]) < 0.02;
}

bool isGoodVertex(size_t ivtx) {
    if (cms2.vtxs_isFake()[ivtx]) return false;
    if (cms2.vtxs_ndof()[ivtx] < 4.) return false;
    if (cms2.vtxs_position()[ivtx].Rho() > 2.0) return false;
    if (fabs(cms2.vtxs_position()[ivtx].Z()) > 24.0) return false;
    return true;
}

unsigned int nGoodVertex() {
  unsigned int nVtx = 0;
  for ( unsigned int i = 0; i < cms2.vtxs_sumpt().size(); ++i ){
    // if (cms2.vtxs_isFake()[i]) continue;
    if (!isGoodVertex(i)) continue;
    nVtx++;
  }
  return nVtx;
}

double dzPV(const LorentzVector& vtx, const LorentzVector& p4, const LorentzVector& pv){
  return (vtx.z()-pv.z()) - ((vtx.x()-pv.x())*p4.x()+(vtx.y()-pv.y())*p4.y())/p4.pt() * p4.z()/p4.pt();
}

bool ww_eld0PV(unsigned int index){
  if ( cms2.vtxs_sumpt().empty() ) return false;
  double sumPtMax = -1;
  int iMax = -1;
  for ( unsigned int i = 0; i < cms2.vtxs_sumpt().size(); ++i ){
    // if (cms2.vtxs_isFake()[i]) continue;
    if (!isGoodVertex(i)) continue;
    if ( cms2.vtxs_sumpt().at(i) > sumPtMax ){
      iMax = i;
      sumPtMax = cms2.vtxs_sumpt().at(i);
    }
  }
  if (iMax<0) return false;
  double dxyPV = cms2.els_d0()[index]-
    cms2.vtxs_position()[iMax].x()*sin(cms2.els_trk_p4()[index].phi())+
    cms2.vtxs_position()[iMax].y()*cos(cms2.els_trk_p4()[index].phi());
  // double dzPV = cms2.els_z0corr()[index]-cms2.vtxs_position()[iMax].z();
  double dzpv = dzPV(cms2.els_vertex_p4()[index], cms2.els_trk_p4()[index], cms2.vtxs_position()[iMax]);
  return fabs(dxyPV) < 0.02 && fabs(dzpv)<1.0;
}

double ww_elIsoVal(unsigned int index){
  float sum = cms2.els_tkIso().at(index);
  if ( fabs(cms2.els_etaSC()[index]) < 1.479 )
    // if ( fabs(cms2.els_p4().at(index).eta()) < 1.479)
    sum += max(0., (cms2.els_ecalIso().at(index) -1.));
  else 
    sum += cms2.els_ecalIso().at(index);
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
  if (cms2.mus_type().at(index) == 8) return false; // not STA
  return true;
}
bool ww_mud0(unsigned int index){
  return fabs(cms2.mus_d0corr()[index]) < 0.02;
}
bool ww_mud0PV(unsigned int index){
  if ( cms2.vtxs_sumpt().empty() ) return false;
  double sumPtMax = -1;
  int iMax = -1;
  for ( unsigned int i = 0; i < cms2.vtxs_sumpt().size(); ++i ){
    if (cms2.vtxs_isFake()[i]) continue;
    if (!isGoodVertex(i)) continue;
    if ( cms2.vtxs_sumpt().at(i) > sumPtMax ){
      iMax = i;
      sumPtMax = cms2.vtxs_sumpt().at(i);
    }
  }
  if (iMax<0) return false;
  double dxyPV = cms2.mus_d0()[index]-
    cms2.vtxs_position()[iMax].x()*sin(cms2.mus_trk_p4()[index].phi())+
    cms2.vtxs_position()[iMax].y()*cos(cms2.mus_trk_p4()[index].phi());
  // double dzpv = cms2.mus_z0corr()[index]-cms2.vtxs_position()[iMax].z();
  double dzpv = dzPV(cms2.mus_vertex_p4()[index], cms2.mus_trk_p4()[index], cms2.vtxs_position()[iMax]);
  return fabs(dxyPV) < 0.02 && fabs(dzpv)<1.0;
}
bool ww_muId(unsigned int index){
  if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) return false; //glb fit chisq
  if (((cms2.mus_type().at(index)) & (1<<1)) == 0)    return false; // global muon
  if (((cms2.mus_type().at(index)) & (1<<2)) == 0)    return false; // tracker muon
  if (cms2.mus_validHits().at(index) < 11)            return false; // # of tracker hits
  if (cms2.mus_gfit_validSTAHits().at(index)==0 ) return false;
  if (cms2.mus_ptErr().at(index)/cms2.mus_p4().at(index).pt()>0.1) return false;
  if (cms2.trks_valid_pixelhits().at(cms2.mus_trkidx().at(index))==0) return false;
  if (cms2.mus_nmatches().at(index)<2) return false;
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
unsigned int numberOfSoftMuons(int i_hyp, bool nonisolated,
			       const std::vector<LorentzVector>& vetojets)
{
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
    if ( nonisolated && ww_muIsoVal(imu)<0.1 && cms2.mus_p4()[imu].pt()>20 ) continue;
    bool skip = false;
    for ( std::vector<LorentzVector>::const_iterator ijet = vetojets.begin();
	  ijet != vetojets.end(); ++ijet )
      if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(*ijet,cms2.mus_p4()[imu])) < 0.3 ) skip=true;
    if ( skip ) continue;
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
bool passedTrigger(TString trigName) {
  if ( find(cms2.hlt_trigNames().begin(), cms2.hlt_trigNames().end(), trigName)
       == cms2.hlt_trigNames().end() ) return false;
  return cms2.passHLTTrigger(trigName);
}

bool passedTriggerRequirements(HypTypeInNtuples type) {
  if ( passedTrigger("HLT_Mu9") ) return true;
  if ( passedTrigger("HLT_Mu15_v1") ) return true;
  if ( passedTrigger("HLT_Ele10_LW_L1R") ) return true;
  if ( passedTrigger("HLT_Ele15_LW_L1R") ) return true;
  if ( passedTrigger("HLT_Ele15_SW_L1R") ) return true;
  if ( passedTrigger("HLT_Ele15_SW_CaloEleId_L1R") ) return true;
  if ( passedTrigger("HLT_Ele17_SW_CaloEleId_L1R") ) return true;
  if ( passedTrigger("HLT_Ele17_SW_TightEleId_L1R") ) return true;
  if ( passedTrigger("HLT_Ele17_SW_TighterEleIdIsol_L1R_v2") ) return true;
  if ( passedTrigger("HLT_Ele17_SW_TighterEleIdIsol_L1R_v3") ) return true;
  return false;
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
getJets(WWJetType type, int i_hyp, double etThreshold, double etaMax, bool sortJets, bool btag)
{
     std::vector<LorentzVector> jets;
     const double vetoCone = 0.3;
     double jec = 1.0;
     
     switch ( type ){
     case jptJet:
        for ( unsigned int i=0; i < cms2.jpts_p4().size(); ++i) {
	  if(applyJEC)
	    jec = jetCorrection(cms2.jpts_p4()[i], jet_corrector_jpt);
	 if ( cms2.jpts_p4()[i].pt() * jec < etThreshold ) continue;
	 if ( btag && !defaultBTag(type,i) ) continue;
	 if ( TMath::Abs(cms2.jpts_p4()[i].eta()) > etaMax ) continue;
	 if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.jpts_p4()[i])) < vetoCone ||
	      TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.jpts_p4()[i])) < vetoCone ) continue;
	 jets.push_back(cms2.jpts_p4()[i] * jec);
       }
       break;
     case pfJet:
       for ( unsigned int i=0; i < cms2.pfjets_p4().size(); ++i) {
	  if(applyJEC)
	    jec = jetCorrection(cms2.pfjets_p4()[i], jet_corrector_pf);
	  if ( cms2.pfjets_p4()[i].pt() * jec < etThreshold ) continue;
	  if ( btag && !defaultBTag(type,i) ) continue;
	  if ( TMath::Abs(cms2.pfjets_p4()[i].eta()) > etaMax ) continue;
	  if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.pfjets_p4()[i])) < vetoCone ||
	       TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.pfjets_p4()[i])) < vetoCone ) continue;
	  jets.push_back(cms2.pfjets_p4()[i] * jec);
       }
       break;
     case GenJet:
       for ( unsigned int i=0; i < cms2.genjets_p4().size(); ++i) {
	 if ( cms2.genjets_p4()[i].pt() < etThreshold ) continue;
	 if ( btag && !defaultBTag(type,i) ) continue;
	 if ( TMath::Abs(cms2.genjets_p4()[i].eta()) > etaMax ) continue;
	 if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.genjets_p4()[i])) < vetoCone ||
	      TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.genjets_p4()[i])) < vetoCone ) continue;
	 jets.push_back(cms2.genjets_p4()[i]);
       }
       break;
     case CaloJet:
       for ( unsigned int i=0; i < cms2.jets_pat_jet_p4().size(); ++i) {
	 if ( cms2.jets_pat_jet_p4()[i].pt() < etThreshold ) continue; // note that this is already corrected
	 if ( btag && !defaultBTag(type,i) ) continue;
	 if ( TMath::Abs(cms2.jets_pat_jet_p4()[i].eta()) > etaMax ) continue;
	 if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.jets_pat_jet_p4()[i])) < vetoCone ||
	      TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.jets_pat_jet_p4()[i])) < vetoCone ) continue;
	 jets.push_back(cms2.jets_pat_jet_p4()[i]);
       }
       break;
     case TrkJet:
       for ( unsigned int i=0; i < cms2.trkjets_p4().size(); ++i) {
	 if(applyJEC)
	   jec = jetCorrection(cms2.trkjets_p4()[i], jet_corrector_trk);
	 if ( cms2.trkjets_p4()[i].pt() < etThreshold ) continue;
	 if ( btag && !defaultBTag(type,i) ) continue;
	 if ( TMath::Abs(cms2.trkjets_p4()[i].eta()) > etaMax ) continue;
	 if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.trkjets_p4()[i])) < vetoCone ||
	      TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.trkjets_p4()[i])) < vetoCone ) continue;
	 jets.push_back(cms2.trkjets_p4()[i] * jec);
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

double BTag(WWJetType type, unsigned int iJet){
  // find a jet in the jets list
  // that matches the current type jet
  LorentzVector jetP4;
  switch ( type ) {
  case jptJet:
    jetP4 = cms2.jpts_p4().at(iJet);
    break;
  case CaloJet:
    return cms2.jets_trackCountingHighEffBJetTag()[iJet];
    break;
  case pfJet:
    return cms2.pfjets_trackCountingHighEffBJetTag()[iJet];
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
  return cms2.jets_trackCountingHighEffBJetTag().at(refJet);
}

bool defaultBTag(WWJetType type, unsigned int iJet){
  return BTag(type,iJet)>2.1;
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
TH2F* hmuFRfakable[4];     // Fake rate study: rate of fakable objects
TH2F* hFakableRateSingleElectron; // Fake rate study: rate of fakable objects
TH2F* hFakableRateSingleMuon; // Fake rate study: rate of fakable objects
TH2F* hFinalRateSingleElectron;   // Fake rate study: rate of final objects
TH2F* hFinalRateSingleMuon;   // Fake rate study: rate of final objects
TH1F* hIsoSingleMuon;             // isolation background 
TH1F* hIsoSingleElectron;         // isolation background 
TH2F* hmaxBtagVsJPTPt;   // btag vs energy distribution for the most energetic jet
TH2F* htoptagz[4];
TH2F* htoptag[4];

// Histograms without the full selections
TH1F* hmaxJPTPt[4];         // energy distribution for the most energetic jet |eta|<5
TH1F* hmaxCaloJetPt[4];     // energy distribution for the most energetic jet |eta|<5
TH1F* hmaxTrkJetPt[4];      // energy distribution for the most energetic jet |eta|<5
TH1F* hmaxCaloTrkJetPt[4];  // energy distribution for the most energetic jet |eta|<5
TH1F* hmaxCaloTrkJet2Et[4]; // energy distribution for the most energetic jet |eta|<5
TH1F* hmaxGenJetPt[4];      // energy distribution for the most energetic jet |eta|<5
TH1F* hmaxPFJetPt[4];      // energy distribution for the most energetic jet |eta|<5

TH1F* hmaxJPTPtHCal[4];         // energy distribution for the most energetic jet |eta|<3
TH1F* hmaxCaloJetPtHCal[4];     // energy distribution for the most energetic jet |eta|<3
TH1F* hmaxTrkJetPtHCal[4];      // energy distribution for the most energetic jet |eta|<3
TH1F* hmaxCaloTrkJetPtHCal[4];  // energy distribution for the most energetic jet |eta|<3
TH1F* hmaxCaloTrkJet2EtHCal[4]; // energy distribution for the most energetic jet |eta|<3
TH1F* hmaxGenJetPtHCal[4];      // energy distribution for the most energetic jet |eta|<3
TH1F* hmaxPFJetPtHCal[4];      // energy distribution for the most energetic jet |eta|<3

TH1F* hmaxJPTPtHF[4];         // energy distribution for the most energetic jet 3<|eta|<5
TH1F* hmaxCaloJetPtHF[4];     // energy distribution for the most energetic jet 3<|eta|<5
TH1F* hmaxTrkJetPtHF[4];      // energy distribution for the most energetic jet 3<|eta|<5
TH1F* hmaxCaloTrkJetPtHF[4];  // energy distribution for the most energetic jet 3<|eta|<5
TH1F* hmaxCaloTrkJet2EtHF[4]; // energy distribution for the most energetic jet 3<|eta|<5
TH1F* hmaxGenJetPtHF[4];      // energy distribution for the most energetic jet 3<|eta|<5
TH1F* hmaxPFJetPtHF[4];      // energy distribution for the most energetic jet 3<|eta|<5

// validation plots
TH1F* hmetVal[4];             // MET validation after ll selection/jetveto
TH1F* hmetProjVal[4];         // Projected MET after ll selection
TH1F* hmaxPFJetPtVal[4];      // leading PFJet Pt after ll selection
TH1F* hdilMassVal[4];         // diLepton Mass after ll selection

// For the DY Estimation

TH1F* hmetInDYEst[4];         // MET in Z window for the DY Estimation
TH1F* hmetOutDYEst[4];        // MET outside Z window for the DY Estimation
TH1F* hdilMassWithMetDYEst[4];// Dilepton mass with MET requirement for DY estimation
TH1F* hdilMassNoMetDYEst[4];  // Dilepton mass without MET requirement for DY estimation

//
// Not cleaned area
//

TH2F* helFRfakable_fakerate[4];     // Fake rate study: rate of fakable objects using the fakerate.cc
TH2F* hmuFRfakable_fakerate[4];     // Fake rate study: rate of fakable objects using the fakerate.cc
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

TH1F* hCentralBquarkEtaAfterVeto;   
TH1F* hForwardBquarkEtaAfterVeto;   

TH2F* hPFJetResponse[4]; // absolute response to Z pt
TH2F* hPFJetRelResponse[4]; // Relative response to Z pt

TH2F* hPFJetResponseWithZero[4]; // absolute response to Z pt with 0 if no btb jet is found
TH2F* hPFJetRelResponseWithZero[4]; // Relative response to Z pt with 0 if no btb jet is found

TH2F* hPFJetResidualResponse; // absolute response to Z pt
TH1F* hnGoodVertex;      // number of good vertexes after baseline selections

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


void find_most_energetic_jets(int i_hyp, double weight, bool realData, double etaMin, double etaMax)
{
  HypothesisType type = getHypothesisType(cms2.hyp_type()[i_hyp]);
  // fill the jet region
  jetregion jet = ALLJET;
  if( etaMin <= 0.0 && etaMax <= 3.0) jet = HCAL;
  if( etaMin >= 3.0 && etaMax <= 5.0) jet = HF;
  
  double vetoCone = 0.3;
  
  if(!realData)
    {
      double genJetMax(0.);
      find_leading_genjet(i_hyp, etaMin, etaMax, vetoCone, genJetMax);
      switch (jet) {
      case ALLJET :
	hmaxGenJetPt[type]->Fill(genJetMax, weight);
	hmaxGenJetPt[3]->Fill(genJetMax, weight);
	break;
      case HCAL :
	hmaxGenJetPtHCal[type]->Fill(genJetMax, weight);
	hmaxGenJetPtHCal[3]->Fill(genJetMax, weight);
	break;
      case HF :
	hmaxGenJetPtHF[type]->Fill(genJetMax, weight);
	hmaxGenJetPtHF[3]->Fill(genJetMax, weight);
	break;
      default:
	break;
      }
    }

  // JPT
  double jptMax(0.);
  int jptMaxIndex(-1);
  find_leading_jptjet(i_hyp, etaMin, etaMax, vetoCone, jptMax, jptMaxIndex);
  if (jptMaxIndex >= 0)
    hmaxBtagVsJPTPt->Fill(jptMax, BTag(jptJet,jptMaxIndex), weight);
  else
    hmaxBtagVsJPTPt->Fill(jptMax, 0.0, weight);

  // Calo 
  double caloJetMax(0.);
  find_leading_calojet(i_hyp, etaMin, etaMax, vetoCone, caloJetMax);
  // TrkJet
  double trkJetMax(0.);
  find_leading_trkjet(i_hyp, etaMin, etaMax, vetoCone, trkJetMax);
  // PFJet 
  double pfJetMax(0.);
  find_leading_pfjet(i_hyp, etaMin, etaMax, vetoCone, pfJetMax);
 
  // Fill the jet Pt histograms
  switch (jet ) {
  case ALLJET:
    hmaxJPTPt[type]->Fill(jptMax, weight);
    hmaxCaloJetPt[type]->Fill(caloJetMax, weight);
    hmaxTrkJetPt[type]->Fill(trkJetMax, weight);
    hmaxCaloTrkJetPt[type]->Fill((caloJetMax+trkJetMax)/2, weight);
    hmaxCaloTrkJet2Et[type]->Fill(std::max(caloJetMax,trkJetMax), weight);
    hmaxPFJetPt[type]->Fill(pfJetMax, weight);
    hmaxJPTPt[3]->Fill(jptMax, weight);
    hmaxCaloJetPt[3]->Fill(caloJetMax, weight);
    hmaxTrkJetPt[3]->Fill(trkJetMax, weight);
    hmaxCaloTrkJetPt[3]->Fill((caloJetMax+trkJetMax)/2, weight);
    hmaxCaloTrkJet2Et[3]->Fill(std::max(caloJetMax,trkJetMax), weight);
    hmaxPFJetPt[3]->Fill(pfJetMax, weight);
    break;
  case HCAL:
    hmaxJPTPtHCal[type]->Fill(jptMax, weight);
    hmaxCaloJetPtHCal[type]->Fill(caloJetMax, weight);
    hmaxTrkJetPtHCal[type]->Fill(trkJetMax, weight);
    hmaxCaloTrkJetPtHCal[type]->Fill((caloJetMax+trkJetMax)/2, weight);
    hmaxCaloTrkJet2EtHCal[type]->Fill(std::max(caloJetMax,trkJetMax), weight);
    hmaxPFJetPtHCal[type]->Fill(pfJetMax, weight);
    hmaxJPTPtHCal[3]->Fill(jptMax, weight);
    hmaxCaloJetPtHCal[3]->Fill(caloJetMax, weight);
    hmaxTrkJetPtHCal[3]->Fill(trkJetMax, weight);
    hmaxCaloTrkJetPtHCal[3]->Fill((caloJetMax+trkJetMax)/2, weight);
    hmaxCaloTrkJet2EtHCal[3]->Fill(std::max(caloJetMax,trkJetMax), weight);
    hmaxPFJetPtHCal[3]->Fill(pfJetMax, weight);
    break;
  case HF:
    hmaxJPTPtHF[type]->Fill(jptMax, weight);
    hmaxCaloJetPtHF[type]->Fill(caloJetMax, weight);
    hmaxTrkJetPtHF[type]->Fill(trkJetMax, weight);
    hmaxCaloTrkJetPtHF[type]->Fill((caloJetMax+trkJetMax)/2, weight);
    hmaxCaloTrkJet2EtHF[type]->Fill(std::max(caloJetMax,trkJetMax), weight);
    hmaxPFJetPtHF[type]->Fill(pfJetMax, weight);
    hmaxJPTPtHF[3]->Fill(jptMax, weight);
    hmaxCaloJetPtHF[3]->Fill(caloJetMax, weight);
    hmaxTrkJetPtHF[3]->Fill(trkJetMax, weight);
    hmaxCaloTrkJetPtHF[3]->Fill((caloJetMax+trkJetMax)/2, weight);
    hmaxCaloTrkJet2EtHF[3]->Fill(std::max(caloJetMax,trkJetMax), weight);
    hmaxPFJetPtHF[3]->Fill(pfJetMax, weight);
     break;
  default:
    break;
  }
}

// ===
// This is duplicating the getJets, FIXME
// ===

void find_leading_genjet(int i_hyp, double etaMin, double etaMax, double vetoCone, double & genJetMax) 
{
  for ( unsigned int i=0; i < cms2.genjets_p4().size(); ++i) {
    if ( cms2.genjets_p4()[i].pt() < genJetMax ) continue;
    if ( TMath::Abs(cms2.genjets_p4()[i].eta()) > etaMax ) continue;
    if ( TMath::Abs(cms2.genjets_p4()[i].eta()) < etaMin ) continue;
    if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.genjets_p4()[i])) < vetoCone ||
	 TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.genjets_p4()[i])) < vetoCone ) continue;
    genJetMax = cms2.genjets_p4()[i].pt();
  }
}

void find_leading_jptjet(int i_hyp, double etaMin, double etaMax, double vetoCone, double & jptMax, int &jptMaxIndex)
{
  double jec = 1.0; 
  for ( unsigned int i=0; i < cms2.jpts_p4().size(); ++i) {
     if ( TMath::Abs(cms2.jpts_p4()[i].eta()) > etaMax ) continue;
     if ( TMath::Abs(cms2.jpts_p4()[i].eta()) < etaMin ) continue;
     if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.jpts_p4()[i])) < vetoCone ||
	  TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.jpts_p4()[i])) < vetoCone ) continue;
     if(applyJEC)
       jec = jetCorrection(cms2.jpts_p4()[i], jet_corrector_jpt);
     //jec =  cms2.jpts_cor()[i]; // use the one in CMS2 ntuple
     if ( cms2.jpts_p4()[i].pt() * jec < jptMax ) continue;
     jptMax = cms2.jpts_p4()[i].pt() * jec;
     jptMaxIndex = i;
   }
}

void find_leading_calojet(int i_hyp, double etaMin, double etaMax, double vetoCone, double & caloJetMax)
{
  double jec = 1.0;
  for ( unsigned int i=0; i < cms2.jets_pat_jet_uncorp4().size(); ++i) {
    if ( TMath::Abs(cms2.jets_pat_jet_uncorp4()[i].eta()) > etaMax ) continue;
    if ( TMath::Abs(cms2.jets_pat_jet_uncorp4()[i].eta()) < etaMin ) continue;
    if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.jets_pat_jet_uncorp4()[i])) < vetoCone ||
	 TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.jets_pat_jet_uncorp4()[i])) < vetoCone ) continue;   
    if(applyJEC)
      jec = jetCorrection(cms2.jets_pat_jet_uncorp4()[i], jet_corrector_calo);
    //jec = 1.0/cms2.jets_pat_noCorrF()[i]; 
    if ( cms2.jets_pat_jet_uncorp4()[i].pt() * jec < caloJetMax ) continue;
    caloJetMax = cms2.jets_pat_jet_uncorp4()[i].pt() * jec;
    //if ( cms2.jets_pat_jet_p4()[i].pt() < caloJetMax ) continue;
    //caloJetMax = cms2.jets_pat_jet_p4()[i].pt();
  }
}

void find_leading_trkjet(int i_hyp, double etaMin, double etaMax, double vetoCone, double & trkJetMax)
{
  double jec = 1.0;
  for ( unsigned int i=0; i < cms2.trkjets_p4().size(); ++i) {
    if ( TMath::Abs(cms2.trkjets_p4()[i].eta()) > etaMax ) continue;
    if ( TMath::Abs(cms2.trkjets_p4()[i].eta()) < etaMin ) continue;
    if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.trkjets_p4()[i])) < vetoCone ||
	 TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.trkjets_p4()[i])) < vetoCone ) continue;
    if(applyJEC)
      jec = jetCorrection(cms2.trkjets_p4()[i], jet_corrector_trk);
    //jec = cms2.trkjets_cor()[i];
    if ( cms2.trkjets_p4()[i].pt() * jec < trkJetMax ) continue;
    trkJetMax = cms2.trkjets_p4()[i].pt() * jec;
  }
}

void find_leading_pfjet(int i_hyp, double etaMin, double etaMax, double vetoCone, double & pfJetMax)
{
  double jec = 1.0;
  for ( unsigned int i=0; i < cms2.pfjets_p4().size(); ++i) {
    if ( TMath::Abs(cms2.pfjets_p4()[i].eta()) > etaMax ) continue;
    if ( TMath::Abs(cms2.pfjets_p4()[i].eta()) < etaMin ) continue;
    if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.pfjets_p4()[i])) < vetoCone ||
	 TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.pfjets_p4()[i])) < vetoCone ) continue;
    if(applyJEC)
      jec = jetCorrection(cms2.pfjets_p4()[i], jet_corrector_pf);
    //jec = cms2.pfjets_cor()[i];
    if ( cms2.pfjets_p4()[i].pt() * jec < pfJetMax ) continue;
    pfJetMax = cms2.pfjets_p4()[i].pt() * jec;
  }
}

void fill_val_plots(int i_hyp, cuts_t cut_passed, double weight)
{
  HypothesisType type = getHypothesisType(cms2.hyp_type()[i_hyp]);
  // Fill MET related plots
  if(CheckCuts( (1<<PASS_JETVETO) | (1<<PASS_LL_FINAL) | (1<<PASS_LT_FINAL)  , cut_passed)) {
    hmetVal[type] -> Fill(metValue(), weight);
    hmetVal[3] -> Fill(metValue(), weight);
    hmetProjVal[type] -> Fill(projectedMet(i_hyp), weight);
    hmetProjVal[3] -> Fill(projectedMet(i_hyp), weight);
  }
  // Fill Leading Jet Pt
  if(CheckCuts(  (1<<PASS_LL_FINAL) | (1<<PASS_LT_FINAL)  , cut_passed)) {
    double pfJetMax(0.);
    find_leading_pfjet(i_hyp, 0.0, 0.5, 0.3, pfJetMax);
    hmaxPFJetPtVal[type]->Fill(pfJetMax, weight);
    hmaxPFJetPtVal[3]->Fill(pfJetMax, weight);
    hdilMassVal[type]->Fill(cms2.hyp_p4()[i_hyp].mass(), weight);
    hdilMassVal[3]->Fill(cms2.hyp_p4()[i_hyp].mass(), weight);
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
  for ( unsigned int i=0; i < cms2.mus_p4().size(); ++i ){
    if ( cms2.mus_p4().at(i).pt() < 20 ) continue;
    if ( fakableMuon(i) ) 
      hFakableRateSingleMuon->Fill(cms2.mus_p4().at(i).pt(), fabs(cms2.mus_p4().at(i).eta()) );
    if ( goodMuonIsolated(i) ) 
      hFinalRateSingleMuon->Fill(cms2.mus_p4().at(i).pt(), fabs(cms2.mus_p4().at(i).eta()) );  
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
					   bool passedLTMuFakableRequirements, 
					   bool passedLLMuFakableRequirements,
					   bool passedLTFinalRequirements,
					   bool passedLLFinalRequirements)
{
    HypothesisType type = getHypothesisType(cms2.hyp_type()[i_hyp]);

    int nbins_eta = helFRfakable_fakerate[3]->GetNbinsX();
    float max_eta = helFRfakable_fakerate[3]->GetXaxis()->GetBinLowEdge(nbins_eta+1);
    int nbins_pt  = helFRfakable_fakerate[3]->GetNbinsY();
    float max_pt  = helFRfakable_fakerate[3]->GetYaxis()->GetBinLowEdge(nbins_pt+1);

    if ( passedLTElFakableRequirements && !passedLTFinalRequirements && passedLLFinalRequirements ){
        helFRfakable[type]->Fill(cms2.hyp_lt_p4().at(i_hyp).pt(),
                fabs(cms2.hyp_lt_p4().at(i_hyp).eta()),
                weight);
        helFRfakable[3]->Fill(cms2.hyp_lt_p4().at(i_hyp).pt(),
                fabs(cms2.hyp_lt_p4().at(i_hyp).eta()),
                weight);

        // overflow goes into last bin
        float eta = fabs(cms2.hyp_lt_p4().at(i_hyp).eta());
        if (eta > max_eta) eta = max_eta-.1;
        float pt = cms2.hyp_lt_p4().at(i_hyp).pt();
        if (pt > max_pt) pt = max_pt-.1;

        helFRfakable_fakerate[type]->Fill(eta,pt,weight);
        helFRfakable_fakerate[3]->Fill(eta,pt,weight);
    }
    if ( passedLLElFakableRequirements && passedLTFinalRequirements && !passedLLFinalRequirements ){
        helFRfakable[type]->Fill(cms2.hyp_ll_p4().at(i_hyp).pt(),
                fabs(cms2.hyp_ll_p4().at(i_hyp).eta()),
                weight);
        helFRfakable[3]->Fill(cms2.hyp_ll_p4().at(i_hyp).pt(),
                fabs(cms2.hyp_ll_p4().at(i_hyp).eta()),
                weight);

        // overflow goes into last bin
        float eta = fabs(cms2.hyp_ll_p4().at(i_hyp).eta());
        if (eta > max_eta) eta = max_eta-.1;
        float pt = cms2.hyp_ll_p4().at(i_hyp).pt();
        if (pt > max_pt) pt = max_pt-.1;

        helFRfakable_fakerate[type]->Fill(eta,pt,weight);
        helFRfakable_fakerate[3]->Fill(eta,pt,weight);
    }
    if ( passedLTMuFakableRequirements && !passedLTFinalRequirements && passedLLFinalRequirements ){
        hmuFRfakable[type]->Fill(cms2.hyp_lt_p4().at(i_hyp).pt(),
                fabs(cms2.hyp_lt_p4().at(i_hyp).eta()),
                weight);
        hmuFRfakable[3]->Fill(cms2.hyp_lt_p4().at(i_hyp).pt(),
                fabs(cms2.hyp_lt_p4().at(i_hyp).eta()),
                weight);

        // overflow goes into last bin
        float eta = fabs(cms2.hyp_lt_p4().at(i_hyp).eta());
        if (eta > max_eta) eta = max_eta-.1;
        float pt = cms2.hyp_lt_p4().at(i_hyp).pt();
        if (pt > max_pt) pt = max_pt-.1;

        hmuFRfakable_fakerate[type]->Fill(eta,pt,weight);
        hmuFRfakable_fakerate[3]->Fill(eta,pt,weight);
    }
    if ( passedLLMuFakableRequirements && passedLTFinalRequirements && !passedLLFinalRequirements ){
        hmuFRfakable[type]->Fill(cms2.hyp_ll_p4().at(i_hyp).pt(),
                fabs(cms2.hyp_ll_p4().at(i_hyp).eta()),
                weight);
        hmuFRfakable[3]->Fill(cms2.hyp_ll_p4().at(i_hyp).pt(),
                fabs(cms2.hyp_ll_p4().at(i_hyp).eta()),
                weight);

        // overflow goes into last bin
        float eta = fabs(cms2.hyp_ll_p4().at(i_hyp).eta());
        if (eta > max_eta) eta = max_eta-.1;
        float pt = cms2.hyp_ll_p4().at(i_hyp).pt();
        if (pt > max_pt) pt = max_pt-.1;

        hmuFRfakable_fakerate[type]->Fill(eta,pt,weight);
        hmuFRfakable_fakerate[3]->Fill(eta,pt,weight);
    }
}

bool 
toptag(WWJetType type, int i_hyp, double minPt,
       std::vector<LorentzVector> ignoreJets=std::vector<LorentzVector>())
{
     std::vector<LorentzVector> jets;
     const double vetoCone    = 0.3;

     switch ( type ){
     case pfJet:
       for ( unsigned int i=0; i < cms2.pfjets_p4().size(); ++i) {
	 if ( cms2.pfjets_p4()[i].pt() < minPt ) continue;
	 bool ignoreJet = false;
	 for ( std::vector<LorentzVector>::const_iterator ijet = ignoreJets.begin();
	       ijet != ignoreJets.end(); ++ijet )
	   if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(*ijet,cms2.pfjets_p4()[i])) < vetoCone ) ignoreJet=true;
	 if ( ignoreJet ) continue;
	 if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.pfjets_p4()[i])) < vetoCone ||
	      TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.pfjets_p4()[i])) < vetoCone ) continue;
	 if ( defaultBTag(type,i) ) return true;
       }
       break;
     case CaloJet:
       for ( unsigned int i=0; i < cms2.jets_p4().size(); ++i) {
	 if ( cms2.jets_p4()[i].pt() < minPt ) continue;
	 bool ignoreJet = false;
	 for ( std::vector<LorentzVector>::const_iterator ijet = ignoreJets.begin();
	       ijet != ignoreJets.end(); ++ijet )
	   if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(*ijet,cms2.jets_p4()[i])) < vetoCone ) ignoreJet=true;
	 if ( ignoreJet ) continue;
	 if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.jets_p4()[i])) < vetoCone ||
	      TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.jets_p4()[i])) < vetoCone ) continue;
// 	 if ( defaultBTag(type,i) && ignoreJets.size()==1 ){
// 	   cout << "b-tagged jet pt: " << cms2.jets_p4()[i].pt() << " \teta: " << cms2.jets_p4()[i].eta() <<
// 	     " \tphi: " << cms2.jets_p4()[i].phi() << endl;
// 	 }
	 if ( defaultBTag(type,i) ) return true;
       }
       break;
     default:
       std::cout << "ERROR: not supported jet type is requested: " << type << " FixIt!" << std::endl;
     }
     return false;
}

bool hypo (int i_hyp, double weight, RooDataSet* dataset, bool zStudy, bool realData) 
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
  // float weight = cms2.evt_scale1fb() * kFactor;

  if (!realData) monitor.nEvtProcessed = cms2.evt_nEvts();
  // monitor.count(cms2, type, "Total number before cuts");
  
  // if ( cms2.hyp_FVFit_prob()[i_hyp] < 0.005 ) return;
  // monitor.count(cms2, type, "after vertex cut");
  
  if ( realData && ! passedTriggerRequirements( hypType(i_hyp) ) )return false;
  monitor.count(cms2, type, "after trigger requirements");

  // Require same sign
  if ( cms2.hyp_lt_id()[i_hyp] * cms2.hyp_ll_id()[i_hyp] > 0 ) return false;

  // Baseline cuts
  if (abs(cms2.hyp_lt_id()[i_hyp]) == 13 && !ww_muBase(cms2.hyp_lt_index()[i_hyp]) ) return false;
  if (abs(cms2.hyp_ll_id()[i_hyp]) == 13 && !ww_muBase(cms2.hyp_ll_index()[i_hyp]) ) return false;
  if (abs(cms2.hyp_lt_id()[i_hyp]) == 11 && !ww_elBase(cms2.hyp_lt_index()[i_hyp]) ) return false;
  if (abs(cms2.hyp_ll_id()[i_hyp]) == 11 && !ww_elBase(cms2.hyp_ll_index()[i_hyp]) ) return false;
  
  monitor.count(cms2,type,"after previous + baseline cuts");
 
  // TEMPORARY
  if (gSystem->Getenv("Sync")) // Synchronization info
  {
    if (nGoodVertex()<1) return false;

    monitor.count(cms2,type,"after previous + primary verex");

    if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && !ww_mud0PV(cms2.hyp_lt_index()[i_hyp]) ) return false;
    if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && !ww_mud0PV(cms2.hyp_ll_index()[i_hyp]) ) return false;
    if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11 && !ww_eld0PV(cms2.hyp_lt_index()[i_hyp]) ) return false;
    if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11 && !ww_eld0PV(cms2.hyp_ll_index()[i_hyp]) ) return false;

    monitor.count(cms2,type,"after previous + d0");

    if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && ww_muIsoVal(cms2.hyp_lt_index()[i_hyp])>0.15 ) return false;
    if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && ww_muIsoVal(cms2.hyp_ll_index()[i_hyp])>0.15 ) return false;
    if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11 && ww_elIsoVal(cms2.hyp_lt_index()[i_hyp])>0.1 ) return false;
    if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11 && ww_elIsoVal(cms2.hyp_ll_index()[i_hyp])>0.1 ) return false;

    monitor.count(cms2,type,"after previous + iso");
    
    if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && !ww_muId(cms2.hyp_lt_index()[i_hyp]) ) return false;
    if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && !ww_muId(cms2.hyp_ll_index()[i_hyp]) ) return false;
    if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11 && 
	! (electronId_VBTF(cms2.hyp_lt_index()[i_hyp], VBTF_35X_80) & (1<<ELEID_ID))  ) return false;
    if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11 && 
	! (electronId_VBTF(cms2.hyp_ll_index()[i_hyp], VBTF_35X_80) & (1<<ELEID_ID))  ) return false;

    monitor.count(cms2,type,"after previous + lepton id");

    if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11 && 
	(fabs(cms2.els_conv_dist().at(cms2.hyp_lt_index()[i_hyp])) < 0.02 &&
	 fabs(cms2.els_conv_dcot().at(cms2.hyp_lt_index()[i_hyp])) < 0.02 )) return false;
    if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11 && 
	(fabs(cms2.els_conv_dist().at(cms2.hyp_ll_index()[i_hyp])) < 0.02 &&
	 fabs(cms2.els_conv_dcot().at(cms2.hyp_ll_index()[i_hyp])) < 0.02 )) return false;

    if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11 && 
	cms2.els_exp_innerlayers39X().at(cms2.hyp_lt_index()[i_hyp]) != 0) return false;
    if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11 && 
	cms2.els_exp_innerlayers39X().at(cms2.hyp_ll_index()[i_hyp]) != 0) return false;

    monitor.count(cms2,type,"after previous + conv rejection");

    if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && !goodMuonIsolated(cms2.hyp_lt_index()[i_hyp]) ) return false;
    if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && !goodMuonIsolated(cms2.hyp_ll_index()[i_hyp]) ) return false;
    if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11 && !goodElectronIsolated(cms2.hyp_lt_index()[i_hyp]) ) return false;
    if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11 && !goodElectronIsolated(cms2.hyp_ll_index()[i_hyp]) ) return false;

    monitor.count(cms2,type,"after previous + lepton id/iso");
    if ( metValue()<20 ) return false;
    monitor.count(cms2,type,"after previous + met>20");
    
    if (cms2.hyp_p4().at(i_hyp).mass2()<0 || 
	cms2.hyp_p4()[i_hyp].mass() < 12) return false;
    monitor.count(cms2,type,"after previous + mll cuts");
    
    if ( type == EE || type == MM) {
      if (inZmassWindow(cms2.hyp_p4()[i_hyp].mass())) return false;
    }
    monitor.count(cms2,type,"after previous + z mass cuts in EE/MM");
    
    if (!passedMetRequirements(i_hyp)) return false;
    monitor.count(cms2,type,"after previous + Full MET cuts: ");
    
    if ( numberOfJets(i_hyp)>0 ) return false;
    monitor.count(cms2,type,"after previous + JetVeto cuts: ");

    if (numberOfSoftMuons(i_hyp,true)>0) return false;
    monitor.count(cms2,type,"after previous + soft muon veto: ");

    if (numberOfExtraLeptons(i_hyp,10)>0) return false;
    monitor.count(cms2,type,"after previous + extra lepton veto: ");
  } // end of Synchronization info

  if (nGoodVertex()<1) return false;
  hnGoodVertex -> Fill(nGoodVertex(), weight);
  
  if (cms2.hyp_p4().at(i_hyp).mass2()<0 || 
      cms2.hyp_p4()[i_hyp].mass() < 12) return false;

  // check electron isolation and id (no selection at this point)
  checkIsolation(i_hyp, weight);
  
  cuts_t cuts_passed = 0;
  
  // == Z mass veto using hyp_leptons for ee and mumu final states
  if ( type == EE || type == MM) {
    if ( !zStudy && !inZmassWindow(cms2.hyp_p4()[i_hyp].mass()))    cuts_passed |= (1<<PASS_ZSEL);
    if ( zStudy  && inZmassWindow(cms2.hyp_p4()[i_hyp].mass()))     cuts_passed |= (1<<PASS_ZSEL);
  }
  if (type == EM)     cuts_passed |= (1<<PASS_ZSEL);

  // Z veto using additional leptons in the event
  // if (additionalZveto()) return;
 
  // == MET
  if (zStudy) cuts_passed |= (1<<PASS_MET);
  if (!zStudy && passedMetRequirements(i_hyp)) { 
    cuts_passed |= (1<<PASS_MET);
  }
  
  // == letpon ID and Isolation
  bool passedLTFinalRequirements = true;
  bool passedLLFinalRequirements = true;
  bool passedLTElFakableRequirements = true;
  bool passedLLElFakableRequirements = true;
  bool passedLTMuFakableRequirements = true;
  bool passedLLMuFakableRequirements = true;

  // Muon quality cuts, including isolation
  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13){
    if ( !goodMuonIsolated(cms2.hyp_lt_index()[i_hyp]) ) passedLTFinalRequirements = false;
    if ( !fakableMuon(cms2.hyp_lt_index()[i_hyp]) ) passedLTMuFakableRequirements = false;
    passedLTElFakableRequirements = false;
  }
  if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13){
    if ( !goodMuonIsolated(cms2.hyp_ll_index()[i_hyp]) ) passedLLFinalRequirements = false;
    if ( !fakableMuon(cms2.hyp_ll_index()[i_hyp]) ) passedLLMuFakableRequirements = false;
    passedLLElFakableRequirements = false;
  } 
  // Electron quality cuts, including isolation
  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11){
    if ( !goodElectronIsolated(cms2.hyp_lt_index()[i_hyp]) ) passedLTFinalRequirements = false;
    if ( !fakableElectron(cms2.hyp_lt_index()[i_hyp]) ) passedLTElFakableRequirements = false;
    passedLTMuFakableRequirements = false;
  }
  if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11){
    if ( !goodElectronIsolated(cms2.hyp_ll_index()[i_hyp]) ) passedLLFinalRequirements = false;
    if ( !fakableElectron(cms2.hyp_ll_index()[i_hyp]) ) passedLLElFakableRequirements = false;
    passedLLMuFakableRequirements = false;
  }

  if ( passedLTFinalRequirements ) cuts_passed |= (1<<PASS_LT_FINAL);
  if ( passedLLFinalRequirements ) cuts_passed |= (1<<PASS_LL_FINAL);

  // == Jet-veto
  const std::vector<LorentzVector>& vetojets(getDefaultJets(i_hyp, false));
  unsigned int nJets = vetojets.size();
  if (nJets==0) 
    cuts_passed |= (1<<PASS_JETVETO);
  if (nJets==1){
    cuts_passed |= (1<<PASS_1JET);
    if ( getDefaultJets(i_hyp,true).size()==1 ){
      cuts_passed |= (1<<PASS_1BJET);
    }
  }
  // trkjet veto
  // if ( !passTrkJetVeto(i_hyp) ) return;
  // == ExtraVeto Muons
  int countmus = numberOfSoftMuons(i_hyp,true);
  int nExtraVetoMuons = numberOfSoftMuons(i_hyp,true,vetojets);
  if ( nExtraVetoMuons == 0) 
    cuts_passed |=   (1<<PASS_SOFTMUVETO);
  if ( numberOfExtraLeptons(i_hyp,10) == 0) 
    cuts_passed |=   (1<<PASS_EXTRALEPTONVETO);
  if ( ! toptag(CaloJet,i_hyp,0,vetojets) )
    cuts_passed |=   (1<<PASS_TOPVETO);

  // -------------------------------------------------------------------//
  // Finished checking the cuts, fill histograms before the final sel   //
  // -------------------------------------------------------------------//

  if (CheckCutsNM1(pass_all, (1<<PASS_LT_FINAL)|(1<<PASS_LL_FINAL), cuts_passed) ) {
    if(dataset)
      getIsolationSidebandsAfterSelections(i_hyp, weight, dataset, passedLTFinalRequirements && passedLLFinalRequirements);
    countFakableObjectsAfterAllSelections(i_hyp, weight, 
					  passedLTElFakableRequirements, passedLLElFakableRequirements, 
					  passedLTMuFakableRequirements, passedLLMuFakableRequirements, 
					  passedLTFinalRequirements, passedLLFinalRequirements);
  }
    
  // Jet-veto effciency studies

  if (CheckCutsNM1(pass_all, (1<<PASS_JETVETO) | (1<<PASS_SOFTMUVETO) | (1<<PASS_EXTRALEPTONVETO), cuts_passed)) {
    find_most_energetic_jets(i_hyp,weight,realData,0.0,5.0);
    find_most_energetic_jets(i_hyp,weight,realData,0.0,3.0);
    find_most_energetic_jets(i_hyp,weight,realData,3.0,5.0);
    if(zStudy) 
      getJetResponseFromZBalance(i_hyp, weight, realData, 0.0, 5.0);
  }
  
  // Make some validation plots 
  if(CheckCuts(  (1<<PASS_LL_FINAL) | (1<<PASS_LT_FINAL), cuts_passed)) {
    fill_val_plots(i_hyp, cuts_passed, weight);
  }
  
  // Do Drell Yan Estimation
  if(CheckCutsNM1(pass_all, (1<<PASS_ZSEL) | (1<<PASS_MET) , cuts_passed))
    fill_dyest_histograms(i_hyp, weight);

  // top background related histograms
  if ( CheckCutsNM1( pass_all, 
		     (1<<PASS_SOFTMUVETO) | (1<<PASS_TOPVETO) | (1<<PASS_JETVETO) | (1<< PASS_ZSEL) | (1<<PASS_MET), 
		     cuts_passed ) ){
    // 2D hist for muon tag counting
    hextramuonsvsnjet[type]->Fill(countmus, nJets, weight);
    hextramuonsvsnjet[3]->Fill(countmus, nJets, weight);
    float tag = 1; // not tagged
    if ( ! CheckCuts(1<<PASS_SOFTMUVETO, cuts_passed) ){
      if ( ! CheckCuts(1<<PASS_TOPVETO, cuts_passed) )
	tag = 4; // both
      else
	tag = 2; // soft muon
    } else {
      if ( ! CheckCuts(1<<PASS_TOPVETO, cuts_passed) )
	tag = 3; // b-tagged
    }
    float evtType = 1; // passed jet veto
    if ( ! CheckCuts(1<<PASS_JETVETO, cuts_passed) ){
    if ( CheckCuts(1<<PASS_1JET,cuts_passed) ){
      if ( CheckCuts(1<<PASS_1BJET,cuts_passed) ){
	  evtType = 3; // have 1 btagged jet
	} else { 
	  evtType = 2; // have 1 non-btagged jet
	}
      } else {
	evtType = 4; // more than 1 jet
      }
    }
    bool zevent = (type == EE || type == MM);
    if ( CheckCuts((1<< PASS_ZSEL), cuts_passed) ) zevent = false;
    if ( zevent ) {
      htoptagz[type]->Fill(tag-0.5,evtType-0.5,weight);
      htoptagz[3]->Fill(tag-0.5,evtType-0.5,weight);
    }
    if ( CheckCuts((1<<PASS_ZSEL)|(1<<PASS_MET), cuts_passed) ){
      htoptag[type]->Fill(tag-0.5,evtType-0.5,weight);
      htoptag[3]->Fill(tag-0.5,evtType-0.5,weight);
    }
  }

  if ( !realData && CheckCutsNM1( pass_all, (1<<PASS_SOFTMUVETO), cuts_passed ) ) {
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
  
  // make the final selections
  if(! CheckCuts(pass_all, cuts_passed)) return false;
  monitor.count(cms2,type,"after all cuts (including soft and extra lepton veto)");
  
  // if ( toptag(jetType(),i_hyp,0) ) return false;
  // if (nExtraVetoMuons) return false;
  // if ( toptag(CaloJet,i_hyp,0) ) return false;
  // monitor.count(cms2,type,"after all cuts + top tagging");

  // -------------------------------------------------------------------//
  // If we made it to here, we passed all cuts and we are ready to fill // 
  // histograms after the full selection                                //
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
     return true;
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

void AddIsoSignalControlSample( int i_hyp, double weight, RooDataSet* dataset, bool realData) {
  if ( !dataset ) return;
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

void initializeHistograms(const char *prefix, bool qcdBackground){
  hypos_total          = new TH1F(Form("%s_hypos_total",prefix),"Total number of hypothesis counts",4,0,4);
  hypos_total_weighted = new TH1F(Form("%s_hypos_total_weighted",prefix),"Total number of hypotheses (weighted)",4,0,4);
  hypos_total_weighted->Sumw2();

  for (unsigned int i=0; i<4; ++i){
    hypos_total->GetXaxis()->SetBinLabel(i+1,HypothesisTypeName(i));
    hypos_total_weighted->GetXaxis()->SetBinLabel(i+1,HypothesisTypeName(i));
  }
  
  const Double_t ptbins[4] = {20,30,80,200};
  // Same binning as fake rate histograms
  const Double_t ptbins_fakerate[6] = {10.,15.,20.,25.,30.,35.};
  const Double_t etabins_fakerate[5] = {0.0,1.0,1.479,2.0,2.5};

  
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
    helFRfakable_fakerate[i] = new TH2F(Form("%s_helFRfakable_fakerate_%s", prefix,HypothesisTypeName(i)), "FR study: rate of fakable objects",4,etabins_fakerate,5,ptbins_fakerate);
    helFRfakable_fakerate[i]->Sumw2();
    hmuFRfakable[i]     = new TH2F(Form("%s_hmuFRfakable_%s",    prefix,HypothesisTypeName(i)), "FR study: rate of fakable objects", 3,ptbins,2,0,3.0);
    hmuFRfakable[i]->Sumw2();
    hmuFRfakable_fakerate[i] = new TH2F(Form("%s_hmuFRfakable_fakerate_%s", prefix,HypothesisTypeName(i)), "FR study: rate of fakable objects",4,etabins_fakerate,5,ptbins_fakerate);
    hmuFRfakable_fakerate[i]->Sumw2();

    const Double_t jetEtbins[12] = {0,10,15,20,25,30,40,50,60,80,100,200};
    hmaxGenJetPt[i]  = new TH1F(Form("%s_hmaxGenJetPt_%s", prefix,HypothesisTypeName(i)), "most energetic jet (|eta|<5) Et (GenJet) JetVeto NM1", 11,jetEtbins);
    hmaxJPTPt[i]  = new TH1F(Form("%s_hmaxJPTPt_%s", prefix,HypothesisTypeName(i)), "most energetic jet (|eta|<5) Et (JPT) JetVeto NM1", 11,jetEtbins);
    hmaxCaloJetPt[i]  = new TH1F(Form("%s_hmaxCaloJetPt_%s", prefix,HypothesisTypeName(i)), "most energetic jet (|eta|<5) Et (CaloJet) JetVeto NM1", 11,jetEtbins);
    hmaxTrkJetPt[i]  = new TH1F(Form("%s_hmaxTrkJetPt_%s", prefix,HypothesisTypeName(i)), "most energetic jet (|eta|<5) Et (TrkJet) JetVeto NM1", 11,jetEtbins);
    hmaxCaloTrkJetPt[i]  = new TH1F(Form("%s_hmaxCaloTrkJetPt_%s", prefix,HypothesisTypeName(i)), "most energetic jet (|eta|<5) Et (average of Calo + Trk Jets) JetVeto NM1", 11,jetEtbins);
    hmaxCaloTrkJet2Et[i]  = new TH1F(Form("%s_hmaxCaloTrkJet2Et_%s", prefix,HypothesisTypeName(i)), "most energetic jet (|eta|<5) Et (Max of Calo + Trk Jets) JetVeto NM1", 11,jetEtbins);
    hmaxPFJetPt[i]  = new TH1F(Form("%s_hmaxPFJetPt_%s", prefix,HypothesisTypeName(i)), "most energetic jet (|eta|<5) Et (PFJet) JetVeto NM1", 11,jetEtbins);

    hmaxGenJetPtHCal[i]  = new TH1F(Form("%s_hmaxGenJetPtHCal_%s", prefix,HypothesisTypeName(i)), "most energetic jet (|eta|<3) Et (GenJet) JetVeto NM1", 11,jetEtbins);
    hmaxJPTPtHCal[i]  = new TH1F(Form("%s_hmaxJPTPtHCal_%s", prefix,HypothesisTypeName(i)), "most energetic jet (|eta|<3) Et (JPT) JetVeto NM1", 11,jetEtbins);
    hmaxCaloJetPtHCal[i]  = new TH1F(Form("%s_hmaxCaloJetPtHCal_%s", prefix,HypothesisTypeName(i)), "most energetic jet (|eta|<3) Et (CaloJet) JetVeto NM1", 11,jetEtbins);
    hmaxTrkJetPtHCal[i]  = new TH1F(Form("%s_hmaxTrkJetPtHCal_%s", prefix,HypothesisTypeName(i)), "most energetic jet (|eta|<3) Et (TrkJet) JetVeto NM1", 11,jetEtbins);
    hmaxCaloTrkJetPtHCal[i]  = new TH1F(Form("%s_hmaxCaloTrkJetPtHCal_%s", prefix,HypothesisTypeName(i)), "most energetic jet (|eta|<3) Et (average of Calo + Trk Jets) JetVeto NM1", 11,jetEtbins);
    hmaxCaloTrkJet2EtHCal[i]  = new TH1F(Form("%s_hmaxCaloTrkJet2EtHCal_%s", prefix,HypothesisTypeName(i)), "most energetic jet (|eta|<3) Et (Max of Calo + Trk Jets) JetVeto NM1", 11,jetEtbins);
    hmaxPFJetPtHCal[i]  = new TH1F(Form("%s_hmaxPFJetPtHCal_%s", prefix,HypothesisTypeName(i)), "most energetic jet (|eta|<3) Et (PFJet) JetVeto NM1", 11,jetEtbins);

    hmaxGenJetPtHF[i]  = new TH1F(Form("%s_hmaxGenJetPtHF_%s", prefix,HypothesisTypeName(i)), "most energetic jet (|eta|<3) Et (GenJet) JetVeto NM1", 11,jetEtbins);
    hmaxJPTPtHF[i]  = new TH1F(Form("%s_hmaxJPTPtHF_%s", prefix,HypothesisTypeName(i)), "most energetic jet (|eta|<3) Et (JPT) JetVeto NM1", 11,jetEtbins);
    hmaxCaloJetPtHF[i]  = new TH1F(Form("%s_hmaxCaloJetPtHF_%s", prefix,HypothesisTypeName(i)), "most energetic jet (|eta|<3) Et (CaloJet) JetVeto NM1", 11,jetEtbins);
    hmaxTrkJetPtHF[i]  = new TH1F(Form("%s_hmaxTrkJetPtHF_%s", prefix,HypothesisTypeName(i)), "most energetic jet (|eta|<3) Et (TrkJet) JetVeto NM1", 11,jetEtbins);
    hmaxCaloTrkJetPtHF[i]  = new TH1F(Form("%s_hmaxCaloTrkJetPtHF_%s", prefix,HypothesisTypeName(i)), "most energetic jet (|eta|<3) Et (average of Calo + Trk Jets) JetVeto NM1", 11,jetEtbins);
    hmaxCaloTrkJet2EtHF[i]  = new TH1F(Form("%s_hmaxCaloTrkJet2EtHF_%s", prefix,HypothesisTypeName(i)), "most energetic jet (|eta|<3) Et (Max of Calo + Trk Jets) JetVeto NM1", 11,jetEtbins);
    hmaxPFJetPtHF[i]  = new TH1F(Form("%s_hmaxPFJetPtHF_%s", prefix,HypothesisTypeName(i)), "most energetic jet (|eta|<3) Et (PFJet) JetVeto NM1", 11,jetEtbins);
   
    hPFJetResponse[i] = new TH2F(Form("%s_hPFJetResponse_%s",prefix,HypothesisTypeName(i) ), Form("%s - PFJet response to Z pt - %s",prefix,HypothesisTypeName(i) ), 11, jetEtbins, 11, jetEtbins);
    hPFJetRelResponse[i] = new TH2F(Form("%s_hPFJetRelResponse_%s",prefix,HypothesisTypeName(i) ), Form("%s - PFJet relative response to Z pt - %s",prefix,HypothesisTypeName(i) ), 20, 0, 100, 20, 0, 2);

    hPFJetResponseWithZero[i] = new TH2F(Form("%s_hPFJetResponseWithZero_%s",prefix,HypothesisTypeName(i) ), Form("%s - PFJet response to Z pt - %s",prefix,HypothesisTypeName(i) ), 11, jetEtbins, 11, jetEtbins);
    hPFJetRelResponseWithZero[i] = new TH2F(Form("%s_hPFJetRelResponseWithZero_%s",prefix,HypothesisTypeName(i) ), Form("%s - PFJet relative response to Z pt - %s",prefix,HypothesisTypeName(i) ), 20, 0, 100, 20, 0, 2);

    hmetVal[i]           = new TH1F(Form("%s_hmetVal_%s",      prefix,HypothesisTypeName(i)), "MET after ll selection/jetveto", 20,0.,100.);
    hmetVal[i]->GetXaxis()->SetTitle("tcMET [GeV]");
    hmetVal[i]->GetYaxis()->SetTitle("Events/(5 GeV)");
    hmetProjVal[i]       = new TH1F(Form("%s_hmetProjVal_%s",      prefix,HypothesisTypeName(i)), "Projected MET after ll selection/jetveto", 20,0.,100.);
    hmetProjVal[i]->GetXaxis()->SetTitle("Projected tcMET [GeV]");
    hmetProjVal[i]->GetYaxis()->SetTitle("Events/(5 GeV)");
    hmaxPFJetPtVal[i]    = new TH1F(Form("%s_hmaxPFJetPtVal_%s",   prefix,HypothesisTypeName(i)), "Leading PFJet Pt after ll selection", 20, 0., 100.);
    hmaxPFJetPtVal[i]->GetXaxis()->SetTitle("Leading PF Jet Pt [GeV]");
    hmaxPFJetPtVal[i]->GetYaxis()->SetTitle("Events/(5 GeV)");
    hdilMassVal[i]   = new TH1F(Form("%s_hdilMassVal_%s",  prefix,HypothesisTypeName(i)), "Di-lepton mass after ll selection", 40, 0., 200.);
    hdilMassVal[i]->GetXaxis()->SetTitle("Dilepton mass [GeV/c^{2}]");
    hdilMassVal[i]->GetYaxis()->SetTitle("Events/(5 GeV)");
        
    hdilMassWithMetDYEst[i] = new TH1F(Form("%s_hdilMassWithMetDYEst_%s",  prefix,HypothesisTypeName(i)), "Di-lepton mass with MET for DY Estimation", 40, 0., 200.);
    hdilMassNoMetDYEst[i] = new TH1F(Form("%s_hdilMassNoMetDYEst_%s",  prefix,HypothesisTypeName(i)), "Di-lepton mass without MET for DY Estimation", 40, 0., 200.);
    hmetInDYEst[i] = new TH1F(Form("%s_hmetInDYEst_%s",  prefix,HypothesisTypeName(i)), "MET in Z mass for DY Estimation", 40, 0., 200.);
    hmetOutDYEst[i] = new TH1F(Form("%s_hmetOutDYEst_%s",  prefix,HypothesisTypeName(i)), "MET outside Z mass for DY Estimation", 40, 0., 200.);


    
    htoptagz[i] = new TH2F(Form("%s_htoptagz_%s",prefix,HypothesisTypeName(i)),
			   "Top tagging on Z-sample",4,0,4,4,0,4);
    htoptagz[i]->Sumw2();
    htoptag[i]  = new TH2F(Form("%s_htoptag_%s",prefix,HypothesisTypeName(i)),
			   "Top tagging for final selection",4,0,4,4,0,4);
    htoptag[i]->Sumw2();



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

    hmaxGenJetPt[i]->Sumw2();
    hmaxJPTPt[i]->Sumw2();
    hmaxCaloJetPt[i]->Sumw2();
    hmaxTrkJetPt[i]->Sumw2();
    hmaxPFJetPt[i]->Sumw2();
    hmaxCaloTrkJetPt[i]->Sumw2();
    hmaxCaloTrkJet2Et[i]->Sumw2();

    hmaxGenJetPtHCal[i]->Sumw2();
    hmaxJPTPtHCal[i]->Sumw2();
    hmaxCaloJetPtHCal[i]->Sumw2();
    hmaxTrkJetPtHCal[i]->Sumw2();
    hmaxPFJetPtHCal[i]->Sumw2();
    hmaxCaloTrkJetPtHCal[i]->Sumw2();
    hmaxCaloTrkJet2EtHCal[i]->Sumw2();

    hmaxGenJetPtHF[i]->Sumw2();
    hmaxJPTPtHF[i]->Sumw2();
    hmaxCaloJetPtHF[i]->Sumw2();
    hmaxTrkJetPtHF[i]->Sumw2();
    hmaxPFJetPtHF[i]->Sumw2();
    hmaxCaloTrkJetPtHF[i]->Sumw2();
    hmaxCaloTrkJet2EtHF[i]->Sumw2();

    hPFJetResponse[i]->Sumw2();
    hPFJetRelResponse[i]->Sumw2();
    hPFJetResponseWithZero[i]->Sumw2();
    hPFJetRelResponseWithZero[i]->Sumw2();

    hmetVal[i]->Sumw2();
    hmetProjVal[i]->Sumw2();
    hmaxPFJetPtVal[i]->Sumw2();
    hdilMassVal[i]->Sumw2();

    hdilMassWithMetDYEst[i]->Sumw2();
    hdilMassNoMetDYEst[i]->Sumw2();
    hmetInDYEst[i]->Sumw2();
    hmetOutDYEst[i]->Sumw2();
   
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
  hmaxBtagVsJPTPt = new TH2F(Form("%s_hmaxBtagVsJPTPt",prefix),   Form("%s - most energetic jet Et (JPT) vs b-tagger output",prefix), 20, 0., 100, 10, 0, 1);
  hmaxBtagVsJPTPt->Sumw2();
  hCentralBquarkEtaAfterVeto = new TH1F(Form("%s_centralBQuarkEtaAfterVeto",prefix), Form("%s - central b quark eta distribution after jet veto",prefix), 20, 0, 10);
  hForwardBquarkEtaAfterVeto = new TH1F(Form("%s_forwardBQuarkEtaAfterVeto",prefix), Form("%s - forward b quark eta distribution after jet veto",prefix), 20, 0, 10);

  if (qcdBackground) {
    hFakableRateSingleElectron = new TH2F(Form("%s_hFakableRateSingleElectron",prefix), "FR study: rate of fakable objects", 3,ptbins,2,0,3.0);
    hFakableRateSingleElectron->Sumw2();
    hFinalRateSingleElectron   = new TH2F(Form("%s_hFinalRateSingleElectron",  prefix), "FR study: rate of final objects", 3,ptbins,2,0,3.0);
    hFinalRateSingleElectron->Sumw2();
    hFakableRateSingleMuon = new TH2F(Form("%s_hFakableRateSingleMuon",prefix), "FR study: rate of fakable objects", 3,ptbins,2,0,3.0);
    hFakableRateSingleMuon->Sumw2();
    hFinalRateSingleMuon   = new TH2F(Form("%s_hFinalRateSingleMuon",  prefix), "FR study: rate of final objects", 3,ptbins,2,0,3.0);
    hFinalRateSingleMuon->Sumw2();
    hIsoSingleMuon             = new TH1F(Form("%s_hIsoSingleMuon",prefix),             "Muon isolation distribution for goodMuonWithoutIsolation", 100,0,10);
    hIsoSingleMuon->Sumw2();
    hIsoSingleElectron         = new TH1F(Form("%s_hIsoSingleElectron",prefix),         "Electron isolation distribution for goodElectronWithoutIsolation", 100,0,10);
    hIsoSingleElectron->Sumw2();
  }

  hPFJetResidualResponse = new TH2F(Form("%s_hPFJetResidualResponse",prefix), Form("%s - PFJet Residual Correction",prefix), 20, 10, 50, 20, 0, 2);
  hPFJetResidualResponse->Sumw2();
  hnGoodVertex      = new TH1F(Form("%s_hnGoodVertex",     prefix), "Number of good vertexes after baseline selections" , 50,0.,50.);	
  hnGoodVertex -> Sumw2();
}



RooDataSet* ScanChain( TChain* chain, 
		       enum Sample sample, 
		       double integratedLumi, // in unit of pb^-1, if negative the weight is 1.
		       double xsec,           // in unit of pb, if negative take it from evt_xsec_excl*evt_kfactor
		       int nProcessedEvents,  // if negative, take it from evt_nEvts
		       bool identifyEvents, 
		       bool qcdBackground,
		       bool zStudy,
		       bool realData,
		       TString cms2_json_file)
{
  const char *prefix = SampleName(sample);
  if ( chain->GetListOfFiles()->GetEntries()==0 ){
    printf("\nERROR: chain is empty for sample: %s\n\n",prefix);
    assert(0);
  }
  // chain->SetParallelUnzip(kTRUE);
  // gErrorIgnoreLevel = 3000; // suppress warnings about missing dictionaries 
  unsigned int nEventsChain = chain->GetEntries();  // number of entries in chain --> number of events from all files
  gErrorIgnoreLevel = -1;
  unsigned int nEventsTotal = 0;

 // declare and create array of histograms
  ofstream selectedEvents(Form("%s.list",prefix));
  std::map<unsigned int, std::set<unsigned int> > runList;
  RooDataSet* dataset = MakeNewDataset(prefix);

  initializeHistograms(prefix,qcdBackground);

  // clear list of duplicates
  already_seen.clear();
  int duplicates_total_n = 0;
  double duplicates_total_weight = 0;
  int nFailedIdentification = 0;
  int nFilteredOut = 0;

  int i_permille_old = 0;
  
  try {
    jetcorr_filenames_jpt.clear();
    jetcorr_filenames_jpt.push_back("files/Spring10_L2Relative_AK5JPT.txt");
    jetcorr_filenames_jpt.push_back("files/Spring10_L3Absolute_AK5JPT.txt");
    if(realData)
      jetcorr_filenames_jpt.push_back("files/Spring10DataV2_L2L3Residual_AK5JPT.txt");
    jet_corrector_jpt= makeJetCorrector(jetcorr_filenames_jpt);
    
    jetcorr_filenames_pf.clear();
    jetcorr_filenames_pf.push_back("files/Spring10_L2Relative_AK5PF.txt");
    jetcorr_filenames_pf.push_back("files/Spring10_L3Absolute_AK5PF.txt");
    if(realData)
      jetcorr_filenames_pf.push_back("files/Spring10DataV2_L2L3Residual_AK5PF.txt");
    jet_corrector_pf= makeJetCorrector(jetcorr_filenames_pf);
    
    jetcorr_filenames_calo.clear();
    jetcorr_filenames_calo.push_back("files/Spring10_L2Relative_AK5Calo.txt");
    jetcorr_filenames_calo.push_back("files/Spring10_L3Absolute_AK5Calo.txt");
    if(realData)
      jetcorr_filenames_calo.push_back("files/Spring10DataV2_L2L3Residual_AK5Calo.txt");
    jet_corrector_calo= makeJetCorrector(jetcorr_filenames_calo);
    
    jetcorr_filenames_trk.clear();
    jetcorr_filenames_trk.push_back("files/Spring10_L2Relative_AK5TRK.txt");
    jetcorr_filenames_trk.push_back("files/Spring10_L3Absolute_AK5TRK.txt");
    jet_corrector_trk= makeJetCorrector(jetcorr_filenames_trk);
  } catch (...){
    cout << "\nFailed to setup correctors needed to get Jet Enetry Scale. Abort\n" << endl;
    assert(0);
  }
  if(realData) {
    if(cms2_json_file=="") {  
      std::cout<<"\n WARNING: Running on the real Data, but not JSON file!"<<std::endl;
      // assert(0);
    } else {
      set_goodrun_file_json(cms2_json_file);
    }
  }
  // file loop
  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  monitor.counters.clear();
  while (TChainElement *currentFile = (TChainElement*)fileIter.Next()) {
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
      // Select the good runs from the json file
      if(realData && cms2_json_file!="") {
	if( !goodrun(cms2.evt_run(), cms2.evt_lumiBlock()) ) continue;
      }
      runList[cms2.evt_run()].insert(cms2.evt_lumiBlock());
      double weight = 1.0;
      if ( !realData && integratedLumi>0 ){
	double mcweight = cms2.genps_weight() > 0.0 ? 1.0 : -1.0;
	   weight = integratedLumi * mcweight * (xsec>0?xsec:cms2.evt_xsec_excl()*cms2.evt_kfactor()) /
	     (nProcessedEvents>0?nProcessedEvents:cms2.evt_nEvts());
	 }       
	 ++nEventsTotal;
	 if (qcdBackground) {
	   // get fake rates
	   extractFakeRateSingleLepton();
	   // isolation
	   extractIsoSingleLepton();
	 }

	 if (cms2.trks_d0().size() == 0) continue;  // needed to get rid of back Monte Carlo events in CMSSW_2_X analysis
	 if (cms2.hyp_type().size() == 0) continue; // skip events without hypothesis
	 EventIdentifier id = { cms2.evt_run(), cms2.evt_event(), cms2.evt_lumiBlock(), cms2.trks_d0()[0], cms2.hyp_lt_p4()[0].pt(), cms2.hyp_lt_p4()[0].eta(), cms2.hyp_lt_p4()[0].phi() };
	 if (is_duplicate(id)) {
	   duplicates_total_n++;
	   if(!realData) duplicates_total_weight += cms2.evt_scale1fb();
	   //cout << "Duplicate event found. Run: " << cms2.evt_run() << ", Event:" << cms2.evt_event() << ", Lumi: " << cms2.evt_lumiBlock() << endl;
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
	 bool goodEvent = false;
	 // find the best candidate with m(ll) closest to the Z mass
	 unsigned int i_hyp_bestZ = bestZHyp();
	 for( unsigned int i_hyp = 0; i_hyp < nHyps; ++i_hyp ) {
	   if(cms2.hyp_p4().at(i_hyp).mass2() < 0 ) break;
	   if(zStudy && (i_hyp != i_hyp_bestZ)) continue;
	   if (hypo(i_hyp, weight, dataset, zStudy, realData))goodEvent=true;
	   AddIsoSignalControlSample(i_hyp, weight, dataset, realData);
	 }
	 if (goodEvent) selectedEvents << cms2.evt_run() << " " <<
	   cms2.evt_lumiBlock() << " " << cms2.evt_event() << endl;
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
  selectedEvents.close();
  if (realData){
      ofstream json("processed.json");
      ofstream json2("processed_detailed.json");
      json << "{";
      json2 << "{";
      bool firstRun(true);
      for ( std::map<unsigned int, std::set<unsigned int> >::const_iterator run = runList.begin();
	    run != runList.end(); ++ run ){
	if ( ! firstRun ) { 
	  json << " ,";
	  json2 << " ,";
	}
	firstRun = false;
	json << '"' << run->first << '"' << ": [[1,9999]]";
	const std::set<unsigned int>& lumis(run->second);
	unsigned int firstLumi = 0;
	unsigned int nLumis = 0;
	json2 << '"' << run->first << '"' << ": [";
	unsigned int nRanges = 0;
	for ( std::set<unsigned int>::const_iterator lumi=lumis.begin();
	      lumi!=lumis.end(); ++lumi )
	  { 
	    if ( firstLumi == 0 ){
	      firstLumi = *lumi;
	      nLumis = 1;
	      continue;
	    }
	    if ( firstLumi+nLumis == *lumi ){
	      nLumis++;
	      continue;
	    }
	    if (nRanges>0) json2 << " ,";
	    json2 << "[" << firstLumi << ", " << firstLumi+nLumis-1 << "]";
	    ++nRanges;
	    firstLumi = *lumi;
	    nLumis = 1;
	  }
	if ( firstLumi != 0 ) {
	  if (nRanges>0) json2 << " ,";
	  json2 << "[" << firstLumi << ", " << firstLumi+nLumis-1 << "]";
	}
	json2 << "]";
      }
      json << "}";
      json2 << "}";
      json << endl;
      json2 << endl;
      json.close();
      json2.close();
  }
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
		    int nProcessedEvents,
		    RooDataSet* output_dataset, 
		    Color_t color, 
		    bool identifyEvents,
		    bool qcdBackground,
		    bool zStudy,
		    bool realData,
		    TString cms2_json_file)
{
  std::vector<string> vec;
  vec.push_back(file_pattern);
  ProcessSample(vec,sample,integratedLumi,xsec,nProcessedEvents,output_dataset,color,identifyEvents,qcdBackground,zStudy,realData,cms2_json_file);
}

void ProcessSample( std::vector<std::string> file_patterns, 
		    Sample sample, 
		    double integratedLumi,
		    double xsec,
		    int nProcessedEvents,
		    RooDataSet* output_dataset, 
		    Color_t color, 
		    bool identifyEvents,
		    bool qcdBackground,
		    bool zStudy,
		    bool realData,
		    TString cms2_json_file)
{
  TChain *tchain = new TChain("Events");
  for ( std::vector<std::string>::const_iterator pattern = file_patterns.begin();
	pattern != file_patterns.end(); ++pattern )
    tchain->Add(pattern->c_str());
  if ( gSystem->Getenv("SkimSamples") ){
    std::cout << "Skimming " << SampleName(sample) << " with skim name " << 
      gSystem->Getenv("SkimName") << " .." << std::endl;
    if ( gSystem->Getenv("Merge") )
      SkimChain(tchain,true);
    else 
      SkimChain(tchain,false);
  } else {
    std::cout << "Processing " << SampleName(sample) << ".." << std::endl;
    RooDataSet* data = ScanChain(tchain,sample,integratedLumi,xsec,nProcessedEvents,identifyEvents,qcdBackground,zStudy,realData,cms2_json_file);
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

struct FileEntry{
  std::string inputDir;
  std::string outputDir;
  std::string fileName;
  unsigned int nIn;
  unsigned int nOut;
  bool operator<(const FileEntry& rhs) const { return inputDir < rhs.inputDir; }
};

bool passedSkimSelection2()
{
  unsigned int nHyps = cms2.hyp_type().size();
  for( unsigned int i_hyp = 0; i_hyp < nHyps; ++i_hyp ) {
    const double min1 = 20;
    const double min2 = 20;
    if ( !( cms2.hyp_lt_p4().at(i_hyp).pt() > min1 && cms2.hyp_ll_p4().at(i_hyp).pt() > min2 ) &&
	 !( cms2.hyp_lt_p4().at(i_hyp).pt() > min2 && cms2.hyp_ll_p4().at(i_hyp).pt() > min1 ) ) continue;

    // charge
    if ( cms2.hyp_lt_id()[i_hyp] * cms2.hyp_ll_id()[i_hyp] > 0 ) continue;
    
    if (nGoodVertex()<1) return false;

    if (cms2.hyp_p4().at(i_hyp).mass2()<0 || 
	cms2.hyp_p4()[i_hyp].mass() < 12) return false;
    if (inZmassWindow(cms2.hyp_p4()[i_hyp].mass())) return false;
    if (!passedMetRequirements(i_hyp)) return false;
    if (numberOfJets(i_hyp)>0) return false;
    if (numberOfSoftMuons(i_hyp,true)>0) return false;
    if (numberOfExtraLeptons(i_hyp,10)>0) return false;
    if ( toptag(CaloJet,i_hyp,0) ) return false;

    // id & iso
    if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && !goodMuonIsolated(cms2.hyp_lt_index()[i_hyp]) ) continue;
    if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && !goodMuonIsolated(cms2.hyp_ll_index()[i_hyp]) ) continue;
    if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11 && !goodElectronIsolated(cms2.hyp_lt_index()[i_hyp]) ) continue;
    if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11 && !goodElectronIsolated(cms2.hyp_ll_index()[i_hyp]) ) continue;

    return true;
  }
  return false;
}

void SkimChain(TChain* chain,bool mergeFiles){
  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  TPRegexp fileNameMatchPattern("(.*?)/([^/]*)$");
  TPRegexp dataDirMatchPattern("^data");
  if ( gSystem->Getenv("DataDir") ) dataDirMatchPattern = TPRegexp(gSystem->Getenv("DataDir"));
  unsigned int nEventsTotal = 0;
  unsigned int nEventsSelectedTotal = 0;
  std::vector<FileEntry> files;
  while (TChainElement *currentFile = (TChainElement*)fileIter.Next()) {
    FileEntry file;
    // cout << "file: " << currentFile->GetTitle() << endl;
    TObjArray* matches = fileNameMatchPattern.MatchS(currentFile->GetTitle());
    assert(matches);
    // cout << "matches->GetLast(): " << matches->GetLast() << endl;
    assert(matches && matches->GetLast()==2);
    TString inputFileName(currentFile->GetTitle());
    TString fileName(((TObjString*)matches->At(2))->GetString());
    file.fileName = fileName.Data();
    TString outputDirectory(((TObjString*)matches->At(1))->GetString());
    file.inputDir = outputDirectory.Data();
    if ( gSystem->Getenv("SkimDir") ) 
      dataDirMatchPattern.Substitute(outputDirectory,gSystem->Getenv("SkimDir"));
    assert(gSystem->Getenv("SkimName"));
    TString skimname = gSystem->Getenv("SkimName");
    outputDirectory += "/" + skimname;
    file.outputDir = outputDirectory.Data();
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
    // tree->SetBranchStatus("EventAuxiliary",0);
    TTree *newtree = tree->CloneTree(0);
    newtree->SetDirectory(output);
    
    cms2.Init(newtree);
    cms2.Init(tree);
    
    // Event Loop
    const unsigned int nEvents = tree->GetEntries();
    unsigned int nEventsSelected = 0;
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
    file.nIn = nEvents;
    file.nOut = nEventsSelected;
    files.push_back(file);
    nEventsSelectedTotal += nEventsSelected;
    output->cd();
    newtree->Write();
    output->Close();
    input->Close();
  }
  cout << Form("Processed events: %u, \tselected: %u\n",nEventsTotal,nEventsSelectedTotal) << endl;
  if (mergeFiles){
    // split files in sets with the same directory name
    sort(files.begin(),files.end());
    std::vector<std::vector<FileEntry> > sets;
    std::string previousDirectory;
    for(std::vector<FileEntry>::const_iterator file = files.begin();
	file != files.end(); ++ file){
      if ( file->outputDir != previousDirectory ){
	sets.push_back(std::vector<FileEntry>());
      }
      previousDirectory = file->outputDir;
      sets.back().push_back(*file);
    }
    // merge sets
    for ( std::vector<std::vector<FileEntry> >::const_iterator aSet = sets.begin();
	  aSet != sets.end(); ++aSet )
      {
	if (aSet->empty()) continue;
	std::ofstream log((aSet->front().outputDir+"/merged.log").c_str());
	if (aSet->size()==1){
	  // simply rename the output file
	  log << aSet->front().inputDir + "/" + aSet->front().fileName <<
	    " \t" << aSet->front().nIn << " \t" << aSet->front().nOut << std::endl;
	  assert(::system(Form("mv -v %s/%s %s/merged.root", aSet->front().outputDir.c_str(),
			       aSet->front().fileName.c_str(), aSet->front().outputDir.c_str()))==0);
	  log.close();
	  continue;
	}
	/*
	std::string mergeCommand(Form("hadd -f %s/merged.root",aSet->front().outputDir.c_str()));
	std::string rmCommand("rm");
	for ( std::vector<FileEntry>::const_iterator file = aSet->begin(); 
	      file != aSet->end(); ++file ){
	  log << file->inputDir + "/" + file->fileName <<
	    " \t" << file->nIn << " \t" << file->nOut << std::endl;
	  mergeCommand += " " + file->outputDir + "/" + file->fileName;
	  rmCommand += " " + file->outputDir + "/" + file->fileName;
	}
	cout << "Merging:\n" << mergeCommand << endl;
	log.close();
	assert(::system(mergeCommand.c_str())==0);
	cout << "Merging is done, remove intermediate files." << endl;
	assert(::system(rmCommand.c_str())==0);
	*/
	cout << "Merging " << aSet->front().outputDir << endl;
	std::string rmCommand("rm");
	TChain* aChain = new TChain("Events");
	for ( std::vector<FileEntry>::const_iterator file = aSet->begin(); 
	      file != aSet->end(); ++file ){
	  log << file->inputDir + "/" + file->fileName <<
	    " \t" << file->nIn << " \t" << file->nOut << std::endl;
	  aChain->Add((file->outputDir + "/" + file->fileName).c_str());
	  rmCommand += " " + file->outputDir + "/" + file->fileName;
	}
	aChain->Merge((aSet->front().outputDir + "/merged.root").c_str());
	cout << "Merging is done, remove intermediate files.\n" << endl;
	assert(::system(rmCommand.c_str())==0);
      }
  }
}

bool passedSkimSelection()
{
  unsigned int nHyps = cms2.hyp_type().size();
  for( unsigned int i_hyp = 0; i_hyp < nHyps; ++i_hyp ) {
    const double min1 = 20;
    const double min2 = 10;
    if ( !( cms2.hyp_lt_p4().at(i_hyp).pt() > min1 && cms2.hyp_ll_p4().at(i_hyp).pt() > min2 ) &&
	 !( cms2.hyp_lt_p4().at(i_hyp).pt() > min2 && cms2.hyp_ll_p4().at(i_hyp).pt() > min1 ) ) continue;

    
    // charge
    if ( cms2.hyp_lt_id()[i_hyp] * cms2.hyp_ll_id()[i_hyp] > 0 ) continue;
    
    // met
    if ( cms2.evt_tcmet() < 20 &&
	 cms2.evt_pfmet() < 20 ) continue;
    
    /*
    // id & iso
    if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && !goodMuonIsolated(cms2.hyp_lt_index()[i_hyp]) ) continue;
    if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && !goodMuonIsolated(cms2.hyp_ll_index()[i_hyp]) ) continue;
    if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11 && !goodElectronIsolated(cms2.hyp_lt_index()[i_hyp]) ) continue;
    if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11 && !goodElectronIsolated(cms2.hyp_ll_index()[i_hyp]) ) continue;
    */

    return true;
  }
  return false;
}


bool CheckCutsNM1(cuts_t apply, cuts_t remove, cuts_t passed)
{           
  if ((passed & (apply & (~remove))) == (apply & (~remove))) return true;
  return false;
}   

bool CheckCuts(cuts_t apply, cuts_t passed)
{
  if ((apply & passed) == apply) return true;
  return false;
}

unsigned int bestZHyp()
{
  unsigned int nHyps = cms2.hyp_type().size();
  if(nHyps<2) return 0;
  unsigned int i_hyp_bestZ = 0;
  double minDev = 20;
  for( unsigned int i_hyp = 0; i_hyp < nHyps; ++i_hyp ) {
    if(fabs(cms2.hyp_p4().at(i_hyp).mass2()-91.1876) < minDev) {
      minDev = fabs(cms2.hyp_p4().at(i_hyp).mass2()-91.1876);
      i_hyp_bestZ = i_hyp;
    }
  }
  return i_hyp_bestZ;
}

void getJetResponseFromZBalance(int i_hyp,  double weight, bool realData, double etaMin, double etaMax)
{ 
  HypothesisType type = getHypothesisType(cms2.hyp_type()[i_hyp]);
  double vetoCone = 0.3;
  double jec = 1.0;
  double pfJetMax(0.);
  int pfJetMaxIndex (-1);
  
  for ( unsigned int i=0; i < cms2.pfjets_p4().size(); ++i) {
    if ( TMath::Abs(cms2.pfjets_p4()[i].eta()) > etaMax ) continue;
    if ( TMath::Abs(cms2.pfjets_p4()[i].eta()) < etaMin ) continue;
    if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.pfjets_p4()[i])) < vetoCone ||
	 TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.pfjets_p4()[i])) < vetoCone ) continue;
    if(applyJEC)
      jec = jetCorrection(cms2.pfjets_p4()[i], jet_corrector_pf);
    if ( cms2.pfjets_p4()[i].pt() * jec < pfJetMax ) continue;
    pfJetMax = cms2.pfjets_p4()[i].pt() * jec;
    pfJetMaxIndex = i;
    
    if(realData) 
      hPFJetResidualResponse->Fill( cms2.pfjets_p4()[i].pt()*cms2.pfjets_cor()[i], jetCorrection(cms2.pfjets_p4()[i], jet_corrector_pf)/ cms2.pfjets_cor()[i], weight);
  }

  if(pfJetMaxIndex == -1) return;

  // apply the Z+1jet selection
  bool zPlusOneJet = true;
  double secJetMax = 10;
 
  // require the leading jet to be back-to-back with the Z
  if ( TMath::ACos(-1.0* TMath::Cos(cms2.hyp_p4()[i_hyp].phi() - cms2.pfjets_p4()[pfJetMaxIndex].phi())) > 0.2 ) zPlusOneJet = false;
  // require no second jet above secJetMax
  for ( int i=0; i < (int) cms2.pfjets_p4().size(); ++i) {
    if(i == pfJetMaxIndex) continue;
    jec = 1.0;
    if(applyJEC)
      jec = jetCorrection(cms2.pfjets_p4()[i], jet_corrector_pf);
    if ( TMath::Abs(cms2.pfjets_p4()[i].eta()) > etaMax ) continue;
    if ( TMath::Abs(cms2.pfjets_p4()[i].eta()) < etaMin ) continue;
    if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.pfjets_p4()[i])) < vetoCone ||
	 TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.pfjets_p4()[i])) < vetoCone ) continue;
    if (cms2.pfjets_p4()[i].pt() * jec > TMath::Max(secJetMax, cms2.pfjets_p4()[i].pt()*jec*0.2) ) {
      zPlusOneJet = false;
    }
  }

  if (zPlusOneJet) {
    jec = 1.0;
    if(applyJEC)
      jec =  jetCorrection(cms2.pfjets_p4()[pfJetMaxIndex], jet_corrector_pf);
    
    hPFJetResponse[type]->Fill( cms2.hyp_p4()[i_hyp].pt(), cms2.pfjets_p4()[pfJetMaxIndex].pt()*jec, weight);
    hPFJetResponse[3]->Fill( cms2.hyp_p4()[i_hyp].pt(), cms2.pfjets_p4()[pfJetMaxIndex].pt()*jec, weight);
    
    hPFJetRelResponse[type]->Fill( cms2.hyp_p4()[i_hyp].pt(), cms2.pfjets_p4()[pfJetMaxIndex].pt()*jec/cms2.hyp_p4()[i_hyp].pt(), weight);
    hPFJetRelResponse[3]->Fill( cms2.hyp_p4()[i_hyp].pt(), cms2.pfjets_p4()[pfJetMaxIndex].pt()*jec/cms2.hyp_p4()[i_hyp].pt(), weight);
    
    hPFJetResponseWithZero[type]->Fill( cms2.hyp_p4()[i_hyp].pt(), cms2.pfjets_p4()[pfJetMaxIndex].pt()*jec, weight);
    hPFJetResponseWithZero[3]->Fill( cms2.hyp_p4()[i_hyp].pt(), cms2.pfjets_p4()[pfJetMaxIndex].pt()*jec, weight);
    
    hPFJetRelResponseWithZero[type]->Fill( cms2.hyp_p4()[i_hyp].pt(), cms2.pfjets_p4()[pfJetMaxIndex].pt()*jec/cms2.hyp_p4()[i_hyp].pt(), weight);
    hPFJetRelResponseWithZero[3]->Fill( cms2.hyp_p4()[i_hyp].pt(), cms2.pfjets_p4()[pfJetMaxIndex].pt()*jec/cms2.hyp_p4()[i_hyp].pt(), weight);
  }
  
  else {
    hPFJetResponseWithZero[type]->Fill( cms2.hyp_p4()[i_hyp].pt(), 0.0, weight);
    hPFJetResponseWithZero[3]->Fill( cms2.hyp_p4()[i_hyp].pt(), 0.0, weight);
    
    hPFJetRelResponseWithZero[type]->Fill( cms2.hyp_p4()[i_hyp].pt(), 0.0, weight);
    hPFJetRelResponseWithZero[3]->Fill( cms2.hyp_p4()[i_hyp].pt(), 0.0, weight);
  }
 }
				
void fill_dyest_histograms(int i_hyp, float weight)
{
  HypothesisType type = getHypothesisType(cms2.hyp_type()[i_hyp]);
  
  // fill the mass histogram
  float mass = cms2.hyp_p4()[i_hyp].mass(); 
  hdilMassNoMetDYEst[type] -> Fill(mass, weight); 
  hdilMassNoMetDYEst[3] -> Fill(mass, weight); 
  
  if (passedMetRequirements (i_hyp) )
    hdilMassWithMetDYEst[type] -> Fill(mass, weight); 
  
  // fill the met histograms for "in" and "out" regions
  if (inZmassWindow(mass)) {
    hmetInDYEst[type] -> Fill(projectedMet(i_hyp), weight); 
    hmetInDYEst[3] -> Fill(projectedMet(i_hyp), weight); 
  }
  else {
    hmetOutDYEst[type] -> Fill(projectedMet(i_hyp), weight); 
    hmetOutDYEst[3] -> Fill(projectedMet(i_hyp), weight); 
  }
}
