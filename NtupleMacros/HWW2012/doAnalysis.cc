const char* config_info = "SmurfV6 selection (Baseline;Tight+Loose;MET20); 42X"; //Skim1
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
#include "Math/VectorUtil.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TPRegexp.h"
#include "monitor.h"
#include "../Tools/goodrun.h"
#include "../Tools/ElectronIDMVA.h"
#include "TTreeCache.h"
using namespace std;

#ifndef __CINT__
#include "CORE/CMS2.h"
#include "CORE/electronSelections.h"
#include "CORE/muonSelections.h"
#include "CORE/jetSelections.h"
#include "CORE/metSelections.h"
#include "jetcorr/FactorizedJetCorrector.h"
#endif

enum jetregion { HCAL, HF, ALLJET};

enum hyp_selection {
  PASSED_BaseLine                = 1UL<<0,
  PASSED_Charge                  = 1UL<<1,
  PASSED_ZVETO                   = 1UL<<2,
  PASSED_ZControlSampleVeryTight = 1UL<<3,   // within Z mass window +/- 5GeV
  PASSED_ZControlSampleTight     = 1UL<<4,   // within Z mass window +/- 10GeV
  PASSED_ZControlSampleLoose     = 1UL<<5,   // within Z mass window +/- 20GeV
  PASSED_MET                     = 1UL<<6,
  PASSED_LT_FINAL                = 1UL<<7,
  PASSED_LT_FO_MU1               = 1UL<<8,
  PASSED_LT_FO_MU2               = 1UL<<9,
  PASSED_LT_FO_ELEV1             = 1UL<<10,
  PASSED_LT_FO_ELEV2             = 1UL<<11,
  PASSED_LT_FO_ELEV3             = 1UL<<12,
  PASSED_LT_FO_ELEV4             = 1UL<<13,
  PASSED_LL_FINAL                = 1UL<<14,
  PASSED_LL_FO_MU1               = 1UL<<15,
  PASSED_LL_FO_MU2               = 1UL<<16,
  PASSED_LL_FO_ELEV1             = 1UL<<17,
  PASSED_LL_FO_ELEV2             = 1UL<<18,
  PASSED_LL_FO_ELEV3             = 1UL<<19,
  PASSED_LL_FO_ELEV4             = 1UL<<20,
  PASSED_JETVETO                 = 1UL<<21,
  PASSED_TopControlSample        = 1UL<<22,  // 2 or more jets
  PASSED_1BJET                   = 1UL<<23,
  PASSED_SOFTMUVETO_NotInJets    = 1UL<<24,
  PASSED_SOFTMUVETO              = 1UL<<25,
  PASSED_EXTRALEPTONVETO         = 1UL<<26,  
  PASSED_TOPVETO_NotInJets       = 1UL<<27,  // exclude jets over threshold from top tagging
  PASSED_TOPVETO                 = 1UL<<28,
  PASSED_Skim1                   = 1UL<<29,  // one fakable object + one final; full met
  PASSED_Trigger                 = 1UL<<30,
  PASSED_Skim3                   = 1UL<<31   // one fakable object + one final
};

// DEFAULT
//wwcuts_t pass_all = PASSED_BaseLine | PASSED_Charge | PASSED_ZVETO | PASSED_MET | PASSED_JETVETO | PASSED_LT_FINAL | PASSED_LL_FINAL | PASSED_SOFTMUVETO | PASSED_EXTRALEPTONVETO | PASSED_TOPVETO;

// wwcuts_t pass_all = PASSED_Skim1;
wwcuts_t pass_all = PASSED_Skim3; // Baseline, Tight+Fakeable, no MET requirement

//wwcuts_t pass_all = PASSED_BaseLine | PASSED_Charge | PASSED_ZVETO | PASSED_MET | PASSED_LT_FINAL | PASSED_LL_FINAL | PASSED_TopControlSample;
// wwcuts_t pass_all = PASSED_BaseLine;

// wwcuts_t pass_all = PASSED_BaseLine | PASSED_Charge | PASSED_LT_FINAL | PASSED_LL_FINAL | PASSED_ZControlSampleTight;

bool applyJEC = true;
bool applyFastJetCorrection = false;
bool lockToCoreSelectors = false;
bool selectBestCandidate = true; // select only one hypothesis per event with the two most energetic leptons
bool useLHeleId = false;
bool useMVAeleId = true;
bool doDYNNLOw = true;
const unsigned int prescale = 1; // DON'T USE ANYTHING BUT 1, unless you know what you are doing

std::vector<std::string> jetcorr_filenames_pfL1FastJetL2L3;
FactorizedJetCorrector *jet_corrector_pfL1FastJetL2L3;
wwcuts_t cuts_passed = 0;

//
// Key analysis method implementation
//

TH1D* HiggsPtKFactor = 0;
ElectronIDMVA* electronIdMVA = 0;
vector<TH2D*>     fDYNNLOKFactorHists;           //vector of hist for Drell-Yan NNLO Kfactor

bool goodElectronTMVA(unsigned int i) {
  //cout << "electronIdMVA.MVAValue=" << electronIdMVA->MVAValue(i, 0) << endl;
  float pt = cms2.els_p4().at(i).pt();
  float etaSC = cms2.els_etaSC().at(i);
  //preselection
  if (fabs(etaSC)<1.479) {
    if (cms2.els_sigmaIEtaIEta().at(i)>0.01 || 
	fabs(cms2.els_dEtaIn().at(i))>0.007 ||
	fabs(cms2.els_dPhiIn().at(i))>0.15 ||
	cms2.els_hOverE().at(i)>0.12 ||
	cms2.els_tkIso().at(i)/pt>0.2 ||
	TMath::Max(cms2.els_ecalIso().at(i) - 1.0, 0.0)/pt>0.20 ||
	cms2.els_hcalIso().at(i)/pt>0.20 ) return 0;
  } else {
    if (cms2.els_sigmaIEtaIEta().at(i)>0.03 || 
	fabs(cms2.els_dEtaIn().at(i))>0.009 ||
	fabs(cms2.els_dPhiIn().at(i))>0.10 ||
	cms2.els_hOverE().at(i)>0.10 ||
	cms2.els_tkIso().at(i)/pt>0.2 ||
	TMath::Max(cms2.els_ecalIso().at(i) - 1.0, 0.0)/pt>0.20 ||
	cms2.els_hcalIso().at(i)/pt>0.20 ) return 0;
  }
  //selection
  if (pt<20){
    if (fabs(etaSC)<1. && electronIdMVA->MVAValue(i, 0)>0.139) return 1;
    else if (fabs(etaSC)>=1. && fabs(etaSC)<1.479 && electronIdMVA->MVAValue(i, 0)>0.525) return 1;
    else if (fabs(etaSC)>=1.479 && fabs(etaSC)<2.5 && electronIdMVA->MVAValue(i, 0)>0.543) return 1;
  } else {
    if (fabs(etaSC)<1. && electronIdMVA->MVAValue(i, 0)>0.947) return 1;
    else if (fabs(etaSC)>=1. && fabs(etaSC)<1.479 && electronIdMVA->MVAValue(i, 0)>0.950) return 1;
    else if (fabs(etaSC)>=1.479 && fabs(etaSC)<2.5 && electronIdMVA->MVAValue(i, 0)>0.884) return 1;
  }
  return 0;
}

bool goodElectronWithoutIsolation(unsigned int i){
  return ww_elBase(i) && ww_elId(i) && ww_eld0PV(i) && ww_eldZPV(i);
}

bool goodElectronIsolated(unsigned int i){
  bool ptcut = cms2.els_p4().at(i).pt() >= 10.0;
  bool core = ptcut && pass_electronSelection( i, electronSelection_smurfV6);
  bool internal = ww_elBase(i) && ww_elId(i) && ww_eld0PV(i) && ww_eldZPV(i) && ww_elIso(i);
  assert(!lockToCoreSelectors || core==internal);
  return internal;
}

bool fakableElectron(unsigned int i, EleFOTypes type){
  if ( cms2.els_p4().at(i).pt() < 10.0 ) return false;
  switch (type){
  case EleFOV1: return pass_electronSelection( i, electronSelectionFO_el_smurf_v1);
  case EleFOV2: return pass_electronSelection( i, electronSelectionFO_el_smurf_v2);
  case EleFOV3: return pass_electronSelection( i, electronSelectionFO_el_smurf_v3);
  case EleFOV4: return pass_electronSelection( i, electronSelectionFO_el_smurf_v4);
  }
  return false;
}

bool goodMuonWithoutIsolation(unsigned int i){
  return ww_muBase(i) && ww_mud0PV(i) && ww_mudZPV(i) && ww_muId(i);
}

bool goodMuonIsolated(unsigned int i){
  bool ptcut = cms2.mus_p4().at(i).pt() >= 10.0;
  bool core = ptcut && muonId(i, NominalSmurfV6);
  bool internal = ww_muBase(i) && ww_mud0PV(i) && ww_mudZPV(i) && ww_muId(i) && ww_muIso(i); 
  assert(!lockToCoreSelectors || core==internal);
  return internal;
}

bool fakableMuon(unsigned int i, MuFOTypes type){
  if ( cms2.mus_p4().at(i).pt() < 10.0 ) return false;
  switch (type){
  case MuFOV1: return muonId(i, muonSelectionFO_mu_smurf_10);
  case MuFOV2: return muonId(i, muonSelectionFO_mu_smurf_04);
  }
  return false;
}

WWJetType jetType(){
  return pfJet;
  // return jptJet;
}

std::vector<JetPair> getDefaultJets(unsigned int i_hyp, bool btagged=false){
  return getJets(jetType(), i_hyp, 30, 5.0, false, btagged); // V1
}

unsigned int numberOfJets(unsigned int i_hyp){
  return getDefaultJets(i_hyp, false).size(); 
}

double metValue(){    return cms2.evt_pfmet(); }
double metPhiValue(){ return cms2.evt_pfmetPhi(); }
double sumetValue(){    return cms2.evt_pfsumet(); }

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

int primaryVertex(){
  //  double sumPtMax = -1;
  //   int iMax = -1;
  //   for ( unsigned int i = 0; i < cms2.vtxs_sumpt().size(); ++i ){
  //     // if (cms2.vtxs_isFake()[i]) continue;
  //     if (!isGoodVertex(i)) continue;
  //     if ( cms2.vtxs_sumpt().at(i) > sumPtMax ){
  //       iMax = i;
  //       sumPtMax = cms2.vtxs_sumpt().at(i);
  //     }
  //   }
  //   if (iMax<0) return false;
  return 0;
}

bool passedMetRequirements(unsigned int i_hyp){
  // if ( cms2.hyp_p4().at(i_hyp).mass()>130 ) return true;
  HypothesisType type = getHypothesisTypeNew(i_hyp);
  // std::vector<LorentzVector> jets = getDefaultJets(i_hyp);
  metStruct trkMET = trackerMET(i_hyp,0.1); //,&jets);
  double pMet = std::min(projectedMet(i_hyp, metValue(), metPhiValue()),
			 projectedMet(i_hyp, trkMET.met, trkMET.metphi));
  // if ( type == EM && cms2.hyp_p4().at(i_hyp).mass()>90 ) return true;
  if ( pMet < 20 ) return false;
  if (type == EE || type == MM) {
    double threshold = 37+nGoodVertex()/2.0;
    // double dmass = fabs(cms2.hyp_p4()[i_hyp].mass()-91);
    // if ( metValue() < 45 ) return false;
    if ( pMet < threshold ) return false;
  }
  return true;
}


//
// Electron Id
//

bool ww_elBase(unsigned int index){
  if (cms2.els_p4().at(index).pt() < 10.0) return false;
  if (fabs(cms2.els_p4().at(index).eta()) > 2.5) return false;
  return true;
}
bool ww_elId(unsigned int index){
  // if( fabs(cms2.els_conv_dist().at(index)) < 0.02 &&
  //     fabs(cms2.els_conv_dcot().at(index)) < 0.02) return false;
  // if (! (electronId_VBTF(index, VBTF_35X_80) & (1<<ELEID_ID)) ) return false;
  // if (! (electronId_VBTF(index, VBTF_35X_70) & (1<<ELEID_ID)) ) return false;
  // if (! (electronId_CIC(index, 4, CIC_SUPERTIGHT) & (1<<ELEID_ID)) ) return false;

  if (useLHeleId) {
    if (cms2.els_p4().at(index).pt()>20 && (passLikelihoodId_v2(index,cms2.els_lh().at(index),0) & (1<<ELEID_ID))!=(1<<ELEID_ID) ) return false; 
    if (cms2.els_p4().at(index).pt()<20 && (passLikelihoodId_v2(index,cms2.els_lh().at(index),0) & (1<<ELEID_ID))!=(1<<ELEID_ID) ) return false;
  }
  if (useMVAeleId){
    if (!goodElectronTMVA(index)) return false;
  } else {
    if (!pass_electronSelection(index, electronSelection_smurfV3_id, false, false) ) return false;
  }
  
  // MIT conversion
  if ( isFromConversionMIT(index) ) return false;

  // conversion rejection - hit based
  if ( cms2.els_exp_innerlayers().at(index) > 0 ) return false;
  // MIT conversion
  // if (! pass_electronSelection(index, (1ll<<ELENOTCONV_MIT), false, false) ) return false;
  // if ( cms2.els_exp_innerlayers39X().at(index) > 0 ) return false;
  //  int ctfIndex = cms2.els_trkidx().at(index);
  // if ( ctfIndex >=0 && 
  //     cms2.els_charge().at(index)!=cms2.trks_charge().at(ctfIndex) ) return false;
  // if ( !electronId_smurf_v2(index) ) return false;
  
  return true;
}
 
bool ww_eld0(unsigned int index){
  return fabs(cms2.els_d0corr()[index]) < 0.02;
}


double dzPV(const LorentzVector& vtx, const LorentzVector& p4, const LorentzVector& pv){
  return (vtx.z()-pv.z()) - ((vtx.x()-pv.x())*p4.x()+(vtx.y()-pv.y())*p4.y())/p4.pt() * p4.z()/p4.pt();
}

bool ww_eld0PV(unsigned int index){
  int vtxIndex = primaryVertex();
  if (vtxIndex<0) return false;
  double dxyPV = cms2.els_d0()[index]-
    cms2.vtxs_position()[vtxIndex].x()*sin(cms2.els_trk_p4()[index].phi())+
    cms2.vtxs_position()[vtxIndex].y()*cos(cms2.els_trk_p4()[index].phi());
  return fabs(dxyPV) < 0.02;
}

bool ww_eldZPV(unsigned int index){
  int vtxIndex = primaryVertex();
  if (vtxIndex<0) return false;
  // double dzPV = cms2.els_z0corr()[index]-cms2.vtxs_position()[iMax].z();
  double dzpv = dzPV(cms2.els_vertex_p4()[index], cms2.els_trk_p4()[index], cms2.vtxs_position()[vtxIndex]);
  return fabs(dzpv)<0.1;
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
  return pass_electronSelection(index, electronSelection_smurfV5_iso);
  //return ww_elIsoVal(index)<0.1;
}

//
// Muon Id
//

bool ww_muBase(unsigned int index){
  if (cms2.mus_p4().at(index).pt() < 10.0) return false;
  if (fabs(cms2.mus_p4().at(index).eta()) > 2.4) return false;
  if (cms2.mus_type().at(index) == 8) return false; // not STA
  return true;
}
bool ww_mud0(unsigned int index){
  return fabs(cms2.mus_d0corr()[index]) < 0.02;
}
double ww_mud0ValuePV(unsigned int index){
  int vtxIndex = primaryVertex();
  if (vtxIndex<0) return 9999;
  double dxyPV = cms2.mus_d0()[index]-
    cms2.vtxs_position()[vtxIndex].x()*sin(cms2.mus_trk_p4()[index].phi())+
    cms2.vtxs_position()[vtxIndex].y()*cos(cms2.mus_trk_p4()[index].phi());
  return fabs(dxyPV);
}

bool ww_mud0PV(unsigned int index){
  if ( cms2.mus_p4().at(index).pt() < 20. ) return ww_mud0ValuePV(index) < 0.01;
  return ww_mud0ValuePV(index) < 0.02;
}
bool ww_mudZPV(unsigned int index, float cut){
  int vtxIndex = primaryVertex();
  if (vtxIndex<0) return false;
  // double dzpv = cms2.mus_z0corr()[index]-cms2.vtxs_position()[iMax].z();
  double dzpv = dzPV(cms2.mus_vertex_p4()[index], cms2.mus_trk_p4()[index], cms2.vtxs_position()[vtxIndex]);
  return fabs(dzpv)<cut;
}
bool ww_muId(unsigned int index){ 
  if (((cms2.mus_type().at(index)) & (1<<2)) == 0)    return false; // tracker muon
  if (cms2.mus_validHits().at(index) < 11)            return false; // # of tracker hits
  // if (cms2.trks_nlayers().at(cms2.mus_trkidx().at(index))<=8) return false;
  if (cms2.mus_ptErr().at(index)/cms2.mus_p4().at(index).pt()>0.1) return false;
  if (cms2.trks_valid_pixelhits().at(cms2.mus_trkidx().at(index))==0) return false;
  if (cms2.mus_trkKink().at(index) > 20.) return false; //kink finder//newcuts
  // if (!isPFMuon(index))return false;
  // global muon
  bool goodMuonGlobalMuon = false;
  if (((cms2.mus_type().at(index)) & (1<<1)) == (1<<1)){
    goodMuonGlobalMuon = true;
    if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) goodMuonGlobalMuon = false; //glb fit chisq
    if (cms2.mus_gfit_validSTAHits().at(index)==0 ) goodMuonGlobalMuon = false;
    if (cms2.mus_nmatches().at(index)<2) goodMuonGlobalMuon = false;
  }
  return goodMuonGlobalMuon || 
    cms2.mus_pid_TMLastStationTight().at(index) == 1; // TM id
}

double ww_muIsoVal(unsigned int index){
  double sum =  cms2.mus_iso03_sumPt().at(index) +
    cms2.mus_iso03_emEt().at(index)  +
    cms2.mus_iso03_hadEt().at(index);
  double pt  = cms2.mus_p4().at(index).pt();
  return sum/pt;
}
bool ww_muIso(unsigned int index){
  if (cms2.mus_p4().at(index).pt()>20) {
    if (TMath::Abs(cms2.mus_p4()[index].eta())<1.479) 
      return muonIsoValuePF(index,0,0.3) < 0.13;
    else 
      return muonIsoValuePF(index,0,0.3) < 0.09;
  } else {
    if (TMath::Abs(cms2.mus_p4()[index].eta())<1.479) 
      return muonIsoValuePF(index,0,0.3) < 0.06;
    else 
      return muonIsoValuePF(index,0,0.3) < 0.05;
  }
//   if ( cms2.mus_p4().at(index).pt() < 20. )
//     return ww_muIsoVal(index)<0.1;
//   return ww_muIsoVal(index)<0.15;
}
unsigned int numberOfSoftMuons(int i_hyp, bool nonisolated,
			       const std::vector<JetPair>& vetojets)
{
  unsigned int nMuons = 0;
  for (int imu=0; imu < int(cms2.mus_charge().size()); ++imu) {
    // quality cuts
    if (  ((cms2.mus_goodmask()[imu]) & (1<<19)) == 0 ) continue; // TMLastStationAngTight
    if ( cms2.mus_p4()[imu].pt() < 3 ) continue;
    if ( ww_mud0ValuePV(imu) > 0.2) continue;
    if ( ! ww_mudZPV(imu,0.2) ) continue;//newcuts, was 0.1
    if ( cms2.mus_validHits()[imu] < 11) continue;
    if ( TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && cms2.hyp_lt_index()[i_hyp] == imu ) continue;
    if ( TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && cms2.hyp_ll_index()[i_hyp] == imu ) continue;
    if ( nonisolated && ww_muIsoVal(imu)<0.1 && cms2.mus_p4()[imu].pt()>20 ) continue;
    bool skip = false;
    for ( std::vector<JetPair>::const_iterator ijet = vetojets.begin();
	  ijet != vetojets.end(); ++ijet )
      if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(ijet->first,cms2.mus_p4()[imu])) < 0.3 ) skip=true;
    if ( skip ) continue;
    ++nMuons;
  }
  return nMuons;
}

typedef std::pair<bool, unsigned int> LeptonPair; // first bool(true-muon, false-electron) , second index

std::vector<LeptonPair> getExtraLeptons(int i_hyp, double minPt){
  std::vector<LeptonPair> leptons;
  for (int i=0; i < int(cms2.mus_charge().size()); ++i) {
    if ( cms2.mus_p4()[i].pt() < minPt ) continue;
    if ( TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && cms2.hyp_lt_index()[i_hyp] == i ) continue;
    if ( TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && cms2.hyp_ll_index()[i_hyp] == i ) continue;
    if ( ! (ww_mud0PV(i) && ww_muId(i) && ww_muIso(i)&&
	    fabs(cms2.mus_p4().at(i).eta()) <2.4) ) continue;
    leptons.push_back(LeptonPair(true,i));
  }
  for (int i=0; i < int(cms2.els_charge().size()); ++i) {
    if ( cms2.els_p4()[i].pt() < minPt ) continue;
    if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.els_p4().at(i)) <0.1) ) continue;
    if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.els_p4().at(i)) <0.1) ) continue;
    if ( !(ww_elId(i) && ww_eld0PV(i) && ww_elIso(i) && 
	   fabs(cms2.els_p4().at(i).eta()) < 2.5) ) continue;
    leptons.push_back(LeptonPair(false,i));
  }
  return leptons;
}

unsigned int numberOfExtraLeptons(int i_hyp, double minPt){
  return getExtraLeptons(i_hyp, minPt).size();
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
  // return true; // no trigger requirements
  // return cms2.filter_ele10mu10IsoId_passed();
  if ( passedTrigger("HLT_Mu17_Ele8_CaloIdL_v1") || 
       passedTrigger("HLT_Mu17_Ele8_CaloIdL_v2") ||
       passedTrigger("HLT_Mu17_Ele8_CaloIdL_v3") ||
       passedTrigger("HLT_Mu17_Ele8_CaloIdL_v4") ||
       passedTrigger("HLT_Mu17_Ele8_CaloIdL_v5") ||
       passedTrigger("HLT_Mu17_Ele8_CaloIdL_v6") ||
       passedTrigger("HLT_Mu17_Ele8_CaloIdL_v8") ||
       passedTrigger("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v1") ||
       passedTrigger("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v3") ||
       passedTrigger("HLT_Mu8_Ele17_CaloIdL_v1") ||
       passedTrigger("HLT_Mu8_Ele17_CaloIdL_v2") ||
       passedTrigger("HLT_Mu8_Ele17_CaloIdL_v3") ||
       passedTrigger("HLT_Mu8_Ele17_CaloIdL_v4") ||
       passedTrigger("HLT_Mu8_Ele17_CaloIdL_v5") || 
       passedTrigger("HLT_Mu8_Ele17_CaloIdL_v6") ||
       passedTrigger("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v1") ||
       passedTrigger("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v3") 
      ) return true;
  if ( passedTrigger("HLT_DoubleMu7_v1") ||
       passedTrigger("HLT_DoubleMu7_v2") ||
       passedTrigger("HLT_Mu13_Mu8_v2") ||
       passedTrigger("HLT_Mu13_Mu8_v3") ||
       passedTrigger("HLT_Mu13_Mu8_v4") ||
       passedTrigger("HLT_Mu13_Mu8_v6") ) return true;
  if ( passedTrigger("HLT_Mu15_v2") ||
       passedTrigger("HLT_Mu24_v1") ||
       passedTrigger("HLT_Mu24_v2") ||
       passedTrigger("HLT_Mu30_v1") ||
       passedTrigger("HLT_Mu30_v2") ||
       passedTrigger("HLT_Mu30_v3") ||
       passedTrigger("HLT_Mu30_v4") ||
       passedTrigger("HLT_Mu30_v5") ||
       passedTrigger("HLT_Mu30_v7") ||
       passedTrigger("HLT_IsoMu17_v5") ||
       passedTrigger("HLT_IsoMu17_v6") ||
       passedTrigger("HLT_IsoMu17_v8") ||
       passedTrigger("HLT_IsoMu17_v9") ||
       passedTrigger("HLT_IsoMu17_v10") ||
       passedTrigger("HLT_IsoMu17_v11") ||
       passedTrigger("HLT_IsoMu17_eta2p1_v1") ||
       passedTrigger("HLT_IsoMu20_eta2p1_v1") ||
       passedTrigger("HLT_IsoMu24_v1") ||
       passedTrigger("HLT_IsoMu24_v2") ||
       passedTrigger("HLT_IsoMu24_v4") ||
       passedTrigger("HLT_IsoMu24_v5") ||
       passedTrigger("HLT_IsoMu24_v6") ||
       passedTrigger("HLT_IsoMu24_v7") ||
       passedTrigger("HLT_IsoMu24_v8") ) return true;
  if ( passedTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1") ||
       passedTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2") ||
       passedTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3") ||
       passedTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v4") ||
       passedTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v5") ||
       passedTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v6") ||
       passedTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5") ||
       passedTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6") ||
       passedTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7") ) return true;
  if ( passedTrigger("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1") ||
       passedTrigger("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2") ||
       passedTrigger("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3") ||
       passedTrigger("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1") ||
       passedTrigger("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2") ||
       passedTrigger("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3") ||
       passedTrigger("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v4") ||
       passedTrigger("HLT_Ele42_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1") ||
       passedTrigger("HLT_Ele52_CaloIdVT_TrkIdT_v1") ||
       passedTrigger("HLT_Ele52_CaloIdVT_TrkIdT_v2") ||
       passedTrigger("HLT_Ele52_CaloIdVT_TrkIdT_v3") ||
       passedTrigger("HLT_Ele65_CaloIdVT_TrkIdT_v3") ) return true;

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

double projectedMet(unsigned int i_hyp, double met, double phi)
{
  double DeltaPhi = nearestDeltaPhi(phi,i_hyp);
  if (DeltaPhi < TMath::Pi()/2) return met*TMath::Sin(DeltaPhi);
  return met;
}

HypTypeInNtuples hypType(unsigned int i_hyp){
  HypTypeInNtuples type = HypTypeInNtuples(cms2.hyp_type().at(i_hyp));
  return type;
}

//
// Jets
//

Bool_t comparePt(JetPair lv1, JetPair lv2) {
   return lv1.first.pt() > lv2.first.pt();
}

std::vector<JetPair> 
getJets(WWJetType type, int i_hyp, double etThreshold, double etaMax, bool sortJets, bool btag)
{
  std::vector<JetPair> jets;
  const double vetoCone = 0.3;
  
  switch ( type ){
  case jptJet:
    for ( unsigned int i=0; i < cms2.jpts_p4().size(); ++i) {
      double jec = 1.0;
      if ( cms2.jpts_p4()[i].pt() * jec < etThreshold ) continue;
      if ( btag && !defaultBTag(type,i) ) continue;
      if ( TMath::Abs(cms2.jpts_p4()[i].eta()) > etaMax ) continue;
      if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.jpts_p4()[i])) < vetoCone ||
	   TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.jpts_p4()[i])) < vetoCone ) continue;
      jets.push_back(JetPair(cms2.jpts_p4()[i] * jec,i));
    }
    break;
  case pfJet:
    for ( unsigned int i=0; i < cms2.pfjets_p4().size(); ++i) {
      double jec = 1.0;
      // cout << cms2.evt_event() << " \traw pt: " << cms2.pfjets_p4().at(i).pt() << endl;
      if(applyJEC){
	jet_corrector_pfL1FastJetL2L3->setRho(cms2.evt_ww_rho_vor());
	jet_corrector_pfL1FastJetL2L3->setJetA(cms2.pfjets_area().at(i));
	jet_corrector_pfL1FastJetL2L3->setJetPt(cms2.pfjets_p4()[i].pt());
	jet_corrector_pfL1FastJetL2L3->setJetEta(cms2.pfjets_p4()[i].eta());
	double corr = jet_corrector_pfL1FastJetL2L3->getCorrection();
	jec *= corr;
	// cout << " \tL1FastJetL2L3 corr: " << corr << " \tpt: " << cms2.pfjets_p4().at(i).pt()*corr << endl;
	
	// jec *= jetCorrection(cms2.pfjets_p4()[i], jet_corrector_pf);
	// cout << " \tL2L3 corr: " << jetCorrection(cms2.pfjets_p4()[i], jet_corrector_pf) << endl;
      }
//       if(applyFastJetCorrection){
// 	jec *= (1-cms2.evt_rho()*cms2.pfjets_area().at(i)/cms2.pfjets_p4().at(i).pt()); //*cms2.pfjets_cor().at(i);// It's full L1Fast*L2*L3
// 	cout << " \tL1 corr: " << (1-cms2.evt_rho()*cms2.pfjets_area().at(i)/cms2.pfjets_p4().at(i).pt()) << " \t" << 
// 	  (1-cms2.evt_rho()*cms2.pfjets_area().at(i)/cms2.pfjets_p4().at(i).pt())*cms2.pfjets_p4().at(i).pt() << 
// 	  " \t" << cms2.evt_rho() << " \t" << cms2.pfjets_area().at(i) << endl;
//       }
      if ( cms2.pfjets_p4()[i].pt() * jec < etThreshold ) continue;
      if ( btag && !defaultBTag(type,i) ) continue;
      if ( TMath::Abs(cms2.pfjets_p4()[i].eta()) > etaMax ) continue;
      if (btag && cms2.pfjets_p4()[i].pt() * jec < 30. && fabs(jetDz(i,0))>2.) continue;//newcuts
      if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.pfjets_p4()[i])) < vetoCone ||
	   TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.pfjets_p4()[i])) < vetoCone ) continue;
      // cout << " \tpassed all cuts" << endl;
      jets.push_back(JetPair(cms2.pfjets_p4()[i] * jec,i));
    }
    break;
  case GenJet:
    for ( unsigned int i=0; i < cms2.genjets_p4().size(); ++i) {
      if ( cms2.genjets_p4()[i].pt() < etThreshold ) continue;
      if ( btag && !defaultBTag(type,i) ) continue;
      if ( TMath::Abs(cms2.genjets_p4()[i].eta()) > etaMax ) continue;
      if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.genjets_p4()[i])) < vetoCone ||
	   TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.genjets_p4()[i])) < vetoCone ) continue;
      jets.push_back(JetPair(cms2.genjets_p4()[i],i));
    }
    break;
    //      case CaloJet:
    //        for ( unsigned int i=0; i < cms2.jets_pat_jet_p4().size(); ++i) {
    // 	 if ( cms2.jets_pat_jet_p4()[i].pt() < etThreshold ) continue; // note that this is already corrected
    // 	 if ( btag && !defaultBTag(type,i) ) continue;
    // 	 if ( TMath::Abs(cms2.jets_pat_jet_p4()[i].eta()) > etaMax ) continue;
    // 	 if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.jets_pat_jet_p4()[i])) < vetoCone ||
    // 	      TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.jets_pat_jet_p4()[i])) < vetoCone ) continue;
    // 	 jets.push_back(cms2.jets_pat_jet_p4()[i]);
    //        }
    //        break;
  case TrkJet:
    for ( unsigned int i=0; i < cms2.trkjets_p4().size(); ++i) {
      double jec = 1.0;
      if ( cms2.trkjets_p4()[i].pt() < etThreshold ) continue;
      if ( btag && !defaultBTag(type,i) ) continue;
      if ( TMath::Abs(cms2.trkjets_p4()[i].eta()) > etaMax ) continue;
      if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.trkjets_p4()[i])) < vetoCone ||
	   TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.trkjets_p4()[i])) < vetoCone ) continue;
      jets.push_back(JetPair(cms2.trkjets_p4()[i] * jec,i));
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

bool inZmassWindow(float mass, double delta=15.0){
  // return ( mass > 76. && mass < 106. );
  return fabs(mass - 91.1876) < delta;
}

double BTag(LorentzVector jetP4){
  int refJet = -1;
  for ( unsigned int i=0; i < cms2.pfjets_p4().size(); ++i) {
    if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(jetP4,cms2.pfjets_p4()[i])) > 0.3 ) continue;
    refJet = i;
  }
  if (refJet == -1){
    // std::cout << "Warning: failed to find a matching jet for b-tagging." << std::endl; 
    return 0.0;
  }
  return cms2.pfjets_trackCountingHighEffBJetTag().at(refJet);
}

double BTag(WWJetType type, unsigned int iJet){
  switch ( type ) {
  case jptJet:
    return BTag(cms2.jpts_p4().at(iJet));
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
  return 0;
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


hypo_monitor monitor(true);


bool toptag(WWJetType type, int i_hyp, double minPt,
	    std::vector<JetPair> ignoreJets=std::vector<JetPair>())
{
  //std::vector<LorentzVector> jets;
  const double vetoCone    = 0.3;

  switch ( type ){
  case pfJet:
    for ( unsigned int i=0; i < cms2.pfjets_p4().size(); ++i) {
      if ( cms2.pfjets_p4()[i].pt() < minPt ) continue;
      bool ignoreJet = false;
      for ( std::vector<JetPair>::const_iterator ijet = ignoreJets.begin();
	    ijet != ignoreJets.end(); ++ijet )
	if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(ijet->first,cms2.pfjets_p4()[i])) < vetoCone ) ignoreJet=true;
      if ( ignoreJet ) continue;
      if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.pfjets_p4()[i])) < vetoCone ||
	   TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.pfjets_p4()[i])) < vetoCone ) continue;
      if ( !defaultBTag(type,i) ) continue;
      // dZ cut
      if (fabs(jetDz(i,0))>2) continue;
      return true;
    }
    break;
  case CaloJet:
    for ( unsigned int i=0; i < cms2.jets_p4().size(); ++i) {
      if ( cms2.jets_p4()[i].pt() < minPt ) continue;
      bool ignoreJet = false;
      for ( std::vector<JetPair>::const_iterator ijet = ignoreJets.begin();
	    ijet != ignoreJets.end(); ++ijet )
	if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(ijet->first,cms2.jets_p4()[i])) < vetoCone ) ignoreJet=true;
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

bool hypoSync(int i_hyp, double weight, bool realData) 
{

  HypothesisType type = getHypothesisTypeNew(i_hyp);

  if (nGoodVertex()<1) return false;
  
  monitor.count(cms2,type,"basic selection",weight);
  
  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && !ww_muId(cms2.hyp_lt_index()[i_hyp]) ) return false;
  if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && !ww_muId(cms2.hyp_ll_index()[i_hyp]) ) return false;
  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11) {
    int index = cms2.hyp_lt_index()[i_hyp];
    if (useLHeleId) {
      if (cms2.els_p4().at(index).pt()>20 && (passLikelihoodId_v2(index,cms2.els_lh().at(index),0) & (1<<ELEID_ID))!=(1<<ELEID_ID) ) return false; 
      if (cms2.els_p4().at(index).pt()<20 && (passLikelihoodId_v2(index,cms2.els_lh().at(index),0) & (1<<ELEID_ID))!=(1<<ELEID_ID) ) return false;
    }
    if (useMVAeleId){
      if (!goodElectronTMVA(index)) return false;
    } else {
      if (!pass_electronSelection(index, electronSelection_smurfV3_id, false, false) ) return false;
    }
  } 
  if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11) {
    int index = cms2.hyp_ll_index()[i_hyp];
    if (useLHeleId) {
      if (cms2.els_p4().at(index).pt()>20 && (passLikelihoodId_v2(index,cms2.els_lh().at(index),0) & (1<<ELEID_ID))!=(1<<ELEID_ID) ) return false; 
      if (cms2.els_p4().at(index).pt()<20 && (passLikelihoodId_v2(index,cms2.els_lh().at(index),0) & (1<<ELEID_ID))!=(1<<ELEID_ID) ) return false;
    }
    if (useMVAeleId){
      if (!goodElectronTMVA(index)) return false;
    } else {
      if (!pass_electronSelection(index, electronSelection_smurfV3_id, false, false) ) return false;
    }
  } 

  monitor.count(cms2,type,"lepton id",weight);

  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && !ww_muIso(cms2.hyp_lt_index()[i_hyp]) ) return false;
  if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && !ww_muIso(cms2.hyp_ll_index()[i_hyp]) ) return false;
  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11 && !ww_elIso(cms2.hyp_lt_index()[i_hyp]) ) return false;
  if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11 && !ww_elIso(cms2.hyp_ll_index()[i_hyp]) ) return false;

  monitor.count(cms2,type,"iso",weight);
  
  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11 && 
      cms2.els_exp_innerlayers().at(cms2.hyp_lt_index()[i_hyp]) != 0) return false;
  if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11 && 
      cms2.els_exp_innerlayers().at(cms2.hyp_ll_index()[i_hyp]) != 0) return false;

  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11 && 
      isFromConversionMIT(cms2.hyp_lt_index()[i_hyp]))
      // ! pass_electronSelection(cms2.hyp_lt_index()[i_hyp], (1ll<<ELENOTCONV_MIT), false, false) )
    return false;
  if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11 && 
      isFromConversionMIT(cms2.hyp_ll_index()[i_hyp]))
      // ! pass_electronSelection(cms2.hyp_ll_index()[i_hyp], (1ll<<ELENOTCONV_MIT), false, false) )
    return false;
  
  monitor.count(cms2,type,"conv rejection",weight);
  
  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && !ww_mud0PV(cms2.hyp_lt_index()[i_hyp]) ) return false;
  if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && !ww_mud0PV(cms2.hyp_ll_index()[i_hyp]) ) return false;
  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11 && !ww_eld0PV(cms2.hyp_lt_index()[i_hyp]) ) return false;
  if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11 && !ww_eld0PV(cms2.hyp_ll_index()[i_hyp]) ) return false;
  
  monitor.count(cms2,type,"d0",weight);
  
  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && !ww_mudZPV(cms2.hyp_lt_index()[i_hyp]) ) return false;
  if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && !ww_mudZPV(cms2.hyp_ll_index()[i_hyp]) ) return false;
  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11 && !ww_eldZPV(cms2.hyp_lt_index()[i_hyp]) ) return false;
  if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11 && !ww_eldZPV(cms2.hyp_ll_index()[i_hyp]) ) return false;
  
  monitor.count(cms2,type,"dZ",weight);

//   if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11 && 
//       (fabs(cms2.els_conv_dist().at(cms2.hyp_lt_index()[i_hyp])) < 0.02 &&
//        fabs(cms2.els_conv_dcot().at(cms2.hyp_lt_index()[i_hyp])) < 0.02 )) return false;
//   if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11 && 
//       (fabs(cms2.els_conv_dist().at(cms2.hyp_ll_index()[i_hyp])) < 0.02 &&
//        fabs(cms2.els_conv_dcot().at(cms2.hyp_ll_index()[i_hyp])) < 0.02 )) return false;
  
  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && !goodMuonIsolated(cms2.hyp_lt_index()[i_hyp]) ) return false;
  if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && !goodMuonIsolated(cms2.hyp_ll_index()[i_hyp]) ) return false;
  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11 && !goodElectronIsolated(cms2.hyp_lt_index()[i_hyp]) ) return false;
  if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11 && !goodElectronIsolated(cms2.hyp_ll_index()[i_hyp]) ) return false;
  
  monitor.count(cms2,type,"lepton id/iso",weight);
  // if ( std::min(metValue(),double(trackerMET(i_hyp,0.2).met))<20 ) return false;
  if ( metValue()<20 ) return false;
  monitor.count(cms2,type,"met>20",weight);
  
  if ( std::max(cms2.hyp_lt_p4().at(i_hyp).pt(),cms2.hyp_ll_p4().at(i_hyp).pt())<20 ) return false;
  if (cms2.hyp_p4().at(i_hyp).mass2()<0 || 
      cms2.hyp_p4()[i_hyp].mass() < 12) return false;
  monitor.count(cms2,type,"tight_pt>20 && mll>12",weight);
  
  if ( type == EE || type == MM) {
    if (inZmassWindow(cms2.hyp_p4()[i_hyp].mass())) return false;
  }
  monitor.count(cms2,type,"z veto",weight);
  
  if (!passedMetRequirements(i_hyp)) return false;
  monitor.count(cms2,type,"Full MET cuts",weight);

  if ( numberOfJets(i_hyp)>0 ) return false;
  monitor.count(cms2,type,"JetVeto cuts",weight);

  std::vector<JetPair> sortedJets = getJets(jetType(), i_hyp, 15.0, 5.0, true, false);
  if (cms2.hyp_type().at(i_hyp)!=1 && cms2.hyp_type().at(i_hyp)!=2 && sortedJets.size()>0) {
    if (fabs(ROOT::Math::VectorUtil::DeltaPhi(sortedJets[0].first,cms2.hyp_p4().at(i_hyp)))>=165.*TMath::Pi()/180.) return false;
  }
  monitor.count(cms2,type,"Jet DPhi",weight);

  if (numberOfSoftMuons(i_hyp,true)>0) return false;
  monitor.count(cms2,type,"soft muon veto",weight);

  if (numberOfExtraLeptons(i_hyp,10)>0) return false;
  monitor.count(cms2,type,"extra lepton veto",weight);

  if ( toptag(jetType(),i_hyp,10/*,vetojets*/) )  return false;//newcuts, was 7
  monitor.count(cms2,type,"top tag",weight);

  if ( TMath::Max(cms2.hyp_ll_p4().at(i_hyp).pt(),cms2.hyp_lt_p4().at(i_hyp).pt())<=30. )  return false;
  monitor.count(cms2,type,"ptMax",weight);

  if ( TMath::Min(cms2.hyp_ll_p4().at(i_hyp).pt(),cms2.hyp_lt_p4().at(i_hyp).pt())<=25. )  return false;
  monitor.count(cms2,type,"ptMin",weight);

  if ( cms2.hyp_p4().at(i_hyp).mass()>=50. )  return false;
  monitor.count(cms2,type,"mll",weight);

  if ( mt(cms2.hyp_p4().at(i_hyp).pt(),cms2.evt_pfmet(),acos(cos(cms2.hyp_p4().at(i_hyp).phi()-metPhiValue())))<=90. || mt(cms2.hyp_p4().at(i_hyp).pt(),cms2.evt_pfmet(),acos(cos(cms2.hyp_p4().at(i_hyp).phi()-metPhiValue())))>=160. )  return false;
  monitor.count(cms2,type,"mT",weight);

  if ( fabs(ROOT::Math::VectorUtil::DeltaPhi(cms2.hyp_ll_p4().at(i_hyp),cms2.hyp_lt_p4().at(i_hyp)))>=60.*TMath::Pi()/180. )  return false;
  monitor.count(cms2,type,"dphill",weight);
  
  if ( cms2.hyp_p4().at(i_hyp).pt() < 45 ) return false;
  if ( (type==EE || type==MM) && cms2.hyp_lt_p4().at(i_hyp).pt() < 15 ) return false;
  if ( (type==EE || type==MM) && cms2.hyp_ll_p4().at(i_hyp).pt() < 15 ) return false;
  monitor.count(cms2,type,"dilep.pt()>45 && lep2.pt()>15",weight);
  return true;
} // end of Synchronization info

bool hypo (int i_hyp, double weight, bool realData) 
{
  cuts_passed = 0;
  if ( std::max(cms2.hyp_lt_p4().at(i_hyp).pt(),cms2.hyp_ll_p4().at(i_hyp).pt())<20 ) return false;
  if ( std::min(cms2.hyp_lt_p4().at(i_hyp).pt(),cms2.hyp_ll_p4().at(i_hyp).pt())<10 ) return false;

  /*
  unsigned int nGenLeptons = 0;
  for ( unsigned int i=0; i<cms2.genps_id().size(); ++i)
    if ( abs(cms2.genps_id().at(i)) == 11 || abs(cms2.genps_id().at(i)) == 13 )
      nGenLeptons++;
  if ( nGenLeptons < 2 ) return;
  */
  // if (cms2.evt_event()!=101838) return;
  HypothesisType type = getHypothesisTypeNew(i_hyp);

  // The event weight including the kFactor (scaled to 1 fb-1)
  // float weight = cms2.evt_scale1fb() * kFactor;

  if (!realData) monitor.nEvtProcessed = cms2.evt_nEvts();
  // monitor.count(cms2, type, "Total number before cuts");
  
  // if ( cms2.hyp_FVFit_prob()[i_hyp] < 0.005 ) return;
  // monitor.count(cms2, type, "after vertex cut");
  
  if (cms2.hyp_p4().at(i_hyp).mass2()<0) return false;
  if (!isGoodVertex(primaryVertex())) return false;
  // if ( realData && ! passedTriggerRequirements( hypType(i_hyp) ) )return false;
  // if ( realData && passedTriggerRequirements( hypType(i_hyp) ) )return false;
  monitor.count(cms2, type, "trigger requirements",weight);
  if (nGoodVertex()<1) return false;

  cuts_passed = PASSED_BaseLine;
  if ( realData && passedTriggerRequirements( hypType(i_hyp) ) )  cuts_passed |= PASSED_Trigger;
  
  // Baseline cuts
  if (abs(cms2.hyp_lt_id()[i_hyp]) == 13 && !ww_muBase(cms2.hyp_lt_index()[i_hyp]) ) cuts_passed &= ~PASSED_BaseLine;
  if (abs(cms2.hyp_ll_id()[i_hyp]) == 13 && !ww_muBase(cms2.hyp_ll_index()[i_hyp]) ) cuts_passed &= ~PASSED_BaseLine;
  if (abs(cms2.hyp_lt_id()[i_hyp]) == 11 && !ww_elBase(cms2.hyp_lt_index()[i_hyp]) ) cuts_passed &= ~PASSED_BaseLine;
  if (abs(cms2.hyp_ll_id()[i_hyp]) == 11 && !ww_elBase(cms2.hyp_ll_index()[i_hyp]) ) cuts_passed &= ~PASSED_BaseLine;
  
  // Require opposite sign
  if ( cms2.hyp_lt_id()[i_hyp] * cms2.hyp_ll_id()[i_hyp] < 0 ) cuts_passed |= PASSED_Charge;

  // monitor.count(cms2,type,"baseline cuts",weight);
 
  if (gSystem->Getenv("Sync")) // Synchronization info
    {
      if ( !CheckCuts( PASSED_BaseLine|PASSED_Charge, cuts_passed ) ) return false;
      if (!hypoSync(i_hyp,weight,realData)) return false;
    }

  if ( std::max(cms2.hyp_lt_p4().at(i_hyp).pt(),cms2.hyp_ll_p4().at(i_hyp).pt())<20 ) cuts_passed &= ~PASSED_BaseLine;
  
  if ( cms2.hyp_p4()[i_hyp].mass() < 12) cuts_passed &= ~PASSED_BaseLine;

  // check electron isolation and id (no selection at this point)
  // checkIsolation(i_hyp, weight);

  // == Z mass veto using hyp_leptons for ee and mumu final states
  if ( type == EE || type == MM) {
    if (!inZmassWindow(cms2.hyp_p4()[i_hyp].mass()))    cuts_passed |= PASSED_ZVETO;
  }
  if (type == EM || type == ME)     cuts_passed |= PASSED_ZVETO;
  if ( inZmassWindow(cms2.hyp_p4()[i_hyp].mass(),5))  cuts_passed |= PASSED_ZControlSampleVeryTight;
  if ( inZmassWindow(cms2.hyp_p4()[i_hyp].mass(),10))  cuts_passed |= PASSED_ZControlSampleTight;
  if ( inZmassWindow(cms2.hyp_p4()[i_hyp].mass(),20))  cuts_passed |= PASSED_ZControlSampleLoose;
    
  // Z veto using additional leptons in the event
  // if (additionalZveto()) return;

  // == MET
  if ( passedMetRequirements(i_hyp) ) cuts_passed |= PASSED_MET;
  
  // == letpon ID and Isolation
  // Muon quality cuts, including isolation
  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13){
    unsigned int index = cms2.hyp_lt_index()[i_hyp];
    if ( goodMuonIsolated(index) )   cuts_passed |= PASSED_LT_FINAL;
    if ( fakableMuon(index,MuFOV1) ) cuts_passed |= PASSED_LT_FO_MU1;
    if ( fakableMuon(index,MuFOV2) ) cuts_passed |= PASSED_LT_FO_MU2;
  }
  if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13){
    unsigned int index = cms2.hyp_ll_index()[i_hyp];
    if ( goodMuonIsolated(index) )   cuts_passed |= PASSED_LL_FINAL;
    if ( fakableMuon(index,MuFOV1) ) cuts_passed |= PASSED_LL_FO_MU1;
    if ( fakableMuon(index,MuFOV2) ) cuts_passed |= PASSED_LL_FO_MU2;
  } 
  // Electron quality cuts, including isolation
  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11){
    unsigned int index = cms2.hyp_lt_index()[i_hyp];
    if ( goodElectronIsolated(index) ) cuts_passed |= PASSED_LT_FINAL;
    if ( fakableElectron(index,EleFOV1) )   cuts_passed |= PASSED_LT_FO_ELEV1;
    if ( fakableElectron(index,EleFOV2) )   cuts_passed |= PASSED_LT_FO_ELEV2;
    if ( fakableElectron(index,EleFOV3) )   cuts_passed |= PASSED_LT_FO_ELEV3;
    if ( fakableElectron(index,EleFOV4) )   cuts_passed |= PASSED_LT_FO_ELEV4;
  }
  if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11){
    unsigned int index = cms2.hyp_ll_index()[i_hyp];
    if ( goodElectronIsolated(index) ) cuts_passed |= PASSED_LL_FINAL;
    if ( fakableElectron(index,EleFOV1) )   cuts_passed |= PASSED_LL_FO_ELEV1;
    if ( fakableElectron(index,EleFOV2) )   cuts_passed |= PASSED_LL_FO_ELEV2;
    if ( fakableElectron(index,EleFOV3) )   cuts_passed |= PASSED_LL_FO_ELEV3;
    if ( fakableElectron(index,EleFOV4) )   cuts_passed |= PASSED_LL_FO_ELEV4;
  }

  // if ( passedLTFinalRequirements || passedLLFinalRequirements ) cuts_passed |= (1<<PASS_PROBE);

  // == Jet-veto
  const std::vector<JetPair>& vetojets(getDefaultJets(i_hyp, false));
  unsigned int nJets = vetojets.size();
  if (nJets==0) cuts_passed |= PASSED_JETVETO;
  if (nJets>1 && 
      getDefaultJets(i_hyp,true).size()>0 )  
    cuts_passed |= PASSED_TopControlSample;
  if (nJets==1){
    if ( getDefaultJets(i_hyp,true).size()==1 ){
      cuts_passed |= PASSED_1BJET;
    }
  }
  // trkjet veto
  // if ( !passTrkJetVeto(i_hyp) ) return;
  // == ExtraVeto Muons
  int countmus = numberOfSoftMuons(i_hyp,true);
  int nExtraVetoMuons = numberOfSoftMuons(i_hyp,true,vetojets);
  if ( countmus == 0)        cuts_passed |=   PASSED_SOFTMUVETO;
  if ( nExtraVetoMuons == 0) cuts_passed |=   PASSED_SOFTMUVETO_NotInJets;
  if ( numberOfExtraLeptons(i_hyp,10) == 0) cuts_passed |= PASSED_EXTRALEPTONVETO;
  if ( ! toptag(jetType(),i_hyp,10,vetojets) ) cuts_passed |= PASSED_TOPVETO_NotInJets;//newcuts, was 7
  if ( ! toptag(jetType(),i_hyp,10) )          cuts_passed |= PASSED_TOPVETO;//newcuts, was 7

  if( CheckCuts(  PASSED_MET | PASSED_BaseLine , cuts_passed)) {
    if ( CheckCuts( PASSED_LT_FINAL | PASSED_LL_FO_MU1,   cuts_passed ) ||
	 CheckCuts( PASSED_LT_FINAL | PASSED_LL_FO_MU2,   cuts_passed ) ||
	 CheckCuts( PASSED_LT_FINAL | PASSED_LL_FO_ELEV1, cuts_passed ) ||
	 CheckCuts( PASSED_LT_FINAL | PASSED_LL_FO_ELEV2, cuts_passed ) ||
	 CheckCuts( PASSED_LT_FINAL | PASSED_LL_FO_ELEV3, cuts_passed ) ||
	 CheckCuts( PASSED_LT_FINAL | PASSED_LL_FO_ELEV4, cuts_passed ) ||
	 CheckCuts( PASSED_LL_FINAL | PASSED_LT_FO_MU1,   cuts_passed ) ||
	 CheckCuts( PASSED_LL_FINAL | PASSED_LT_FO_MU2,   cuts_passed ) ||
	 CheckCuts( PASSED_LL_FINAL | PASSED_LT_FO_ELEV1, cuts_passed ) ||
	 CheckCuts( PASSED_LL_FINAL | PASSED_LT_FO_ELEV2, cuts_passed ) ||
	 CheckCuts( PASSED_LL_FINAL | PASSED_LT_FO_ELEV3, cuts_passed ) ||
	 CheckCuts( PASSED_LL_FINAL | PASSED_LT_FO_ELEV4, cuts_passed ) )
      cuts_passed |= PASSED_Skim1;
  }
  if( CheckCuts(  PASSED_BaseLine, cuts_passed)) {
    if ( CheckCuts( PASSED_LT_FINAL | PASSED_LL_FO_MU1,   cuts_passed ) ||
	 CheckCuts( PASSED_LT_FINAL | PASSED_LL_FO_MU2,   cuts_passed ) ||
	 CheckCuts( PASSED_LT_FINAL | PASSED_LL_FO_ELEV1, cuts_passed ) ||
	 CheckCuts( PASSED_LT_FINAL | PASSED_LL_FO_ELEV2, cuts_passed ) ||
	 CheckCuts( PASSED_LT_FINAL | PASSED_LL_FO_ELEV3, cuts_passed ) ||
	 CheckCuts( PASSED_LT_FINAL | PASSED_LL_FO_ELEV4, cuts_passed ) ||
	 CheckCuts( PASSED_LL_FINAL | PASSED_LT_FO_MU1,   cuts_passed ) ||
	 CheckCuts( PASSED_LL_FINAL | PASSED_LT_FO_MU2,   cuts_passed ) ||
	 CheckCuts( PASSED_LL_FINAL | PASSED_LT_FO_ELEV1, cuts_passed ) ||
	 CheckCuts( PASSED_LL_FINAL | PASSED_LT_FO_ELEV2, cuts_passed ) ||
	 CheckCuts( PASSED_LL_FINAL | PASSED_LT_FO_ELEV3, cuts_passed ) ||
	 CheckCuts( PASSED_LL_FINAL | PASSED_LT_FO_ELEV4, cuts_passed ) )
      cuts_passed |= PASSED_Skim3;
  }
  
  monitor.count(cms2,type,"all cuts (including soft and extra lepton veto)",weight);
  return true;

}//end of void hypo

double mt(double pt1, double pt2, double dphi){
  return 2*sqrt(pt1*pt2)*fabs(sin(dphi/2));
}

void FillSmurfNtuple(SmurfTree& tree, unsigned int i_hyp, 
		     double weight, enum SmurfTree::DataType sample, wwcuts_t cuts_passed){
  tree.InitVariables();
  tree.run_         = cms2.evt_run();
  tree.event_       = cms2.evt_event();
  tree.lumi_        = cms2.evt_lumiBlock();
  tree.nvtx_        = nGoodVertex();
  tree.scale1fb_    = weight;
  tree.met_         = metValue();
  tree.sumet_       = sumetValue();
  tree.metPhi_      = metPhiValue();
  metStruct trkMET  = trackerMET(i_hyp,0.1);
  tree.trackMet_    = trkMET.met;
  tree.trackMetPhi_ = trkMET.metphi;
  tree.pTrackMet_   = projectedMet(i_hyp, trkMET.met, trkMET.metphi);

  bool ltIsFirst = true;
  if ( cms2.hyp_lt_p4().at(i_hyp).pt()<cms2.hyp_ll_p4().at(i_hyp).pt() ) ltIsFirst = false;
  tree.type_ = SmurfTree::Type(cms2.hyp_type().at(i_hyp));
  if ( tree.type_ == SmurfTree::em || tree.type_ == SmurfTree::me ){
    if ( ltIsFirst )
      tree.type_ = abs(cms2.hyp_lt_id().at(i_hyp))==11 ? SmurfTree::em : SmurfTree::me;
    else
      tree.type_ = abs(cms2.hyp_lt_id().at(i_hyp))==11 ? SmurfTree::me : SmurfTree::em;
  }

  tree.lep1_  = ltIsFirst ? cms2.hyp_lt_p4().at(i_hyp) : cms2.hyp_ll_p4().at(i_hyp);
  tree.lep2_  = ltIsFirst ? cms2.hyp_ll_p4().at(i_hyp) : cms2.hyp_lt_p4().at(i_hyp);
  tree.lq1_   = ltIsFirst ? cms2.hyp_lt_charge().at(i_hyp) : cms2.hyp_ll_charge().at(i_hyp);
  tree.lq2_   = ltIsFirst ? cms2.hyp_ll_charge().at(i_hyp) : cms2.hyp_lt_charge().at(i_hyp);
  tree.lid1_  = ltIsFirst ? cms2.hyp_lt_id().at(i_hyp) : cms2.hyp_ll_id().at(i_hyp);
  tree.lid2_  = ltIsFirst ? cms2.hyp_ll_id().at(i_hyp) : cms2.hyp_lt_id().at(i_hyp);
  const std::vector<JetPair>& jets = getJets(jetType(), i_hyp, 0, 5.0, true, false); //GC apply sorting
  if (jets.size()>0){
    tree.jet1_ = jets.at(0).first;
    tree.jet1Btag_ = BTag(jetType(),jets.at(0).second);
  }
  if (jets.size()>1){
    tree.jet2_ = jets.at(1).first;
    tree.jet2Btag_ = BTag(jetType(),jets.at(1).second);
  }
  if (jets.size()>2){
    tree.jet3_ = jets.at(2).first;
    tree.jet3Btag_ = BTag(jetType(),jets.at(2).second);
  }
  tree.njets_ = numberOfJets(i_hyp);
  tree.dilep_ = cms2.hyp_p4().at(i_hyp);
  tree.pmet_ = projectedMet(i_hyp, metValue(), metPhiValue());
  tree.dPhi_ = fabs(ROOT::Math::VectorUtil::DeltaPhi(cms2.hyp_lt_p4().at(i_hyp),cms2.hyp_ll_p4().at(i_hyp)));
  tree.dR_   = ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4().at(i_hyp),cms2.hyp_ll_p4().at(i_hyp));
  if (jets.size()>0) {
    tree.dPhiLep1Jet1_ = fabs(ROOT::Math::VectorUtil::DeltaPhi(cms2.hyp_lt_p4().at(i_hyp),jets.at(0).first));
    tree.dRLep1Jet1_   = ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4().at(i_hyp),jets.at(0).first);
    tree.dPhiLep2Jet1_ = fabs(ROOT::Math::VectorUtil::DeltaPhi(cms2.hyp_ll_p4().at(i_hyp),jets.at(0).first));
    tree.dRLep2Jet1_   = ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4().at(i_hyp),jets.at(0).first);
    tree.dPhiDiLepJet1_= fabs(ROOT::Math::VectorUtil::DeltaPhi(cms2.hyp_p4().at(i_hyp),jets.at(0).first));
    if (!ltIsFirst){
      std::swap(tree.dPhiLep1Jet1_,tree.dPhiLep2Jet1_);
      std::swap(tree.dRLep1Jet1_,tree.dRLep2Jet1_);
    }
  }
  tree.dPhiDiLepMET_ = acos(cos(cms2.hyp_p4().at(i_hyp).phi()-metPhiValue()));
  tree.dPhiLep1MET_ = acos(cos(tree.lep1_.phi()-metPhiValue()));
  tree.dPhiLep2MET_ = acos(cos(tree.lep2_.phi()-metPhiValue()));

  tree.mt_ = mt(tree.dilep_.pt(),tree.met_,tree.dPhiDiLepMET_);
  tree.mt1_ = mt(tree.lep1_.pt(),tree.met_,tree.dPhiLep1MET_);
  tree.mt2_ = mt(tree.lep2_.pt(),tree.met_,tree.dPhiLep2MET_);

  if (sample!=SmurfTree::data){
    tree.genmet_ = cms2.gen_met();
    tree.genmetPhi_ = cms2.gen_metPhi();
    tree.lep1McId_ = ltIsFirst ? cms2.hyp_lt_mc_id().at(i_hyp) : cms2.hyp_ll_mc_id().at(i_hyp);
    tree.lep2McId_ = ltIsFirst ? cms2.hyp_ll_mc_id().at(i_hyp) : cms2.hyp_lt_mc_id().at(i_hyp);
    tree.lep1MotherMcId_  = ltIsFirst ? cms2.hyp_lt_mc_motherid().at(i_hyp) : cms2.hyp_ll_mc_motherid().at(i_hyp);// GC
    tree.lep2MotherMcId_  = ltIsFirst ? cms2.hyp_ll_mc_motherid().at(i_hyp) : cms2.hyp_lt_mc_motherid().at(i_hyp);// GC
    if (jetType()==pfJet) { //GC fixme in case they are not pfjets 
      if (jets.size()>0){
	tree.jet1McId_ = cms2.pfjets_mc_id().at(jets.at(0).second);
      }
      if (jets.size()>1){
	tree.jet2McId_ = cms2.pfjets_mc_id().at(jets.at(1).second);
      }
      if (jets.size()>2){
	tree.jet3McId_ = cms2.pfjets_mc_id().at(jets.at(2).second);
      }
    }
    tree.processId_ = cms2.genps_signalProcessID();
    tree.higgsPt_ = -9999.;
    tree.sfWeightHPt_ = 1.;
    if (tree.processId_==10010 || tree.processId_==10001) {
      float pt = getHiggsPt();
      tree.higgsPt_ = pt;
      if (tree.processId_==10010) {
	TString ss = tree.name(sample);//get higgs mass from sample name
	ss.Remove(0,3);//get rid of hww
	tree.sfWeightHPt_ = getHiggsPtWeight(pt,ss.Atoi());
      }
    }

    if (doDYNNLOw) {
      if (sample==SmurfTree::dyee||sample==SmurfTree::dymm||sample==SmurfTree::dytt) {
	if (isDYee()||isDYmm()||isDYtt()) {
	  float pt = getZPt();
	  float rap = getZRapidity();
	  float mass = getZMass();
	  float kfact =  getDYNNLOWeight(pt,rap,mass);
	  tree.scale1fb_ *= kfact;
	}
      }
    }

    for (unsigned int nbc=0;nbc<cms2.puInfo_nPUvertices().size();++nbc) {
      if (cms2.puInfo_bunchCrossing().at(nbc)==0) tree.npu_ = cms2.puInfo_nPUvertices().at(nbc);
      else if (cms2.puInfo_bunchCrossing().at(nbc)==-1) tree.npuMinusOne_ = cms2.puInfo_nPUvertices().at(nbc);
      else if (cms2.puInfo_bunchCrossing().at(nbc)==+1) tree.npuPlusOne_ = cms2.puInfo_nPUvertices().at(nbc);
    }

  }

  tree.dstype_ = sample;

  // fill cuts
  tree.cuts_ = 0;
  if (ltIsFirst){
    if ( cuts_passed & PASSED_LT_FINAL )    tree.cuts_ |= SmurfTree::Lep1FullSelection;
    if ( cuts_passed & PASSED_LT_FO_MU1 )   tree.cuts_ |= SmurfTree::Lep1LooseMuV1;
    if ( cuts_passed & PASSED_LT_FO_MU2 )   tree.cuts_ |= SmurfTree::Lep1LooseMuV2;
    if ( cuts_passed & PASSED_LT_FO_ELEV1 ) tree.cuts_ |= SmurfTree::Lep1LooseEleV1;
    if ( cuts_passed & PASSED_LT_FO_ELEV2 ) tree.cuts_ |= SmurfTree::Lep1LooseEleV2;
    if ( cuts_passed & PASSED_LT_FO_ELEV3 ) tree.cuts_ |= SmurfTree::Lep1LooseEleV3;
    if ( cuts_passed & PASSED_LT_FO_ELEV4 ) tree.cuts_ |= SmurfTree::Lep1LooseEleV4;
    if ( cuts_passed & PASSED_LL_FINAL )    tree.cuts_ |= SmurfTree::Lep2FullSelection;
    if ( cuts_passed & PASSED_LL_FO_MU1 )   tree.cuts_ |= SmurfTree::Lep2LooseMuV1;
    if ( cuts_passed & PASSED_LL_FO_MU2 )   tree.cuts_ |= SmurfTree::Lep2LooseMuV2;
    if ( cuts_passed & PASSED_LL_FO_ELEV1 ) tree.cuts_ |= SmurfTree::Lep2LooseEleV1;
    if ( cuts_passed & PASSED_LL_FO_ELEV2 ) tree.cuts_ |= SmurfTree::Lep2LooseEleV2;
    if ( cuts_passed & PASSED_LL_FO_ELEV3 ) tree.cuts_ |= SmurfTree::Lep2LooseEleV3;
    if ( cuts_passed & PASSED_LL_FO_ELEV4 ) tree.cuts_ |= SmurfTree::Lep2LooseEleV4;
  } else {
    if ( cuts_passed & PASSED_LL_FINAL )    tree.cuts_ |= SmurfTree::Lep1FullSelection;
    if ( cuts_passed & PASSED_LL_FO_MU1 )   tree.cuts_ |= SmurfTree::Lep1LooseMuV1;
    if ( cuts_passed & PASSED_LL_FO_MU2 )   tree.cuts_ |= SmurfTree::Lep1LooseMuV2;
    if ( cuts_passed & PASSED_LL_FO_ELEV1 ) tree.cuts_ |= SmurfTree::Lep1LooseEleV1;
    if ( cuts_passed & PASSED_LL_FO_ELEV2 ) tree.cuts_ |= SmurfTree::Lep1LooseEleV2;
    if ( cuts_passed & PASSED_LL_FO_ELEV3 ) tree.cuts_ |= SmurfTree::Lep1LooseEleV3;
    if ( cuts_passed & PASSED_LL_FO_ELEV4 ) tree.cuts_ |= SmurfTree::Lep1LooseEleV4;
    if ( cuts_passed & PASSED_LT_FINAL )    tree.cuts_ |= SmurfTree::Lep2FullSelection;
    if ( cuts_passed & PASSED_LT_FO_MU1 )   tree.cuts_ |= SmurfTree::Lep2LooseMuV1;
    if ( cuts_passed & PASSED_LT_FO_MU2 )   tree.cuts_ |= SmurfTree::Lep2LooseMuV2;
    if ( cuts_passed & PASSED_LT_FO_ELEV1 ) tree.cuts_ |= SmurfTree::Lep2LooseEleV1;
    if ( cuts_passed & PASSED_LT_FO_ELEV2 ) tree.cuts_ |= SmurfTree::Lep2LooseEleV2;
    if ( cuts_passed & PASSED_LT_FO_ELEV3 ) tree.cuts_ |= SmurfTree::Lep2LooseEleV3;
    if ( cuts_passed & PASSED_LT_FO_ELEV4 ) tree.cuts_ |= SmurfTree::Lep2LooseEleV4;
  }
  if ( cuts_passed & PASSED_BaseLine )                tree.cuts_ |= SmurfTree::BaseLine;
  if ( cuts_passed & PASSED_Charge )                  tree.cuts_ |= SmurfTree::ChargeMatch;
  if ( cuts_passed & PASSED_MET )                     tree.cuts_ |= SmurfTree::FullMET;
  if ( cuts_passed & PASSED_ZVETO )                   tree.cuts_ |= SmurfTree::ZVeto;
  if ( cuts_passed & PASSED_Trigger )                 tree.cuts_ |= SmurfTree::Trigger;
  if ( !(cuts_passed & PASSED_TOPVETO) || !(cuts_passed & PASSED_SOFTMUVETO) )           
    tree.cuts_ |= SmurfTree::TopTag;
  else
    tree.cuts_ |= SmurfTree::TopVeto;
  if ( !(cuts_passed & PASSED_TOPVETO_NotInJets) ||
       !(cuts_passed & PASSED_SOFTMUVETO_NotInJets) ) tree.cuts_ |= SmurfTree::TopTagNotInJets;
  if ( cuts_passed & PASSED_1BJET )                   tree.cuts_ |= SmurfTree::OneBJet;
  if ( cuts_passed & PASSED_EXTRALEPTONVETO )         tree.cuts_ |= SmurfTree::ExtraLeptonVeto;

  // find third lepton and fill its information
  tree.nSoftMuons_ = numberOfSoftMuons(i_hyp,true);

  // highest Btag for soft jets (jets below the threshold)
  tree.jetLowBtag_ = -999.;
  for (unsigned int i=tree.njets_; i<jets.size(); ++i)
    if ( tree.jetLowBtag_ < BTag(jetType(),jets.at(i).second) )
      tree.jetLowBtag_ = BTag(jetType(),jets.at(i).second);

  // third lepton
  const std::vector<LeptonPair>& leptons(getExtraLeptons(i_hyp,0));
  for (std::vector<LeptonPair>::const_iterator lep=leptons.begin();
       lep!=leptons.end(); ++lep){
    if (lep->first){
      // muons
      if ( tree.lep3_.pt() < cms2.mus_p4().at(lep->second).pt() ){
	tree.lep3_ = cms2.mus_p4().at(lep->second);
	tree.lq3_   = cms2.mus_charge().at(lep->second);
	tree.lid3_ = tree.lq3_>0?-13:13;
	if (sample!=SmurfTree::data){// GC
	  tree.lep3McId_ = cms2.mus_mc_id().at(lep->second);
	  tree.lep3MotherMcId_  = cms2.mus_mc_motherid().at(lep->second);
	}
      }
    } else {
      if ( tree.lep3_.pt() < cms2.els_p4().at(lep->second).pt() ){
	tree.lep3_ = cms2.els_p4().at(lep->second);
	tree.lq3_   = cms2.els_charge().at(lep->second);
	tree.lid3_ = tree.lq3_>0?-11:11;
	if (sample!=SmurfTree::data){// GC
	  tree.lep3McId_ = cms2.els_mc_id().at(lep->second);
	  tree.lep3MotherMcId_  = cms2.els_mc_motherid().at(lep->second);
	}
      }
    }
  }
  if (tree.lep3_.pt()>0 && jets.size()>0) {
    tree.dPhiLep3Jet1_ = fabs(ROOT::Math::VectorUtil::DeltaPhi(tree.lep3_,jets.at(0).first));
    tree.dRLep3Jet1_   = ROOT::Math::VectorUtil::DeltaR(tree.lep3_,jets.at(0).first);
  }
}

int bestHypothesis(const std::vector<std::pair<unsigned int, unsigned int> >& candidates){
  int best = -1;
  for( unsigned int i = 0; i < candidates.size(); ++i ) {
    unsigned int i_hyp = candidates.at(i).first;
    if (best<0){
      best = i_hyp;
      continue;
    }
    if ( std::max(cms2.hyp_lt_p4().at(i_hyp).pt(), cms2.hyp_ll_p4().at(i_hyp).pt()) >= //GC add = in case the lepton is the same 
	 std::max(cms2.hyp_lt_p4().at(best).pt(),  cms2.hyp_ll_p4().at(best).pt()) &&
	 std::min(cms2.hyp_lt_p4().at(i_hyp).pt(), cms2.hyp_ll_p4().at(i_hyp).pt()) >= 
	 std::min(cms2.hyp_lt_p4().at(best).pt(),  cms2.hyp_ll_p4().at(best).pt()) )
      best = i_hyp;
  }
  return best;
}

void ScanChain( TChain* chain, 
		enum SmurfTree::DataType sample, 
		double integratedLumi, // in unit of pb^-1, if negative the weight is 1.
		double xsec,           // in unit of pb, if negative take it from evt_xsec_excl*evt_kfactor
		int nProcessedEvents,  // if negative, take it from evt_nEvts
		bool identifyEvents, 
		bool realData,
		TString cms2_json_file)
{
  std::string prefix = SmurfTree::name(sample);
  if ( chain->GetListOfFiles()->GetEntries()==0 ){
    printf("\nERROR: chain is empty for sample: %s\n\n",prefix.c_str());
    assert(0);
  }
  if ( gSystem->Getenv("ZSelection") )
    pass_all = PASSED_ZControlSampleTight | PASSED_JETVETO | PASSED_LT_FINAL | PASSED_LL_FINAL | 
      PASSED_SOFTMUVETO | PASSED_EXTRALEPTONVETO | PASSED_TOPVETO;
  
  if ( gSystem->Getenv("ZSelectionVeryTight") )
    pass_all = PASSED_ZControlSampleVeryTight | PASSED_JETVETO | PASSED_LT_FINAL | PASSED_LL_FINAL | 
      PASSED_SOFTMUVETO | PASSED_EXTRALEPTONVETO | PASSED_TOPVETO;

  if ( gSystem->Getenv("TopSelection") )
    pass_all = PASSED_ZVETO | PASSED_TopControlSample | PASSED_MET | 
      PASSED_LT_FINAL | PASSED_LL_FINAL | PASSED_EXTRALEPTONVETO;

  unsigned int nEventsChain = chain->GetEntries();  // number of entries in chain --> number of events from all files
  gErrorIgnoreLevel = -1;
  unsigned int nEventsTotal = 0;
  
  // make smurf ntuples
  gSystem->MakeDirectory("smurf");
  TFile* fSmurf = TFile::Open(Form("smurf/%s.root",prefix.c_str()),"RECREATE");
  assert(fSmurf);
  SmurfTree smurfTree;
  smurfTree.CreateTree();
  smurfTree.tree_->SetDirectory(fSmurf);

  std::map<unsigned int, std::set<unsigned int> > runList;

  // clear list of duplicates
  already_seen.clear();
  int duplicates_total_n = 0;
  double duplicates_total_weight = 0;
  int nFailedIdentification = 0;
  int nFilteredOut = 0;

  int i_permille_old = 0;

  try { 
    jetcorr_filenames_pfL1FastJetL2L3.clear();
    if (realData) {
      jetcorr_filenames_pfL1FastJetL2L3.push_back("files/GR_R_42_V14_AK5PF_L1FastJet.txt");
      jetcorr_filenames_pfL1FastJetL2L3.push_back("files/GR_R_42_V14_AK5PF_L2Relative.txt");
      jetcorr_filenames_pfL1FastJetL2L3.push_back("files/GR_R_42_V14_AK5PF_L3Absolute.txt");
      jetcorr_filenames_pfL1FastJetL2L3.push_back("files/GR_R_42_V19_AK5PF_L2L3Residual.txt");
    } else {
      jetcorr_filenames_pfL1FastJetL2L3.push_back("files/START42_V13_AK5PF_L1FastJet.txt");
      jetcorr_filenames_pfL1FastJetL2L3.push_back("files/START42_V13_AK5PF_L2Relative.txt");
      jetcorr_filenames_pfL1FastJetL2L3.push_back("files/START42_V13_AK5PF_L3Absolute.txt");
    }
    jet_corrector_pfL1FastJetL2L3= makeJetCorrector(jetcorr_filenames_pfL1FastJetL2L3);
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

  if (doDYNNLOw && (sample==SmurfTree::dyee || sample==SmurfTree::dymm || sample==SmurfTree::dytt) ) {
    //TFile *tmpFile = new TFile( "files/DYNNLOKFactor_PowhegToFEWZ.root", "READ");
    //const int nMassBins = 13;
    TFile *tmpFile = new TFile( "files/fewz_powheg_weights_stepwise_2011_fine7.root", "READ");
    const int nMassBins = 41;
    for(int i=0; i<nMassBins; i++){         
      TString hname = TString::Format("weight_%02d",i+1);
      TH2D *tmpHist = (TH2D*)tmpFile->Get(hname);
      tmpHist->SetDirectory(0);
      fDYNNLOKFactorHists.push_back(tmpHist);
    }
    tmpFile->Close();
    delete tmpFile;
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
    TTreeCache::SetLearnEntries(10);
    tree->SetCacheSize(128*1024*1024);

    cms2.Init(tree);  // set branch addresses for TTree tree
    
    TStopwatch t;
    //Event Loop
    unsigned int nEvents = tree->GetEntries();

    for( unsigned int event = 0; event < nEvents; ++event) {
      tree->LoadTree(event);
      cms2.GetEntry(event);  // get entries for Event number event from branches of TTree tree
      if (cms2.evt_event()%prescale!=0) continue;
      //if (cms2.evt_event()!=67523347) continue;
      // if (cms2.evt_event()!=13595393 && cms2.evt_event()!=227560649) continue;
      // Select the good runs from the json file
      if(realData && cms2_json_file!="") {
	if( !goodrun(cms2.evt_run(), cms2.evt_lumiBlock()) ) continue;
      }
      runList[cms2.evt_run()].insert(cms2.evt_lumiBlock());
      double weight = 1.0;
      if ( !realData && integratedLumi>0 ){
	double mcweight = cms2.genps_weight() > 0.0 ? 1.0 : -1.0;
	weight = integratedLumi * mcweight * (xsec>0?xsec:cms2.evt_xsec_excl()*cms2.evt_kfactor()*cms2.evt_filt_eff()) /
	     (nProcessedEvents>0?nProcessedEvents:cms2.evt_nEvts());
	 }       
	 ++nEventsTotal;

	 if (cms2.trks_d0().size() == 0) continue;  // needed to get rid of back Monte Carlo events in CMSSW_2_X analysis
	 if (cms2.hyp_type().size() == 0) continue; // skip events without hypothesis
	 EventIdentifier id(cms2,realData);
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
	 std::vector<std::pair<unsigned int,unsigned int> > candidates;
	 for( unsigned int i_hyp = 0; i_hyp < nHyps; ++i_hyp ) {
	   cuts_passed = 0;
	   if(cms2.hyp_p4().at(i_hyp).mass2() < 0 ) continue;
	   if ( hypo(i_hyp, weight,  realData) ) 
	     candidates.push_back(std::pair<unsigned int,unsigned int>(i_hyp,cuts_passed));
	 }
	 if ( !candidates.empty() ){
	   int best = bestHypothesis(candidates);
	   for ( unsigned int i=0; i<candidates.size(); ++i){
	     unsigned int i_hyp = candidates.at(i).first;
	     if (!selectBestCandidate|| int(i_hyp)==best){
	       FillSmurfNtuple(smurfTree,i_hyp,weight,sample,candidates.at(i).second);
	       smurfTree.tree_->Fill();
	     }
	   }
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
    delete tree;
    delete f;
  }
  //monitor.print();
  //monitor.makeHistograms(prefix.c_str());
  // monitor.printEvents(3);
  if ( nEventsChain != nEventsTotal ) {
    printf("ERROR: number of events from files (%d) is not equal to total number"
	      " of events (%d)\n", nEventsChain, nEventsTotal);
  }
  printf("Total number of skipped events due to bad identification: %d (%0.0f %%)\n",   
	   nFailedIdentification, nFailedIdentification*100.0/(nEventsChain+1e-5));
  printf("Total number of filtered out events: %d (%0.0f %%)\n",   
	   nFilteredOut, nFilteredOut*100.0/(nEventsChain+1e-5));

  fSmurf->cd(); 
  smurfTree.tree_->Write();
  smurfTree.info_.SetTitle(config_info);
  smurfTree.info_.Write();
  fSmurf->Close();

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
}
EventIdentifier::EventIdentifier(CMS2& cms2, bool isData){
  data = isData;
  run = cms2.evt_run();
  event = cms2.evt_event();
  lumi = cms2.evt_lumiBlock();
  trks_d0 = cms2.trks_d0().empty() ? 0 : cms2.trks_d0().front();
}
   
bool EventIdentifier::operator < (const EventIdentifier &other) const
{
     if (run != other.run)  return run < other.run;
     if (event != other.event) return event < other.event;
     if (data) return false;
     if (fabs(trks_d0 - other.trks_d0) > 1e-6 * trks_d0) return trks_d0 < other.trks_d0;
     return false;
}

bool EventIdentifier::operator == (const EventIdentifier &other) const
{
     if (run != other.run) return false;
     if (event != other.event) return false;
     if (data) return true;
     if (fabs(trks_d0 - other.trks_d0) > 1e-6 * trks_d0)  return false;
     return true;
}

// filter events by process
bool filterByProcess( enum SmurfTree::DataType sample ) {
  switch (sample) {
  case SmurfTree::dyee: 
    return isDYee();
  case SmurfTree::dymm:
    return isDYmm();
  case SmurfTree::dytt:
    return isDYtt();
  case SmurfTree::qqww:
    return isWW();
  case SmurfTree::wz:
    return isWZ();
  case SmurfTree::zz:
    return isZZ();
  default:
    return true;
  }
}
  
bool isIdentified( enum SmurfTree::DataType sample ) {
  switch (sample) {
  case SmurfTree::dyee:
  case SmurfTree::dymm:
  case SmurfTree::dytt:
    return getDrellYanType()!=999;
  case SmurfTree::qqww:
  case SmurfTree::wz:
  case SmurfTree::zz:
    return getVVType()!=999;
  default:
    return true;
  }
}

void ProcessSample( std::string file_pattern, 
		    SmurfTree::DataType sample, 
		    double integratedLumi,
		    double xsec,
		    int nProcessedEvents,
		    bool identifyEvents,
		    bool realData,
		    TString cms2_json_file)
{
  std::vector<string> vec;
  vec.push_back(file_pattern);
  ProcessSample(vec,sample,integratedLumi,xsec,nProcessedEvents,identifyEvents,realData,cms2_json_file);
}

void ProcessSample( std::vector<std::string> file_patterns, 
		    SmurfTree::DataType sample, 
		    double integratedLumi,
		    double xsec,
		    int nProcessedEvents,
		    bool identifyEvents,
		    bool realData,
		    TString cms2_json_file)
{
  electronIdMVA = new ElectronIDMVA();
  electronIdMVA->Initialize("BDTG method", 2,
			   "./files/Subdet0LowPt_WithIPInfo_BDTG.weights.xml",
			   "./files/Subdet1LowPt_WithIPInfo_BDTG.weights.xml",
			   "./files/Subdet2LowPt_WithIPInfo_BDTG.weights.xml",
			   "./files/Subdet0HighPt_WithIPInfo_BDTG.weights.xml",
			   "./files/Subdet1HighPt_WithIPInfo_BDTG.weights.xml",
			   "./files/Subdet2HighPt_WithIPInfo_BDTG.weights.xml");

  TChain *tchain = new TChain("Events");
  for ( std::vector<std::string>::const_iterator pattern = file_patterns.begin();
	pattern != file_patterns.end(); ++pattern )
    tchain->Add(pattern->c_str());

  std::cout << "Processing " << SmurfTree::name(sample) << ".." << std::endl;
  ScanChain(tchain,sample,integratedLumi,xsec,nProcessedEvents,identifyEvents,realData,cms2_json_file);

  delete electronIdMVA;
  delete tchain;
}

bool CheckCutsNM1(wwcuts_t apply, wwcuts_t remove, wwcuts_t passed)
{           
  if ((passed & (apply & (~remove))) == (apply & (~remove))) return true;
  return false;
}   

bool CheckCuts(wwcuts_t apply, wwcuts_t passed)
{
  if ((apply & passed) == apply) return true;
  return false;
}

float getHiggsPt() {
  for (unsigned int i=0; i<cms2.genps_id().size(); ++i) {
    if (cms2.genps_status().at(i) == 3 && cms2.genps_id().at(i) == 25) {
      return cms2.genps_p4().at(i).pt();
    }
  }
  return -1.;
}

float getHiggsPtWeight(float pt, int higgsMass) {
  if (pt<0) return 1.;
  if (HiggsPtKFactor!=0 && TString(HiggsPtKFactor->GetName()).Contains(Form("%i",higgsMass))) 
    return HiggsPtKFactor->GetBinContent( HiggsPtKFactor->GetXaxis()->FindFixBin(pt));
  TFile *fHiggsPtKFactorFile = TFile::Open("./files/ggHWW_KFactors_PowhegToHQT.root");
  assert(fHiggsPtKFactorFile);
  char kfactorHistName[100];
  sprintf(kfactorHistName, "KFactor_PowhegToHQT_mH%d", higgsMass);
  HiggsPtKFactor = (TH1D*)(fHiggsPtKFactorFile->Get(kfactorHistName));
  if (HiggsPtKFactor) HiggsPtKFactor->SetDirectory(0);
  assert(HiggsPtKFactor);
  fHiggsPtKFactorFile->Close();
  delete fHiggsPtKFactorFile;
  //cout << "Using new hist for higgs pT k-factors: " << HiggsPtKFactor->GetName() << endl;
  return HiggsPtKFactor->GetBinContent( HiggsPtKFactor->GetXaxis()->FindFixBin(pt));
}

float getZPt() {
  for (unsigned int i = 0; i < cms2.genps_id().size(); ++i) {
    if ( cms2.genps_id().at(i) == 23 ){
      return cms2.genps_p4().at(i).pt();
    }
  }
  return -999.;
}

float getZRapidity() {
  for (unsigned int i = 0; i < cms2.genps_id().size(); ++i) {
    if ( cms2.genps_id().at(i) == 23 ){
      return cms2.genps_p4().at(i).Rapidity();
    }
  }
  return -999.;
}

float getZMass() {
  for (unsigned int i = 0; i < cms2.genps_id().size(); ++i) {
    if ( cms2.genps_id().at(i) == 23 ){
      return cms2.genps_p4().at(i).mass();
    }
  }
  return -999.;
}

float getDYNNLOWeight(float pt, float rap, float mass) {

  float DYNNLOKFactor = 1.;
  //Find Mass Bin
  const UInt_t nMassBins = 41;
  const Double_t massBinLimits[nMassBins+1] = {15,   20,  25,  30,  35,  40,  45,  50,  55,  60, 
					       64,   68,  72,  76,  81,  86,  91,  96, 101, 106, 
					       110, 115, 120, 126, 133, 141, 150, 160, 171, 185, 
					       200, 220, 243, 273, 320, 380, 440, 510, 600, 1000, 
					       1500}; //41 bins
  //const UInt_t nMassBins = 13;
  //const Double_t massBinLimits[nMassBins+1] = {0,20,30,40,50,60,76,86,96,106,120,150,200,600}; // 13 bins
  Int_t massBinIndex = -1 ;
  for(UInt_t binIndex=0; binIndex < nMassBins; ++binIndex){
    if( mass >= massBinLimits[binIndex] && mass < massBinLimits[binIndex+1]) {
      massBinIndex = binIndex;
      break;
    }
  }
  //Found the mass bin
  if (massBinIndex >= 0 && massBinIndex < Int_t(nMassBins) ) {
    UInt_t ptBin = fDYNNLOKFactorHists[massBinIndex]->GetXaxis()->FindFixBin(pt);
    UInt_t yBin = fDYNNLOKFactorHists[massBinIndex]->GetYaxis()->FindFixBin(rap);

    if(Int_t(ptBin) == fDYNNLOKFactorHists[massBinIndex]->GetNbinsX() + 1)
      ptBin = fDYNNLOKFactorHists[massBinIndex]->GetNbinsX();
    if(ptBin == 0)
      ptBin = 1;
    if(Int_t(yBin) == fDYNNLOKFactorHists[massBinIndex]->GetNbinsY() + 1)
      yBin = fDYNNLOKFactorHists[massBinIndex]->GetNbinsY();
    if(yBin == 0)
      yBin = 1;
    DYNNLOKFactor = fDYNNLOKFactorHists[massBinIndex]->GetBinContent( ptBin, yBin);
  } 
  return DYNNLOKFactor;

}
