const char* config_info = "SmurfV8 selection (Baseline;Tight+Loose;MET20); 42X"; //Skim1
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
//#include "../Tools/ElectronIDMVA.h"
#include "../Tools/MuonIDMVA.h"
#include "../Tools/EGammaMvaEleEstimator.h"
#include "../Tools/MuonMVAEstimator.h"
#include "../Tools/MuonEffectiveArea.h"
#include "TTreeCache.h"

#include "../HWW2012CORE/analysisObjects.h"
#include "../HWW2012CORE/analysisTools.h"
#include "../HWW2012CORE/analysisSelections.h"
#include "../HWW2012CORE/DYMVA.h"

using namespace std;
using namespace tas;

#ifndef __CINT__
#include "CORE/CMS2.h"
#include "CORE/electronSelections.h"
#include "CORE/muonSelections.h"
#include "CORE/jetSelections.h"
#include "CORE/metSelections.h"
#include "jetcorr/FactorizedJetCorrector.h"
#endif

// DEFAULT
//wwcuts_t pass_all = PASSED_BaseLine | PASSED_Charge | PASSED_ZVETO | PASSED_MET | PASSED_JETVETO | PASSED_LT_FINAL | PASSED_LL_FINAL | PASSED_SOFTMUVETO | PASSED_EXTRALEPTONVETO | PASSED_TOPVETO;

// wwcuts_t pass_all = PASSED_Skim1;
wwcuts_t pass_all = PASSED_Skim3; // Baseline, Tight+Fakeable, no MET requirement

//wwcuts_t pass_all = PASSED_BaseLine | PASSED_Charge | PASSED_ZVETO | PASSED_MET | PASSED_LT_FINAL | PASSED_LL_FINAL | PASSED_TopControlSample;
//wwcuts_t pass_all = PASSED_BaseLine;

// wwcuts_t pass_all = PASSED_BaseLine | PASSED_Charge | PASSED_LT_FINAL | PASSED_LL_FINAL | PASSED_ZControlSampleTight;

bool applyJEC = true;
bool lockToCoreSelectors = false;
bool applyFastJetCorrection = false;
bool selectBestCandidate = true; // select only one hypothesis per event with the two most energetic leptons
bool useLHeleId = false;
int useMVAeleId = 1;//zero means off, otherwise it's the mva version
bool useMVAmuId = false;
bool doDYNNLOw = false;
const unsigned int prescale = 1; // DON'T USE ANYTHING BUT 1, unless you know what you are doing

std::vector<std::string> jetcorr_filenames_pfL1FastJetL2L3;
FactorizedJetCorrector *jet_corrector_pfL1FastJetL2L3;
wwcuts_t cuts_passed = 0;

//
// Key analysis method implementation
//

TH1D* HiggsPtKFactor = 0;
//ElectronIDMVA* electronIdMVA = 0;
MuonIDMVA* muonIdMVA = 0;
EGammaMvaEleEstimator* egammaMvaEleEstimator = 0;
MuonMVAEstimator* muonMVAEstimator = 0;
vector<TH2D*>     fDYNNLOKFactorHists;           //vector of hist for Drell-Yan NNLO Kfactor

hypo_monitor monitor(true);


bool hypoSync(int i_hyp, double weight, bool realData) 
{

  HypothesisType type = getHypothesisTypeNew(i_hyp);

  if (nGoodVertex()<1) return false;
  
  monitor.count(cms2,type,"basic selection",weight);

  if (cms2.hyp_lt_id()[i_hyp]*cms2.hyp_ll_id()[i_hyp]>0) return false;
  if ( std::min(cms2.hyp_lt_p4().at(i_hyp).pt(),cms2.hyp_ll_p4().at(i_hyp).pt())<10 ) return false;
  if ( std::max(cms2.hyp_lt_p4().at(i_hyp).pt(),cms2.hyp_ll_p4().at(i_hyp).pt())<20 ) return false;
  
  if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && std::fabs(cms2.mus_p4().at(cms2.hyp_ll_index()[i_hyp]).eta()) >= 2.4) return false;
  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && std::fabs(cms2.mus_p4().at(cms2.hyp_lt_index()[i_hyp]).eta()) >= 2.4) return false;
  if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11 && std::fabs(cms2.els_p4().at(cms2.hyp_ll_index()[i_hyp]).eta()) >= 2.5) return false;
  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11 && std::fabs(cms2.els_p4().at(cms2.hyp_lt_index()[i_hyp]).eta()) >= 2.5) return false;
  monitor.count(cms2,type,"charge + pt 20/10 + eta",weight);


  // 
  std::vector<Int_t> nullMu; // null identified muons 
  std::vector<Int_t> nullEle; // null identified electrons 
  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && !ww_muIso(cms2.hyp_lt_index()[i_hyp], muonMVAEstimator, nullMu, nullEle) ) return false;  
  if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && !ww_muIso(cms2.hyp_ll_index()[i_hyp], muonMVAEstimator, nullMu, nullEle) ) return false; 
  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11 && !ww_elIso(cms2.hyp_lt_index()[i_hyp]) ) return false;
  if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11 && !ww_elIso(cms2.hyp_ll_index()[i_hyp]) ) return false;
  monitor.count(cms2,type,"lepton iso",weight);


  //
  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11 && cms2.els_exp_innerlayers().at(cms2.hyp_lt_index()[i_hyp]) != 0) return false;
  if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11 && cms2.els_exp_innerlayers().at(cms2.hyp_ll_index()[i_hyp]) != 0) return false;
  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11 && isFromConversionMIT(cms2.hyp_lt_index()[i_hyp])) return false;
  if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11 && isFromConversionMIT(cms2.hyp_ll_index()[i_hyp])) return false;
  monitor.count(cms2,type,"conv rejection",weight);


  //
  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && !ww_mud0PV(cms2.hyp_lt_index()[i_hyp]) ) return false;
  if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && !ww_mud0PV(cms2.hyp_ll_index()[i_hyp]) ) return false;
  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11 && !ww_eld0PV(cms2.hyp_lt_index()[i_hyp]) ) return false;
  if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11 && !ww_eld0PV(cms2.hyp_ll_index()[i_hyp]) ) return false;
  monitor.count(cms2,type,"d0",weight);


  //
  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && !ww_mudZPV(cms2.hyp_lt_index()[i_hyp]) ) return false;
  if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && !ww_mudZPV(cms2.hyp_ll_index()[i_hyp]) ) return false;
  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11 && !ww_eldZPV(cms2.hyp_lt_index()[i_hyp]) ) return false;
  if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11 && !ww_eldZPV(cms2.hyp_ll_index()[i_hyp]) ) return false;
  monitor.count(cms2,type,"dZ",weight);
 
 
  //
  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && !ww_muId(cms2.hyp_lt_index()[i_hyp],useMVAmuId,muonIdMVA) ) return false;
  if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && !ww_muId(cms2.hyp_ll_index()[i_hyp],useMVAmuId,muonIdMVA) ) return false;
  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11) {
    int index = cms2.hyp_lt_index()[i_hyp];   
    if (useLHeleId) {
      if (cms2.els_p4().at(index).pt()>20 && (passLikelihoodId_v2(index,cms2.els_lh().at(index),0) & (1<<ELEID_ID))!=(1<<ELEID_ID) ) return false; 
      if (cms2.els_p4().at(index).pt()<20 && (passLikelihoodId_v2(index,cms2.els_lh().at(index),0) & (1<<ELEID_ID))!=(1<<ELEID_ID) ) return false;
    }
    if (useMVAeleId>0){
      if (!goodElectronTMVA(egammaMvaEleEstimator, useMVAeleId, index) ) return false;
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
    if (useMVAeleId>0){
      if (!goodElectronTMVA(egammaMvaEleEstimator, useMVAeleId, index)) return false;
    } else {
      if (!pass_electronSelection(index, electronSelection_smurfV3_id, false, false) ) return false;
    }
  } 
  monitor.count(cms2,type,"lepton id",weight);


  //
  if ( TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && 
  	   !goodMuonIsolated(cms2.hyp_lt_index()[i_hyp], lockToCoreSelectors, useMVAmuId, muonIdMVA, muonMVAEstimator,  nullMu, nullEle) ) return false;  
  if ( TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && 
       !goodMuonIsolated(cms2.hyp_ll_index()[i_hyp], lockToCoreSelectors, useMVAmuId, muonIdMVA, muonMVAEstimator,  nullMu, nullEle) ) return false;   
  if ( TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11 && 
  	   !goodElectronIsolated(cms2.hyp_lt_index()[i_hyp], useLHeleId, useMVAeleId, egammaMvaEleEstimator, lockToCoreSelectors) ) return false;
  if ( TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11 && 
  	   !goodElectronIsolated(cms2.hyp_ll_index()[i_hyp], useLHeleId, useMVAeleId, egammaMvaEleEstimator, lockToCoreSelectors) ) return false;
  monitor.count(cms2,type,"lepton id/iso pt20/10",weight);


  //
  if ( metValue()<20 ) return false;
  monitor.count(cms2,type,"met>20",weight);


  //
  if (cms2.hyp_p4().at(i_hyp).mass2()<0 || cms2.hyp_p4()[i_hyp].mass() < 12) return false;
  monitor.count(cms2,type,"mll>12 (20->12 for SF)",weight);


  //
  if ( type == EE || type == MM) {
    if (inZmassWindow(cms2.hyp_p4()[i_hyp].mass())) return false;
  }
  monitor.count(cms2,type,"z veto",weight);


  // 
  if ( minmet(i_hyp) < 20 )  return false; 
  monitor.count(cms2,type,"minMet > 20 GeV",weight);


  //
  int njets = numberOfJets(i_hyp, applyJEC, jet_corrector_pfL1FastJetL2L3); 
  std::vector<JetPair> sortedJets = getJets(jetType(), i_hyp, 15.0, 4.7, applyJEC, jet_corrector_pfL1FastJetL2L3, true, false); // new cut
  if (cms2.hyp_type().at(i_hyp)!=1 && cms2.hyp_type().at(i_hyp)!=2 && sortedJets.size()>0) {
 	  if (njets>= 2 && fabs(ROOT::Math::VectorUtil::DeltaPhi(sortedJets[0].first+sortedJets[1].first,cms2.hyp_p4().at(i_hyp)))>=165.*TMath::Pi()/180.) return false;
  } 
  monitor.count(cms2,type,"Jet DPhi only for 2jet bin",weight);
  

  //
  if (numberOfSoftMuons(i_hyp,true)>0) return false;
  monitor.count(cms2,type,"soft muon veto",weight);


  //
  if (numberOfExtraLeptons(i_hyp,10, useLHeleId, useMVAeleId, egammaMvaEleEstimator, useMVAmuId, muonIdMVA,muonMVAEstimator,  nullMu, nullEle)>0) return false;
  monitor.count(cms2,type,"extra lepton veto",weight);


  //
  if ( toptag(jetType(),i_hyp, 10, jet_corrector_pfL1FastJetL2L3) )  return false;
  monitor.count(cms2,type,"top tag",weight);


  //
  if ( cms2.hyp_p4().at(i_hyp).pt() < 45 ) return false;
  monitor.count(cms2,type,"dilep.pt()>45",weight);
 

  // 
  if ( !passedMetRequirements(i_hyp,jet_corrector_pfL1FastJetL2L3) ) return false;
  monitor.count(cms2,type,"Full MET cuts",weight);


  //
  if ( njets==0 )
    monitor.count(cms2,type,"0 jet bin",weight);
  
  //
  if ( njets==1 )
    monitor.count(cms2,type,"1 jet bin",weight);

  //
  if ( njets==2 || njets==3 ) {
    if (fabs(sortedJets[0].first.eta())>=4.7 || fabs(sortedJets[1].first.eta())>=4.7) return false;
    if (njets==3 && sortedJets[2].first.pt()>30 && ((sortedJets[0].first.eta()-sortedJets[2].first.eta() > 0 && sortedJets[1].first.eta()-sortedJets[2].first.eta() < 0) ||
					(sortedJets[1].first.eta()-sortedJets[2].first.eta() > 0 && sortedJets[0].first.eta()-sortedJets[2].first.eta() < 0)) ) return false;
    monitor.count(cms2,type,"2 jet bin",weight);
  }
  

  //
  if ( njets!=0 ) return false;
  monitor.count(cms2,type,"JetVeto cuts",weight);


  //
  if ( njets!=0 ) return false;
  if ( TMath::Max(cms2.hyp_ll_p4().at(i_hyp).pt(),cms2.hyp_lt_p4().at(i_hyp).pt())<=30. )  return false;
  monitor.count(cms2,type,"ptMax",weight);


  //
  if ( njets!=0 ) return false;
  if ( TMath::Min(cms2.hyp_ll_p4().at(i_hyp).pt(),cms2.hyp_lt_p4().at(i_hyp).pt())<=25. )  return false;
  monitor.count(cms2,type,"ptMin",weight);


  //
  if ( njets!=0 ) return false;
  if ( cms2.hyp_p4().at(i_hyp).mass()>=50. )  return false;
  monitor.count(cms2,type,"mll",weight);


  //
  if ( njets!=0 ) return false;
  if ( mt(cms2.hyp_p4().at(i_hyp).pt(),cms2.evt_pfmet(),acos(cos(cms2.hyp_p4().at(i_hyp).phi()-metPhiValue())))<80.)  return false;
  monitor.count(cms2,type,"mT",weight);


  //
  if ( njets!=0 ) return false;
  if ( fabs(ROOT::Math::VectorUtil::DeltaPhi(cms2.hyp_ll_p4().at(i_hyp),cms2.hyp_lt_p4().at(i_hyp)))>=60.*TMath::Pi()/180. )  return false;
  monitor.count(cms2,type,"dphill",weight);

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
  monitor.count(cms2, type, "Total number before cuts");
  
  // if ( cms2.hyp_FVFit_prob()[i_hyp] < 0.005 ) return;
  // monitor.count(cms2, type, "after vertex cut");
  
  if (cms2.hyp_p4().at(i_hyp).mass2()<0) return false;
  if (!isGoodVertex(primaryVertex())) return false;
  // if ( realData && ! passedTriggerRequirements( hypType(i_hyp) ) )return false;
  // if ( realData && passedTriggerRequirements( hypType(i_hyp) ) )return false;
  monitor.count(cms2, type, "trigger requirements",weight);
  if (nGoodVertex()<1) return false;

  cuts_passed = PASSED_BaseLine;
  if ( realData && passedTriggerRequirements() )  cuts_passed |= PASSED_Trigger;
 
  // Baseline cuts
  if (abs(cms2.hyp_lt_id()[i_hyp]) == 13 && !ww_muBase(cms2.hyp_lt_index()[i_hyp]) ) cuts_passed &= ~PASSED_BaseLine;
  if (abs(cms2.hyp_ll_id()[i_hyp]) == 13 && !ww_muBase(cms2.hyp_ll_index()[i_hyp]) ) cuts_passed &= ~PASSED_BaseLine;
  if (abs(cms2.hyp_lt_id()[i_hyp]) == 11 && !ww_elBase(cms2.hyp_lt_index()[i_hyp]) ) cuts_passed &= ~PASSED_BaseLine;
  if (abs(cms2.hyp_ll_id()[i_hyp]) == 11 && !ww_elBase(cms2.hyp_ll_index()[i_hyp]) ) cuts_passed &= ~PASSED_BaseLine;
  
  // Require opposite sign
  if ( cms2.hyp_lt_id()[i_hyp] * cms2.hyp_ll_id()[i_hyp] < 0 ) cuts_passed |= PASSED_Charge;

  //monitor.count(cms2,type,"baseline cuts",weight);
 
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
  if ( passedMetRequirements(i_hyp, jet_corrector_pfL1FastJetL2L3 ) ) cuts_passed |= PASSED_MET;

  std::vector<Int_t> nullMu; // null identified muons 
  std::vector<Int_t> nullEle; // null identified electrons  
  
  // == letpon ID and Isolation
  // Muon quality cuts, including isolation
  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13){
    unsigned int index = cms2.hyp_lt_index()[i_hyp];
    if ( goodMuonIsolated(index, lockToCoreSelectors, useMVAmuId, muonIdMVA, muonMVAEstimator,  nullMu, nullEle) )   cuts_passed |= PASSED_LT_FINAL;
//    if ( fakableMuon(index,MuFOV1, muonMVAEstimator,  nullMu, nullEle) && !goodMuonIsolated(index, lockToCoreSelectors, useMVAmuId, muonIdMVA, muonMVAEstimator,  nullMu, nullEle)) cuts_passed |= PASSED_LT_FO_MU1;
    if ( fakableMuon(index,MuFOV2, muonMVAEstimator,  nullMu, nullEle) && !goodMuonIsolated(index, lockToCoreSelectors, useMVAmuId, muonIdMVA, muonMVAEstimator,  nullMu, nullEle)) cuts_passed |= PASSED_LT_FO_MU2;
  }
  if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13){
    unsigned int index = cms2.hyp_ll_index()[i_hyp];
    if ( goodMuonIsolated(index, lockToCoreSelectors, useMVAmuId, muonIdMVA, muonMVAEstimator,  nullMu, nullEle) )   cuts_passed |= PASSED_LL_FINAL;
//    if ( fakableMuon(index, MuFOV1, muonMVAEstimator,  nullMu, nullEle) && !goodMuonIsolated(index, lockToCoreSelectors, useMVAmuId, muonIdMVA, muonMVAEstimator,  nullMu, nullEle)) cuts_passed |= PASSED_LL_FO_MU1;
    if ( fakableMuon(index, MuFOV2, muonMVAEstimator,  nullMu, nullEle) && !goodMuonIsolated(index, lockToCoreSelectors, useMVAmuId, muonIdMVA, muonMVAEstimator,  nullMu, nullEle)) cuts_passed |= PASSED_LL_FO_MU2;
  } 
  // Electron quality cuts, including isolation
  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11){
    unsigned int index = cms2.hyp_lt_index()[i_hyp];
    if ( goodElectronIsolated(index, useLHeleId, useMVAeleId, egammaMvaEleEstimator, lockToCoreSelectors) ) cuts_passed |= PASSED_LT_FINAL;
//    if ( fakableElectron(index,EleFOV1) && !goodElectronIsolated(index, useLHeleId, useMVAeleId, egammaMvaEleEstimator, lockToCoreSelectors) )   cuts_passed |= PASSED_LT_FO_ELEV1;
//    if ( fakableElectron(index,EleFOV2) && !goodElectronIsolated(index, useLHeleId, useMVAeleId, egammaMvaEleEstimator, lockToCoreSelectors))   cuts_passed |= PASSED_LT_FO_ELEV2;
//    if ( fakableElectron(index,EleFOV3) && !goodElectronIsolated(index, useLHeleId, useMVAeleId, egammaMvaEleEstimator, lockToCoreSelectors))   cuts_passed |= PASSED_LT_FO_ELEV3;
    if ( fakableElectron(index,EleFOV4) && !goodElectronIsolated(index, useLHeleId, useMVAeleId, egammaMvaEleEstimator, lockToCoreSelectors))   cuts_passed |= PASSED_LT_FO_ELEV4;
  }
  if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11){
    unsigned int index = cms2.hyp_ll_index()[i_hyp];
    if ( goodElectronIsolated(index, useLHeleId, useMVAeleId, egammaMvaEleEstimator, lockToCoreSelectors) ) cuts_passed |= PASSED_LL_FINAL;
//    if ( fakableElectron(index,EleFOV1) && !goodElectronIsolated(index, useLHeleId, useMVAeleId, egammaMvaEleEstimator, lockToCoreSelectors))   cuts_passed |= PASSED_LL_FO_ELEV1;
//    if ( fakableElectron(index,EleFOV2) && !goodElectronIsolated(index, useLHeleId, useMVAeleId, egammaMvaEleEstimator, lockToCoreSelectors))   cuts_passed |= PASSED_LL_FO_ELEV2;
//    if ( fakableElectron(index,EleFOV3) && !goodElectronIsolated(index, useLHeleId, useMVAeleId, egammaMvaEleEstimator, lockToCoreSelectors))   cuts_passed |= PASSED_LL_FO_ELEV3;
    if ( fakableElectron(index,EleFOV4) && !goodElectronIsolated(index, useLHeleId, useMVAeleId, egammaMvaEleEstimator, lockToCoreSelectors) )   cuts_passed |= PASSED_LL_FO_ELEV4;
  }

  // if ( passedLTFinalRequirements || passedLLFinalRequirements ) cuts_passed |= (1<<PASS_PROBE);

  // == Jet-veto
  const std::vector<JetPair>& vetojets(getDefaultJets(i_hyp, applyJEC, jet_corrector_pfL1FastJetL2L3, false));
  unsigned int nJets = vetojets.size();
  if (nJets==0) cuts_passed |= PASSED_JETVETO;
  if (nJets>1 && 
      getDefaultJets(i_hyp, applyJEC, jet_corrector_pfL1FastJetL2L3, true).size()>0 )  
    cuts_passed |= PASSED_TopControlSample;
  if (nJets==1){
    if ( getDefaultJets(i_hyp, applyJEC, jet_corrector_pfL1FastJetL2L3, true).size()==1 ){
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
  if ( numberOfExtraLeptons(i_hyp,10, useLHeleId, useMVAeleId, egammaMvaEleEstimator, useMVAmuId, muonIdMVA,muonMVAEstimator,  nullMu, nullEle) == 0) cuts_passed |= PASSED_EXTRALEPTONVETO;
  if ( ! toptag(jetType(),i_hyp,10,jet_corrector_pfL1FastJetL2L3,vetojets) ) cuts_passed |= PASSED_TOPVETO_NotInJets;//newcuts, was 7
  if ( ! toptag(jetType(),i_hyp,10,jet_corrector_pfL1FastJetL2L3) )          cuts_passed |= PASSED_TOPVETO;//newcuts, was 7

  if( CheckCuts(  PASSED_MET | PASSED_BaseLine , cuts_passed)) {
    if ( CheckCuts( PASSED_LT_FINAL | PASSED_LL_FINAL,   cuts_passed ) ||
     CheckCuts( PASSED_LT_FINAL | PASSED_LL_FO_MU1,   cuts_passed ) ||
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
    if ( CheckCuts( PASSED_LT_FINAL | PASSED_LL_FINAL,   cuts_passed ) ||
     CheckCuts( PASSED_LT_FINAL | PASSED_LL_FO_MU1,   cuts_passed ) ||
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

  if(! CheckCuts(pass_all, cuts_passed)) return false;

  monitor.count(cms2,type,"all cuts (including soft and extra lepton veto)",weight);
  return true;

}//end of void hypo

void FillSmurfNtuple(SmurfTree& tree, unsigned int i_hyp, 
		     double weight, enum SmurfTree::DataType sample, wwcuts_t cuts_passed){


  tree.InitVariables();
  tree.run_         = cms2.evt_run();
  tree.event_       = cms2.evt_event();
  tree.lumi_        = cms2.evt_lumiBlock();
  tree.nvtx_        = nGoodVertex();
  tree.scale1fb_    = weight;
  tree.met_         = metValue();
  tree.metSig_      = cms2.evt_pfmetSignificance();
  tree.sumet_       = sumetValue();
  tree.metPhi_      = metPhiValue();
//  metStruct trkMET  = trackerMET(i_hyp,0.1);
//  tree.trackMet_    = trkMET.met;
//  tree.trackMetPhi_ = trkMET.metphi;
//  tree.pTrackMet_   = projectedMet(i_hyp, trkMET.met, trkMET.metphi);
  tree.trackMet_    = cms2.trk_met()[i_hyp];
  tree.trackMetPhi_ = cms2.trk_metPhi()[i_hyp];
  tree.pTrackMet_   = projectedMet(i_hyp, tree.trackMet_, tree.trackMetPhi_);

  // PDF stuffs
  if (sample!=SmurfTree::data && sample!=SmurfTree::dyttDataDriven){
	  tree.Q_    	= cms2.pdfinfo_scale();
	  tree.id1_   	= cms2.pdfinfo_id1();
	  tree.x1_    	= cms2.pdfinfo_x1(); 
	  tree.pdf1_    = cms2.pdfinfo_pdf1();
	  tree.id2_    	= cms2.pdfinfo_id2(); 
	  tree.x2_    	= cms2.pdfinfo_x2(); 
	  tree.pdf2_    = cms2.pdfinfo_pdf2(); 
  } 

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

  // SC eta for electrons and  mu p4 eta for muons
  if (ltIsFirst) {
    tree.lep1DetEta_ =  abs(cms2.hyp_lt_id().at(i_hyp))==11 ?
						cms2.els_etaSC().at(cms2.hyp_lt_index().at(i_hyp)) :  cms2.mus_p4().at(cms2.hyp_lt_index().at(i_hyp)).eta(); 
    tree.lep2DetEta_ =  abs(cms2.hyp_ll_id().at(i_hyp))==11 ?
						cms2.els_etaSC().at(cms2.hyp_ll_index().at(i_hyp)) :  cms2.mus_p4().at(cms2.hyp_ll_index().at(i_hyp)).eta(); 
  } else {
    tree.lep1DetEta_ =  abs(cms2.hyp_ll_id().at(i_hyp))==11 ?
						cms2.els_etaSC().at(cms2.hyp_ll_index().at(i_hyp)) :  cms2.mus_p4().at(cms2.hyp_ll_index().at(i_hyp)).eta(); 
    tree.lep2DetEta_ =  abs(cms2.hyp_lt_id().at(i_hyp))==11 ?
						cms2.els_etaSC().at(cms2.hyp_lt_index().at(i_hyp)) :  cms2.mus_p4().at(cms2.hyp_lt_index().at(i_hyp)).eta(); 
  }

  std::vector<Int_t> nullMu; // null identified muons 
  std::vector<Int_t> nullEle; // null identified electrons  
  
  if (ltIsFirst) {
    tree.lmva1_ = abs(cms2.hyp_lt_id().at(i_hyp))==11 ? 
      egammaMvaEleEstimator->mvaValue(cms2.hyp_lt_index().at(i_hyp), false) : 
      muonMVAEstimator->mvaValueIso( cms2.hyp_lt_index().at(i_hyp), cms2.evt_ww_rho(), MuonEffectiveArea::kMuEAFall11MC, nullEle, nullMu, false );
    tree.lmva2_ = abs(cms2.hyp_ll_id().at(i_hyp))==11 ? 
      egammaMvaEleEstimator->mvaValue(cms2.hyp_ll_index().at(i_hyp), false) : 
      muonMVAEstimator->mvaValueIso( cms2.hyp_ll_index().at(i_hyp), cms2.evt_ww_rho(), MuonEffectiveArea::kMuEAFall11MC, nullEle, nullMu, false );
  } else {
    tree.lmva1_ = abs(cms2.hyp_ll_id().at(i_hyp))==11 ? 
      egammaMvaEleEstimator->mvaValue(cms2.hyp_ll_index().at(i_hyp), false) : 
      muonMVAEstimator->mvaValueIso( cms2.hyp_ll_index().at(i_hyp), cms2.evt_ww_rho(), MuonEffectiveArea::kMuEAFall11MC, nullEle, nullMu, false );
    tree.lmva2_ = abs(cms2.hyp_lt_id().at(i_hyp))==11 ? 
      egammaMvaEleEstimator->mvaValue(cms2.hyp_lt_index().at(i_hyp), false) : 
      muonMVAEstimator->mvaValueIso( cms2.hyp_lt_index().at(i_hyp), cms2.evt_ww_rho(), MuonEffectiveArea::kMuEAFall11MC, nullEle, nullMu, false );
  }


  const std::vector<JetPair>& jets = getJets(jetType(), i_hyp, 0, 4.7, applyJEC, jet_corrector_pfL1FastJetL2L3, true, false);
  if (jets.size()>0){
    tree.jet1_ = jets.at(0).first;
    tree.jet1Btag_ = BTag(jetType(),jets.at(0).second, jets.at(0).first.pt());
    tree.jet1ProbBtag_ = cms2.pfjets_jetBProbabilityBJetTag()[jets.at(0).second];
    //tree.jet1Dz_ = jetDz(0,0); // FIXME : commented out for embedded sample
  }
  if (jets.size()>1){
    tree.jet2_ = jets.at(1).first;
    tree.jet2Btag_ = BTag(jetType(),jets.at(1).second, jets.at(1).first.pt());
    tree.jet2ProbBtag_ = cms2.pfjets_jetBProbabilityBJetTag()[jets.at(1).second];
  }
  if (jets.size()>2){
    tree.jet3_ = jets.at(2).first;
    tree.jet3Btag_ = BTag(jetType(),jets.at(2).second, jets.at(2).first.pt());
    tree.jet3ProbBtag_ = cms2.pfjets_jetBProbabilityBJetTag()[jets.at(2).second];
  }
  tree.njets_ = numberOfJets(i_hyp, applyJEC, jet_corrector_pfL1FastJetL2L3);
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

  if (sample!=SmurfTree::data && sample!=SmurfTree::dyttDataDriven){
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
	tree.sfWeightHPt_ = getHiggsPtWeight(pt,ss.Atoi(), HiggsPtKFactor);
      }
	}

	if (doDYNNLOw) {
		if (sample==SmurfTree::dyee||sample==SmurfTree::dymm||sample==SmurfTree::dytt) {
			if (isDYee()||isDYmm()||isDYtt()) {
				float pt = getZPt();
				float rap = getZRapidity();
				float mass = getZMass();
				float kfact =  getDYNNLOWeight(pt,rap,mass, fDYNNLOKFactorHists);
				tree.scale1fb_ *= kfact;
			}
		}
	}
   
    // true nPU
    for (unsigned int nbc=0;nbc<cms2.puInfo_trueNumInteractions().size();++nbc) {
      if (cms2.puInfo_bunchCrossing().at(nbc)==0) tree.npu_ = cms2.puInfo_trueNumInteractions().at(nbc);
      else if (cms2.puInfo_bunchCrossing().at(nbc)==-1) tree.npuMinusOne_ = cms2.puInfo_trueNumInteractions().at(nbc);
      else if (cms2.puInfo_bunchCrossing().at(nbc)==+1) tree.npuPlusOne_ = cms2.puInfo_trueNumInteractions().at(nbc);
//    for (unsigned int nbc=0;nbc<cms2.puInfo_nPUvertices().size();++nbc) {
//      if (cms2.puInfo_bunchCrossing().at(nbc)==0) tree.npu_ = cms2.puInfo_nPUvertices().at(nbc);
//      else if (cms2.puInfo_bunchCrossing().at(nbc)==-1) tree.npuMinusOne_ = cms2.puInfo_nPUvertices().at(nbc);
//      else if (cms2.puInfo_bunchCrossing().at(nbc)==+1) tree.npuPlusOne_ = cms2.puInfo_nPUvertices().at(nbc);
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
    if ( tree.jetLowBtag_ < BTag(jetType(),jets.at(i).second, jets.at(i).first.pt()) )
      tree.jetLowBtag_ = BTag(jetType(),jets.at(i).second, jets.at(i).first.pt());

  // third lepton
  const std::vector<LeptonPair>& leptons(getExtraLeptons(i_hyp,0, useLHeleId, useMVAeleId, egammaMvaEleEstimator, useMVAmuId, muonIdMVA,muonMVAEstimator,  nullMu, nullEle));
  for (std::vector<LeptonPair>::const_iterator lep=leptons.begin();
		  lep!=leptons.end(); ++lep){
	  if (lep->first){
		  // muons
		  if ( tree.lep3_.pt() < cms2.mus_p4().at(lep->second).pt() ){
			  tree.lep3_ = cms2.mus_p4().at(lep->second);
			  tree.lq3_   = cms2.mus_charge().at(lep->second);
			  tree.lid3_ = tree.lq3_>0?-13:13;
			  if (sample!=SmurfTree::data && sample!=SmurfTree::dyttDataDriven){// GC
				  tree.lep3McId_ = cms2.mus_mc_id().at(lep->second);
				  tree.lep3MotherMcId_  = cms2.mus_mc_motherid().at(lep->second);
			  }
		  }
	  } else {
		  if ( tree.lep3_.pt() < cms2.els_p4().at(lep->second).pt() ){
			  tree.lep3_ = cms2.els_p4().at(lep->second);
			  tree.lq3_   = cms2.els_charge().at(lep->second);
			  tree.lid3_ = tree.lq3_>0?-11:11;
			  if (sample!=SmurfTree::data && sample!=SmurfTree::dyttDataDriven){// GC
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
		TString cms2_json_file,
		int beginrun, 
		int endrun)
{

  std::string prefix = SmurfTree::name(sample);
  if ( chain->GetListOfFiles()->GetEntries()==0 ){
    printf("\nERROR: chain is empty for sample: %s\n\n",prefix.c_str());
    assert(0);
  }

    TObjArray *listOfFiles = chain->GetListOfFiles();
    TIter fileIter(listOfFiles);
    TChainElement *file0 = ((TChainElement*)fileIter.Next());
    if (prefix=="ttbar" || prefix=="dymm" || prefix=="dyee"){
        if (TString(file0->GetTitle()).Contains("madgraph")) prefix+="_mg";
    }
    if (prefix=="tw"){
        if (TString(file0->GetTitle()).Contains("DS")) prefix+="_ds";
    }
        if (prefix=="qqww"){
            if (TString(file0->GetTitle()).Contains("mcatnlo")) prefix="ww_mcnlo";
            if (TString(file0->GetTitle()).Contains("scaleup")) prefix+="_up";
            if (TString(file0->GetTitle()).Contains("scaledown")) prefix+="_down";
        }
    if (prefix=="wz" || prefix=="zz"){
        if (TString(file0->GetTitle()).Contains("pythia")) prefix+="_py";
    }
    if (prefix=="wgstar"){
        prefix="wg3l";
    }
    if (prefix=="dyttDataDriven"){
        prefix="data-emb-tau123";
    }
    fileIter.Reset();

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
  gSystem->MakeDirectory(Form("smurf", beginrun, endrun));
  TFile* fSmurf = TFile::Open(Form("smurf/%s_%d_%d.root", prefix.c_str(), beginrun, endrun),"RECREATE");
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
      jetcorr_filenames_pfL1FastJetL2L3.push_back("files/GR_P_V42_AN3_L1FastJet_AK5PF.txt");
      jetcorr_filenames_pfL1FastJetL2L3.push_back("files/GR_P_V42_AN3_L2Relative_AK5PF.txt");
      jetcorr_filenames_pfL1FastJetL2L3.push_back("files/GR_P_V42_AN3_L3Absolute_AK5PF.txt");
      jetcorr_filenames_pfL1FastJetL2L3.push_back("files/GR_P_V42_AN3_L2L3Residual_AK5PF.txt");
    } else {
      jetcorr_filenames_pfL1FastJetL2L3.push_back("files/START53_V15_L1FastJet_AK5PF.txt");
      jetcorr_filenames_pfL1FastJetL2L3.push_back("files/START53_V15_L2Relative_AK5PF.txt");
      jetcorr_filenames_pfL1FastJetL2L3.push_back("files/START53_V15_L3Absolute_AK5PF.txt");
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
    } else if (cms2_json_file.Contains(".json")) {
      set_goodrun_file_json(cms2_json_file);
    } else set_goodrun_file(cms2_json_file);
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
	
          if (int(cms2.evt_run())<beginrun || int(cms2.evt_run())>=endrun) continue;
	  if (cms2.evt_event()%prescale!=0) continue;
	  //if (cms2.evt_event()!=20494961 && cms2.evt_event()!=787812949 && cms2.evt_event()!=73349194 && cms2.evt_event()!=68905080) continue;
	  //cout << "event: " << cms2.evt_event() << endl;

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
//	   printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
//	   "\033[0m\033[32m <---\033[0m\015", i_permille/10.);
//	   fflush(stdout);
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
  monitor.print();
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
      ofstream json(Form("processedjson/processed_%d_%d.json", beginrun, endrun));
      ofstream json2(Form("processedjson/processed_detailed_%d_%d.json", beginrun, endrun));
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

void ProcessSample( std::string file_pattern, 
		    SmurfTree::DataType sample, 
		    double integratedLumi,
		    double xsec,
		    int nProcessedEvents,
		    bool identifyEvents,
		    bool realData,
		    TString cms2_json_file,
			int beginrun, 
			int endrun)
{
  std::vector<string> vec;
  vec.push_back(file_pattern);
  ProcessSample(vec,sample,integratedLumi,xsec,nProcessedEvents,identifyEvents,realData,cms2_json_file,beginrun,endrun);
}

void ProcessSample( std::vector<std::string> file_patterns, 
		    SmurfTree::DataType sample, 
		    double integratedLumi,
		    double xsec,
		    int nProcessedEvents,
		    bool identifyEvents,
		    bool realData,
		    TString cms2_json_file,
			int beginrun, 
			int endrun)
{
/*
  if (useMVAeleId>0) {
    electronIdMVA = new ElectronIDMVA();
    if (useMVAeleId==2) {
      electronIdMVA->Initialize("BDTG method", 2,
				"files/Subdet0LowPt_WithIPInfo_BDTG.weights.xml",
				"files/Subdet1LowPt_WithIPInfo_BDTG.weights.xml",
				"files/Subdet2LowPt_WithIPInfo_BDTG.weights.xml",
				"files/Subdet0HighPt_WithIPInfo_BDTG.weights.xml",
				"files/Subdet1HighPt_WithIPInfo_BDTG.weights.xml",
				"files/Subdet2HighPt_WithIPInfo_BDTG.weights.xml");
    } else if (useMVAeleId==3){
      electronIdMVA->Initialize("BDTG method", 3,
				"files/Subdet0LowPt_IDIsoCombined_BDTG.weights.xml",
				"files/Subdet1LowPt_IDIsoCombined_BDTG.weights.xml",
				"files/Subdet2LowPt_IDIsoCombined_BDTG.weights.xml",
				"files/Subdet0HighPt_IDIsoCombined_BDTG.weights.xml",
				"files/Subdet1HighPt_IDIsoCombined_BDTG.weights.xml",
				"files/Subdet2HighPt_IDIsoCombined_BDTG.weights.xml");
    } else {
      cout << "error: electron id mva version not supported" << endl;
      return;
    }
  }
*/
  if (useMVAmuId) {
    muonIdMVA = new MuonIDMVA();
    muonIdMVA->Initialize("BDTG method", 1,
			  "files/BarrelPtBin0_IDIsoCombined_BDTG.weights.xml",
			  "files/EndcapPtBin0_IDIsoCombined_BDTG.weights.xml",
			  "files/BarrelPtBin1_IDIsoCombined_BDTG.weights.xml",
			  "files/EndcapPtBin1_IDIsoCombined_BDTG.weights.xml",
			  "files/BarrelPtBin2_IDIsoCombined_BDTG.weights.xml",
			  "files/EndcapPtBin2_IDIsoCombined_BDTG.weights.xml");
  }

  // --------------- EGamma Id MVA  --------------------------
  vector<std::string> egammaweights;
  egammaweights.push_back("files/Electrons_BDTG_TrigV0_Cat1.weights.xml"); 
  egammaweights.push_back("files/Electrons_BDTG_TrigV0_Cat2.weights.xml"); 
  egammaweights.push_back("files/Electrons_BDTG_TrigV0_Cat3.weights.xml"); 
  egammaweights.push_back("files/Electrons_BDTG_TrigV0_Cat4.weights.xml"); 
  egammaweights.push_back("files/Electrons_BDTG_TrigV0_Cat5.weights.xml"); 
  egammaweights.push_back("files/Electrons_BDTG_TrigV0_Cat6.weights.xml"); 
  egammaMvaEleEstimator = new EGammaMvaEleEstimator();
  egammaMvaEleEstimator->initialize("BDT", EGammaMvaEleEstimator::kTrig, true, egammaweights );
  // --------------- EGamma Id MVA  --------------------------

  // --------------- Muon RingIso MVA  --------------------------
  vector<std::string> muonisoweights;
  muonisoweights.push_back("files/MuonIsoMVA_sixie-BarrelPt5To10_V0_BDTG.weights.xml");
  muonisoweights.push_back("files/MuonIsoMVA_sixie-EndcapPt5To10_V0_BDTG.weights.xml");
  muonisoweights.push_back("files/MuonIsoMVA_sixie-BarrelPt10ToInf_V0_BDTG.weights.xml");
  muonisoweights.push_back("files/MuonIsoMVA_sixie-EndcapPt10ToInf_V0_BDTG.weights.xml");
  muonisoweights.push_back("files/MuonIsoMVA_sixie-Tracker_V0_BDTG.weights.xml");
  muonisoweights.push_back("files/MuonIsoMVA_sixie-Global_V0_BDTG.weights.xml");
  muonMVAEstimator = new MuonMVAEstimator();
  muonMVAEstimator->initialize( "MuonIso_BDTG_IsoRings", MuonMVAEstimator::kIsoRings, true, muonisoweights );
  // --------------- Muon RingIso MVA  --------------------------

  TChain *tchain = new TChain("Events");
  for ( std::vector<std::string>::const_iterator pattern = file_patterns.begin();
	pattern != file_patterns.end(); ++pattern )
    tchain->Add(pattern->c_str());

  std::cout << "Number of events to process : " << tchain->GetEntries() << endl;

  std::cout << "Processing " << SmurfTree::name(sample) << ".." << std::endl;
  ScanChain(tchain,sample,integratedLumi,xsec,nProcessedEvents,identifyEvents,realData,cms2_json_file,beginrun,endrun);

//  if (useMVAeleId>0) delete electronIdMVA;
  if (useMVAmuId)    delete muonIdMVA;
  delete tchain;
  delete muonMVAEstimator;
  delete egammaMvaEleEstimator;
}

