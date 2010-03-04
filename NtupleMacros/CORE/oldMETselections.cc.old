//===========================================================
//
// Various selection functions are kept here
//
//============================================================
#include <assert.h>
#include <algorithm>
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TDatabasePDG.h"
#include "selections.h"

// CMS2 includes
#include "CMS2.h"
#include "utilities.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

//--------------------------------------------
// 'simple' MET selection
//--------------------------------------------
bool metSimple (float threshold, const TVector3& corr) {
  TVector3 hyp_met;
  hyp_met.SetPtEtaPhi(cms2.evt_tcmet(), 0, cms2.evt_tcmetPhi());
  hyp_met += corr;
  if (hyp_met.Pt() < threshold) return false;
  return true;
}

//--------------------------------------------
// pT(ll)/MET ratio selection
//--------------------------------------------
bool metBalance (int i_hyp, const TVector3& corr) {
  TVector3 hyp_met;
  hyp_met.SetPtEtaPhi(cms2.evt_tcmet(), 0, cms2.evt_tcmetPhi());
  hyp_met += corr;
  if(hyp_met.Pt()/cms2.hyp_p4()[i_hyp].pt() < 0.6 &&
      acos(cos(hyp_met.Phi()-cms2.hyp_p4()[i_hyp].phi() - 3.1416)) < 0.25 ) return false;
  return true;
}

//--------------------------------------------
// pMET selection
//--------------------------------------------
bool metProjected (int i_hyp, const TVector3& corr) {
  TVector3 hyp_met;
  hyp_met.SetPtEtaPhi(cms2.evt_tcmet(), 0, cms2.evt_tcmetPhi());
  hyp_met += corr;
  double metspec = MetSpecial(hyp_met.Pt(), hyp_met.Phi(), i_hyp);
  if ( metspec < 20 ) return false;
  return true;
}

//--------------------------------------------
// PASS5 MET
//--------------------------------------------
bool pass5Met (int i_hyp, const TVector3& corr) {
  // for e-e and mu-mu
  if (cms2.hyp_type()[i_hyp] == 0 || cms2.hyp_type()[i_hyp] == 3) {
    if(!metSimple(45.0, corr)) return false;
    if(!metBalance(i_hyp, corr))      return false;
    if(!metProjected(i_hyp, corr))    return false;
  }
  // for e-mu and mu-e
  if (cms2.hyp_type()[i_hyp] == 1 || cms2.hyp_type()[i_hyp] == 2) {
    if(!metSimple(20.0, corr)) return false;
    if(!metProjected(i_hyp, corr))    return false;
  }
  return true;
}

//--------------------------------------------
// Pass 2 MET selection
//--------------------------------------------
bool pass2Met (int i_hyp, const TVector3& corr) {
  // for e-e and mu-mu
  TVector3 hyp_met;
  hyp_met.SetPtEtaPhi(cms2.evt_tcmet(), 0, cms2.evt_tcmetPhi());
  hyp_met += corr;
  if (cms2.hyp_type()[i_hyp] == 0 || cms2.hyp_type()[i_hyp] == 3) {
    if (hyp_met.Pt() < 30) return false;
    //    if ( TMath::Abs(hyp_p4[i_hyp]->mass()-90.0)<10.0) return false;
    if( hyp_met.Pt()/cms2.hyp_p4()[i_hyp].pt()<0.6 && 
	acos(cos(hyp_met.Phi()-cms2.hyp_p4()[i_hyp].phi()-3.1416))<0.25 ) return false;
  }
  // for e-mu and mu-e
  if (cms2.hyp_type()[i_hyp] == 1 || cms2.hyp_type()[i_hyp] == 2) {
    if (hyp_met.Pt() < 20) return false;
  }
  return true;
}

double MetSpecial(double Met, double MetPhi, int i_hyp)
{
  //Warning, this was designed to work in a dilepton environment - NOT a trilepton  
  double DeltaPhi = nearestDeltaPhi(MetPhi,i_hyp);

  if (DeltaPhi < TMath::Pi()/2) return Met*TMath::Sin(DeltaPhi);
  else return Met;

  return -1.0;
}//END MetSpecial calculation  

double MetSpecialTrilep(double Met, double MetPhi, int i_hyp)
{
  //Warning, this was designed to work in a trilepton environment - NOT a dilepton  
  double DeltaPhi = nearestDeltaPhiTrilep(MetPhi,i_hyp);

  if (DeltaPhi < TMath::Pi()/2) return Met*TMath::Sin(DeltaPhi);
  else return Met;

  return -1.0;
}//END MetSpecial calculation  

//--------------------------------------------
// Pass 4 MET selection
// Use MetSpecial from CDF for now
//--------------------------------------------
bool pass4Met(int i_hyp, const TVector3& corr) {
  TVector3 hyp_met;
  hyp_met.SetPtEtaPhi(cms2.evt_tcmet(), 0, cms2.evt_tcmetPhi());
  hyp_met += corr;
  double metspec = MetSpecial(hyp_met.Pt(), hyp_met.Phi(), i_hyp);
  if (cms2.hyp_type()[i_hyp] == 0 || cms2.hyp_type()[i_hyp] == 3) {
    if ( metspec < 20 ) return false;
    //if ( metspec < 20 && hyp_p4->mass() < 90 ) return false;
    if ( hyp_met.Pt() < 45 ) return false;
  }
  else if (cms2.hyp_type()[i_hyp] == 1 || cms2.hyp_type()[i_hyp] == 2) {
    //if ( metspec < 20 && hyp_p4->mass() < 90 ) return false;
    if ( metspec < 20 ) return false;
  }
  return true;
}

//-------------------------------------------
// plain SUMET cut for reducing SM contribution 
// in WW selection to 10%
//-------------------------------------------
bool sumEt10(double sumEt) {
  //  std::cout<<"cutting on sumEt10: "<<sumEt<<std::endl;
  if ( sumEt < 225 ) return false; 
  return true;
}

//-------------------------------------------
// plain SUMET cut for reducing SM contribution 
// in WW selection to 1%
//-------------------------------------------
bool sumEt1(double sumEt) {
  if ( sumEt < 525 ) return false; // good for ~1%? currently "optimizing"
  return true;
}

//-------------------------------------------
// plain MET cut for reducing SM contribution 
// in WW selection to 10%
//-------------------------------------------
bool met10(int i_hyp, const TVector3& corr) {
  TVector3 hyp_met;
  hyp_met.SetPtEtaPhi(cms2.evt_metMuonJESCorr(), 0, cms2.evt_metMuonJESCorrPhi());
  hyp_met += corr;
  if ( hyp_met.Pt() < 110 ) return false;
  return true;
} 
  
//-------------------------------------------
// plain MET cut for reducing SM contribution 
// in WW selection to 1%
//-------------------------------------------
bool met1(int i_hyp, const TVector3& corr) {
  TVector3 hyp_met;
  hyp_met.SetPtEtaPhi(cms2.evt_metMuonJESCorr(), 0, cms2.evt_metMuonJESCorrPhi());
  hyp_met += corr;
  if ( hyp_met.Pt() < 175 ) return false;
  return true;
}

// hyp-dependent met variables: make sure that the MET is corrected for the muon used in the hypothesis
// in case it's not: return corrected value
double met_pat_metCor_hyp(unsigned int hypIdx){
  if (cms2.hyp_type()[hypIdx] ==3) return cms2.met_pat_metCor();
  double lmetx = cms2.met_pat_metCor()*cos(cms2.met_pat_metPhiCor());
  double lmety = cms2.met_pat_metCor()*sin(cms2.met_pat_metPhiCor());

  unsigned int i_lt = cms2.hyp_lt_index()[hypIdx];
  unsigned int i_ll = cms2.hyp_ll_index()[hypIdx];
  if (abs(cms2.hyp_lt_id()[hypIdx])==13 && cms2.mus_met_flag()[i_lt] == 0){
    lmetx+= - cms2.mus_met_deltax()[i_lt] - cms2.mus_p4()[i_lt].x();
    lmety+= - cms2.mus_met_deltay()[i_lt] - cms2.mus_p4()[i_lt].y();
  }
  if (abs(cms2.hyp_ll_id()[hypIdx])==13 && cms2.mus_met_flag()[i_ll] == 0){
    lmetx+= - cms2.mus_met_deltax()[i_ll] - cms2.mus_p4()[i_ll].x();
    lmety+= - cms2.mus_met_deltay()[i_ll] - cms2.mus_p4()[i_ll].y();
  }
  return sqrt(lmetx*lmetx+lmety*lmety);
}

double met_pat_metPhiCor_hyp(unsigned int hypIdx){
  if (cms2.hyp_type()[hypIdx] ==3) return cms2.met_pat_metPhiCor();
  double lmetx = cms2.met_pat_metCor()*cos(cms2.met_pat_metPhiCor());
  double lmety = cms2.met_pat_metCor()*sin(cms2.met_pat_metPhiCor());

  unsigned int i_lt = cms2.hyp_lt_index()[hypIdx];
  unsigned int i_ll = cms2.hyp_ll_index()[hypIdx];
  if (abs(cms2.hyp_lt_id()[hypIdx])==13 && cms2.mus_met_flag()[i_lt] == 0){
    lmetx+= - cms2.mus_met_deltax()[i_lt] - cms2.mus_p4()[i_lt].x();
    lmety+= - cms2.mus_met_deltay()[i_lt] - cms2.mus_p4()[i_lt].y();
  }
  if (abs(cms2.hyp_ll_id()[hypIdx])==13 && cms2.mus_met_flag()[i_ll] == 0){
    lmetx+= - cms2.mus_met_deltax()[i_ll] - cms2.mus_p4()[i_ll].x();
    lmety+= - cms2.mus_met_deltay()[i_ll] - cms2.mus_p4()[i_lt].y();
  }
  return atan2(lmety,lmetx);
}

// do the same with tcmet now
double evt_tcmet_hyp(unsigned int hypIdx){
  if (cms2.hyp_type()[hypIdx] ==3) return cms2.evt_tcmet();
  double lmetx = cms2.evt_tcmet()*cos(cms2.evt_tcmetPhi());
  double lmety = cms2.evt_tcmet()*sin(cms2.evt_tcmetPhi());

  unsigned int i_lt = cms2.hyp_lt_index()[hypIdx];
  unsigned int i_ll = cms2.hyp_ll_index()[hypIdx];
  if (abs(cms2.hyp_lt_id()[hypIdx])==13){
    if(cms2.mus_tcmet_flag()[i_lt] == 0){
      lmetx+= - cms2.mus_met_deltax()[i_lt] - cms2.mus_p4()[i_lt].x();
      lmety+= - cms2.mus_met_deltay()[i_lt] - cms2.mus_p4()[i_lt].y();
    } else if (cms2.mus_tcmet_flag()[i_lt] == 4){
      lmetx+= - cms2.mus_tcmet_deltax()[i_lt] - cms2.mus_met_deltax()[i_lt] - cms2.mus_p4()[i_lt].x(); 
      lmety+= - cms2.mus_tcmet_deltay()[i_lt] - cms2.mus_met_deltay()[i_lt] - cms2.mus_p4()[i_lt].y(); 
    }
  }
  if (abs(cms2.hyp_ll_id()[hypIdx])==13){
    if(cms2.mus_tcmet_flag()[i_ll] == 0){ 
      lmetx+= - cms2.mus_met_deltax()[i_ll] - cms2.mus_p4()[i_ll].x(); 
      lmety+= - cms2.mus_met_deltay()[i_ll] - cms2.mus_p4()[i_ll].y(); 
    } else if (cms2.mus_tcmet_flag()[i_ll] == 4){ 
      lmetx+= - cms2.mus_tcmet_deltax()[i_ll] - cms2.mus_met_deltax()[i_ll] - cms2.mus_p4()[i_ll].x();  
      lmety+= - cms2.mus_tcmet_deltay()[i_ll] - cms2.mus_met_deltay()[i_ll] - cms2.mus_p4()[i_ll].y();  
    } 
  }
  return sqrt(lmetx*lmetx+lmety*lmety);
}

double evt_tcmetPhi_hyp(unsigned int hypIdx){
  if (cms2.hyp_type()[hypIdx] ==3) return cms2.evt_tcmetPhi();
  double lmetx = cms2.evt_tcmet()*cos(cms2.evt_tcmetPhi());
  double lmety = cms2.evt_tcmet()*sin(cms2.evt_tcmetPhi());

  unsigned int i_lt = cms2.hyp_lt_index()[hypIdx];
  unsigned int i_ll = cms2.hyp_ll_index()[hypIdx];
  if (abs(cms2.hyp_lt_id()[hypIdx])==13){
    if(cms2.mus_tcmet_flag()[i_lt] == 0){
      lmetx+= - cms2.mus_met_deltax()[i_lt] - cms2.mus_p4()[i_lt].x();
      lmety+= - cms2.mus_met_deltay()[i_lt] - cms2.mus_p4()[i_lt].y();
    } else if (cms2.mus_tcmet_flag()[i_lt] == 4){
      lmetx+= - cms2.mus_tcmet_deltax()[i_lt] - cms2.mus_met_deltax()[i_lt] - cms2.mus_p4()[i_lt].x(); 
      lmety+= - cms2.mus_tcmet_deltay()[i_lt] - cms2.mus_met_deltay()[i_lt] - cms2.mus_p4()[i_lt].y(); 
    }
  }
  if (abs(cms2.hyp_ll_id()[hypIdx])==13){
    if(cms2.mus_tcmet_flag()[i_ll] == 0){ 
      lmetx+= - cms2.mus_met_deltax()[i_ll] - cms2.mus_p4()[i_ll].x(); 
      lmety+= - cms2.mus_met_deltay()[i_ll] - cms2.mus_p4()[i_ll].y(); 
    } else if (cms2.mus_tcmet_flag()[i_ll] == 4){ 
      lmetx+= - cms2.mus_tcmet_deltax()[i_ll] - cms2.mus_met_deltax()[i_ll] - cms2.mus_p4()[i_ll].x();  
      lmety+= - cms2.mus_tcmet_deltay()[i_ll] - cms2.mus_met_deltay()[i_ll] - cms2.mus_p4()[i_ll].y();  
    } 
  }
  return atan2(lmety,lmetx);
}

// met cut for ttbar dilepton analysis...
// includes a boolean to switch to tcmet
bool passMet_OF20_SF30(int hypIdx, bool useTcMet) {
  float mymet;
  if (useTcMet) {
    mymet = evt_tcmet_hyp(hypIdx);
  } else {
    mymet = met_pat_metCor_hyp(hypIdx);
  }
  if  (cms2.hyp_type().at(hypIdx) == 0 || cms2.hyp_type().at(hypIdx) == 3) {
    if (mymet < 30) return false;
  }
  
  if (cms2.hyp_type().at(hypIdx) == 1 || cms2.hyp_type().at(hypIdx) == 2) {
    if (mymet < 20) return false;
  }
  return true;
}  

//  ***** The following two functions should be deprecated *********************
//  ***** and substituted by the preceeding one            *********************
// event-level pat-met: emu met >20, mm,em met>30
bool passPatMet_OF20_SF30(float metx, float mety, int hypIdx){
  float mymet = sqrt(metx*metx + mety*mety);
  if  (cms2.hyp_type().at(hypIdx) == 0 || cms2.hyp_type().at(hypIdx) == 3) {
    if (mymet < 30) return false;
  }
  
  if (cms2.hyp_type().at(hypIdx) == 1 || cms2.hyp_type().at(hypIdx) == 2) {
    if (mymet < 20) return false;
  }
  return true;
}

// event-level pat-met: emu met >20, mm,em met>30
bool passPatMet_OF20_SF30(int hypIdx){
  return passPatMet_OF20_SF30(met_pat_metCor_hyp(hypIdx)*cos(met_pat_metPhiCor_hyp(hypIdx)), 
			      met_pat_metCor_hyp(hypIdx)*sin(met_pat_metPhiCor_hyp(hypIdx)),
			      hypIdx);
}

//**************************************************************************
// met cut for ttbar dilepton analysis...
// includes a boolean to switch to tcmet
bool passMet_OF30_SF50(int hypIdx, bool useTcMet) {
  float mymet;
  if (useTcMet) {
    mymet = evt_tcmet_hyp(hypIdx);
  } else {
    mymet = met_pat_metCor_hyp(hypIdx);
  }
  if  (cms2.hyp_type().at(hypIdx) == 0 || cms2.hyp_type().at(hypIdx) == 3) {
    if (mymet < 50) return false;
  }
  
  if (cms2.hyp_type().at(hypIdx) == 1 || cms2.hyp_type().at(hypIdx) == 2) {
    if (mymet < 30) return false;
  }
  return true;
}  

//  ***** The following two functions should be deprecated *********************
//  ***** and substituted by the preceeding one            *********************
// event-level pat-met: emu met >20, mm,em met>30
bool passPatMet_OF30_SF50(float metx, float mety, int hypIdx){
  float mymet = sqrt(metx*metx + mety*mety);
  if  (cms2.hyp_type().at(hypIdx) == 0 || cms2.hyp_type().at(hypIdx) == 3) {
    if (mymet < 50) return false;
  }
  
  if (cms2.hyp_type().at(hypIdx) == 1 || cms2.hyp_type().at(hypIdx) == 2) {
    if (mymet < 30) return false;
  }
  return true;
}
// event-level pat-met: emu met >30, mm,em met>50
bool passPatMet_OF30_SF50(int hypIdx){
  return passPatMet_OF30_SF50(met_pat_metCor_hyp(hypIdx)*cos(met_pat_metPhiCor_hyp(hypIdx)), 
			      met_pat_metCor_hyp(hypIdx)*sin(met_pat_metPhiCor_hyp(hypIdx)),
			      hypIdx);
}

bool passMetVJets09(float value, bool useTcMet) {
  float mymet;
  if (useTcMet) {
    mymet = cms2.evt_tcmet();
  } else {
    mymet = cms2.met_pat_metCor();
  }
  if (mymet <= value) return false;
  return true;
}
