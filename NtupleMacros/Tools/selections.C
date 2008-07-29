//===========================================================
//
// Various selection functions are kept here
//
//============================================================
#include "Math/LorentzVector.h"
#include "TMath.h"
// #include "selections.h"
//#include <vector>
//#include "CMS1.h"
#include "TDatabasePDG.h"

//----------------------------------------------------------------
// A ridicolusly simple function, but since the Z veto is used 
// in two places, might as well centralize it to keep consistency
//----------------------------------------------------------------
bool inZmassWindow (float mass) {
  if (mass > 76. && mass < 106.) return true;
  return false;
}
//----------------------------------------------------------------
// true electron
//----------------------------------------------------------------
// bool trueElectron(int index) {
//   if ( abs(els_mc_id[index]) != 11 ) return false;
//   if ( abs(els_mc_motherid[index]) == 23 || abs(els_mc_motherid[index]) == 24) return true;
//   return false;
// }
//----------------------------------------------------------------
// true muon
//----------------------------------------------------------------
// bool trueMuon(int index) {
//    if ( abs(mus_mc_id[index]) != 13 ) return false;
//    if ( abs(mus_mc_motherid[index]) == 23 || abs(mus_mc_motherid[index]) == 24) return true;
//   return false;
// }
//----------------------------------------------------------------
// Electron ID without isolation
//----------------------------------------------------------------
bool goodElectronWithoutIsolation(int index) {
  if ( els_tightId[index]     !=  1) return false;
  if ( els_closestMuon[index] != -1) return false;
  if ( abs(els_d0[index]) > 0.025)   return false;
  return true;
}
//----------------------------------------------------------------
// Muon ID without isolation
//---------------------------------------------------------------
bool goodMuonWithoutIsolation(int index) {
  if (mus_gfit_chi2[index]/mus_gfit_ndof[index] > 5.) return false;
  if (abs(mus_d0[index])   > 0.25) return false;
  if (mus_validHits[index] < 7)    return false;
  return true;
}
//-----------------------------------------------------------
// Electron Isolation
//-----------------------------------------------------------
bool passElectronIsolation(int index) {
  double sum = els_tkIso[index];
  double pt  = els_p4[index].pt();
   if ( pt/(pt+sum) < 0.92) return false;
  return true;  
} 
//-----------------------------------------------------------
// Muon Isolation
//-----------------------------------------------------------
bool passMuonIsolation(int index) {
  double sum =  mus_iso03_sumPt[index] +  
                mus_iso03_emEt[index]  +
                mus_iso03_hadEt[index];
  double pt  = mus_p4[index].pt(); 
  if ( pt/(pt+sum) < 0.92) return false;
  return true;  
}
//--------------------------------------------
// Muon ID with isolation
//--------------------------------------------
bool goodMuonIsolated(int index) {
  if (!goodMuonWithoutIsolation(index)) return false;
  if (!passMuonIsolation(index))       return false;
  return true;
}
//--------------------------------------------
// Electron ID with isolation
//--------------------------------------------
bool goodElectronIsolated(int index) {
  if (!goodElectronWithoutIsolation(index)) return false;
  if (!passElectronIsolation(index))       return false;
  return true;
}
