//===========================================================
//
// Various selection functions are kept here
//
//============================================================
#include <assert.h>
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TDatabasePDG.h"
#include "selections.h"

// CMS2 includes
#include "CMS2.h"
#include "utilities.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

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
bool trueElectron(int index) {
  if ( TMath::Abs(cms2.els_mc_id()[index]) != 11 ) return false;
//   if ( TMath::Abs(cms2.els_mc_motherid()[index]) == 23 || TMath::Abs(cms2.els_mc_motherid()[index]) == 24) return true;
  return true;
}
// //----------------------------------------------------------------
// // true muon
// //----------------------------------------------------------------
bool trueMuon(int index) {
  if ( TMath::Abs(cms2.mus_mc_id()[index]) != 13 ) return false;
//    if ( TMath::Abs(cms2.mus_mc_motherid()[index]) == 23 || TMath::Abs(cms2.mus_mc_motherid()[index]) == 24) return true;
  return true;
}
//----------------------------------------------------------------
// Electron ID without isolation
//----------------------------------------------------------------
bool goodElectronWithoutIsolation(int index) {
  if ( cms2.els_tightId22XMinMatteo().at(index)     !=  1) return false;
  if ( cms2.els_closestMuon().at(index) != -1) return false;
  if ( TMath::Abs(cms2.els_d0corr().at(index)) > 0.025)   return false;
  return true;
}

//----------------------------------------------------------------
// MC truth tag - Muon from W
//----------------------------------------------------------------
bool trueMuonFromW_WJets(int index) {

  // added requirement for muon to have pt>20
    double pt_cut = 20.;
  //  double pt_cut = 0.;
  
  // try to find out if there is a muon
  // from W. Currently valid ONLY on WJets!! 
  if ( TMath::Abs(cms2.mus_mc_id()[index]) == 13 && ( TMath::Abs(cms2.mus_mc_motherid()[index]) == 24) && cms2.mus_p4().at(index).pt() > pt_cut ) {
    //    std::cout<<"best match - part id: "<<(cms2.mus_mc_id()[index])<<" mother: "<<(cms2.mus_mc_motherid()[index])<<std::endl;
    return true;
  }
  // need additional loop over status 3 particles
  // if it contains a W and an muon from W, we claim
  // that the muon is from W - valid for WJets
  unsigned int nGen = cms2.genps_id().size();
  for( unsigned int iGen = 0; iGen < nGen; ++iGen){
    // just in case the mother link does work :)
    if ( cms2.genps_status()[iGen] == 3 && TMath::Abs(cms2.genps_id()[iGen]) == 13 && TMath::Abs(cms2.genps_id_mother()[iGen]) == 24 && cms2.genps_p4()[iGen].pt() > pt_cut ) {
      //    std::cout<<"2nd best match part id: "<<cms2.genps_id()[iGen]<<" mother: "<<cms2.genps_id_mother()[iGen]<<std::endl;
      return true;
    }
    // also return true if there is a status 3 muon
    if ( cms2.genps_status()[iGen] == 3 && TMath::Abs(cms2.genps_id()[iGen]) == 13  && cms2.genps_p4()[iGen].pt() > pt_cut ) {
      //    std::cout<<"3rd best match part id: "<<cms2.genps_id()[iGen]<<" mother: NOT CHECKED"<<std::endl;
      return true;
    }
    //      if ( TMath::Abs(cms2.genps_id_mother()[iGen]) == 23 ) return true;
    //      if ( TMath::Abs(cms2.genps_id_mother()[iGen]) == 24 ) return true;
  }
  
  return false;
}


//----------------------------------------------------------------
// MC truth tag - Electron from W
//----------------------------------------------------------------
bool trueElectronFromW_WJets(int index) {
  // try to find out if there is a electron
  // from W. Currently valid ONLY on WJets!! 
  if ( TMath::Abs(cms2.els_mc_id()[index]) == 11 && ( TMath::Abs(cms2.els_mc_motherid()[index]) == 24) ) {
    //    std::cout<<"best match - part id: "<<(cms2.els_mc_id()[index])<<" mother: "<<(cms2.els_mc_motherid()[index])<<std::endl;
    return true;
  }
  // need additional loop over status 3 particles
  // if it contains a W and an electron from W, we claim
  // that the electron is from W - valid for WJets
  unsigned int nGen = cms2.genps_id().size();
  for( unsigned int iGen = 0; iGen < nGen; ++iGen){
    // just in case the mother link does work :)
    if ( cms2.genps_status()[iGen] == 3 && TMath::Abs(cms2.genps_id()[iGen]) == 11 && TMath::Abs(cms2.genps_id_mother()[iGen]) == 24 ) {
      //      std::cout<<"2nd best match part id: "<<cms2.genps_id()[iGen]<<" mother: "<<cms2.genps_id_mother()[iGen]<<std::endl;
      return true;
    }
    // also return true if there is a status 3 electron
    if ( cms2.genps_status()[iGen] == 3 && TMath::Abs(cms2.genps_id()[iGen]) == 11 ) {
      //      std::cout<<"3rd best match part id: "<<cms2.genps_id()[iGen]<<" mother: NOT CHECKED"<<std::endl;
      return true;
    }
  }
  
  return false;
}


//----------------------------------------------------------------
// Electron ID without isolation or d0 cut
//----------------------------------------------------------------
bool goodElectronWithoutIsolationWithoutd0(int index) {
  if ( cms2.els_tightId22XMinMatteo().at(index)     !=  1) return false;
  if ( cms2.els_closestMuon().at(index) != -1) return false;
  return true;
}
//----------------------------------------------------------------
// loose Electron ID without isolation
//----------------------------------------------------------------
bool goodLooseElectronWithoutIsolation(int index) {
  if ( cms2.els_looseId().at(index)     !=  1) return false;
  if ( cms2.els_closestMuon().at(index) != -1) return false;
  if ( TMath::Abs(cms2.els_d0corr().at(index)) > 0.025)   return false;
  return true;
}
//----------------------------------------------------------------
// Muon ID without isolation
//---------------------------------------------------------------
bool goodMuonWithoutIsolation(int index) {
  if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) > 10.) return false;
  if (TMath::Abs(cms2.mus_d0corr().at(index))   > 0.2) return false;
  if (cms2.mus_validHits().at(index) < 11)    return false;
  if ( (cms2.mus_type().at(index)&0x2)==0 ) return false;
  return true;
}
//-----------------------------------------------------------
// Electron Isolation using pat iso
//-----------------------------------------------------------
double el_rel_iso (int index, bool use_calo_iso)
{
     double sum = cms2.els_pat_trackIso().at(index);
     if (use_calo_iso)
        sum += cms2.els_pat_ecalIso()[index] + cms2.els_pat_hcalIso()[index];
  
     double pt  = cms2.els_p4().at(index).pt();
     return pt / (pt + sum + 1e-5);
}
bool passElectronIsolation(int index, bool use_calo_iso) 
{
     const double cut = 0.92;
     return el_rel_iso(index, use_calo_iso) > cut;
}

bool passElectronIsolationLoose(int index, bool use_calo_iso) 
{
  const double cut = 0.8; // leads to 91 pred, 81 obs
     return el_rel_iso(index, use_calo_iso) > cut;
} 

bool passElectronIsolationLoose2(int index, bool use_calo_iso) 
{
  //     const double cut = 0.85; leads to 107 pred, 81 obs
     const double cut = 0.75;
     return el_rel_iso(index, use_calo_iso) > cut;
} 

//=================================================
// lepton isolation: use id of the lepton do decode
//-------------------------------------------------
bool passLeptonIsolation(int id, int index, bool use_ele_calo_iso){
  if (abs(id)==11) return passElectronIsolation(index, use_ele_calo_iso);
  if (abs(id)==13) return passMuonIsolation(index);
  return false;
}

//-----------------------------------------------------------
// Electron Isolation using ECAL clusters
//-----------------------------------------------------------
double el_rel_iso_1_6 (int index, bool use_calo_iso)
{
     double sum = cms2.els_tkIso().at(index);
     if (use_calo_iso)
	  sum += cms2.els_ecalJuraIso()[index] + cms2.els_hcalConeIso()[index];
     double pt  = cms2.els_p4().at(index).pt();
     return pt / (pt + sum+1e-5);
}

bool passElectronIsolation_1_6 (int index, bool use_calo_iso) 
{
     const double cut = 0.92;
     return el_rel_iso_1_6(index, use_calo_iso) > cut;
}

bool passElectronIsolationLoose_1_6 (int index, bool use_calo_iso) 
{
     const double cut = 0.8; // leads to 91 pred, 81 obs
     return el_rel_iso_1_6(index, use_calo_iso) > cut;
} 

bool passElectronIsolationLoose2_1_6 (int index, bool use_calo_iso) 
{
  //     const double cut = 0.85; leads to 107 pred, 81 obs
     const double cut = 0.75;
     return el_rel_iso_1_6(index, use_calo_iso) > cut;
} 
//-----------------------------------------------------------
// Muon Isolation
//-----------------------------------------------------------
double mu_rel_iso (int index)
{
     double sum =  cms2.mus_iso03_sumPt().at(index) +  
	  cms2.mus_iso03_emEt().at(index)  +
	  cms2.mus_iso03_hadEt().at(index);
     double pt  = cms2.mus_p4().at(index).pt();
     return pt / (pt+sum+1e-5);

}
bool passMuonIsolation(int index) 
{
     return mu_rel_iso(index) > 0.92;
}

bool passMuonIsolationLoose(int index) 
{
  //     const double cut = 0.85; leads to 107 pred, 81 obs
     const double cut = 0.75;
     return mu_rel_iso(index) > cut;
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
bool goodElectronIsolated(int index, bool use_calo_iso) {
  // if (!goodElectronWithoutIsolationWithoutd0(index)) return false;
  if (!goodElectronWithoutIsolation(index)) return false;
  if (!passElectronIsolation(index, use_calo_iso))       return false;
  return true;
}
//--------------------------------------------
// "super-tight" Electron ID (to kill EM fakes)
//--------------------------------------------
bool supertightElectron (int index)
{
     if (TMath::Abs(cms2.els_p4()[index].eta()) > 1.479)
	  return false;
     if (cms2.els_sigmaPhiPhi()[index] > 0.018)
  	  return false;
     if (cms2.els_sigmaEtaEta()[index] > 0.009)
	  return false;
     if (cms2.els_eOverPIn()[index] > 2 || cms2.els_eOverPIn()[index] < 0.75)
	  return false;
     if (cms2.els_charge()[index] * cms2.els_dPhiIn()[index] > 0.04)
	  return false;
     if (cms2.els_charge()[index] * cms2.els_dPhiOut()[index] < -1e-3)
	  return false;
     if (cms2.els_dEtaIn()[index] > 0.0025 || cms2.els_dEtaIn()[index] < -0.0025)
	  return false;
     return true;
}
//--------------------------------------------
// tighter cut on delta phi in (to kill more fakes)
//--------------------------------------------
bool deltaPhiInElectron (int index)
{
     return cms2.els_charge()[index] * cms2.els_dPhiIn()[index] < 0.04;
}

//--------------------------------------------
// Fudge to reject events where the LT muon
// has a badly measured momentum
//--------------------------------------------
bool muonReconstructionCleaning(int i_hyp, float threshold)
{
   // only apply to mm hyp type
   if (cms2.hyp_type()[i_hyp] == 0)
   {
      if (fabs((cms2.hyp_lt_trk_p4()[i_hyp].Pt()
             /cms2.hyp_lt_p4()[i_hyp].Pt()) - 1) > threshold) return false;
   }
   return true;
}

//
// Refactorized MET cuts
//
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
//
// end refactorized MET cuts
//

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

double nearestDeltaPhi(double Phi, int i_hyp)
{
  //WARNING!  This was designed to work in a dilepton environment - NOT a trilepton 
  double tightDPhi = TMath::Min(TMath::Abs(cms2.hyp_lt_p4()[i_hyp].Phi() - Phi), 2*TMath::Pi() - TMath::Abs(cms2.hyp_lt_p4()[i_hyp].Phi() - Phi));
  double looseDPhi = TMath::Min(TMath::Abs(cms2.hyp_ll_p4()[i_hyp].Phi() - Phi), 2*TMath::Pi() - TMath::Abs(cms2.hyp_ll_p4()[i_hyp].Phi() - Phi));

  return TMath::Min(tightDPhi, looseDPhi);

}//END nearest DeltaPhi                                                                                                                                 

double MetSpecial(double Met, double MetPhi, int i_hyp)
{
  //Warning, this was designed to work in a dilepton environment - NOT a trilepton  
  double DeltaPhi = nearestDeltaPhi(MetPhi,i_hyp);

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

//-------------------------------------------------
// Auxiliary function to scan the doc line and 
// identify DY-> ee vs mm vs tt
//-------------------------------------------------
int getDrellYanType() {
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

int getVVType() {
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


//--------------------------------------------
// Booleans for DY
//------------------------------------------
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

//--------------------------------------------
// Booleans for VVjets
//------------------------------------------

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

//--------------------------------------------------------------------
// Veto events if there are two leptons in the 
// event that make the Z mass.  This uses the mus and els
// blocks, ie, it is a veto that can use the 3rd (4th,5th,..)
// lepton in the event.
//
// Both leptons must be 20 GeV, and pass the same cuts as 
// the hypothesis leptons, except that one of them can be non-isolated
//---------------------------------------------------------------------
bool additionalZveto() {

  // true if we want to veto this event
  bool veto=false;

  // first, look for Z->mumu
  for (unsigned int i=0; i < cms2.mus_p4().size(); i++) {
    if (cms2.mus_p4().at(i).pt() < 20.)     continue;
    if (!goodMuonWithoutIsolation(i)) continue;

    for (unsigned int j=i+1; j < cms2.mus_p4().size(); j++) {
      if (cms2.mus_p4().at(j).pt() < 20.) continue;
      if (!goodMuonWithoutIsolation(j)) continue;
      if (cms2.mus_charge().at(i) == cms2.mus_charge().at(j)) continue;

      // At least one of them has to pass isolation
      if (!passMuonIsolation(i) && !passMuonIsolation(j)) continue;

      // Make the invariant mass
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > 
                vec = cms2.mus_p4().at(i) + cms2.mus_p4().at(j);
      if ( inZmassWindow(vec.mass()) ) return true;

    }
  }

  // now, look for Z->ee
  for (unsigned int i=0; i < cms2.els_p4().size(); i++) {
    if (cms2.els_p4().at(i).pt() < 20.)     continue;
    if (!goodElectronWithoutIsolation(i)) continue;

    for (unsigned int j=i+1; j<cms2.els_p4().size(); j++) {
      if (cms2.els_p4().at(j).pt() < 20.) continue;
      if (!goodElectronWithoutIsolation(j)) continue;
      if (cms2.els_charge().at(i) == cms2.els_charge().at(j)) continue;

      // At least one of them has to pass isolation
      if (!passElectronIsolation(i, false) && !passElectronIsolation(j, false)) continue;

      // Make the invariant mass
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > 
                vec = cms2.els_p4().at(i) + cms2.els_p4().at(j);
      if ( inZmassWindow(vec.mass()) ) return true;

    }
  }
  // done
  return veto;
}

//------------------------------------------------------------
// Not a selection function per se, but useful nonetheless:
// dumps the documentation lines for this event
//------------------------------------------------------------
void dumpDocLines() {
  int size = cms2.genps_id().size();
  // Initialize particle database
  TDatabasePDG *pdg = new TDatabasePDG();
  std::cout << "------------------------------------------" << std::endl;
  for (int j=0; j<size; j++) {
    cout << setw(9) << left << pdg->GetParticle(cms2.genps_id().at(j))->GetName() << " "
         << setw(7) << right << setprecision(4) << cms2.genps_p4().at(j).pt() << "  "
         << setw(7) << right << setprecision(4) << cms2.genps_p4().at(j).phi() << "  "
         << setw(10) << right << setprecision(4) << cms2.genps_p4().at(j).eta() << "  "
         << setw(10) << right << setprecision(4) << cms2.genps_p4().at(j).mass() << endl;
  }
  // Clean up
  delete pdg;
}

//----------------------------------------------------------------------
// search and destroy extra muons.  If they're stiff and isolated,
// they're probably from WZ, ZZ or (in emu) DYmm (with a spurious
// electron).  If they're soft or non-isolated, use them to tag top
// events
// ----------------------------------------------------------------------
bool passTriLepVeto (int i_dilep)
{
     double tag_mu_pt = tagMuonPt(i_dilep);
     if (tag_mu_pt < 0) // no mu
	  return true;
     if (tag_mu_pt < 20) // soft
	  return true;
     double tag_mu_iso = tagMuonRelIso(i_dilep);
     if (tag_mu_iso < 0.9) // non-isolated
	  return true;
     // we've found a muon that's stiff and isolated.  Fail the veto.
     return false;
}

//int numberOfExtraMuons(int i_hyp, bool nonisolated = false){
int numberOfExtraMuons(int i_hyp, bool nonisolated){
  unsigned int nMuons = 0;
  for (int imu=0; imu < int(cms2.mus_charge().size()); ++imu) {
    // quality cuts
    if (  ((cms2.mus_goodmask()[imu]) & (1<<14)) == 0 ) continue; // TMLastStationOptimizedLowPtTight
    if ( cms2.mus_p4()[imu].pt() < 3 ) continue;
    if ( TMath::Abs(cms2.mus_d0corr()[imu]) > 0.2) continue;
    if ( cms2.mus_validHits()[imu] < 11) continue;
    if ( TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && cms2.hyp_lt_index()[i_hyp] == imu ) continue;
    if ( TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && cms2.hyp_ll_index()[i_hyp] == imu ) continue;
    if ( nonisolated && mu_rel_iso(imu) > 0.90 && cms2.mus_p4()[imu].pt()>20 ) continue;
    ++nMuons;
  }
  return nMuons;
}

bool passMuonBVeto_1_6 (int i_dilep, bool soft_nonisolated)
{
     if (soft_nonisolated) {
	  double tag_mu_pt = tagMuonPt(i_dilep);
	  if (tag_mu_pt < 0) // no mu
	       return true;
	  if (tag_mu_pt < 20) // soft
	       return false;
	  return true;
     } else {
	  unsigned int mus_in_hyp = 0;
	  if (TMath::Abs(cms2.hyp_lt_id()[i_dilep]) == 13)
	       mus_in_hyp++;
	  if (TMath::Abs(cms2.hyp_ll_id()[i_dilep]) == 13)
	       mus_in_hyp++;
	  return cms2.mus_p4().size() <= mus_in_hyp;
     }
     assert(false);
     return false;
}

// If there is an extra muon in the event, return its index.  (Otherwise -1) 
int tagMuonIdx (int i_dilep)
{
     for (unsigned int i = 0; i < cms2.mus_p4().size(); ++i) {
	  if (TMath::Abs(cms2.hyp_lt_id()[i_dilep]) == 13 &&
	      cms2.hyp_lt_index()[i_dilep] == int(i))
	       continue;
	  if (TMath::Abs(cms2.hyp_ll_id()[i_dilep]) == 13 &&
	      cms2.hyp_ll_index()[i_dilep] == int(i))
	       continue;
	  return i;
     }
     return -1;
}

// return the pt of the first muon that's not lt or ll
double tagMuonPt (int i_dilep)
{
     int idx = tagMuonIdx(i_dilep);
     if (idx == -1)
	  return -1;
     return cms2.mus_p4()[idx].pt();
}

// same for rel iso
double tagMuonRelIso (int i_dilep)
{
     int idx = tagMuonIdx(i_dilep);
     if (idx == -1)
	  return -1;
     return mu_rel_iso(idx);
}

//--------------------------------
//
// Functions related to trkjet veto
// This needs ot be rewritten once we
// have added appropriate variables into ntuple itself.
//
//--------------------------------
int NjetVeto(std::vector<TLorentzVector>& Jet, double min_et) {
  int njets = 0;
  for(unsigned int i=0; i<Jet.size() ; ++i) {
    if ( Jet[i].Perp() >= min_et) {
      njets++;
    }
  }
  return njets;
}

int nTrkJets(int i_hyp){
  std::vector<TLorentzVector> trkjets;
  double jetet = 0;
  double jeteta = 3.0;
  // TrkJets & CaloJet save it after the lepton subtraction

  for ( unsigned int itrkjet=0; itrkjet<cms2.trkjets_p4().size(); ++itrkjet) {
    if ((TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11 && dRbetweenVectors(cms2.hyp_lt_p4()[i_hyp],cms2.trkjets_p4()[itrkjet]) < 0.4)||
	(TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11 && dRbetweenVectors(cms2.hyp_ll_p4()[i_hyp],cms2.trkjets_p4()[itrkjet]) < 0.4)
	) continue;
    TLorentzVector p(cms2.trkjets_p4()[itrkjet].Px(), cms2.trkjets_p4()[itrkjet].Py(), cms2.trkjets_p4()[itrkjet].Pz(), cms2.trkjets_p4()[itrkjet].E());
    if (p.Perp() < jetet) continue;
    if (TMath::Abs(p.Eta()) > jeteta) continue;
    trkjets.push_back(p);
  }

  return NjetVeto(trkjets, 15);
}

// count JPT jets excluding those that are close to the hypothesis leptons
unsigned int nJPTs (int i_hyp)
{
     return nJPTs(i_hyp, 20);
}

std::vector<LorentzVector> JPTs(int i_hyp, double etThreshold)
{
     std::vector<LorentzVector> ret;
     const double etaMax      = 3.0;
     const double vetoCone    = 0.4;
     
     for ( unsigned int i=0; i < cms2.jpts_p4().size(); ++i) {
	  if ( cms2.jpts_p4()[i].Et() < etThreshold ) continue;
	  if ( TMath::Abs(cms2.jpts_p4()[i].eta()) > etaMax ) continue;
	  if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.jpts_p4()[i])) < vetoCone ||
	       TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.jpts_p4()[i])) < vetoCone ) continue;
	  ret.push_back(cms2.jpts_p4()[i]);
     }
     return ret;
}

unsigned int nJPTs(int i_hyp, double etThreshold)
{
     return JPTs(i_hyp, etThreshold).size();
}

bool passTrkJetVeto(int i_hyp)
{
     return nTrkJets(i_hyp) == 0;
}

double reliso_lt (int i_hyp, bool use_calo_iso)
{
     // muons do it one way:
     if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13) {
	  return mu_rel_iso(cms2.hyp_lt_index()[i_hyp]);
     }
     // electrons do it another way:
     if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11) {
	  return el_rel_iso(cms2.hyp_lt_index()[i_hyp], use_calo_iso);
     }
     // mysterions are not well handled:
     return 0;
}

double reliso_ll (int i_hyp, bool use_calo_iso)
{
     // muons do it one way:
     if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13) {
	  return mu_rel_iso(cms2.hyp_ll_index()[i_hyp]);
     }
     // electrons do it another way:
     if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11) {
	  return el_rel_iso(cms2.hyp_ll_index()[i_hyp], use_calo_iso);
     }
     // mysterions are not well handled:
     return 0;
}

bool trueMuonFromW(int index) {

  bool muIsFromW = false;

  if( TMath::Abs(cms2.mus_mc_id()[index]) == 13 && TMath::Abs(cms2.mus_mc_motherid()[index]) == 24 ) muIsFromW = true;

  return muIsFromW;
}

bool trueElectronFromW(int index) {

  bool elIsFromW = false;

  if( TMath::Abs(cms2.els_mc_id()[index]) == 11 && TMath::Abs(cms2.els_mc_motherid()[index]) == 24 ) elIsFromW = true;

  return elIsFromW;
}

int conversionPartner (int i_el)
{
     LorentzVector el_p4 = cms2.els_p4()[i_el];
     double min_dcottheta = 5e-4;
     int idx = -1;
     const double el_cottheta = el_p4.pz() / el_p4.pt();
     for (unsigned int i = 0; i < cms2.trks_trk_p4().size(); ++i) {
	  // if the track is within some ratio of the electron momentum,
	  // reject it
	  if (cms2.trks_trk_p4()[i].pt() > 0.1 * el_p4.pt())
	       continue;
	  const double tk_cottheta = cms2.trks_trk_p4()[i].pz() / 
	       cms2.trks_trk_p4()[i].pt();
	  if (TMath::Abs(tk_cottheta - el_cottheta) < min_dcottheta) {
	       min_dcottheta = TMath::Abs(tk_cottheta - el_cottheta);
	       idx = i;
	  }
     }
     if (min_dcottheta < 5e-4)
	  return idx;
     return -1;
}

double conversionDeltaPhi (int i_conv, int i_el)
{
     if (i_conv == -1)
	  return -1;
     double dphi = ROOT::Math::VectorUtil::DeltaPhi(cms2.trks_trk_p4()[i_conv],
						    cms2.els_p4()[i_el]);
     return TMath::Abs(dphi);
}

/* missing ntuple variables
bool passCaloTrkjetCombo ()
{
     double dr = 999;
     double max_sumpt = 0;
     for (unsigned int i = 0; i < cms2.jets_p4().size(); ++i) {
	  double min_dr = 999;
	  double sumpt = 0;
	  for (unsigned int j = 0; j < cms2.trkjets_p4().size(); ++j) {
	       double dr = ROOT::Math::VectorUtil::DeltaR(cms2.trkjets_p4()[j],
							  cms2.jets_p4()[i]);
	       if (dr < min_dr) {
		    min_dr = dr;
		    sumpt = cms2.trkjets_p4()[j].pt() + 
			 cms2.jets_p4()[i].pt() * cms2.jets_tq_noCorrF()[i];
	       }
	  }
	  if (sumpt > max_sumpt) {
	       max_sumpt = sumpt;
	       dr = min_dr;
	  }
     }
     if (max_sumpt > 15 && dr < 0.5)
	  return false;
     return true;
}
*/

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

bool passTrackIsolation(int index){
  //
  // track isolation
  //
  const double cut = 0.92;
  double pt  = cms2.trks_trk_p4()[index].pt();
  return (pt / (pt + trkIsolation(index))) > cut;
}

int passTrackZVeto(int hyp_index) {
  //
  // build list of trk-isolated tracks with minimal quality cuts for dilepton hypothesis
  // nHits > 7
  // pT >= 1.0
  //
  // combine tracks to Z candidates and check invariant mass window for 3 modi:
  // 1: combine one track with either lt or ll of dilepton hyp, if a cand. is within the z window, reject event (return 1)
  // 2: combine two tracks, if a cand. is within the z window ,reject event (return 2)
  // 3: 1 && 2
  //
  // when combining 2 objects, require:
  // - delta z0 < 0.5 cm.
  // - delta R > 0.1
  // - opposite sign

  double dZCut = 0.5;
  bool mode1 = false;
  bool mode2 = false;

  // store trk indices of ll and lt
  unsigned int llTrkIdx = 1000000;
  unsigned int ltTrkIdx = 1000000;
  if ( TMath::Abs(cms2.hyp_lt_id()[hyp_index]) == 13 ) 
    ltTrkIdx = cms2.mus_trkidx()[cms2.hyp_lt_index()[hyp_index]];
  else if ( TMath::Abs(cms2.hyp_lt_id()[hyp_index]) == 11 )
    ltTrkIdx = cms2.els_trkidx()[cms2.hyp_lt_index()[hyp_index]];
  if ( TMath::Abs(cms2.hyp_ll_id()[hyp_index]) == 13 ) 
    llTrkIdx = cms2.mus_trkidx()[cms2.hyp_ll_index()[hyp_index]];
  else if ( TMath::Abs(cms2.hyp_ll_id()[hyp_index]) == 11 )
    llTrkIdx = cms2.els_trkidx()[cms2.hyp_ll_index()[hyp_index]];

  // form trk-isolated track collection
  std::vector<int> isoTrks;
  for ( unsigned int trk = 0;
	trk < cms2.trks_trk_p4().size();
	++trk ) {
    // exclude lt or ll
    if ( trk == llTrkIdx ) continue;
    if ( trk == ltTrkIdx ) continue;
    // require at least 8 hits
    if ( cms2.trks_validHits()[trk] <= 7 ) continue;
    // require pt >= 1. GeV
    if ( cms2.trks_trk_p4()[trk].pt() < 1.0 ) continue;
    // require trk-isolation
    if ( ! passTrackIsolation(trk) ) continue;
    isoTrks.push_back(trk);
  }

  // loop over all selected trk-isolated tracks
  for ( unsigned int trk1 = 0;
	trk1 < isoTrks.size();
	++trk1 ) {
    // try to combine with ll or lt, if in Z mass window, veto event
    if ( !mode1 ) {
      // require opposite sign
      if ( (cms2.hyp_lt_charge()[hyp_index] * cms2.trks_charge()[trk1]) < 0 ) {
	// require delta z0 < 0.5 cm
	double dZ = TMath::Abs(cms2.hyp_lt_z0()[hyp_index] - cms2.trks_z0()[trk1]);
	if ( dZ < dZCut ) {
	  // require delta R > 0.1
	  double dR = ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[hyp_index], cms2.trks_trk_p4()[trk1]);
	  if ( dR > 0.1 ) {
	    // if inv. mass of candidate within z mass window, veto event
	    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > 
	      vec = cms2.hyp_lt_p4()[hyp_index] + cms2.trks_trk_p4()[trk1];
	    if ( inZmassWindow(vec.mass()) ) {
	      mode1 = true;
	    }
	  }
	}
      }
      // require opposite sign
      if ( (cms2.hyp_ll_charge()[hyp_index] * cms2.trks_charge()[trk1]) < 0 ) {
	// require delta z0 < 0.5 cm
	double dZ = TMath::Abs(cms2.hyp_ll_z0()[hyp_index] - cms2.trks_z0()[trk1]);
	if ( dZ < dZCut ) {
	  // require delta R > 0.1
	  double dR = ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[hyp_index], cms2.trks_trk_p4()[trk1]);
	  if ( dR > 0.1 ) {
	    // if inv. mass of candidate within z mass window, veto event
	    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > 
	      vec = cms2.hyp_ll_p4()[hyp_index] + cms2.trks_trk_p4()[trk1];
	    if ( inZmassWindow(vec.mass()) ) {
	      mode1 = true;
	    }
	  }
	}
      }
      // try to combine one track with another track from the collection of selected trk-isolated tracks
    }
    for ( unsigned int trk2 = trk1+1;
	  trk2 < isoTrks.size();
	  ++trk2 ) {
      if ( !mode2 ) {
	// require opposite sign
	if ( (cms2.trks_charge()[trk2] * cms2.trks_charge()[trk1]) < 0 ) {
	  // require delta z0 < 0.5 cm
	  double dZ = TMath::Abs(cms2.trks_z0()[trk2] - cms2.trks_z0()[trk1]);
	  if ( dZ < dZCut ) {
	    // require delta R > 0.1
	    double dR = ROOT::Math::VectorUtil::DeltaR(cms2.trks_trk_p4()[trk2], cms2.trks_trk_p4()[trk1]);
	    if ( dR > 0.1 ) {
	      // if inv. mass of candidate within z mass window, veto event
	      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > 
		vec = cms2.trks_trk_p4()[trk2] + cms2.trks_trk_p4()[trk1];
	      if ( inZmassWindow(vec.mass()) ) {
		mode2 = true;
	      }
	    }
	  }
	}
      }
    }
  }

  if ( mode1 && mode2 ) return 3;
  if ( mode1 ) return 1;
  if ( mode2 ) return 2;
  
  return -1;
}


// count genp leptons
//-------------------------------------------------- 
// Returns the number of e,mu, and tau in the doc lines 
//----------------------------------------------------- 
void leptonGenpCount(int& nele, int& nmuon, int& ntau) { 
  nele=0; 
  nmuon=0; 
  ntau=0; 
  int size = cms2.genps_id().size(); 
  for (int jj=0; jj<size; jj++) { 
    if (abs(cms2.genps_id().at(jj)) == 11) nele++; 
    if (abs(cms2.genps_id().at(jj)) == 13) nmuon++; 
    if (abs(cms2.genps_id().at(jj)) == 15) ntau++; 
  } 
} 


// TTbar-->dilepton selections for Fall08-based analysis
double muonTrkIsolationPAT(int index){ 
  double sum =  cms2.mus_pat_trackIso().at(index); 
  double pt  = cms2.mus_p4().at(index).pt(); 
  return  pt/(pt+sum); 
} 
double muonCalIsolationPAT(int index){ 
  double sum =  cms2.mus_pat_caloIso().at(index); 
  double pt  = cms2.mus_p4().at(index).pt(); 
  return  pt/(pt+sum); 
} 

// met cut for ttbar dilepton analysis...
// includes a boolean to switch to tcmet
bool passMet_OF20_SF30(int hypIdx, bool useTcMet) {
  float mymet;
  if (useTcMet) {
    mymet = cms2.evt_tcmet();
  } else {
    mymet = cms2.met_pat_metCor();
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
  return passPatMet_OF20_SF30(cms2.met_pat_metCor()*cos(cms2.met_pat_metPhiCor()), 
			      cms2.met_pat_metCor()*sin(cms2.met_pat_metPhiCor()),
			      hypIdx);
}
//**************************************************************************
// met cut for ttbar dilepton analysis...
// includes a boolean to switch to tcmet
bool passMet_OF30_SF50(int hypIdx, bool useTcMet) {
  float mymet;
  if (useTcMet) {
    mymet = cms2.evt_tcmet();
  } else {
    mymet = cms2.met_pat_metCor();
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
  return passPatMet_OF30_SF50(cms2.met_pat_metCor()*cos(cms2.met_pat_metPhiCor()), 
			      cms2.met_pat_metCor()*sin(cms2.met_pat_metPhiCor()),
			      hypIdx);
}

//**************************************************************************

//-----------------------------------------------------------------------------------------------
//New selections for the common TTDil working group
//-----------------------------------------------------------------------------------------------
//
// loose lepton definitions 
//

bool electron20Eta2p4(int index){
  if (cms2.els_p4().at(index).pt() < 20 )   return false;
  if (fabs(cms2.els_p4().at(index).eta()) > 2.4 ) return false;
  
  return true;
}


bool looseElectronSelectionNoIsoTTDil08(int index) {
  if ( ! electron20Eta2p4(index) ) return false;
  if ( cms2.els_looseId().at(index)     !=  1) return false;
  if ( fabs(cms2.els_d0corr().at(index)) > 0.040)   return false;
  if ( cms2.els_closestMuon().at(index) != -1) return false; 

  return true;
}

//
double electronTrkIsolationPAT(int index){ 
  double sum = cms2.els_pat_trackIso().at(index); 
  double pt  = cms2.els_p4().at(index).pt(); 
  return pt/(pt+sum); 
} 
double electronCalIsolationPAT(int index){ 
  double sum = cms2.els_pat_caloIso().at(index); 
  double pt  = cms2.els_p4().at(index).pt(); 
  return pt/(pt+sum); 
} 


// factorize this stuff out
float electronTrkIsolationTTDil08(int index){
  return electronTrkIsolationPAT(index);
}
float electronCalIsolationTTDil08(int index){
  return electronCalIsolationPAT(index);
}

bool looseElectronSelectionTTDil08(int index) {
  if ( ! looseElectronSelectionNoIsoTTDil08(index) ) return false;

  if ( electronTrkIsolationTTDil08(index) < 0.5 ) return false;
  if ( electronCalIsolationTTDil08(index) < 0.5 ) return false;

  return true;
}

bool passElectronIsolationTTDil08(int index){
  if ( electronTrkIsolationTTDil08(index) < 0.9 ) return false;
  if ( electronCalIsolationTTDil08(index) < 0.8 ) return false;

  return true;
}

bool muon20Eta2p4(int index){
  if (cms2.mus_p4().at(index).pt() < 20 ) return false;
  if (fabs(cms2.mus_p4().at(index).eta()) >2.4 ) return false;

  return true;
}

bool looseMuonSelectionNoIsoTTDil08(int index) {

  if (! muon20Eta2p4(index) ) return false;

  if(!(2 & cms2.mus_type().at(index))) return false;
  if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) > 10.) return false;
  //  if (fabs(cms2.mus_d0corr().at(index))   > 0.25) return false;
  if (cms2.mus_validHits().at(index) < 11)    return false;
  
  return true;
}

bool lepton20Eta2p4(int id, int index){
  if (abs(id)==11) return electron20Eta2p4(index);
  if (abs(id)==13) return muon20Eta2p4(index);
  return false;
}

// factorize this stuff out
float muonTrkIsolationTTDil08(int index){
  return muonTrkIsolationPAT(index);
}
float muonCalIsolationTTDil08(int index){
  return muonCalIsolationPAT(index);
}

float leptonTrkIsolationTTDil08(int id, int index){
  if (abs(id) == 11) return electronTrkIsolationTTDil08(index);
  if (abs(id) == 13) return muonTrkIsolationTTDil08(index);
  return -1;
}
float leptonCalIsolationTTDil08(int id, int index){
  if (abs(id) == 11) return electronCalIsolationTTDil08(index);
  if (abs(id) == 13) return muonCalIsolationTTDil08(index);
  return -1;
}

bool looseMuonSelectionTTDil08(int index) {
  if (! looseMuonSelectionNoIsoTTDil08(index) ) return false;
  
  if ( muonTrkIsolationTTDil08(index) < 0.5 ) return false;
  if ( muonCalIsolationTTDil08(index) < 0.5 ) return false;

  return true;
}

bool passMuonIsolationTTDil08(int index) {
  if ( muonTrkIsolationTTDil08(index) < 0.9 ) return false;
  if ( muonCalIsolationTTDil08(index) < 0.9 ) return false;

  return true;
}

//now make it figure out which lepton type it is
bool passLeptonIsolationTTDil08(int id, int index){
  if (abs(id) == 11) return passElectronIsolationTTDil08(index);
  if (abs(id) == 13) return passMuonIsolationTTDil08(index);
  return false;
}

bool looseLeptonSelectionNoIsoTTDil08(int id, int index){
  if (abs(id) == 11) return looseElectronSelectionNoIsoTTDil08(index);  
  if (abs(id) == 13) return looseMuonSelectionNoIsoTTDil08(index);
  return false;
}
bool looseLeptonSelectionTTDil08(int id, int index){
  if (abs(id) == 11) return looseElectronSelectionTTDil08(index);  
  if (abs(id) == 13) return looseMuonSelectionTTDil08(index);
  return false;
}

//-----------------------------------------------------------------------------------------------
//New selections for the common SUSY Dilepton working group
//-----------------------------------------------------------------------------------------------
//
double inv_mu_rel_iso(int index)
{
  double sum =  cms2.mus_iso03_sumPt().at(index) +
    cms2.mus_iso03_emEt().at(index)  +
    cms2.mus_iso03_hadEt().at(index);
  double pt  = cms2.mus_p4().at(index).pt();
  return sum/pt;
}

double inv_el_rel_iso(int index, bool use_calo_iso)
{
  double sum = cms2.els_pat_trackIso().at(index);
  if (use_calo_iso)
    sum += cms2.els_pat_ecalIso()[index] + cms2.els_pat_hcalIso()[index];
  double pt  = cms2.els_p4().at(index).pt();
  return sum/pt;
}

bool passMuonIsolationVJets09(int index)
{
  const double cut = 0.1;
  return inv_mu_rel_iso(index) < cut;
}

bool passElectronIsolationVJets09(int index, bool use_calo_iso)
{
  const double cut = 0.1;
  return inv_el_rel_iso(index, use_calo_iso) < cut;
}

bool passLeptonIsolationVJets09(int id, int index){
  if (abs(id) == 11) return passElectronIsolationVJets09(index, true);
  if (abs(id) == 13) return passMuonIsolationVJets09(index);
  return false;
}

bool looseElectronSelectionVJets09(int index) {
  if (fabs(cms2.els_p4().at(index).eta()) > 2.5) return false;
  //  if ( cms2.els_tightId22XMinMatteo().at(index)     !=  1) return false;
  //  if ( cms2.els_looseId().at(index)     !=  1) return false;
  if ( cms2.els_pat_robustTightId().at(index)   < 0.5 ) return false; // SUSY group 

  if ( fabs(cms2.els_d0corr().at(index)) >= 0.2)   return false; // check if it is 0.2
  //  if ( cms2.els_closestMuon().at(index) != -1) return false;
  return true;
}

bool looseMuonSelectionVJets09(int index) {
  if (fabs(cms2.mus_p4().at(index).eta()) > 2.1) return false;
  if (((cms2.mus_type().at(index)) & (1<<1)) == 0) return false;
  if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) return false; // should < 10
  if (fabs(cms2.mus_d0corr().at(index))   >= 0.2) return false;
  if (cms2.mus_validHits().at(index) < 11)    return false;
  if (cms2.mus_pat_ecalvetoDep().at(index) >= 4) return false; // ECalE < 4
  if (cms2.mus_pat_hcalvetoDep().at(index) >= 6) return false; // HCalE < 6
  return true;
}
bool passLeptonIDVJets09(int id, int index){
  if (abs(id) == 11) return looseElectronSelectionVJets09(index);
  if (abs(id) == 13) return looseMuonSelectionVJets09(index);
  return false;
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

int numberOfExtraMuonsVJets09(int i_hyp){
  unsigned int nMuons = 0;
  for (int imu=0; imu < int(cms2.mus_p4().size()); imu++) {
    // quality cuts
    if ( cms2.mus_p4()[imu].pt() < 20 ) continue;
    if ( fabs(cms2.mus_p4()[imu].eta()) >  2.1 ) continue;
    if ((( cms2.mus_type()[imu]) & (1<<1)) == 0) continue;
    if (cms2.mus_gfit_chi2()[imu]/cms2.mus_gfit_ndof()[imu] >= 10) continue;
    if ( fabs(cms2.mus_d0corr()[imu]) >= 0.2) continue;
    if ( cms2.mus_validHits()[imu] < 11) continue;
    if (cms2.mus_pat_ecalvetoDep()[imu] >= 4) continue; 
    if (cms2.mus_pat_hcalvetoDep()[imu] >= 6) continue; 
    if ( TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && cms2.hyp_lt_index()[i_hyp] == imu ) continue;
    if ( TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && cms2.hyp_ll_index()[i_hyp] == imu ) continue;
    if ( inv_mu_rel_iso(imu) >= 0.1 ) continue;
    nMuons++;
  }
  return nMuons;
}


int numberOfExtraElectronsVJets09(int i_hyp){
  unsigned int nElec = 0;
  for (int iel=0; iel < int(cms2.els_p4().size()); iel++) {
    // quality cuts
    if ( cms2.els_p4()[iel].pt() < 20 ) continue;
    if ( fabs(cms2.els_p4()[iel].eta()) >  2.5 ) continue;
    if ( cms2.els_pat_robustTightId()[iel]   < 0.5  ) continue; // check
    if ( fabs(cms2.els_d0corr()[iel]) >= 0.2)  continue; 
    if ( TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11 && cms2.hyp_lt_index()[i_hyp] == iel ) continue;
    if ( TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11 && cms2.hyp_ll_index()[i_hyp] == iel ) continue;
    if ( inv_el_rel_iso(iel, true) >= 0.1 ) continue;
    nElec++;
  }
  return nElec;
}

//------------------------------------------------------------------------------------
// SUSY dilepton cuts 09 for TAS

bool compareEt(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > lv1, 
                 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > lv2) {
  return lv1.pt() > lv2.pt();
}

// Fakeble Objects

bool isFakeableElSUSY09(int iEl){
  
  Double_t pt = cms2.els_p4()[iEl].Pt();
  Double_t eta = cms2.els_p4()[iEl].Eta();
  
  if( pt < 10)  return false;
  if( fabs( eta ) > 2.4 )  return false;
  //reject if electron is close to muon;
  if ( cms2.els_closestMuon()[iEl] > -1)  return false;
  // Isolation
  if  (inv_el_relsusy_iso(iEl, true) > 0.4) return false;
  if ( conversionElectron(iEl))                             return false;
  if ( isChargeFlip(iEl))                                   return false;
  if(!GoodSusyElectronWithoutIsolationNoD0(iEl) )           return false;
  return true;
}

bool isChargeFlip(int elIndex){
  //true if electron is likely to be a charge flip
  if ((cms2.els_trkidx().at(elIndex) >= 0) && (cms2.els_charge().at(elIndex) != cms2.trks_charge().at(cms2.els_trkidx().at(elIndex)))) return true;

  return false;
}

bool isNumElSUSY09(int iEl) {
  
  Double_t pt = cms2.els_p4()[iEl].Pt();
  Double_t eta = cms2.els_p4()[iEl].Eta();
  if( pt < 10) return false;
  if( fabs( eta ) > 2.4 ) return false;
  //reject if electron is close to muon;
  if ( cms2.els_closestMuon()[iEl] > -1)  return false;
  if ( cms2.els_hOverE()[iEl]   > 0.2 )              return false;
  if (! GoodSusyElectronWithoutIsolation(iEl)) return false; 
  if (! GoodSusyElectronWithIsolation(iEl, true)) return false;
  if ( conversionElectron(iEl))                             return false;
  if ( isChargeFlip(iEl))                                   return false;
  return true;
}

bool isFakeableMuSUSY09(int iMu) {

  Double_t pt = cms2.mus_p4()[iMu].Pt();
  Double_t eta = cms2.mus_p4()[iMu].Eta();
  if (!((cms2.mus_type().at(iMu)) & (1<<1)) ) return false; // global muon
  if (!((cms2.mus_type().at(iMu)) & (1<<2)) ) return false; // tracker muon
  if( pt < 10)  return false;
  if( fabs( eta ) > 2.4 ) return false;
  if( cms2.mus_gfit_chi2()[iMu]/cms2.mus_gfit_ndof()[iMu] > 20) return false;
  if (inv_mu_relsusy_iso(iMu) > 0.4 ) return false;
  if (fabs(cms2.mus_d0corr().at(iMu))   >= 0.02) return false;
  return true;


}

bool isNumMuSUSY09(int iMu) {

  Double_t pt = cms2.mus_p4()[iMu].Pt();
  Double_t eta = cms2.mus_p4()[iMu].Eta();
  
  if (pt < 10)  return false; 
  if (fabs( eta ) > 2.4 ) return false; 
  if (!GoodSusyMuonWithoutIsolation(iMu)) return false;
  if (!GoodSusyMuonWithIsolation(iMu)) return false;
  return true;
}

bool GoodSusyElectronWithoutIsolation(int index) { 
  if ( cms2.els_tightId22XMinMatteo().at(index)     !=  1) return false; 
  if ( fabs(cms2.els_d0corr().at(index)) >= 0.02)   return false; 
  if ( cms2.els_closestMuon().at(index) != -1) return false; 
  if ( TMath::Abs(cms2.els_p4()[index].eta()) > 2.4) return false;
  return true; 
} 
 
bool GoodSusyMuonWithoutIsolation(int index) { 
  if (((cms2.mus_type().at(index)) & (1<<1)) == 0) return false; // global muon
  if (((cms2.mus_type().at(index)) & (1<<2)) == 0) return false; // tracker muon
  if (cms2.mus_validHits().at(index) < 11)    return false; 
  if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) return false; 
  if (fabs(cms2.mus_d0corr().at(index))   >= 0.02) return false; 
  if (cms2.mus_pat_ecalvetoDep().at(index) >= 4) return false; // ECalE < 4 
  if (cms2.mus_pat_hcalvetoDep().at(index) >= 6) return false; // HCalE < 6 
  if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.4) return false;
  return true; 
} 

double inv_mu_relsusy_iso(int index)
{
  double sum =  cms2.mus_iso03_sumPt().at(index) +
    cms2.mus_iso03_emEt().at(index)  +
    cms2.mus_iso03_hadEt().at(index);
  double pt  = cms2.mus_p4().at(index).pt();
  return sum/max(pt,20.);
}

double inv_el_relsusy_iso(int index, bool use_calo_iso)
{
  double sum = cms2.els_pat_trackIso().at(index);
  if (use_calo_iso)
     sum += max(0., (cms2.els_pat_ecalIso().at(index) -2.)) + cms2.els_pat_hcalIso().at(index);
  double pt  = cms2.els_p4().at(index).pt();
  return sum/max(pt,20.);
}

bool GoodSusyMuonWithIsolation(int index)
{
  const double cut = 0.1;
  return inv_mu_relsusy_iso(index) < cut;
}

bool GoodSusyElectronWithIsolation(int index, bool use_calo_iso)
{
  const double cut = 0.1;
  return inv_el_relsusy_iso(index, use_calo_iso) < cut;
}

bool GoodSusyLeptonIsolation(int id, int index){
  if (abs(id) == 11) return GoodSusyElectronWithIsolation(index, true);
  if (abs(id) == 13) return GoodSusyMuonWithIsolation(index);
  return false;
}

bool PassSusyElectronIsolationLoose(int index, bool use_calo_iso)
{
  const double cut = 0.4; 
  return inv_el_relsusy_iso(index, use_calo_iso) < cut;
}

bool PassSusyElectronIsolation(int index, bool use_calo_iso)
{
  const double cut = 0.1;
  return inv_el_relsusy_iso(index, use_calo_iso) < cut;
}


bool GoodSusyLeptonID(int id, int index){ 
  if (abs(id) == 11) return GoodSusyElectronWithoutIsolation(index); 
  if (abs(id) == 13) return GoodSusyMuonWithoutIsolation(index); 
  return false; 
}    

bool GoodSusyTrigger(int dilType){
  bool hlt_ele15_lw_l1r = cms2.passHLTTrigger("HLT_Ele15_SW_L1R");
  bool hltMu9           = cms2.passHLTTrigger("HLT_Mu9");
  // bool hltdiMu3         = cms2.passHLTTrigger("HLT_DoubleMu3");
  //  bool hltdiEle10       = cms2.passHLTTrigger("HLT_DoubleEle10_SWL1R");
  // bool hltdiEle10       = cms2.passHLTTrigger("HLT_DoubleEle5_SW_L1R");

  if (dilType == 0 && ! (hltMu9) ) return false;
  if ((dilType == 1 || dilType == 2) && ! (hltMu9 || hlt_ele15_lw_l1r)) return false;
  if (dilType == 3 && ! hlt_ele15_lw_l1r) return false;

  return true;
}

bool GoodSusyElectronWithoutIsolationNoD0(int index) { 
  if ( cms2.els_tightId22XMinMatteo().at(index)     !=  1) return false; 
  if ( cms2.els_closestMuon().at(index) != -1) return false; 
  if ( TMath::Abs(cms2.els_p4()[index].eta()) > 2.4) return false;
  return true; 
} 

int numberOfExtraElectronsSUSY(int i_hyp){ 
  unsigned int nElec = 0; 
  for (int iel=0; iel < int(cms2.els_p4().size()); iel++) { 
    if ( cms2.els_p4()[iel].pt() < 10 ) continue; 
    if ( fabs(cms2.els_p4()[iel].eta()) > 2.4 ) continue; 
    if (! GoodSusyElectronWithoutIsolation(iel)) continue; 
    if (! GoodSusyElectronWithIsolation(iel, true)) continue;
    if ( TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11 && cms2.hyp_lt_index()[i_hyp] == iel ) continue; 
    if ( TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11 && cms2.hyp_ll_index()[i_hyp] == iel ) continue; 
    nElec++; 
  } 
  return nElec; 
} 

int numberOfExtraMuonsSUSY(int i_hyp){ 
  unsigned int nMuons = 0; 
  for (int imu=0; imu < int(cms2.mus_p4().size()); imu++) { 
    if ( cms2.mus_p4()[imu].pt() < 10 ) continue; 
    if ( fabs(cms2.mus_p4()[imu].eta()) > 2.4 ) continue; 
    if (!GoodSusyMuonWithoutIsolation(imu)) continue;
    if (!GoodSusyMuonWithIsolation(imu)) continue; 
    if ( TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && cms2.hyp_lt_index()[i_hyp] == imu ) continue; 
    if ( TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && cms2.hyp_ll_index()[i_hyp] == imu ) continue; 
    nMuons++; 
  } 
  return nMuons; 
} 

vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > getCaloJets(int i_hyp) {
  vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > calo_jets;
  calo_jets.clear();
  
  for (unsigned int jj=0; jj < cms2.jets_p4().size(); ++jj) {
    if ((dRbetweenVectors(cms2.hyp_lt_p4()[i_hyp],cms2.jets_p4()[jj]) < 0.4)||
	(dRbetweenVectors(cms2.hyp_ll_p4()[i_hyp],cms2.jets_p4()[jj]) < 0.4)
	) continue;
    if (cms2.jets_p4()[jj].pt() < 30) continue;
    if (fabs(cms2.jets_p4()[jj].Eta()) > 2.4) continue;
    //    if (cms2.jets_emFrac()[jj] < 0.1) continue;
    calo_jets.push_back(cms2.jets_p4()[jj]);
  }
  
  if (calo_jets.size() > 1) {
    sort(calo_jets.begin(), calo_jets.end(),  compareEt);
  }
  return calo_jets;
}


vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > getJPTJets(int i_hyp) {
  vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > jpt_jets;
  jpt_jets.clear();
  
  for (unsigned int jj=0; jj < cms2.jpts_p4().size(); ++jj) {
    if ((dRbetweenVectors(cms2.hyp_lt_p4()[i_hyp],cms2.jpts_p4()[jj]) < 0.4)||
	(dRbetweenVectors(cms2.hyp_ll_p4()[i_hyp],cms2.jpts_p4()[jj]) < 0.4)
	) continue;
    if (cms2.jpts_p4()[jj].Et() < 30) continue;
    if (fabs(cms2.jpts_p4()[jj].Eta()) > 2.4) continue;
//    cout << cms2.jpts_emFrac()[jj] << endl;
//    if (cms2.jpts_emFrac()[jj] < 0.1) continue;
    jpt_jets.push_back(cms2.jpts_p4()[jj]);
  }
  
  if (jpt_jets.size() > 1) {
    sort(jpt_jets.begin(), jpt_jets.end(),  compareEt);
  }
  return jpt_jets;
}

int ttbarconstituents(int i_hyp){ 

  // Catagories WW = 1, WO = 2, and OO = 3

  bool isTrueLepton_ll = false;
  bool isTrueLepton_lt = false;

  isTrueLepton_ll = ( (abs(cms2.hyp_ll_id()[i_hyp]) == abs(cms2.hyp_ll_mc_id()[i_hyp]) &&
		       abs(cms2.hyp_ll_mc_motherid()[i_hyp]) < 50 //I wish I could match to W or Z explicitely, not in MGraph
		       )
		      || (cms2.hyp_ll_mc_id()[i_hyp]==22 && 
			  TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.hyp_ll_mc_p4()[i_hyp])) <0.05
			  && abs(cms2.hyp_ll_id()[i_hyp]) == abs(cms2.hyp_ll_mc_motherid()[i_hyp])
			  )
		      );

  isTrueLepton_lt = ( (abs(cms2.hyp_lt_id()[i_hyp]) == abs(cms2.hyp_lt_mc_id()[i_hyp]) &&
		       abs(cms2.hyp_lt_mc_motherid()[i_hyp]) < 50 //I wish I could match to W or Z explicitely, not in MGraph
		       )
		      || (cms2.hyp_lt_mc_id()[i_hyp]==22 && 
			  TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.hyp_lt_mc_p4()[i_hyp])) <0.05
			  && abs(cms2.hyp_lt_id()[i_hyp]) == abs(cms2.hyp_lt_mc_motherid()[i_hyp])
			  )
		      );

  bool ttbarlep = false;
  bool ttbarother = false;

  //  if (genpCountPDGId(11,13,15) == 1) ttbarother = true;
  if (genpCountPDGId(11,13,15) == 2) ttbarlep = true;
  

  bool isrealW_ll = false;
  bool isrealW_lt = false;     

  { 
    int els_mo = 0; 
    int mus_mo = 0; 
    int els_id = 0;
    int mus_id = 0;
    if (abs(cms2.hyp_ll_id()[i_hyp]) == 11) els_mo = abs(cms2.els_mc3_motherid()[cms2.hyp_ll_index()[i_hyp]]); 
    if (abs(cms2.hyp_ll_id()[i_hyp]) == 13) mus_mo = abs(cms2.mus_mc3_motherid()[cms2.hyp_ll_index()[i_hyp]]); 
    if (abs(cms2.hyp_ll_id()[i_hyp]) == 11) els_id = abs(cms2.els_mc3_id()[cms2.hyp_ll_index()[i_hyp]]);
    if (abs(cms2.hyp_ll_id()[i_hyp]) == 13) mus_id = abs(cms2.mus_mc3_id()[cms2.hyp_ll_index()[i_hyp]]);
    if ( (abs(cms2.hyp_ll_mc_motherid()[i_hyp]) == 24) || (els_mo == 24) || (mus_mo == 24) || (els_id == 24) || (mus_id == 24)) isrealW_ll = true; 
  } 
 
  { 
    int els_mo = 0;  
    int mus_mo = 0;  
    int els_id = 0;
    int mus_id = 0;
    if (abs(cms2.hyp_lt_id()[i_hyp]) == 11) els_mo = abs(cms2.els_mc3_motherid()[cms2.hyp_lt_index()[i_hyp]]); 
    if (abs(cms2.hyp_lt_id()[i_hyp]) == 13) mus_mo = abs(cms2.mus_mc3_motherid()[cms2.hyp_lt_index()[i_hyp]]); 
    if (abs(cms2.hyp_lt_id()[i_hyp]) == 11) els_id = abs(cms2.els_mc3_id()[cms2.hyp_lt_index()[i_hyp]]);
    if (abs(cms2.hyp_lt_id()[i_hyp]) == 13) mus_id = abs(cms2.mus_mc3_id()[cms2.hyp_lt_index()[i_hyp]]);
    if ((abs(cms2.hyp_lt_mc_motherid()[i_hyp]) == 24) || (els_mo == 24) || (mus_mo == 24) || (els_id == 24) || (mus_id == 24)) isrealW_lt = true; 
  } 

  bool isLepll = false;
  bool isLeplt = false;

  if (isTrueLepton_ll || isrealW_ll) isLepll = true;
  if (isTrueLepton_lt || isrealW_lt) isLeplt = true;
  
  if (isLepll && isLeplt && ttbarlep) {
    return 1;
  }  else if (isLepll || isLeplt )  {
    return 2;
  } else {
    return 3;
  }

}

bool additionalZvetoSUSY09(int i_hyp) {

  // true if we want to veto this event
  bool veto=false;

  // first, look for Z->mumu
  for (unsigned int i=0; i < cms2.mus_p4().size(); i++) {
    bool hypLep1 = false;
    if (cms2.mus_p4().at(i).pt() < 10.)     continue;
    if (!GoodSusyMuonWithoutIsolation(i)) continue;
    if ( TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && cms2.hyp_lt_index()[i_hyp] == i ) hypLep1 = true;
    if ( TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && cms2.hyp_ll_index()[i_hyp] == i ) hypLep1 = true;

    for (unsigned int j=i+1; j < cms2.mus_p4().size(); j++) {
      bool hypLep2 = false;
      if (cms2.mus_p4().at(j).pt() < 10.) continue;
      if (!GoodSusyMuonWithoutIsolation(j)) continue;
      if (cms2.mus_charge().at(i) == cms2.mus_charge().at(j)) continue;
      if ( TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && cms2.hyp_lt_index()[i_hyp] == j ) hypLep2 = true;
      if ( TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && cms2.hyp_ll_index()[i_hyp] == j ) hypLep2 = true;
      // At least one of them has to pass isolation
      if (!GoodSusyMuonWithIsolation(i) && !GoodSusyMuonWithIsolation(j)) continue;
      if ( hypLep1 && hypLep2 ) continue;
      if ( !hypLep1 && !hypLep2 ) continue;
      // Make the invariant mass
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >
                vec = cms2.mus_p4().at(i) + cms2.mus_p4().at(j);
      if ( inZmassWindow(vec.mass()) ) return true;

    }
  }

  // now, look for Z->ee
  for (unsigned int i=0; i < cms2.els_p4().size(); i++) {
    bool hypLep1 = false;
    if (cms2.els_p4().at(i).pt() < 10.)     continue;
    if (! GoodSusyElectronWithoutIsolation(i)) continue;
    if ( TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11 && cms2.hyp_lt_index()[i_hyp] == i ) hypLep1 = true;
    if ( TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11 && cms2.hyp_ll_index()[i_hyp] == i ) hypLep1 = true;

    for (unsigned int j=i+1; j<cms2.els_p4().size(); j++) {
      bool hypLep2 = false;
      if (cms2.els_p4().at(j).pt() < 10.) continue;
      if (! GoodSusyElectronWithoutIsolation(j)) continue;
      if (cms2.els_charge().at(i) == cms2.els_charge().at(j)) continue;
      // At least one of them has to pass isolation
      if (! GoodSusyElectronWithIsolation(i, true) && ! GoodSusyElectronWithIsolation(j, true)) continue;
      if ( TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11 && cms2.hyp_lt_index()[i_hyp] == j ) hypLep2 = true;
      if ( TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11 && cms2.hyp_ll_index()[i_hyp] == j ) hypLep2 = true;
      if ( hypLep1 && hypLep2 ) continue;
      if ( !hypLep1 && !hypLep2 ) continue;
      // Make the invariant mass
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >
                vec = cms2.els_p4().at(i) + cms2.els_p4().at(j);
      if ( inZmassWindow(vec.mass()) ) return true;

    }
  }
  // done
  return veto;
}



//--------------------------------------------------------------------
// Veto events if there are two leptons in the 
// event that make the Z mass.  This uses the mus and els
// blocks, ie, it is a veto that can use the 3rd (4th,5th,..)
// lepton in the event.
//
// Both leptons must be 20 GeV, and pass the same cuts as 
// the hypothesis leptons, except that one of them can be non-isolated
//---------------------------------------------------------------------
bool additionalZvetoTTDil08() {

  // true if we want to veto this event
  bool veto=false;

  // first, look for Z->mumu
  for (unsigned int i=0; i<cms2.mus_p4().size(); i++) {
    if (cms2.mus_p4().at(i).pt() < 20.)     continue;
    if (!looseMuonSelectionNoIsoTTDil08(i)) continue;

    for (unsigned int j=i+1; j<cms2.mus_p4().size(); j++) {
      if (cms2.mus_p4().at(j).pt() < 20.) continue;
      if (!looseMuonSelectionNoIsoTTDil08(j)) continue;

      if (cms2.mus_charge().at(i) == cms2.mus_charge().at(j)) continue;

      // At least one of them has to pass isolation
      if (!passMuonIsolationTTDil08(i) && !passMuonIsolationTTDil08(j)) continue;

      // Make the invariant mass
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > 
	vec = cms2.mus_p4().at(i) + cms2.mus_p4().at(j);
      if ( inZmassWindow(vec.mass()) ) return true;

    }
  }

  // now, look for Z->ee
  for (unsigned int i=0; i<cms2.evt_nels(); i++) {
    if (cms2.els_p4().at(i).pt() < 20.)     continue;
    if (!looseElectronSelectionNoIsoTTDil08(i)) continue;

    for (unsigned int j=i+1; j<cms2.evt_nels(); j++) {
      if (cms2.els_p4().at(j).pt() < 20.) continue;
      if (!looseElectronSelectionNoIsoTTDil08(j)) continue;

      if (cms2.els_charge().at(i) == cms2.els_charge().at(j)) continue;

      // At least one of them has to pass isolation
      if (!passElectronIsolationTTDil08(i) && !passElectronIsolationTTDil08(j)) continue;

      // Make the invariant mass
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > 
	vec = cms2.els_p4().at(i) + cms2.els_p4().at(j);
      if ( inZmassWindow(vec.mass()) ) return true;

    }
  }
  // done
  return veto;
}

// 
// true if there is a muon not in the hypothesis
bool haveExtraMuon(int hypIdx){
  bool result = false;

  int nHMus = 0;
  if (abs(cms2.hyp_lt_id().at(hypIdx))==13){
    nHMus++;
  }
  if (abs(cms2.hyp_ll_id().at(hypIdx))==13){
    nHMus++;
  }

  int nEvtMus = cms2.mus_p4().size();
  result = (nEvtMus - nHMus) > 0;

  return result;
}

// true if there is a muon not in the hypothesis
bool haveExtraMuon5(int hypIdx){
  double minPtCut = 5;
  bool result = false;

  int nHMus = 0;
  if (abs(cms2.hyp_lt_id().at(hypIdx))==13){
    if (cms2.hyp_lt_p4().at(hypIdx).pt() > minPtCut ) nHMus++;
  }
  if (abs(cms2.hyp_ll_id().at(hypIdx))==13){
    if (cms2.hyp_ll_p4().at(hypIdx).pt() > minPtCut) nHMus++;
  }

  int nEvtMus = 0;
  for (unsigned int iMu = 0; iMu < cms2.mus_p4().size(); ++iMu){
    if (cms2.mus_p4().at(iMu).pt() > minPtCut) nEvtMus++;
  }
  result = (nEvtMus - nHMus) > 0;

  return result;
}



// Trigger-related selections

//Bits passing for hyp type
bool passTriggersMu9orLisoE15(int dilType) {
  //old bit based method
  //bool hlt_ele15_lw_l1r = ((cms2.evt_HLT2() & (1<<(50-32))) != 0);
  //bool hltMu9 = ((cms2.evt_HLT3() & (1<<(82-64))) != 0);
  
  //TString method
  bool hlt_ele15_lw_l1r = cms2.passHLTTrigger("HLT_Ele15_SW_L1R");
  bool hltMu9           = cms2.passHLTTrigger("HLT_Mu9");
  
  if (dilType == 0 && ! (hltMu9) ) return false;
  if ((dilType == 1 || dilType == 2) && ! (hltMu9 || hlt_ele15_lw_l1r)) return false;
  if (dilType == 3 && ! hlt_ele15_lw_l1r) return false;     

  return true;
}


bool passTriggersTTDil08JanTrial(int dilType) {
  //trigger selections used in AN09/047 (at least as of v4 on 04-10-09
  bool hlt_Mu15_L1Mu7 = cms2.passHLTTrigger("HLT_Mu15_L1Mu7");
  bool hlt_DoubleMu3 = cms2.passHLTTrigger("HLT_DoubleMu3");
  bool hlt_IsoEle10_Mu10_L1R = cms2.passHLTTrigger("HLT_IsoEle10_Mu10_L1R");
  bool passMuMutriggers = (hlt_Mu15_L1Mu7 ||  hlt_DoubleMu3 ||  hlt_IsoEle10_Mu10_L1R);
  
  bool hlt_IsoEle18_L1R = cms2.passHLTTrigger("HLT_IsoEle18_L1R");
  bool hlt_DoubleIsoEle12_L1R = cms2.passHLTTrigger("HLT_DoubleIsoEle12_L1R");
  bool passEEtriggers = (hlt_IsoEle18_L1R || hlt_DoubleIsoEle12_L1R );

  bool passEMutriggers = (hlt_Mu15_L1Mu7 || hlt_IsoEle18_L1R || hlt_IsoEle10_Mu10_L1R || hlt_DoubleMu3);

  if (dilType == 0 && ! passMuMutriggers ) return false;
  if ((dilType == 1 || dilType == 2)  && ! passEMutriggers ) return false;
  if (dilType == 3 && ! passEEtriggers ) return false; 
  return true;
}


int genpCountPDGId(int id0, int id1, int id2){ 
  int count = 0; 
  int size = cms2.genps_id().size(); 
  for (int jj=0; jj<size; jj++) { 
    if (abs(cms2.genps_id().at(jj)) == id0) count++; 
    if (abs(cms2.genps_id().at(jj)) == id1) count++; 
    if (abs(cms2.genps_id().at(jj)) == id2) count++; 
  } 
  return count; 
} 


int genpDileptonType(){
  //0 mumu; 1 emu; 2 ee
  
  unsigned int nmus = 0;
  unsigned int nels = 0;
  int size = cms2.genps_id().size();
  for (int jj=0; jj<size; jj++) {
    if (abs(cms2.genps_id().at(jj)) == 11) nels++;
    if (abs(cms2.genps_id().at(jj)) == 13) nmus++;
  }

  if ((nels + nmus) != 2){
    return -1;
  }

  int dilType = -1;
  if (nmus == 2) dilType = 0;
  if (nels == 2) dilType = 2;
  if (nels == 1 && nmus == 1) dilType = 1;
  return dilType;
}



int eventDilIndexByWeightTTDil08(const std::vector<unsigned int>& goodHyps, int& strasbourgDilType, bool printDebug){
  int result = -1;
  unsigned int nGoodHyps = goodHyps.size();
  if ( nGoodHyps == 0 ) return result;

  float maxWeight = -1;
  unsigned int maxWeightIndex = 9999;
  
  for (unsigned int hypIdxL=0; hypIdxL < nGoodHyps; ++hypIdxL){
    unsigned int hypIdx = goodHyps[hypIdxL];
    float hypWeight_lt = 0;
    float hypWeight_ll = 0;
    float hypWeight_iso = 0;
    float hypWeight = 0;
    unsigned int i_lt = cms2.hyp_lt_index().at(hypIdx);
    unsigned int i_ll = cms2.hyp_ll_index().at(hypIdx);

    int id_lt = cms2.hyp_lt_id().at(hypIdx);
    int id_ll = cms2.hyp_ll_id().at(hypIdx);

    float isoTk_lt = leptonTrkIsolationTTDil08(id_lt, i_lt);
    float isoTk_ll = leptonTrkIsolationTTDil08(id_ll, i_ll);

    float isoCal_lt = leptonCalIsolationTTDil08(id_lt, i_lt);
    float isoCal_ll = leptonCalIsolationTTDil08(id_ll, i_ll);

    //ad-hoc selection of weights
    if (abs(id_lt) == 11){
      //I want to select "trk & cal"-isolated ones
      hypWeight_iso += (isoTk_lt*isoCal_lt - 0.25); //shift by 0.25 to be positive-definite
      if (cms2.els_tightId().at(i_lt)) hypWeight_lt += 0.2;
    }
    if (abs(id_lt) == 13){
      //I want to select "trk & cal"-isolated ones	    
      hypWeight_iso += (isoTk_lt*isoCal_lt - 0.25);//shift by 0.25 to be positive-definite
      hypWeight_lt += 0.4;
    }
    if (abs(id_ll) == 11){
      //I want to select "trk & cal"-isolated ones
      hypWeight_iso *= (isoTk_ll*isoCal_ll - 0.25); //shift by 0.25 to be positive-definite
      if (cms2.els_tightId().at(i_ll)) hypWeight_ll += 0.2;
    }
    if (abs(id_ll) == 13){
      //I want to select "trk & cal"-isolated ones
      hypWeight_iso *= (isoTk_ll*isoCal_ll - 0.25); //shift by 0.25 to be positive-definite
      hypWeight_ll += 0.4;
    }
    float pt_lt = cms2.hyp_lt_p4().at(hypIdx).pt();
    float pt_ll = cms2.hyp_ll_p4().at(hypIdx).pt();
    hypWeight_lt += (1. - 20./pt_lt*20./pt_lt);
    hypWeight_ll += (1. - 20./pt_ll*20./pt_ll);

    hypWeight = hypWeight_ll*hypWeight_lt*hypWeight_iso; //again, desire to have both good

    if (hypWeight > maxWeight){
      maxWeight = hypWeight;
      maxWeightIndex = hypIdx;
    }
  }


  //Now let's implement the Strasbourg-type disambiguation/dispatch
  //ee
  {
    std::vector<unsigned int> looseEls(0);
    std::vector<unsigned int> looseMus(0);
    for (unsigned int iEl =0; iEl < cms2.els_p4().size(); ++iEl){
      if (looseElectronSelectionTTDil08(iEl)){
	looseEls.push_back(iEl);
      }
    }
    for (unsigned int iMu =0; iMu < cms2.mus_p4().size(); ++iMu){
      if (looseMuonSelectionTTDil08(iMu)){
	looseMus.push_back(iMu);
      }
    }
	    
    bool pass_elec = false;
    if (looseEls.size()>1){
      if (cms2.els_charge().at(looseEls[0]) != cms2.els_charge().at(looseEls[1])){
	pass_elec = true;
      }
      if (looseMus.size()>0 && cms2.mus_p4().at(looseMus[0]).pt() > cms2.els_p4().at(looseEls[1]).pt()) pass_elec = false;
      if (looseMus.size()>0 && 
	  ( ( muonTrkIsolationTTDil08(looseMus[0]) > electronTrkIsolationTTDil08(looseEls[0]) 
	      && cms2.mus_charge().at(looseMus[0]) != cms2.els_charge().at(looseEls[0]) )
	    || ( muonTrkIsolationTTDil08(looseMus[0]) > electronTrkIsolationTTDil08(looseEls[1])
		 && cms2.mus_charge().at(looseMus[0]) != cms2.els_charge().at(looseEls[0]))
	    )
	  ) pass_elec = false; 
    }
    bool pass_muon = false;
    if (looseMus.size()>1){
      for (unsigned int iMu=0; iMu < looseMus.size(); ++iMu){
	for (unsigned int jMu=iMu+1; jMu < looseMus.size(); ++jMu){
	  if (cms2.mus_charge().at(looseMus[iMu]) != cms2.mus_charge().at(looseMus[jMu])) pass_muon = true;
	}
      }
      if (looseEls.size()>0 && cms2.els_p4().at(looseEls[0]).pt() > cms2.mus_p4().at(looseMus[1]).pt()
	  && cms2.mus_charge().at(looseMus[1]) != cms2.els_charge().at(looseEls[0])) pass_muon = false;
    }
    bool pass_elecmuon = false;
    if (looseMus.size() > 0 && looseEls.size() > 0){
      if (! pass_elec && ! pass_muon ){
	if (cms2.mus_charge().at(looseMus[0]) != cms2.els_charge().at(looseEls[0])) pass_elecmuon = true;
	if (! pass_elecmuon && looseEls.size()>1){
	  if (cms2.mus_charge().at(looseMus[0]) != cms2.els_charge().at(looseEls[0])) pass_elecmuon = true;
	}
      }
    }

    unsigned int passStatus = 0;
    if (pass_muon) passStatus++;
    if (pass_elecmuon) passStatus++;
    if (pass_elec) passStatus++;
    if (passStatus > 1) std::cout<<"ERROR: inconsistent assignment"<<std::endl;
    if (passStatus == 1){
      if (pass_muon) strasbourgDilType = 0;
      if (pass_elecmuon) strasbourgDilType = 1;
      if (pass_elec) strasbourgDilType = 2;
    }
  }

  if (printDebug){
    int genpDilType = genpDileptonType();
    if (genpDilType>=0 ){ std::cout<<"Dil type "<<genpDilType<<std::endl;
      if (nGoodHyps > 1){
	int maxWeightType = cms2.hyp_type().at(maxWeightIndex);
	if ((maxWeightType == 0 && genpDilType == 0)
	    || ( (maxWeightType == 1 || maxWeightType == 2) && genpDilType == 1)
	    || (maxWeightType == 3 && genpDilType == 2)){
	  std::cout<<"Dil type "<<genpDilType<<" ; Strasbourg dil type "<<strasbourgDilType 
		   <<" assigned correctly by maxWeight method";
	  std::cout<<" out of"; for(unsigned int iih=0;iih<nGoodHyps;++iih)std::cout<<" "<<cms2.hyp_type().at(goodHyps[iih]);
	  std::cout<<std::endl;
	} else {
	  std::cout<<"Dil type "<<genpDilType<<" ; Strasbourg dil type "<<strasbourgDilType 
		   <<" assigned incorrectly by maxWeight method";
	  std::cout<<" out of"; for(unsigned int iih=0;iih<nGoodHyps;++iih)std::cout<<" "<<cms2.hyp_type().at(goodHyps[iih]);
	  std::cout<<std::endl;	    
	}
      }
    }
  }

  result = maxWeightIndex;
  return result;
}


//-------------------------------------------------------------
// TTDil08 Fake Rate selections (FO, numerator selections)
//-------------------------------------------------------------


//Is the electron a Numerator electron? 
//Same selections as the TTDil selections for electrons
bool isNumElTTDil08(int iEl) {
  
  Double_t pt = cms2.els_p4()[iEl].Pt();
  Double_t eta = cms2.els_p4()[iEl].Eta();
  
   if( pt < 20)
    return false;
  if( fabs( eta ) > 2.4 )  
    return false;
  //reject if electron is close to muon;
  if( cms2.els_closestMuon()[iEl] > -1) 
    return false;

  //now the numerator cuts
  //add in track isolation
  if( pt/(pt + cms2.els_pat_trackIso()[iEl]) < 0.9)
    return false;
  //add in calo isolation
  if( pt/(pt + cms2.els_pat_caloIso()[iEl]) < 0.8)
    return false;
  //corrected d0 cut
  if( fabs(cms2.els_d0corr()[iEl]) > 0.04) 
    return false;
  //electron id cut
  if( ! cms2.els_looseId()[iEl] )
    return false;

  return true;
}
//------------------------------------------------------------
// is it a FO electron?
//------------------------------------------------------------
bool isFakeableElTTDil08(int iEl) {
  
  Double_t pt = cms2.els_p4()[iEl].Pt();
  Double_t eta = cms2.els_p4()[iEl].Eta();
  
  //base selections:
  //make up some loose electron selection for now
  if( pt < 20)
    return false;
  if( fabs( eta ) > 2.4 )  
    return false;
  //reject if electron is close to muon;
  if( cms2.els_closestMuon()[iEl] > -1) 
    return false;
        
  //check if the electron passes the FO
  //object selections
  //add in track isolation
  if( pt/(pt + cms2.els_pat_trackIso()[iEl]) < 0.7)   //0.7 
    return false;
  //add in calo isolation
  if( pt/(pt + cms2.els_pat_caloIso()[iEl]) < 0.6)   //0.6 
    return false;
        
  return true;
}

//------------------------------------------------------------
//is numerator muon?
//------------------------------------------------------------

bool isNumMuTTDil08(int iMu) {

  Double_t pt = cms2.mus_p4()[iMu].Pt();
  Double_t eta = cms2.mus_p4()[iMu].Eta();

  if(!(2 & cms2.mus_type()[iMu]))
    return false;

  if( pt < 20)
    return false;
  if( fabs( eta ) > 2.4 )  
    return false;

  //now the numerator cuts
  if( cms2.mus_gfit_chi2()[iMu]/cms2.mus_gfit_ndof()[iMu] > 10)
    return false;
  //add in track isolation
  if( pt/(pt + cms2.mus_pat_trackIso()[iMu]) < 0.9)
    return false;
  //add in calo isolation
  if( pt/(pt + cms2.mus_pat_caloIso()[iMu]) < 0.9)
    return false;

  //nHits cut
  if( cms2.mus_validHits()[iMu] < 11 )
    return false;

  return true;
}
                            
//------------------------------------------------------------
//   is FO Mu
//------------------------------------------------------------

bool isFakeableMuTTDil08(int iMu) {

  Double_t pt = cms2.mus_p4()[iMu].Pt();
  Double_t eta = cms2.mus_p4()[iMu].Eta();
  
  //base selections:
  //loose muon selection
  //only globalMuons
  if(!(2 & cms2.mus_type()[iMu]))
    return false;
  
  if( pt < 20)
    return false;
  if( fabs( eta ) > 2.4 )  
    return false;
  if( cms2.mus_gfit_chi2()[iMu]/cms2.mus_gfit_ndof()[iMu] > 20)
    return false;
  //check if the muon passes the FO
  //object selections
  //add in track isolation
  if( pt/(pt + cms2.mus_pat_trackIso()[iMu]) < 0.7) 
    return false;
  //add in calo isolation
  if( pt/(pt + cms2.mus_pat_caloIso()[iMu]) < 0.7) 
    return false;

  return true;
}

//------------------------------------------------------------

bool trueGammaFromMuon(int electron) {
  // true gamma reconstructed as electron 
  // gamma coming from muon
  if(TMath::Abs(cms2.els_mc_id()[electron]) == 22 && TMath::Abs(cms2.els_mc_motherid()[electron]) == 13) { // ok, abs of photon makes no sense ;)
    //    std::cout<<"Gamma from muon event - r: " << cms2.evt_run() << " e: " << cms2.evt_event() << " l: " << cms2.evt_lumiBlock() << std::endl;
    return true;
  }
//   if( cms2.els_mc_motherid()[electron] == 22 ) {
//     std::cout<<"Electron with gamma mother - r: " << cms2.evt_run() << " e: " << cms2.evt_event() << " l: " << cms2.evt_lumiBlock() << std::endl;
//   }
//   if( cms2.els_mc_id()[electron] == 22 ) {
//     std::cout<<"***"<<std::endl;
//     std::cout<<"Electron which is a really a gamma - r: " << cms2.evt_run() << " e: " << cms2.evt_event() << " l: " << cms2.evt_lumiBlock() << std::endl;
//     std::cout<<"el mc id: "<<cms2.els_mc_id()[electron]<<" el mother id "<< cms2.els_mc_motherid()[electron]<<std::endl;
//     std::cout<<"***"<<std::endl;
//   }

  return false;
}
//----------------------------------------------------
// Utility function to calculate dR between vectors
//-----------------------------------------------------
double dRBetweenVectors(LorentzVector v1, LorentzVector v2) {
    double deta = v1.eta() - v2.eta();
    double dphi = fabs(v1.phi() - v2.phi());
    if (dphi > TMath::Pi()) dphi = TMath::TwoPi() - dphi;
    return sqrt(deta*deta+dphi*dphi);
}

bool conversionElectron(int electron) {
  // true if electron is a conversion electron
  if( fabs(cms2.els_conv_dist()[electron]) < 0.02 && fabs(cms2.els_conv_dcot()[electron]) < 0.02)
    return true;

  return false;
}
