
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
// #include "selections.h"
//#include <vector>
//#include "CMS1.h"
#include "TDatabasePDG.h"
#include "selections.h"
#ifdef TOOLSLIB
#include "CMS2.h"
#include "matchTools.h"
#include "utilities.h"
#else
#include "../Tools/matchTools.C"
#endif

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
// bool trueElectron(int index) {
//   if ( abs(cms2.els_mc_id()[index]) != 11 ) return false;
//   if ( abs(cms2.els_mc_motherid()[index]) == 23 || abs(cms2.els_mc_motherid()[index]) == 24) return true;
//   return false;
// }
// //----------------------------------------------------------------
// // true muon
// //----------------------------------------------------------------
// bool trueMuon(int index) {
//    if ( abs(cms2.mus_mc_id()[index]) != 13 ) return false;
//    if ( abs(cms2.mus_mc_motherid()[index]) == 23 || abs(cms2.mus_mc_motherid()[index]) == 24) return true;
//   return false;
// }
//----------------------------------------------------------------
// Electron ID without isolation
//----------------------------------------------------------------
bool goodElectronWithoutIsolation(int index) {
  if ( cms2.els_tightId().at(index)     !=  1) return false;
  if ( cms2.els_closestMuon().at(index) != -1) return false;
  if ( abs(cms2.els_d0().at(index)) > 0.025)   return false;
  return true;
}
//----------------------------------------------------------------
// Electron ID without isolation or d0 cut
//----------------------------------------------------------------
bool goodElectronWithoutIsolationWithoutd0(int index) {
  if ( cms2.els_tightId().at(index)     !=  1) return false;
  if ( cms2.els_closestMuon().at(index) != -1) return false;
  return true;
}
//----------------------------------------------------------------
// loose Electron ID without isolation
//----------------------------------------------------------------
bool goodLooseElectronWithoutIsolation(int index) {
  if ( cms2.els_looseId().at(index)     !=  1) return false;
  if ( cms2.els_closestMuon().at(index) != -1) return false;
  if ( abs(cms2.els_d0().at(index)) > 0.025)   return false;
  return true;
}
//----------------------------------------------------------------
// Muon ID without isolation
//---------------------------------------------------------------
bool goodMuonWithoutIsolation(int index) {
  if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) > 5.) return false;
  if (abs(cms2.mus_d0().at(index))   > 0.25) return false;
  if (cms2.mus_validHits().at(index) < 7)    return false;
  return true;
}
//-----------------------------------------------------------
// Electron Isolation
//-----------------------------------------------------------
double el_rel_iso (int index, bool use_calo_iso)
{
     double sum = cms2.els_tkIso().at(index);
     if (use_calo_iso)
	  sum += cms2.els_ecalJuraIso()[index] + cms2.els_hcalConeIso()[index];
     double pt  = cms2.els_p4().at(index).pt();
     return pt / (pt + sum);
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
//-----------------------------------------------------------
// Muon Isolation
//-----------------------------------------------------------
double mu_rel_iso (int index)
{
     double sum =  cms2.mus_iso03_sumPt().at(index) +  
	  cms2.mus_iso03_emEt().at(index)  +
	  cms2.mus_iso03_hadEt().at(index);
     double pt  = cms2.mus_p4().at(index).pt();
     return pt / (pt+sum);

}
bool passMuonIsolation(int index) 
{
     return mu_rel_iso(index) > 0.92;
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
  if (!goodElectronWithoutIsolation(index)) return false;
  if (!passElectronIsolation(index, use_calo_iso))       return false;
  return true;
}
//--------------------------------------------
// "super-tight" Electron ID (to kill EM fakes)
//--------------------------------------------
bool supertightElectron (int index)
{
     if (fabs(cms2.els_p4()[index].eta()) > 1.479)
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
// Pass 2 MET selection
//--------------------------------------------
bool pass2Met (int i_hyp, const TVector3& corr) {
  // for e-e and mu-mu
  TVector3 hyp_met;
  hyp_met.SetPtEtaPhi(cms2.hyp_met()[i_hyp], 0, cms2.hyp_metPhi()[i_hyp]);
  hyp_met += corr;
  if (cms2.hyp_type()[i_hyp] == 0 || cms2.hyp_type()[i_hyp] == 3) {
    if (hyp_met.Pt() < 30) return false;
    //    if ( fabs(hyp_p4[i_hyp]->mass()-90.0)<10.0) return false;
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
  hyp_met.SetPtEtaPhi(cms2.hyp_met()[i_hyp], 0, cms2.hyp_metPhi()[i_hyp]);
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
  if ( sumEt < 550 ) return false; // good for ~1%? currently "optimizing"
  return true;
}

//-------------------------------------------------
// Auxiliary function to scan the doc line and 
// identify DY-> ee vs mm vs tt
//-------------------------------------------------
int getDrellYanType() {
  bool foundZ;
  int size = cms2.genps_id().size();
  for (int jj=0; jj<size; jj++) {
    if (cms2.genps_id().at(jj) == 23) {
      foundZ = true;
      if (jj+3 > size) {
	std::cout << 
	  "Found Z but not enough room in doc lines for leptons?" << std::endl;
        return 999;
      }
      if (abs(cms2.genps_id().at(jj+1)) == 11) return 0;  //DY->ee
      if (abs(cms2.genps_id().at(jj+1)) == 13) return 1;  //DY->mm
      if (abs(cms2.genps_id().at(jj+1)) == 15) return 2;  //DY->tautau
    }
  }
  std::cout << "Does not look like a DY event" << std::endl;
  return 999;
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

bool passMuonBVeto (int i_dilep, bool soft_nonisolated)
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
	  if (abs(cms2.hyp_lt_id()[i_dilep]) == 13)
	       mus_in_hyp++;
	  if (abs(cms2.hyp_ll_id()[i_dilep]) == 13)
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
	  if (abs(cms2.hyp_lt_id()[i_dilep]) == 13 &&
	      cms2.hyp_lt_index()[i_dilep] == int(i))
	       continue;
	  if (abs(cms2.hyp_ll_id()[i_dilep]) == 13 &&
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
    if ((abs(cms2.hyp_lt_id()[i_hyp]) == 11 && dRbetweenVectors(cms2.hyp_lt_p4()[i_hyp],cms2.trkjets_p4()[itrkjet]) < 0.4)||
	(abs(cms2.hyp_ll_id()[i_hyp]) == 11 && dRbetweenVectors(cms2.hyp_ll_p4()[i_hyp],cms2.trkjets_p4()[itrkjet]) < 0.4)
	) continue;
    TLorentzVector p(cms2.trkjets_p4()[itrkjet].Px(), cms2.trkjets_p4()[itrkjet].Py(), cms2.trkjets_p4()[itrkjet].Pz(), cms2.trkjets_p4()[itrkjet].E());
    if (p.Perp() < jetet) continue;
    if (fabs(p.Eta()) > jeteta) continue;
    trkjets.push_back(p);
  }

  return NjetVeto(trkjets, 15);
}

bool passTrkJetVeto(int i_hyp)
{
     return nTrkJets(i_hyp) == 0;
}

double reliso_lt (int i_hyp, bool use_calo_iso)
{
     // muons do it one way:
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 13) {
	  return mu_rel_iso(cms2.hyp_lt_index()[i_hyp]);
     }
     // electrons do it another way:
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11) {
	  return el_rel_iso(cms2.hyp_lt_index()[i_hyp], use_calo_iso);
     }
     // mysterions are not well handled:
     return 0;
}

double reliso_ll (int i_hyp, bool use_calo_iso)
{
     // muons do it one way:
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 13) {
	  return mu_rel_iso(cms2.hyp_ll_index()[i_hyp]);
     }
     // electrons do it another way:
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 11) {
	  return el_rel_iso(cms2.hyp_ll_index()[i_hyp], use_calo_iso);
     }
     // mysterions are not well handled:
     return 0;
}

bool trueMuonFromW(int index) {

  bool muIsFromW = false;

  if( abs(cms2.mus_mc_id()[index]) == 13 && abs(cms2.mus_mc_motherid()[index]) == 24 ) muIsFromW = true;

  return muIsFromW;
}

bool isFakeDenominatorElectron(int index) {
  //
  // returns true if input fulfills certain cuts
  //

  // cut definition
  float pt_cut        		= 15.;
  float eta_cut       		= 2.5;
  float hOverE_cut    		= 0.2;

  bool result = true;

  if ( cms2.els_p4()[index].Pt()  < pt_cut )            result = false;
  if ( std::abs(cms2.els_p4()[index].Eta()) > eta_cut ) result = false;
  if ( !passElectronIsolation(index) )          	result = false;
  if ( cms2.els_hOverE()[index]   > hOverE_cut )        result = false;

  return result;

}

bool isFakeNumeratorElectron(int index, int type=0) { 
  //
  // 1=loose, 2=tight
  //
  // returns true if input fulfills certain cuts
  //
  
  // cut definition
  float pt_cut        		= 15;
  float eta_cut       		= 2.5;

  bool result = true;

  if ( cms2.els_p4()[index].Pt()  < pt_cut )                 result = false;
  if ( std::abs(cms2.els_p4()[index].Eta()) > eta_cut )      result = false;
  if ( !passElectronIsolation(index) )          	result = false;
  if ( type == 1 ) {
    // loose
    if ( !goodLooseElectronWithoutIsolation(index) )   result = false;
  } else if ( type == 2 ) {
    // tight
    if ( !goodElectronWithoutIsolation(index) )   result = false;
  } else {
    cout << "WARNING: wrong electron type detected, please select loose (1) or tight (2)" << endl;
  }

  return result;
  
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
	  if (fabs(tk_cottheta - el_cottheta) < min_dcottheta) {
	       min_dcottheta = fabs(tk_cottheta - el_cottheta);
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
     return fabs(dphi);
}

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

//-------------------------------------------
// plain MET cut for reducing SM contribution 
// in WW selection to 10%
//-------------------------------------------
bool met10(int i_hyp, const TVector3& corr) {
  TVector3 hyp_met;
  hyp_met.SetPtEtaPhi(cms2.hyp_met()[i_hyp], 0, cms2.hyp_metPhi()[i_hyp]);
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
  hyp_met.SetPtEtaPhi(cms2.hyp_met()[i_hyp], 0, cms2.hyp_metPhi()[i_hyp]);
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
  // 3: apply both 1 and 2 (return 3)
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
  if ( abs(cms2.hyp_lt_id()[hyp_index]) == 13 ) 
    ltTrkIdx = cms2.mus_trkidx()[cms2.hyp_lt_index()[hyp_index]];
  else if ( abs(cms2.hyp_lt_id()[hyp_index]) == 11 )
    ltTrkIdx = cms2.els_trkidx()[cms2.hyp_lt_index()[hyp_index]];
  if ( abs(cms2.hyp_ll_id()[hyp_index]) == 13 ) 
    ltTrkIdx = cms2.mus_trkidx()[cms2.hyp_ll_index()[hyp_index]];
  else if ( abs(cms2.hyp_ll_id()[hyp_index]) == 11 )
    ltTrkIdx = cms2.els_trkidx()[cms2.hyp_ll_index()[hyp_index]];

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
	double dZ = fabs(cms2.hyp_lt_z0()[hyp_index] - cms2.trks_z0()[trk1]);
	if ( dZ < dZCut ) {
	  // require delta R > 0.1
	  double dR = ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[hyp_index], cms2.trks_trk_p4()[trk1]);
	  if ( dR > 0.1 ) {
	    // if inv. mass of candidate within z mass window, veto event
	    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > 
	      vec = cms2.hyp_lt_p4()[hyp_index] + cms2.trks_trk_p4()[trk1];
	    if ( inZmassWindow(vec.mass()) ) mode1 = true;
	  }
	}
      }
      // require opposite sign
      if ( (cms2.hyp_ll_charge()[hyp_index] * cms2.trks_charge()[trk1]) < 0 ) {
	// require delta z0 < 0.5 cm
	double dZ = fabs(cms2.hyp_ll_z0()[hyp_index] - cms2.trks_z0()[trk1]);
	if ( dZ < dZCut ) {
	  // require delta R > 0.1
	  double dR = ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[hyp_index], cms2.trks_trk_p4()[trk1]);
	  if ( dR > 0.1 ) {
	    // if inv. mass of candidate within z mass window, veto event
	    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > 
	      vec = cms2.hyp_ll_p4()[hyp_index] + cms2.trks_trk_p4()[trk1];
	    if ( inZmassWindow(vec.mass()) ) mode1 = true;
	  }
	}
      }
      // try to combine one track with another track from the collection of selected trk-isolated tracks
    } else if ( !mode2 ) {
      for ( unsigned int trk2 = trk1+1;
	    trk2 < isoTrks.size();
	    ++trk2 ) {
	// require opposite sign
	if ( (cms2.trks_charge()[trk2] * cms2.trks_charge()[trk1]) < 0 ) {
	  // require delta z0 < 0.5 cm
	  double dZ = fabs(cms2.trks_z0()[trk2] - cms2.trks_z0()[trk1]);
	  if ( dZ < dZCut ) {
	    // require delta R > 0.1
	    double dR = ROOT::Math::VectorUtil::DeltaR(cms2.trks_trk_p4()[trk2], cms2.trks_trk_p4()[trk1]);
	    if ( dR > 0.1 ) {
	      // if inv. mass of candidate within z mass window, veto event
	      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > 
		vec = cms2.trks_trk_p4()[trk2] + cms2.trks_trk_p4()[trk1];
	      if ( inZmassWindow(vec.mass()) ) mode2 = true;
	    }
	  }
	}
      }
    }
  }

  if ( mode1 || mode2 ) return 3;
  if ( mode1 ) return 1;
  if ( mode2 ) return 2;
  
  return -1;
}
