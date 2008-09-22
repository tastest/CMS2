//===========================================================
//
// Various selection functions are kept here
//
//============================================================
#include "Math/LorentzVector.h"
#include "TMath.h"
#include "TLorentzVector.h"
// #include "selections.h"
//#include <vector>
//#include "CMS1.h"
#include "TDatabasePDG.h"
#include "../Tools/matchTools.h"

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
bool passElectronIsolation(int index) {
  double sum = cms2.els_tkIso().at(index);
  double pt  = cms2.els_p4().at(index).pt();
   if ( pt/(pt+sum) < 0.92) return false;
  return true;  
} 
//-----------------------------------------------------------
// Muon Isolation
//-----------------------------------------------------------
bool passMuonIsolation(int index) {
  double sum =  cms2.mus_iso03_sumPt().at(index) +  
                cms2.mus_iso03_emEt().at(index)  +
                cms2.mus_iso03_hadEt().at(index);
  double pt  = cms2.mus_p4().at(index).pt();
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
//--------------------------------------------
// Pass 2 MET selection
//--------------------------------------------
bool pass2Met (int i_hyp) {
  // for e-e and mu-mu
  if (cms2.hyp_type()[i_hyp] == 0 || cms2.hyp_type()[i_hyp] == 3) {
    if (cms2.hyp_met()[i_hyp] < 30) return false;
    //    if ( fabs(hyp_p4[i_hyp]->mass()-90.0)<10.0) return false;
    if( cms2.hyp_met()[i_hyp]/cms2.hyp_p4()[i_hyp].pt()<0.6 && 
	acos(cos(cms2.hyp_metPhi()[i_hyp]-cms2.hyp_p4()[i_hyp].phi()-3.1416))<0.25 ) return false;
  }
  // for e-mu and mu-e
  if (cms2.hyp_type()[i_hyp] == 1 || cms2.hyp_type()[i_hyp] == 2) {
    if (cms2.hyp_met()[i_hyp] < 20) return false;
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
bool pass4Met(int i_hyp) {
  double metspec = MetSpecial(cms2.hyp_met()[i_hyp], cms2.hyp_metPhi()[i_hyp], i_hyp);
  if (cms2.hyp_type()[i_hyp] == 0 || cms2.hyp_type()[i_hyp] == 3) {
    if ( metspec < 20 ) return false;
    //if ( metspec < 20 && hyp_p4->mass() < 90 ) return false;
    if ( cms2.hyp_met()[i_hyp] < 45 ) return false;
  }
  else if (cms2.hyp_type()[i_hyp] == 1 || cms2.hyp_type()[i_hyp] == 2) {
    //if ( metspec < 20 && hyp_p4->mass() < 90 ) return false;
    if ( metspec < 20 ) return false;
  }
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
      if (!passElectronIsolation(i) && !passElectronIsolation(j)) continue;

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


bool passTrkJetVeto(int i_hyp){
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

  return (NjetVeto(trkjets, 15) == 0);
}

