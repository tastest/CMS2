 //===========================================================
//
// Various selection functions are kept here
//
//============================================================
#include "Math/LorentzVector.h"
#include "TMath.h"
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
  if ( els_tightId()[index]     !=  1) return false;
  if ( els_closestMuon()[index] != -1) return false;
  if ( abs(els_d0()[index]) > 0.025)   return false;
  return true;
}
//----------------------------------------------------------------
// Muon ID without isolation
//---------------------------------------------------------------
bool goodMuonWithoutIsolation(int index) {
  if (mus_gfit_chi2()[index]/mus_gfit_ndof()[index] > 5.) return false;
  if (abs(mus_d0()[index])   > 0.25) return false;
  if (mus_validHits()[index] < 7)    return false;
  return true;
}
//-----------------------------------------------------------
// Electron Isolation
//-----------------------------------------------------------
bool passElectronIsolation(int index) {
  double sum = els_tkIso()[index];
  double pt  = els_p4()[index].pt();
   if ( pt/(pt+sum) < 0.92) return false;
  return true;  
} 
//-----------------------------------------------------------
// Muon Isolation
//-----------------------------------------------------------
bool passMuonIsolation(int index) {
  double sum =  mus_iso03_sumPt()[index] +  
                mus_iso03_emEt()[index]  +
                mus_iso03_hadEt()[index];
  double pt  = mus_p4()[index].pt(); 
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

//-------------------------------------------------
// Auxiliary function to scan the doc line and 
// identify DY-> ee vs mm vs tt
//-------------------------------------------------
int getDrellYanType() {
  bool foundZ;
  int size = genps_id().size();
  for (int jj=0; jj<size; jj++) {
    if (genps_id()[jj] == 23) {
      foundZ = true;
      if (jj+3 > size) {
	std::cout << 
         "Found Z but not enough room in doc lines for leptons?" << std::endl;
	return 999;
      }
      if (abs(genps_id()[jj+1]) == 11) return 0;  //DY->ee
      if (abs(genps_id()[jj+1]) == 13) return 1;  //DY->mm
      if (abs(genps_id()[jj+1]) == 15) return 2;  //DY->tautau
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
//
// Because isolation is not in the els block, for Z->ee veto we are 
// forced to use the hyp_block for one of them (if we want to put 
// an isolation requirement).  This would be a problem only for 
// Z->ee + 2 other hyp leptons (unlikely)
// We also assume that the leptons in hyp block have already passed isolation
//---------------------------------------------------------------------
bool additionalZveto() {

  // true if we want to veto this event
  bool veto=false;

  // first, look for Z->mumu
  for (unsigned int i=0; i<mus_p4().size(); i++) {
    if (mus_p4()[i].pt() < 20.)     continue;
    if (!goodMuonWithoutIsolation(i)) continue;

    for (unsigned int j=i+1; j<mus_p4().size(); j++) {
      if (mus_p4()[j].pt() < 20.) continue;
      if (!goodMuonWithoutIsolation(j)) continue;
      if (mus_charge()[i] == mus_charge()[j]) continue;

      // At least one of them has to pass isolation
      if (!passMuonIsolation(i) && !passMuonIsolation(j)) continue;

      // Make the invariant mass
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > 
                vec = mus_p4()[i] + mus_p4()[j];
      if ( inZmassWindow(vec.mass()) ) return true;

    }
  }

  // Now look at Z->ee.  
  for (unsigned int i=0; i<evt_nels(); i++) {
    if (els_p4()[i].pt() < 20.)     continue;
    if (!goodElectronWithoutIsolation(i)) continue;

    for (unsigned int j=i+1; j<evt_nels(); j++) {
      if (els_p4()[j].pt() < 20.) continue;
      if (!goodElectronWithoutIsolation(j)) continue;
      if (els_charge()[i] == els_charge()[j]) continue;

      // At least one of them has to pass isolation
      if (!passElectronIsolation(i) && !passElectronIsolation(j)) continue;

      // Make the invariant mass
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > 
                vec = els_p4()[i] + els_p4()[j];
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
  int size = genps_id().size();
  // Initialize particle database
  TDatabasePDG *pdg = new TDatabasePDG();
  std::cout << "------------------------------------------" << std::endl;
  for (int j=0; j<size; j++) {
    cout << setw(9) << left << pdg->GetParticle(genps_id()[j])->GetName() << " "
         << setw(7) << right << setprecision(4) << genps_p4()[j].pt() << "  "
         << setw(7) << right << setprecision(4) << genps_p4()[j].phi() << "  "
         << setw(10) << right << setprecision(4) << genps_p4()[j].eta() << "  "
         << setw(10) << right << setprecision(4) << genps_p4()[j].mass() << endl;
  }
  // Clean up
  delete pdg;
}

void correctMETmuons_crossedE(double& met, double& metPhi, 
			      double muon_pt, double muon_phi,
			      double muon_track_theta, double muon_track_phi,
			      double mu_crossedem_dep, double mu_crossedhad_dep, double mu_crossedho_dep ) {

  // first, account for muon momentum
  double metx =  met*cos(metPhi);
  double mety =  met*sin(metPhi);
  double pt0  =  muon_pt; 
  double phi0 =  muon_phi; 
  metx -= pt0*cos(phi0);
  mety -= pt0*sin(phi0);
  
  
  met = sqrt(metx*metx+mety*mety);
  metPhi = atan2(mety, metx);
   
   double muEx = 0.0;
   double muEy = 0.0;
   
   
   // use muon position at the outer most state of the silicon track if 
   // TrackExtra is available and momentum direction at the origin 
   // otherwise. Both should be fine.
   // NOTICE: MET is built out of towers, which are 5x5 ECAL crystals + one 
   // element of HCAL and HO. Muon energy is reported for individual crossed 
   // elements of all the detectors and 3x3 elements of each time as an 
   // alternative way of energy calculation.
   double theta = muon_track_theta;
   double phi   = muon_track_phi;
	 
   muEx += ( mu_crossedem_dep + mu_crossedhad_dep + mu_crossedho_dep )*sin(theta)*cos( phi );
   muEy += ( mu_crossedem_dep + mu_crossedhad_dep + mu_crossedho_dep )*sin(theta)*sin( phi );
   
   
   metx = met*cos(metPhi) + muEx;
   mety = met*sin(metPhi) + muEy;
   met   = sqrt(metx*metx + mety*mety);
   metPhi = atan2(mety, metx);
}

