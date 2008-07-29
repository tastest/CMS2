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
// Electron ID without isolation
//----------------------------------------------------------------
bool goodElectronWithoutIsolation(int index) {
  if ( els_tightId.at(index)     !=  1) return false;
  if ( els_closestMuon.at(index) != -1) return false;
  if ( abs(els_d0.at(index)) > 0.025)   return false;
  return true;
}
//----------------------------------------------------------------
// Muon ID without isolation
//---------------------------------------------------------------
bool goodMuonWithoutIsolation(int index) {
  if (mus_gfit_chi2.at(index)/mus_gfit_ndof.at(index) > 5.) return false;
  if (abs(mus_d0.at(index))   > 0.25) return false;
  if (mus_validHits.at(index) < 7)    return false;
  return true;
}
//-----------------------------------------------------------
// Electron Isolation
//-----------------------------------------------------------
bool passElectronIsolation(int index) {
  double sum = els_tkIso.at(index);
  double pt  = els_p4.at(index).pt();
   if ( pt/(pt+sum) < 0.92) return false;
  return true;  
} 
//-----------------------------------------------------------
// Muon Isolation
//-----------------------------------------------------------
bool passMuonIsolation(int index) {
  double sum =  mus_iso03_sumPt.at(index) +  
                mus_iso03_emEt.at(index)  +
                mus_iso03_hadEt.at(index);
  double pt  = mus_p4.at(index).pt();
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
  if (hyp_type[i_hyp] == 0 || hyp_type[i_hyp] == 3) {
    if (hyp_met[i_hyp] < 30) return false;
    //    if ( fabs(hyp_p4[i_hyp]->mass()-90.0)<10.0) return false;
    if( hyp_met[i_hyp]/hyp_p4[i_hyp].pt()<0.6 && 
	acos(cos(hyp_metPhi[i_hyp]-hyp_p4[i_hyp].phi()-3.1416))<0.25 ) return false;
  }
  // for e-mu and mu-e
  if (hyp_type[i_hyp] == 1 || hyp_type[i_hyp] == 2) {
    if (hyp_met[i_hyp] < 20) return false;
  }
  return true;
}


#if 0
//-------------------------------------------------
// Auxiliary function to scan the doc line and 
// identify DY-> ee vs mm vs tt
//-------------------------------------------------
int getDrellYanType() {
  bool foundZ;
  int size = genps_id.size();
  for (int jj=0; jj<size; jj++) {
    if (genps_id.at(jj) == 23) {
      foundZ = true;
      if (jj+3 > size) {
	std::cout << 
         "Found Z but not enough room in doc lines for leptons?" << std::endl;
	return 999;
      }
      if (abs(genps_id.at(jj+1)) == 11) return 0;  //DY->ee
      if (abs(genps_id.at(jj+1)) == 13) return 1;  //DY->mm
      if (abs(genps_id.at(jj+1)) == 15) return 2;  //DY->tautau
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
#endif
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
  for (int i=0; i < mus_p4.size(); i++) {
    if (mus_p4.at(i).pt() < 20.)     continue;
    if (!goodMuonWithoutIsolation(i)) continue;

    for (int j=i+1; j < mus_p4.size(); j++) {
      if (mus_p4.at(j).pt() < 20.) continue;
      if (!goodMuonWithoutIsolation(j)) continue;
      if (mus_charge.at(i) == mus_charge.at(j)) continue;

      // At least one of them has to pass isolation
      if (!passMuonIsolation(i) && !passMuonIsolation(j)) continue;

      // Make the invariant mass
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > 
                vec = mus_p4.at(i) + mus_p4.at(j);
      if ( inZmassWindow(vec.mass()) ) return true;

    }
  }

  // now, look for Z->ee
  for (int i=0; i < els_p4.size(); i++) {
    if (els_p4.at(i).pt() < 20.)     continue;
    if (!goodElectronWithoutIsolation(i)) continue;

    for (int j=i+1; j<els_p4.size(); j++) {
      if (els_p4.at(j).pt() < 20.) continue;
      if (!goodElectronWithoutIsolation(j)) continue;
      if (els_charge.at(i) == els_charge.at(j)) continue;

      // At least one of them has to pass isolation
      if (!passElectronIsolation(i) && !passElectronIsolation(j)) continue;

      // Make the invariant mass
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > 
                vec = els_p4.at(i) + els_p4.at(j);
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
  int size = genps_id.size();
  // Initialize particle database
  TDatabasePDG *pdg = new TDatabasePDG();
  std::cout << "------------------------------------------" << std::endl;
  for (int j=0; j<size; j++) {
    cout << setw(9) << left << pdg->GetParticle(genps_id.at(j))->GetName() << " "
         << setw(7) << right << setprecision(4) << genps_p4.at(j).pt() << "  "
         << setw(7) << right << setprecision(4) << genps_p4.at(j).phi() << "  "
         << setw(10) << right << setprecision(4) << genps_p4.at(j).eta() << "  "
         << setw(10) << right << setprecision(4) << genps_p4.at(j).mass() << endl;
  }
  // Clean up
  delete pdg;
}

