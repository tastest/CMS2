//===========================================================
//
// Various selection functions are kept here
//
//============================================================
#include "Math/LorentzVector.h"
#include "TMath.h"
//#include "CMS2_SW2.h"
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
// (Changed Impact parameter cut to 400 micron for startup)
//----------------------------------------------------------------
bool goodElectronWithoutIsolation(int index) {
  if ( els_tightId().at(index)     !=  1) return false;
  if ( els_closestMuon().at(index) != -1) return false;
  if ( fabs(els_d0corr().at(index)) > 0.040)   return false;
  return true;
}
//----------------------------------------------------------------
// Muon ID without isolation
//---------------------------------------------------------------
bool goodMuonWithoutIsolation(int index) {
  if (mus_gfit_chi2().at(index)/mus_gfit_ndof().at(index) > 5.) return false;
  if (fabs(mus_d0corr().at(index))   > 0.25) return false;
  if (mus_validHits().at(index) < 7)    return false;
  return true;
}
//-----------------------------------------------------------
// Electron Isolation
//-----------------------------------------------------------
double electronIsolation(int index){
  double sum = els_tkIso().at(index);
  double pt  = els_p4().at(index).pt();
  return pt/(pt+sum);
}
double electronTrkIsolation(int index){
  double sum = els_tkIso().at(index);
  double pt  = els_p4().at(index).pt();
  return pt/(pt+sum);
}
double electronCalIsolation(int index){
  double sum = els_pat_caloIso().at(index);
  double pt  = els_p4().at(index).pt();
  return pt/(pt+sum);
}

bool passElectronIsolation(int index) {
  return ( electronIsolation(index) >= 0.92 );
} 
bool passElectronAntiIsolation(int index) {
  double iso = electronIsolation(index);
  return ( iso > 0.6 && iso < 0.92 );
}

//-----------------------------------------------------------
// Muon Isolation
//-----------------------------------------------------------
double muonIsolation(int index){
  double sum =  mus_iso03_sumPt().at(index) +  
    mus_iso03_emEt().at(index)  +
    mus_iso03_hadEt().at(index);
  double pt  = mus_p4().at(index).pt();
  return  pt/(pt+sum);
}
double muonTrkIsolation(int index){
  double sum =  mus_iso03_sumPt().at(index);
  double pt  = mus_p4().at(index).pt();
  return  pt/(pt+sum);
}
double muonCalIsolation(int index){
  double sum =  mus_iso03_emEt().at(index)  +
    mus_iso03_hadEt().at(index);
  double pt  = mus_p4().at(index).pt();
  return  pt/(pt+sum);
}

bool passMuonIsolation(int index) {
  return ( muonIsolation(index) >= 0.92 );
}
bool passMuonAntiIsolation(int index) {
  double iso = muonIsolation(index);
  return ( iso > 0.6 && iso < 0.92 );
}
//---------------------------------------
// Anti-isolation test: 0 - both isolated; 
//                      1 - one (e in emu) is anti-isolated
//                      2 - one (mu in emu) is anti-isolated
//                      3 - both are  anti-isolated
//--------------------------------------
bool passDilAntiIsolation(int isoLooseMode, int hypIdx){
  assert (isoLooseMode >= 0 && isoLooseMode < 4);

  if (isoLooseMode == 0){
    // Muon quality cuts, including isolation
    if (abs(hyp_lt_id().at(hypIdx)) == 13 && !passMuonIsolation(hyp_lt_index().at(hypIdx)) ) return false;
    if (abs(hyp_ll_id().at(hypIdx)) == 13 && !passMuonIsolation(hyp_ll_index().at(hypIdx)) ) return false;
    
    // Electron quality cuts, including isolation 
    if (abs(hyp_lt_id().at(hypIdx)) == 11 && !passElectronIsolation(hyp_lt_index().at(hypIdx)) ) return false;
    if (abs(hyp_ll_id().at(hypIdx)) == 11 && !passElectronIsolation(hyp_ll_index().at(hypIdx)) ) return false;    
  } else {
    if (isoLooseMode == 3){
      // Muon quality cuts, including isolation
      if (abs(hyp_lt_id().at(hypIdx)) == 13 && !passMuonAntiIsolation(hyp_lt_index().at(hypIdx)) ) return false;
      if (abs(hyp_ll_id().at(hypIdx)) == 13 && !passMuonAntiIsolation(hyp_ll_index().at(hypIdx)) ) return false;
      
      // Electron quality cuts, including isolation 
      if (abs(hyp_lt_id().at(hypIdx)) == 11 && !passElectronAntiIsolation(hyp_lt_index().at(hypIdx)) ) return false;
      if (abs(hyp_ll_id().at(hypIdx)) == 11 && !passElectronAntiIsolation(hyp_ll_index().at(hypIdx)) ) return false;	  
    } else if (abs(hyp_lt_id().at(hypIdx))!= abs(hyp_ll_id().at(hypIdx))){
      if (isoLooseMode == 2){
	//mu is anti-isolated for emu mode
	// Muon quality cuts, including isolation
	if (abs(hyp_lt_id().at(hypIdx)) == 13 && !passMuonAntiIsolation(hyp_lt_index().at(hypIdx)) ) return false;
	if (abs(hyp_ll_id().at(hypIdx)) == 13 && !passMuonAntiIsolation(hyp_ll_index().at(hypIdx)) ) return false;
	
	// Electron quality cuts, including isolation 
	if (abs(hyp_lt_id().at(hypIdx)) == 11 && !passElectronIsolation(hyp_lt_index().at(hypIdx)) ) return false;
	if (abs(hyp_ll_id().at(hypIdx)) == 11 && !passElectronIsolation(hyp_ll_index().at(hypIdx)) ) return false;	  	  
      } else { //isoLooseMode == 1
	//e is anti-isolated for emu mode
	// Muon quality cuts, including isolation
	if (abs(hyp_lt_id().at(hypIdx)) == 13 && !passMuonIsolation(hyp_lt_index().at(hypIdx)) ) return false;
	if (abs(hyp_ll_id().at(hypIdx)) == 13 && !passMuonIsolation(hyp_ll_index().at(hypIdx)) ) return false;
	
	// Electron quality cuts, including isolation 
	if (abs(hyp_lt_id().at(hypIdx)) == 11 && !passElectronAntiIsolation(hyp_lt_index().at(hypIdx)) ) return false;
	if (abs(hyp_ll_id().at(hypIdx)) == 11 && !passElectronAntiIsolation(hyp_ll_index().at(hypIdx)) ) return false;
      }
    } else {
      //one-only is isolated
      if ( abs(hyp_lt_id().at(hypIdx)) == 11 && 
	   !(
	     ( passElectronAntiIsolation(hyp_lt_index().at(hypIdx)) && passElectronIsolation(hyp_ll_index().at(hypIdx)) )
	     ||
	     ( passElectronAntiIsolation(hyp_ll_index().at(hypIdx)) && passElectronIsolation(hyp_lt_index().at(hypIdx)) )
	     )
	   ) return false;
      if ( abs(hyp_lt_id().at(hypIdx)) == 13 && 
	   !(
	     ( passMuonAntiIsolation(hyp_lt_index().at(hypIdx)) && passMuonIsolation(hyp_ll_index().at(hypIdx)) )
	     ||
	     ( passMuonAntiIsolation(hyp_ll_index().at(hypIdx)) && passMuonIsolation(hyp_lt_index().at(hypIdx)) )
	     )
	   ) return false;
    }
    
  }
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
bool pass2Met(int hypIdx) {
  // for e-e and mu-mu
  if (hyp_type().at(hypIdx) == 0 || hyp_type().at(hypIdx) == 3) {
    if (hyp_met().at(hypIdx) < 30) return false;
    
    if( hyp_met().at(hypIdx)/hyp_p4().at(hypIdx).pt()<0.6 && 
	acos(cos(hyp_metPhi().at(hypIdx)-hyp_p4().at(hypIdx).phi()-3.1416))<0.25 ) return false;
  }
  // for e-mu and mu-e
  if (hyp_type().at(hypIdx) == 1 || hyp_type().at(hypIdx) == 2) {
    if (hyp_met().at(hypIdx) < 20) return false;
  }
  return true;
}

// event-level pat-met: emu met >20, mm,em met>30
bool passPatMet_OF20_SF30(int hypIdx){
  if  (hyp_type().at(hypIdx) == 0 || hyp_type().at(hypIdx) == 3) {
    if (met_pat_metCor() < 30) return false;
  }
  
  if (hyp_type().at(hypIdx) == 1 || hyp_type().at(hypIdx) == 2) {
    if (met_pat_metCor() < 20) return false;
  }
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
    if (genps_id().at(jj) == 23) {
      foundZ = true;
      if (jj+3 > size) {
	std::cout << 
	  "Found Z but not enough room in doc lines for leptons?" << std::endl;
	return 999;
      }
      if (abs(genps_id().at(jj+1)) == 11) return 0;  //DY->ee
      if (abs(genps_id().at(jj+1)) == 13) return 1;  //DY->mm
      if (abs(genps_id().at(jj+1)) == 15) return 2;  //DY->tautau
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

int genpCountPDGId(int id0, int id1=-1, int id2=-1){
  int count = 0;
  int size = genps_id().size();
  for (int jj=0; jj<size; jj++) {
    if (abs(genps_id().at(jj)) == id0) count++;
    if (abs(genps_id().at(jj)) == id1) count++;
    if (abs(genps_id().at(jj)) == id2) count++;
  }
  return count;
}

int genpDileptonType(){
  //0 mumu; 1 emu; 2 ee
  
  uint nmus = 0;
  uint nels = 0;
  int size = genps_id().size();
  for (int jj=0; jj<size; jj++) {
    if (abs(genps_id().at(jj)) == 11) nels++;
    if (abs(genps_id().at(jj)) == 13) nmus++;
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
    cout << setw(9) << left << pdg->GetParticle(genps_id().at(j))->GetName() << " "
         << setw(7) << right << setprecision(4) << genps_p4().at(j).pt() << "  "
         << setw(7) << right << setprecision(4) << genps_p4().at(j).phi() << "  "
         << setw(10) << right << setprecision(4) << genps_p4().at(j).eta() << "  "
         << setw(10) << right << setprecision(4) << genps_p4().at(j).mass() << endl;
  }
  // Clean up
  delete pdg;
}
//--------------------------------------------------
// Returns the number of e,mu, and tau in the doc lines
//-----------------------------------------------------
void leptonCount(int& nele, int& nmuon, int& ntau) {
  nele=0;
  nmuon=0;
  ntau=0;
  int size = genps_id().size();
  for (int jj=0; jj<size; jj++) {
    if (abs(genps_id().at(jj)) == 11) nele++;
    if (abs(genps_id().at(jj)) == 13) nmuon++;
    if (abs(genps_id().at(jj)) == 15) ntau++;
  }
}

// true if there is a muon not in the hypothesis
bool haveExtraMuon(int hypIdx){
  bool result = false;

  int nHMus = 0;
  if (abs(hyp_lt_id().at(hypIdx))==13){
    nHMus++;
  }
  if (abs(hyp_ll_id().at(hypIdx))==13){
    nHMus++;
  }

  int nEvtMus = mus_p4().size();
  result = (nEvtMus - nHMus) > 0;

  return result;
}

// true if there is a muon not in the hypothesis
bool haveExtraMuon5(int hypIdx){
  double minPtCut = 5;
  bool result = false;

  int nHMus = 0;
  if (abs(hyp_lt_id().at(hypIdx))==13){
    if (hyp_lt_p4().at(hypIdx).pt() > minPtCut ) nHMus++;
  }
  if (abs(hyp_ll_id().at(hypIdx))==13){
    if (hyp_ll_p4().at(hypIdx).pt() > minPtCut) nHMus++;
  }

  int nEvtMus = 0;
  for (uint iMu = 0; iMu < mus_p4().size(); ++iMu){
    if (mus_p4().at(iMu).pt() > minPtCut) nEvtMus++;
  }
  result = (nEvtMus - nHMus) > 0;

  return result;
}


//-----------------------------------------------------------------------------------------------
//New selections for the common TTDil working group
//-----------------------------------------------------------------------------------------------
//
// loose lepton definitions 
//

bool electron20Eta2p4(int index){
  if (els_p4().at(index).pt() < 20 )   return false;
  if (fabs(els_p4().at(index).eta()) > 2.4 ) return false;
  
  return true;
}


bool looseElectronSelectionNoIsoTTDil08(int index) {
  if ( ! electron20Eta2p4(index) ) return false;
  if ( els_looseId().at(index)     !=  1) return false;
  if ( fabs(els_d0corr().at(index)) > 0.040)   return false;

  return true;
}


bool looseElectronSelectionTTDil08(int index) {
  if ( ! looseElectronSelectionNoIsoTTDil08(index) ) return false;

  if ( electronTrkIsolation(index) < 0.5 ) return false;
  if ( electronCalIsolation(index) < 0.5 ) return false;

  return true;
}

bool passElectronIsolationTTDil08(int index){
  if ( electronTrkIsolation(index) < 0.9 ) return false;
  if ( electronCalIsolation(index) < 0.8 ) return false;

  return true;
}

bool muon20Eta2p4(int index){
  if (mus_p4().at(index).pt() < 20 ) return false;
  if (fabs(mus_p4().at(index).eta()) >2.4 ) return false;

  return true;
}

bool looseMuonSelectionNoIsoTTDil08(int index) {

  if (! muon20Eta2p4(index) ) return false;

  if(!(2 & mus_type().at(index))) return false;
  if (mus_gfit_chi2().at(index)/mus_gfit_ndof().at(index) > 10.) return false;
  //  if (fabs(mus_d0corr().at(index))   > 0.25) return false;
  if (mus_validHits().at(index) < 11)    return false;
  
  return true;
}

bool looseMuonSelectionTTDil08(int index) {
  if (! looseMuonSelectionNoIsoTTDil08(index) ) return false;
  
  if ( muonTrkIsolation(index) < 0.5 ) return false;
  if ( muonCalIsolation(index) < 0.5 ) return false;

  return true;
}

bool passMuonIsolationTTDil08(int index) {
  if ( muonTrkIsolation(index) < 0.9 ) return false;
  if ( muonCalIsolation(index) < 0.9 ) return false;

  return true;
}

//------------------------------------------------------------------------------------

//--------------------------------------------------------------------
// Veto events if there are two leptons in the 
// event that make the Z mass.  This uses the mus and els
// blocks, ie, it is a veto that can use the 3rd (4th,5th,..)
// lepton in the event.
//
// Both leptons must be 20 GeV, and pass the same cuts as 
// the hypothesis leptons, except that one of them can be non-isolated
//---------------------------------------------------------------------
bool additionalZveto(bool useTTDil08 = false) {

  // true if we want to veto this event
  bool veto=false;

  // first, look for Z->mumu
  for (uint i=0; i<mus_p4().size(); i++) {
    if (mus_p4().at(i).pt() < 20.)     continue;
    if (useTTDil08){
      if (!looseMuonSelectionNoIsoTTDil08(i)) continue;
    } else {
      if (!goodMuonWithoutIsolation(i)) continue;
    }

    for (uint j=i+1; j<mus_p4().size(); j++) {
      if (mus_p4().at(j).pt() < 20.) continue;
      if (useTTDil08){
	if (!looseMuonSelectionNoIsoTTDil08(j)) continue;
      } else {
	if (!goodMuonWithoutIsolation(j)) continue;
      }
      if (mus_charge().at(i) == mus_charge().at(j)) continue;

      // At least one of them has to pass isolation
      if (useTTDil08){
	if (!passMuonIsolationTTDil08(i) && !passMuonIsolationTTDil08(j)) continue;
      } else {
	if (!passMuonIsolation(i) && !passMuonIsolation(j)) continue;
      }

      // Make the invariant mass
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > 
	vec = mus_p4().at(i) + mus_p4().at(j);
      if ( inZmassWindow(vec.mass()) ) return true;

    }
  }

  // now, look for Z->ee
  for (uint i=0; i<evt_nels(); i++) {
    if (els_p4().at(i).pt() < 20.)     continue;
    if (useTTDil08){
      if (!looseElectronSelectionNoIsoTTDil08(i)) continue;
    } else {
      if (!goodElectronWithoutIsolation(i)) continue;
    }

    for (uint j=i+1; j<evt_nels(); j++) {
      if (els_p4().at(j).pt() < 20.) continue;
      if (useTTDil08){
	if (!looseElectronSelectionNoIsoTTDil08(j)) continue;
      } else {
	if (!goodElectronWithoutIsolation(j)) continue;
      }
      if (els_charge().at(i) == els_charge().at(j)) continue;

      // At least one of them has to pass isolation
      if (!passElectronIsolationTTDil08(i) && !passElectronIsolationTTDil08(j)) continue;

      // Make the invariant mass
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > 
	vec = els_p4().at(i) + els_p4().at(j);
      if ( inZmassWindow(vec.mass()) ) return true;

    }
  }
  // done
  return veto;
}


