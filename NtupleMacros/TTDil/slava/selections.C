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
double electronTrkIsolationPAT(int index){
  double sum = els_pat_trackIso().at(index);
  double pt  = els_p4().at(index).pt();
  return pt/(pt+sum);
}
double electronCalIsolationPAT(int index){
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
double muonTrkIsolationPAT(int index){
  double sum =  mus_pat_trackIso().at(index);
  double pt  = mus_p4().at(index).pt();
  return  pt/(pt+sum);
}
double muonCalIsolationPAT(int index){
  double sum =  mus_pat_caloIso().at(index);
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
bool passPatMet_OF20_SF30(float metx, float mety, int hypIdx){
  float mymet = sqrt(metx*metx + mety*mety);
  if  (hyp_type().at(hypIdx) == 0 || hyp_type().at(hypIdx) == 3) {
    if (mymet < 30) return false;
  }
  
  if (hyp_type().at(hypIdx) == 1 || hyp_type().at(hypIdx) == 2) {
    if (mymet < 20) return false;
  }
  return true;
}
// event-level pat-met: emu met >20, mm,em met>30
bool passPatMet_OF20_SF30(int hypIdx){
  return passPatMet_OF20_SF30(met_pat_metCor()*cos(met_pat_metPhiCor()), 
			      met_pat_metCor()*sin(met_pat_metPhiCor()),
			      hypIdx);
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
  if ( els_closestMuon().at(index) != -1) return false; 

  return true;
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


// Trigger-related selections

//Bits passing for hyp type
bool passTriggersMu9orLisoE15(int dilType) {
  bool hltLooseIsoEle15_LW_L1R = ((evt_HLT2() & (1<<(47-32))) != 0);
  bool hltMu9 = ((evt_HLT3() & (1<<(82-64))) != 0);
  
  if (dilType == 0 && ! (hltMu9) ) return false;
  if ((dilType == 1 || dilType == 2) && ! (hltMu9 || hltLooseIsoEle15_LW_L1R)) return false;
  if (dilType == 3 && ! hltLooseIsoEle15_LW_L1R) return false;     

  return true;
}


bool passTriggersTTDil08JanTrial(int dilType) {
  bool hltIsoEle18_L1R = ((evt_HLT2() & (1<<(45-32))) != 0);
  bool hltDoubleIsoEle12_L1R = ((evt_HLT2() & (1<<(54-32))) != 0); 
  bool hltMu15_L1Mu7 = ((evt_HLT3() & (1<<(86-64))) != 0); 
  bool hltDoubleMu3 = ((evt_HLT3() & (1<<(90-64))) != 0);
  bool hltIsoEle10_Mu10_L1R = ((evt_HLT4() & (1<<(126-96))) != 0);
  
  if (dilType == 0 && ! (hltMu15_L1Mu7 || hltDoubleMu3) ) return false;
  if ((dilType == 1 || dilType == 2) 
      && ! (hltIsoEle18_L1R || hltMu15_L1Mu7 || hltIsoEle10_Mu10_L1R)) return false;
  if (dilType == 3 && ! (hltIsoEle18_L1R || hltDoubleIsoEle12_L1R) ) return false; 
  return true;
}


int eventDilIndexByWeightTTDil08(const std::vector<unsigned int>& goodHyps, int& strasbourgDilType, bool printDebug = false){
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
    unsigned int i_lt = hyp_lt_index().at(hypIdx);
    unsigned int i_ll = hyp_ll_index().at(hypIdx);

    int id_lt = hyp_lt_id().at(hypIdx);
    int id_ll = hyp_ll_id().at(hypIdx);

    float isoTk_lt = leptonTrkIsolationTTDil08(id_lt, i_lt);
    float isoTk_ll = leptonTrkIsolationTTDil08(id_ll, i_ll);

    float isoCal_lt = leptonCalIsolationTTDil08(id_lt, i_lt);
    float isoCal_ll = leptonCalIsolationTTDil08(id_ll, i_ll);

    //ad-hoc selection of weights
    if (abs(id_lt) == 11){
      //I want to select "trk & cal"-isolated ones
      hypWeight_iso += (isoTk_lt*isoCal_lt - 0.25); //shift by 0.25 to be positive-definite
      if (els_tightId().at(i_lt)) hypWeight_lt += 0.2;
    }
    if (abs(id_lt) == 13){
      //I want to select "trk & cal"-isolated ones	    
      hypWeight_iso += (isoTk_lt*isoCal_lt - 0.25);//shift by 0.25 to be positive-definite
      hypWeight_lt += 0.4;
    }
    if (abs(id_ll) == 11){
      //I want to select "trk & cal"-isolated ones
      hypWeight_iso *= (isoTk_ll*isoCal_ll - 0.25); //shift by 0.25 to be positive-definite
      if (els_tightId().at(i_ll)) hypWeight_ll += 0.2;
    }
    if (abs(id_ll) == 13){
      //I want to select "trk & cal"-isolated ones
      hypWeight_iso *= (isoTk_ll*isoCal_ll - 0.25); //shift by 0.25 to be positive-definite
      hypWeight_ll += 0.4;
    }
    float pt_lt = hyp_lt_p4().at(hypIdx).pt();
    float pt_ll = hyp_ll_p4().at(hypIdx).pt();
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
    for (unsigned int iEl =0; iEl < els_p4().size(); ++iEl){
      if (looseElectronSelectionTTDil08(iEl)){
	looseEls.push_back(iEl);
      }
    }
    for (unsigned int iMu =0; iMu < mus_p4().size(); ++iMu){
      if (looseMuonSelectionTTDil08(iMu)){
	looseMus.push_back(iMu);
      }
    }
	    
    bool pass_elec = false;
    if (looseEls.size()>1){
      if (els_charge().at(looseEls[0]) != els_charge().at(looseEls[1])){
	pass_elec = true;
      }
      if (looseMus.size()>0 && mus_p4().at(looseMus[0]).pt() > els_p4().at(looseEls[1]).pt()) pass_elec = false;
      if (looseMus.size()>0 && 
	  ( ( muonTrkIsolation(looseMus[0]) > electronTrkIsolation(looseEls[0]) 
	      && mus_charge().at(looseMus[0]) != els_charge().at(looseEls[0]) )
	    || ( muonTrkIsolation(looseMus[0]) > electronTrkIsolation(looseEls[1])
		 && mus_charge().at(looseMus[0]) != els_charge().at(looseEls[0]))
	    )
	  ) pass_elec = false; 
    }
    bool pass_muon = false;
    if (looseMus.size()>1){
      for (unsigned int iMu=0; iMu < looseMus.size(); ++iMu){
	for (unsigned int jMu=iMu+1; jMu < looseMus.size(); ++jMu){
	  if (mus_charge().at(looseMus[iMu]) != mus_charge().at(looseMus[jMu])) pass_muon = true;
	}
      }
      if (looseEls.size()>0 && els_p4().at(looseEls[0]).pt() > mus_p4().at(looseMus[1]).pt()
	  && mus_charge().at(looseMus[1]) != els_charge().at(looseEls[0])) pass_muon = false;
    }
    bool pass_elecmuon = false;
    if (looseMus.size() > 0 && looseEls.size() > 0){
      if (! pass_elec && ! pass_muon ){
	if (mus_charge().at(looseMus[0]) != els_charge().at(looseEls[0])) pass_elecmuon = true;
	if (! pass_elecmuon && looseEls.size()>1){
	  if (mus_charge().at(looseMus[0]) != els_charge().at(looseEls[0])) pass_elecmuon = true;
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
	int maxWeightType = hyp_type().at(maxWeightIndex);
	if ((maxWeightType == 0 && genpDilType == 0)
	    || ( (maxWeightType == 1 || maxWeightType == 2) && genpDilType == 1)
	    || (maxWeightType == 3 && genpDilType == 2)){
	  std::cout<<"Dil type "<<genpDilType<<" ; Strasbourg dil type "<<strasbourgDilType 
		   <<" assigned correctly by maxWeight method";
	  std::cout<<" out of"; for(unsigned int iih=0;iih<nGoodHyps;++iih)std::cout<<" "<<hyp_type().at(goodHyps[iih]);
	  std::cout<<std::endl;
	} else {
	  std::cout<<"Dil type "<<genpDilType<<" ; Strasbourg dil type "<<strasbourgDilType 
		   <<" assigned incorrectly by maxWeight method";
	  std::cout<<" out of"; for(unsigned int iih=0;iih<nGoodHyps;++iih)std::cout<<" "<<hyp_type().at(goodHyps[iih]);
	  std::cout<<std::endl;	    
	}
      }
    }
  }

  result = maxWeightIndex;
  return result;
}

