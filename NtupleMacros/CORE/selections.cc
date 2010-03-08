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

//----------------------------------------------------------------
// Simple function that tells you whether or not a given set of 
// subdetectors were valid
//----------------------------------------------------------------
bool isSubDetectorGood( unsigned int cuts ) {
     return ( ( cms2.evt_detectorStatus() & cuts ) == cuts );
}

//----------------------------------------------------------------
// Simple function that tells you whether or not a track passed 
// a particular quality flag.
//----------------------------------------------------------------
bool isTrackQuality( int index, int cuts ) {
     return ( ( cms2.trks_qualityMask().at(index) & cuts ) == cuts );
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
  if ( cms2.els_egamma_tightId().at(index)     !=  1) return false;
  if ( cms2.els_closestMuon().at(index) != -1) return false;
  if ( TMath::Abs(cms2.els_d0corr().at(index)) > 0.025)   return false;
  return true;
}


//----------------------------------------------------------------
// Electron ID without isolation or d0 cut
//----------------------------------------------------------------
bool goodElectronWithoutIsolationWithoutd0(int index) {
  if ( cms2.els_egamma_tightId().at(index)     !=  1) return false;
  if ( cms2.els_closestMuon().at(index) != -1) return false;
  return true;
}
//----------------------------------------------------------------
// loose Electron ID without isolation
//----------------------------------------------------------------
bool goodLooseElectronWithoutIsolation(int index) {
  if ( cms2.els_egamma_looseId().at(index)     !=  1) return false;
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

//
// electron isolation tk Jura

float recomputeTrackIsolation(int eleIndex, float pTCut = 0.7, float z0Cut = 0.2, 
				float strip = 0.01, float dRIn = 0.015, float dROut = 0.3)
{
 
        float isoSum = 0.0;
        for (size_t i = 0; i < cms2.trks_trk_p4().size(); ++i)
        {
        
                float dEta = cms2.trks_trk_p4()[i].Eta() - cms2.els_trk_p4()[eleIndex].Eta();
                float dPhi = acos(cos(cms2.trks_trk_p4()[i].Phi() - cms2.els_trk_p4()[eleIndex].Phi()));
                float dR = sqrt(dEta*dEta + dPhi*dPhi);
                const float &pT = cms2.trks_trk_p4()[i].Pt();
                
                if (pT < pTCut) continue;
                if (fabs(cms2.trks_z0()[i] - cms2.els_z0()[eleIndex]) > z0Cut) continue;
                if (dR < dRIn) continue;
                if (dR > dROut) continue;
                
                // jurassic isolation
                if (fabs(dEta) > strip) isoSum += pT;

        }

        return isoSum;

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

double nearestDeltaPhiJet(double Phi, int i_hyp) {

  double result = TMath::Pi();

  std::vector<LorentzVector> jets = JPTsTrilep(i_hyp, 20);

  for ( unsigned int jet = 0;
	jet < jets.size();
	++jet ) {
    double delta = TMath::Min(TMath::Abs(jets[jet].Phi() - Phi), 2*TMath::Pi() - TMath::Abs(jets[jet].Phi() - Phi));
    if ( delta < result ) 
      result = delta;
  }

  return result;
}

double nearestDeltaPhi(double Phi, int i_hyp)
{
  //WARNING!  This was designed to work in a dilepton environment - NOT a trilepton 
  double tightDPhi = TMath::Min(TMath::Abs(cms2.hyp_lt_p4()[i_hyp].Phi() - Phi), 2*TMath::Pi() - TMath::Abs(cms2.hyp_lt_p4()[i_hyp].Phi() - Phi));
  double looseDPhi = TMath::Min(TMath::Abs(cms2.hyp_ll_p4()[i_hyp].Phi() - Phi), 2*TMath::Pi() - TMath::Abs(cms2.hyp_ll_p4()[i_hyp].Phi() - Phi));

  return TMath::Min(tightDPhi, looseDPhi);

}//END nearest DeltaPhi                                                                                                                                 

double nearestDeltaPhiTrilep(double Phi, int i_hyp)
{
  //WARNING!  This was designed to work in a trilepton environment - NOT a dilepton 

  LorentzVector first,second,third;
  if (abs(cms2.hyp_trilep_first_type()[i_hyp]) == 1) {
    first = cms2.mus_p4()[cms2.hyp_trilep_first_index()[i_hyp]];
  } else {
    first = cms2.els_p4()[cms2.hyp_trilep_first_index()[i_hyp]];
  }
  if (abs(cms2.hyp_trilep_second_type()[i_hyp]) == 1) {
    second = cms2.mus_p4()[cms2.hyp_trilep_second_index()[i_hyp]];
  } else {
    second = cms2.els_p4()[cms2.hyp_trilep_second_index()[i_hyp]];
  }
  if (abs(cms2.hyp_trilep_third_type()[i_hyp]) == 1) {
    third = cms2.mus_p4()[cms2.hyp_trilep_third_index()[i_hyp]];
  } else {
    third = cms2.els_p4()[cms2.hyp_trilep_third_index()[i_hyp]];
  }
  
  double firstDPhi = TMath::Min(TMath::Abs(first.Phi() - Phi), 2*TMath::Pi() - TMath::Abs(first.Phi() - Phi));
  double secondDPhi = TMath::Min(TMath::Abs(second.Phi() - Phi), 2*TMath::Pi() - TMath::Abs(second.Phi() - Phi));
  double thirdDPhi = TMath::Min(TMath::Abs(third.Phi() - Phi), 2*TMath::Pi() - TMath::Abs(third.Phi() - Phi));

  return TMath::Min(TMath::Min(firstDPhi, secondDPhi),thirdDPhi);

}//END nearest DeltaPhi                                                                                                                                 



bool additionalZveto() {
//--------------------------------------------------------------------
// Veto events if there are two leptons in the 
// event that make the Z mass.  This uses additionalZcounter below...
//---------------------------------------------------------------------
  bool veto = false;
  if (additionalZcounter() > 0) veto = true;
  return veto;
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

std::vector<LorentzVector> JPTsTrilep(int i_hyp, double etThreshold)
{
     std::vector<LorentzVector> ret;
     const double etaMax      = 3.0;
     const double vetoCone    = 0.4;

     LorentzVector first,second,third;
     if (abs(cms2.hyp_trilep_first_type()[i_hyp]) == 1) {
       first = cms2.mus_p4()[cms2.hyp_trilep_first_index()[i_hyp]];
     } else {
       first = cms2.els_p4()[cms2.hyp_trilep_first_index()[i_hyp]];
     }
     if (abs(cms2.hyp_trilep_second_type()[i_hyp]) == 1) {
       second = cms2.mus_p4()[cms2.hyp_trilep_second_index()[i_hyp]];
     } else {
       second = cms2.els_p4()[cms2.hyp_trilep_second_index()[i_hyp]];
     }
     if (abs(cms2.hyp_trilep_third_type()[i_hyp]) == 1) {
       third = cms2.mus_p4()[cms2.hyp_trilep_third_index()[i_hyp]];
     } else {
       third = cms2.els_p4()[cms2.hyp_trilep_third_index()[i_hyp]];
     }
     
     for ( unsigned int i=0; i < cms2.jpts_p4().size(); ++i) {
	  if ( cms2.jpts_p4()[i].Et() < etThreshold ) continue;
	  if ( TMath::Abs(cms2.jpts_p4()[i].eta()) > etaMax ) continue;
	  if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(first,cms2.jpts_p4()[i])) < vetoCone ||
	       TMath::Abs(ROOT::Math::VectorUtil::DeltaR(second,cms2.jpts_p4()[i])) < vetoCone ||
	       TMath::Abs(ROOT::Math::VectorUtil::DeltaR(third,cms2.jpts_p4()[i])) < vetoCone ) continue;
	  ret.push_back(cms2.jpts_p4()[i]);
     }
     return ret;
}

unsigned int nJPTsTrilep(int i_hyp, double etThreshold)
{
     return JPTsTrilep(i_hyp, etThreshold).size();
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

bool passTrackIsolation(int index){
  //
  // track isolation
  //
  const double cut = 0.92;
  double pt  = cms2.trks_trk_p4()[index].pt();
  return (pt / (pt + trkIsolation(index))) > cut;
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
  if ( cms2.els_egamma_looseId().at(index)     !=  1) return false;
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
  double sum = cms2.els_tkIso().at(index);
  if (use_calo_iso)
    sum += cms2.els_ecalIso()[index] + cms2.els_hcalIso()[index];
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
  //  if ( cms2.els_egamma_tightId().at(index)     !=  1) return false;
  //  if ( cms2.els_egamma_looseId().at(index)     !=  1) return false;
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

bool comparePt (const LorentzVector &lv1, 
		const LorentzVector &lv2) 
{
     return lv1.pt() > lv2.pt();
}

bool GoodSusyElectronWithoutIsolation(int index) { 
  if ( cms2.els_egamma_tightId().at(index)     !=  1) return false; 
  if ( fabs(cms2.els_d0corr().at(index)) >= 0.02)   return false; 
  if ( cms2.els_closestMuon().at(index) != -1) return false; 
  if ( TMath::Abs(cms2.els_p4()[index].eta()) > 2.4) return false;
// New
//  if ( conversionElectron(index)) return false;
//  if ( isChargeFlip(index)) return false;

  return true; 
} 

bool GoodSusyElectronWithoutIsolationNoD0(int index) { 
  if ( cms2.els_egamma_tightId().at(index)     !=  1) return false; 
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

bool isNumElSUSY09(int iEl) {
  Double_t pt = cms2.els_p4()[iEl].Pt();
  if( pt < 10) return false;
  if (!GoodSusyElectronWithoutIsolation(iEl)) return false;
  if (!PassSusyElectronIsolation(iEl, true)) return false;
  return true;
}

bool isNumMuSUSY09(int iMu) {
  Double_t pt = cms2.mus_p4()[iMu].Pt();
  if (pt < 10)  return false;
  if (!GoodSusyMuonWithoutIsolation(iMu)) return false;
  if (!PassSusyMuonIsolation(iMu)) return false;

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
  //  const double cut = 0.1;
  return GoodSusyMuonWithoutIsolation(index) && PassSusyMuonIsolation(index);
}

bool PassSusyMuonIsolation(int index)
{
  const double cut = 0.1;
  return inv_mu_relsusy_iso(index) < cut;
}

bool GoodSusyElectronWithIsolation(int index, bool use_calo_iso)
{
  //  const double cut = 0.1;
  return GoodSusyElectronWithoutIsolation(index) && PassSusyElectronIsolation(index, use_calo_iso);
}

bool PassSusyElectronIsolation(int index, bool use_calo_iso)
{
  const double cut = 0.1;
  return inv_el_relsusy_iso(index, use_calo_iso) < cut;
}

bool PassSusyElectronIsolationLoose(int index, bool use_calo_iso)
{
  const double cut = 0.4; //v50_0
  //  const double cut = 0.25; //v50_1
  return inv_el_relsusy_iso(index, use_calo_iso) < cut;
}

bool GoodSusyLeptonWithIsolation(int id, int index){
  if (abs(id) == 11) return GoodSusyElectronWithIsolation(index, true);
  if (abs(id) == 13) return GoodSusyMuonWithIsolation(index);
  return false;
}

bool PassSusyLeptonIsolation(int id, int index){
  if (abs(id) == 11) return PassSusyElectronIsolation(index, true);
  if (abs(id) == 13) return PassSusyMuonIsolation(index);
  return false;
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

int numberOfExtraElectronsSUSY(int i_hyp){ 
  unsigned int nElec = 0; 
  for (int iel=0; iel < int(cms2.els_p4().size()); iel++) { 
    if ( cms2.els_p4()[iel].pt() < 10 ) continue; 
    if (fabs(cms2.els_p4()[iel].eta()) > 2.4 ) continue; 
    if (!GoodSusyElectronWithoutIsolation(iel)) continue; 
    if (!PassSusyElectronIsolation(iel, true)) continue;
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
    if (!PassSusyMuonIsolation(imu)) continue; 
    if ( TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && cms2.hyp_lt_index()[i_hyp] == imu ) continue; 
    if ( TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && cms2.hyp_ll_index()[i_hyp] == imu ) continue; 
    nMuons++; 
  } 
  return nMuons; 
} 
//   jets_p4   
std::vector<LorentzVector> getCaloJets(int i_hyp) {
  std::vector<LorentzVector> calo_jets;
  calo_jets.clear();
  
  for (unsigned int jj=0; jj < cms2.jets_p4().size(); ++jj) {
    if ((dRbetweenVectors(cms2.hyp_lt_p4()[i_hyp],cms2.jets_p4()[jj]) < 0.4)||
	(dRbetweenVectors(cms2.hyp_ll_p4()[i_hyp],cms2.jets_p4()[jj]) < 0.4)
	) continue;
    if (cms2.jets_p4()[jj].pt() < 30) continue;
    if (fabs(cms2.jets_p4()[jj].Eta()) > 2.4) continue;
    //fkw July21 2009 if (cms2.jets_emFrac()[jj] < 0.1) continue;
    calo_jets.push_back(cms2.jets_p4()[jj]);
  }
  
  if (calo_jets.size() > 1) {
       sort(calo_jets.begin(), calo_jets.end(),  comparePt);
  }
  return calo_jets;
}


std::vector<LorentzVector> getJPTJets(int i_hyp) {
  std::vector<LorentzVector> jpt_jets;
  jpt_jets.clear();
  
  for (unsigned int jj=0; jj < cms2.jpts_p4().size(); ++jj) {
    if ((dRbetweenVectors(cms2.hyp_lt_p4()[i_hyp],cms2.jpts_p4()[jj]) < 0.4)||
	(dRbetweenVectors(cms2.hyp_ll_p4()[i_hyp],cms2.jpts_p4()[jj]) < 0.4)
	) continue;
    if (cms2.jpts_p4()[jj].pt() < 30) continue;
    if (fabs(cms2.jpts_p4()[jj].Eta()) > 2.4) continue;
//    if (cms2.jpts_emFrac()[jj] < 0.1) continue;
    jpt_jets.push_back(cms2.jpts_p4()[jj]);
  }
  
  if (jpt_jets.size() > 1) {
    sort(jpt_jets.begin(), jpt_jets.end(),  comparePt);
  }
  return jpt_jets;
}

bool additionalZvetoSUSY09(int i_hyp) {

  // true if we want to veto this event
  bool veto=false;

  // first, look for Z->mumu
  for (unsigned int i=0; i < cms2.mus_p4().size(); i++) {
    bool hypLep1 = false;
    if (cms2.mus_p4().at(i).pt() < 10.)     continue;
    if (!GoodSusyMuonWithoutIsolation(i)) continue;
    if ( TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && cms2.hyp_lt_index()[i_hyp] == int(i) ) hypLep1 = true;
    if ( TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && cms2.hyp_ll_index()[i_hyp] == int(i) ) hypLep1 = true;

    for (unsigned int j=i+1; j < cms2.mus_p4().size(); j++) {
      bool hypLep2 = false;
      if (cms2.mus_p4().at(j).pt() < 10.) continue;
      if (!GoodSusyMuonWithoutIsolation(j)) continue;
      if (cms2.mus_charge().at(i) == cms2.mus_charge().at(j)) continue;
      if ( TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && cms2.hyp_lt_index()[i_hyp] == int(j) ) hypLep2 = true;
      if ( TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && cms2.hyp_ll_index()[i_hyp] == int(j) ) hypLep2 = true;
      // At least one of them has to pass isolation
      if (!PassSusyMuonIsolation(i) && !PassSusyMuonIsolation(j)) continue;
      if ( hypLep1 && hypLep2 ) continue;
      if ( !hypLep1 && !hypLep2 ) continue;
      // Make the invariant mass
      LorentzVector vec = cms2.mus_p4().at(i) + cms2.mus_p4().at(j);
      if ( inZmassWindow(vec.mass()) ) return true;

    }
  }

  // now, look for Z->ee
  for (unsigned int i=0; i < cms2.els_p4().size(); i++) {
    bool hypLep1 = false;
    if (cms2.els_p4().at(i).pt() < 10.)     continue;
    if (! GoodSusyElectronWithoutIsolation(i)) continue;
    if ( TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11 && cms2.hyp_lt_index()[i_hyp] == int(i) ) hypLep1 = true;
    if ( TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11 && cms2.hyp_ll_index()[i_hyp] == int(i) ) hypLep1 = true;
    for (unsigned int j=i+1; j<cms2.els_p4().size(); j++) {
      bool hypLep2 = false;
      if (cms2.els_p4().at(j).pt() < 10.) continue;
      if (! GoodSusyElectronWithoutIsolation(j)) continue;
      if (cms2.els_charge().at(i) == cms2.els_charge().at(j)) continue;
      // At least one of them has to pass isolation
      if (!PassSusyElectronIsolation(i, true) && ! PassSusyElectronIsolation(j, true)) continue;
      if ( TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11 && cms2.hyp_lt_index()[i_hyp] == int(j) ) hypLep2 = true;
      if ( TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11 && cms2.hyp_ll_index()[i_hyp] == int(j) ) hypLep2 = true;
      if ( hypLep1 && hypLep2 ) continue;
      if ( !hypLep1 && !hypLep2 ) continue;
      // Make the invariant mass
      LorentzVector vec = cms2.els_p4().at(i) + cms2.els_p4().at(j);
      if ( inZmassWindow(vec.mass()) ) return true;

    }
  }
  // done
  return veto;
}

// For Fake rates

bool isFakeableElSUSY09(int iEl){
 
  Double_t pt = cms2.els_p4()[iEl].Pt();
  Double_t eta = cms2.els_p4()[iEl].Eta();
 
  if( pt < 10)  return false;
  if( fabs( eta ) > 2.4 )  return false;
// New
  if ( conversionElectron(iEl)) return false;
  if ( isChargeFlip(iEl)) return false;
     
  //reject if electron is close to muon;
  if ( cms2.els_closestMuon()[iEl] > -1)  return false;
  // Isolation
  if  (inv_el_relsusy_iso(iEl, true) > 0.4) return false;
  // H/E
  if ( cms2.els_hOverE()[iEl]   > 0.2 ) return false;

  return true;
}

// Muons
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
      LorentzVector vec = cms2.mus_p4().at(i) + cms2.mus_p4().at(j);
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
      LorentzVector vec = cms2.els_p4().at(i) + cms2.els_p4().at(j);
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

int genpCountPDGId_Pt20h24(int id0, int id1, int id2){ 
  int count = 0; 
  int size = cms2.genps_id().size(); 
  for (int jj=0; jj<size; jj++) { 
    if ( cms2.genps_p4()[jj].pt() < 20  ||  fabs(cms2.genps_p4()[jj].eta())>2.4) continue;
    if (abs(cms2.genps_id()[jj]) == id0) count++; 
    if (abs(cms2.genps_id()[jj]) == id1) count++; 
    if (abs(cms2.genps_id()[jj]) == id2) count++; 
  } 
  return count; 
} 


bool matchesMCTruthDilExtended(unsigned int hypIdx){
  //this better be in the selections.cc
  bool isTrueLepton_ll = false;
  bool isTrueLepton_lt = false;
  isTrueLepton_ll = ( (abs(cms2.hyp_ll_id()[hypIdx]) == abs(cms2.hyp_ll_mc_id()[hypIdx]) &&
		       abs(cms2.hyp_ll_mc_motherid()[hypIdx]) < 50 //I wish I could match to W or Z explicitely, not in MGraph
		       )
		      || (cms2.hyp_ll_mc_id()[hypIdx]==22 && 
			  TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[hypIdx],cms2.hyp_ll_mc_p4()[hypIdx])) <0.05
			  && abs(cms2.hyp_ll_id()[hypIdx]) == abs(cms2.hyp_ll_mc_motherid()[hypIdx])
			  )
		      );
  isTrueLepton_lt = ( (abs(cms2.hyp_lt_id()[hypIdx]) == abs(cms2.hyp_lt_mc_id()[hypIdx]) &&
		       abs(cms2.hyp_lt_mc_motherid()[hypIdx]) < 50 //I wish I could match to W or Z explicitely, not in MGraph
		       )
		      || (cms2.hyp_lt_mc_id()[hypIdx]==22 && 
			  TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[hypIdx],cms2.hyp_lt_mc_p4()[hypIdx])) <0.05
			  && abs(cms2.hyp_lt_id()[hypIdx]) == abs(cms2.hyp_lt_mc_motherid()[hypIdx])
			  )
		      );
  return (isTrueLepton_lt && isTrueLepton_ll);  
}


int eventDilIndexByMaxMass(const std::vector<unsigned int>& goodHyps, bool printDebug){
  int result = -1;
  int strasbourgDilType = -1;
  unsigned int nGoodHyps = goodHyps.size();
  if ( nGoodHyps == 0 ) return result;

  float maxWeight = -1;
  unsigned int maxWeightIndex = 9999;
  
  for (unsigned int hypIdxL=0; hypIdxL < nGoodHyps; ++hypIdxL){
    unsigned int hypIdx = goodHyps[hypIdxL];
    float hypWeight = cms2.hyp_p4()[hypIdx].mass();
    if (hypWeight > maxWeight){
      maxWeight = hypWeight;
      maxWeightIndex = hypIdx;
    }
  }

  if (printDebug){
    int genpDilType = genpDileptonType();
    if (genpDilType>=0 ){ std::cout<<"Dil type "<<genpDilType<<std::endl;
      if (nGoodHyps > 1){
	int maxWeightType = cms2.hyp_type()[maxWeightIndex];
	if ((maxWeightType == 0 && genpDilType == 0)
	    || ( (maxWeightType == 1 || maxWeightType == 2) && genpDilType == 1)
	    || (maxWeightType == 3 && genpDilType == 2)){
	  std::cout<<"Dil type "<<genpDilType<<" ; Strasbourg dil type "<<strasbourgDilType 
		   <<" assigned correctly by maxWeight method";
	  std::cout<<" out of"; for(unsigned int iih=0;iih<nGoodHyps;++iih)std::cout<<" "<<cms2.hyp_type()[goodHyps[iih]];
	  std::cout<<std::endl;
	} else {
	  std::cout<<"Dil type "<<genpDilType<<" ; Strasbourg dil type "<<strasbourgDilType 
		   <<" assigned incorrectly by maxWeight method";
	  std::cout<<" out of"; for(unsigned int iih=0;iih<nGoodHyps;++iih)std::cout<<" "<<cms2.hyp_type()[goodHyps[iih]];
	  std::cout<<std::endl;	    
	}
      }
    } else{
      if (genpCountPDGId(11,13,15) == 2){
	std::cout<<"TauDil type "<<std::endl;
      }
    }
    int nMCTruth = 0;
    for(unsigned int iih=0;iih<nGoodHyps;++iih) if (matchesMCTruthDilExtended(goodHyps[iih])) nMCTruth++;
    std::cout<<"Ne: "<<genpCountPDGId_Pt20h24(11)<<" nmu: "<<genpCountPDGId_Pt20h24(13)<<" ntau: "<<genpCountPDGId_Pt20h24(15)
	     <<" ngood "<<nGoodHyps
	     <<" hyp_typeM: "<<cms2.hyp_type()[maxWeightIndex]<<" matchMC "<<(matchesMCTruthDilExtended(maxWeightIndex)? 1 : 0)
	     <<" nMatches "<<nMCTruth
	     <<" ltid "<< cms2.hyp_lt_id()[maxWeightIndex]
	     <<" ltmcid "<< cms2.hyp_lt_mc_id()[maxWeightIndex]<<" ltmcmid "<< cms2.hyp_lt_mc_motherid()[maxWeightIndex]
	     <<" llid "<< cms2.hyp_ll_id()[maxWeightIndex]
	     <<" llmcid "<< cms2.hyp_ll_mc_id()[maxWeightIndex]<<" llmcmid "<< cms2.hyp_ll_mc_motherid()[maxWeightIndex]
	     <<std::endl;
  }

  result = maxWeightIndex;
  return result;
}

int eventDilIndexByWeightTTDil08(const std::vector<unsigned int>& goodHyps, int& strasbourgDilType, bool printDebug, bool usePtOnlyForWeighting){
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
      if (! usePtOnlyForWeighting && cms2.els_egamma_tightId().at(i_lt)) hypWeight_lt += 0.2;
    }
    if (abs(id_lt) == 13){
      //I want to select "trk & cal"-isolated ones	    
      hypWeight_iso += (isoTk_lt*isoCal_lt - 0.25);//shift by 0.25 to be positive-definite
      if (! usePtOnlyForWeighting) hypWeight_lt += 0.4;
    }
    if (abs(id_ll) == 11){
      //I want to select "trk & cal"-isolated ones
      hypWeight_iso *= (isoTk_ll*isoCal_ll - 0.25); //shift by 0.25 to be positive-definite
      if (! usePtOnlyForWeighting && cms2.els_egamma_tightId().at(i_ll)) hypWeight_ll += 0.2;
    }
    if (abs(id_ll) == 13){
      //I want to select "trk & cal"-isolated ones
      hypWeight_iso *= (isoTk_ll*isoCal_ll - 0.25); //shift by 0.25 to be positive-definite
      if (! usePtOnlyForWeighting) hypWeight_ll += 0.4;
    }
    float pt_lt = cms2.hyp_lt_p4().at(hypIdx).pt();
    float pt_ll = cms2.hyp_ll_p4().at(hypIdx).pt();
    hypWeight_lt += (1. - 20./pt_lt*20./pt_lt);
    hypWeight_ll += (1. - 20./pt_ll*20./pt_ll);
    
    if (usePtOnlyForWeighting){
      hypWeight = hypWeight_ll*hypWeight_lt; //again, desire to have both good
    } else {
      hypWeight = hypWeight_ll*hypWeight_lt*hypWeight_iso; //again, desire to have both good
    }

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
    int nMCTruth = 0;
    for(unsigned int iih=0;iih<nGoodHyps;++iih) if (matchesMCTruthDilExtended(goodHyps[iih])) nMCTruth++;
    std::cout<<"Ne: "<<genpCountPDGId_Pt20h24(11)<<" nmu: "<<genpCountPDGId_Pt20h24(13)<<" ntau: "<<genpCountPDGId_Pt20h24(15)
	     <<" ngood "<<nGoodHyps<<" SBtype "<<strasbourgDilType
	     <<" hyp_typeM: "<<cms2.hyp_type()[maxWeightIndex]<<" matchMC "<<(matchesMCTruthDilExtended(maxWeightIndex)? 1 : 0)
	     <<" nMatches "<<nMCTruth
	     <<" ltid "<< cms2.hyp_lt_id()[maxWeightIndex]
	     <<" ltmcid "<< cms2.hyp_lt_mc_id()[maxWeightIndex]<<" ltmcmid "<< cms2.hyp_lt_mc_motherid()[maxWeightIndex]
	     <<" llid "<< cms2.hyp_ll_id()[maxWeightIndex]
	     <<" llmcid "<< cms2.hyp_ll_mc_id()[maxWeightIndex]<<" llmcmid "<< cms2.hyp_ll_mc_motherid()[maxWeightIndex]
	     <<std::endl;    
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
  if( ! cms2.els_egamma_looseId()[iEl] )
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

//----------------------------------------------------
// Conversions stuff
//----------------------------------------------------

//old cut to find conversions
bool conversionElectron(int electron) {
  // true if electron is a conversion electron
  if( fabs(cms2.els_conv_dist()[electron]) < 0.02 && fabs(cms2.els_conv_dcot()[electron]) < 0.02)
    return true;

  return false;
}

//utility function to get the dist and delta cot theta
std::pair<float, float> getConversionInfo(LorentzVector trk1_p4, 
					  int trk1_q, float trk1_d0, 
					  LorentzVector trk2_p4,
					  int trk2_q, float trk2_d0,
					  float bField) {
  
  
  double tk1Curvature = -0.3*bField*(trk1_q/trk1_p4.pt())/100.;
  double rTk1 = fabs(1./tk1Curvature);
  double xTk1 = (1./tk1Curvature - trk1_d0)*cos(trk1_p4.phi());
  double yTk1 = (1./tk1Curvature - trk1_d0)*sin(trk1_p4.phi());
    
  double tk2Curvature = -0.3*bField*(trk2_q/trk2_p4.pt())/100.;
  double rTk2 = fabs(1./tk2Curvature);
  double xTk2 = (1./tk2Curvature - trk2_d0)*cos(trk2_p4.phi());
  double yTk2 = (1./tk2Curvature - trk2_d0)*sin(trk2_p4.phi());
	 
  double dist = sqrt(pow(xTk1-xTk2, 2) + pow(yTk1-yTk2 , 2));
  dist = dist - (rTk1 + rTk2);

  double dcot = 1/tan(trk1_p4.theta()) - 1/tan(trk2_p4.theta());

  return make_pair(dist, dcot);
  
}

//new conversion loop to determine if the electron is from a conversion or not
bool isconversionElectron09(int elIdx) {

  for(unsigned int tkIdx = 0; tkIdx < cms2.trks_trk_p4().size(); tkIdx++) {
       if(ROOT::Math::VectorUtil::DeltaR(cms2.els_trk_p4()[elIdx], cms2.trks_trk_p4()[tkIdx]) > 0.3)
      continue;
    //skip the electron's track
    if(cms2.els_trkidx()[elIdx] == int(tkIdx) && cms2.els_trkshFrac()[elIdx] > 0.45)
      continue;
    //ship non-opp sign candidates
    if(cms2.trks_charge()[tkIdx] + cms2.els_charge()[elIdx] != 0)
      continue;
    
    std::pair<float, float> temp = getConversionInfo(cms2.els_trk_p4()[elIdx], cms2.els_charge()[elIdx], cms2.els_d0()[elIdx], 
						     cms2.trks_trk_p4()[tkIdx], cms2.trks_charge()[tkIdx], cms2.trks_d0()[tkIdx],
						     cms2.evt_bField());
    
    if(fabs(temp.first) < 0.02 && fabs(temp.second) < 0.02)
      return true;
    
  }//track loop
  
  return false;
  
}
//---------------------------------------------------
//End conversion functions
//---------------------------------------------------




// false if below ptcut, aboveabsEtaCut, below dRCut wrt hypothesis
bool isGoodDilHypJet(unsigned int jetIdx, unsigned int hypIdx, double ptCut, double absEtaCut, double dRCut, bool muJetClean){
  if (cms2.jets_p4()[jetIdx].pt()< ptCut || fabs(cms2.jets_p4()[jetIdx].eta())> absEtaCut) return false;
  double dR_ll = ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[hypIdx],cms2.jets_p4()[jetIdx]);
  double dR_lt = ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[hypIdx],cms2.jets_p4()[jetIdx]);
  
  if (abs(cms2.hyp_ll_id()[hypIdx]) == 11){
    if (dR_ll < dRCut) return false;
  }
  if (abs(cms2.hyp_lt_id()[hypIdx]) == 11){
    if (dR_lt < dRCut) return false;
  }

  if (muJetClean){
    if (abs(cms2.hyp_ll_id()[hypIdx]) == 13){
      if (dR_ll < dRCut) return false;
    }
    if (abs(cms2.hyp_lt_id()[hypIdx]) == 13){
      if (dR_lt < dRCut) return false;
    }
  }

  return true;

}

// false if below ptcut, aboveabsEtaCut, below dRCut wrt hypothesis
bool isGoodDilHypJPTJet(unsigned int jetIdx, unsigned int hypIdx, double ptCut, double absEtaCut, double dRCut){
  if (cms2.jpts_p4()[jetIdx].pt()< ptCut || fabs(cms2.jpts_p4()[jetIdx].eta())> absEtaCut) return false;
  double dR_ll = ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[hypIdx],cms2.jpts_p4()[jetIdx]);
  double dR_lt = ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[hypIdx],cms2.jpts_p4()[jetIdx]);
  
  if (dR_ll < dRCut) return false;
  if (dR_lt < dRCut) return false;

  return true;

}


int findPrimTrilepZ(int i_hyp, double &mass) {
  // find primary Z candidate in trilepton hyp
  //
  // return: 1: first lepton was not used in Z candidate
  //         2: second lepton was not used in Z candidate
  //         3: second lepton was not used in Z candidate
  //         999: no Z candidate could be found
  //
  //         900: error code, something went wrong

  // z mass array, coding:
  // index 0: first, second
  // index 1: first, third
  // index 2: second, third
  double z_mass[3] = {0.,0.,0.};

  // check if first and second lepton form Z candidate
  if ( abs(cms2.hyp_trilep_first_type()[i_hyp]) == abs(cms2.hyp_trilep_second_type()[i_hyp]) ) {
    if ( abs(cms2.hyp_trilep_first_type()[i_hyp]) == 1 ) {
      if ( cms2.mus_charge()[cms2.hyp_trilep_first_index()[i_hyp]] !=  cms2.mus_charge()[cms2.hyp_trilep_second_index()[i_hyp]] ) {
	LorentzVector z = cms2.mus_p4()[cms2.hyp_trilep_first_index()[i_hyp]] + cms2.mus_p4()[cms2.hyp_trilep_second_index()[i_hyp]];
// 	if ( inZmassWindow(z.mass()) )
	  z_mass[0] = z.mass();
      }
    } else {
      if ( cms2.els_charge()[cms2.hyp_trilep_first_index()[i_hyp]] !=  cms2.els_charge()[cms2.hyp_trilep_second_index()[i_hyp]] ) {
	LorentzVector z = cms2.els_p4()[cms2.hyp_trilep_first_index()[i_hyp]] + cms2.els_p4()[cms2.hyp_trilep_second_index()[i_hyp]];
// 	if ( inZmassWindow(z.mass()) )
	  z_mass[0] = z.mass();
      }
    }
  }
  // check if first and third lepton form Z candidate
  if ( abs(cms2.hyp_trilep_first_type()[i_hyp]) == abs(cms2.hyp_trilep_third_type()[i_hyp]) ) {
    if ( abs(cms2.hyp_trilep_first_type()[i_hyp]) == 1 ) {
      if ( cms2.mus_charge()[cms2.hyp_trilep_first_index()[i_hyp]] !=  cms2.mus_charge()[cms2.hyp_trilep_third_index()[i_hyp]] ) {
	LorentzVector z = cms2.mus_p4()[cms2.hyp_trilep_first_index()[i_hyp]] + cms2.mus_p4()[cms2.hyp_trilep_third_index()[i_hyp]];
// 	if ( inZmassWindow(z.mass()) )
	  z_mass[1] = z.mass();
      }
    } else {
      if ( cms2.els_charge()[cms2.hyp_trilep_first_index()[i_hyp]] !=  cms2.els_charge()[cms2.hyp_trilep_third_index()[i_hyp]] ) {
	LorentzVector z = cms2.els_p4()[cms2.hyp_trilep_first_index()[i_hyp]] + cms2.els_p4()[cms2.hyp_trilep_third_index()[i_hyp]];
// 	if ( inZmassWindow(z.mass()) )
	  z_mass[1] = z.mass();
      }
    }
  }
  // check if second and third lepton form Z candidate
  if ( abs(cms2.hyp_trilep_second_type()[i_hyp]) == abs(cms2.hyp_trilep_third_type()[i_hyp]) ) {
    if ( abs(cms2.hyp_trilep_second_type()[i_hyp]) == 1 ) {
      if ( cms2.mus_charge()[cms2.hyp_trilep_second_index()[i_hyp]] !=  cms2.mus_charge()[cms2.hyp_trilep_third_index()[i_hyp]] ) {
	LorentzVector z = cms2.mus_p4()[cms2.hyp_trilep_second_index()[i_hyp]] + cms2.mus_p4()[cms2.hyp_trilep_third_index()[i_hyp]];
// 	if ( inZmassWindow(z.mass()) )
	  z_mass[2] = z.mass();
      }
    } else {
      if ( cms2.els_charge()[cms2.hyp_trilep_second_index()[i_hyp]] !=  cms2.els_charge()[cms2.hyp_trilep_third_index()[i_hyp]] ) {
	LorentzVector z = cms2.els_p4()[cms2.hyp_trilep_second_index()[i_hyp]] + cms2.els_p4()[cms2.hyp_trilep_third_index()[i_hyp]];
// 	if ( inZmassWindow(z.mass()) )
	  z_mass[2] = z.mass();
      }
    }
  }

  // check which combination is nearest to Z mass and return unsused lepton
//   if ( z_mass[0] == 0 && z_mass[1] == 0 && z_mass[2] == 0 )
//     return 999;
  int ret = 0;
  if ( fabs( z_mass[0] - 91) <= fabs( z_mass[1] - 91) && fabs( z_mass[0] - 91) <= fabs( z_mass[2] - 91) ) {
    mass = z_mass[0];
    ret = 3;
  } else if ( fabs( z_mass[1] - 91) <= fabs( z_mass[0] - 91) && fabs( z_mass[1] - 91) <= fabs( z_mass[2] - 91) ) {
    mass = z_mass[1];
    ret = 2;
  } else if ( fabs( z_mass[2] - 91) <= fabs( z_mass[0] - 91) && fabs( z_mass[2] - 91) <= fabs( z_mass[1] - 91) ) {
    mass = z_mass[2];
    ret = 1;
  }

  if ( !inZmassWindow(mass) ) {
    return 999;
  } else {
    return ret;
  }

  return 900;
    
}

bool vetoAddZ(int i_hyp, int unusedLepton, double &mass) {
  // veto event if unused lepton (not used to form primary Z) and a isolated track form a second Z, only for trilepton
  LorentzVector lepton;
  if ( unusedLepton == 1 ) {
    if ( abs(cms2.hyp_trilep_first_type()[i_hyp]) == 1 ) 
      lepton = cms2.mus_p4()[cms2.hyp_trilep_first_index()[i_hyp]];
    else
      lepton = cms2.els_p4()[cms2.hyp_trilep_first_index()[i_hyp]];
  } else if ( unusedLepton == 2 ) {
    if ( abs(cms2.hyp_trilep_second_type()[i_hyp]) == 1 ) 
      lepton = cms2.mus_p4()[cms2.hyp_trilep_second_index()[i_hyp]];
    else
      lepton = cms2.els_p4()[cms2.hyp_trilep_second_index()[i_hyp]];
  } else if ( unusedLepton == 3 ) {
    if ( abs(cms2.hyp_trilep_third_type()[i_hyp]) == 1 ) 
      lepton = cms2.mus_p4()[cms2.hyp_trilep_third_index()[i_hyp]];
    else
      lepton = cms2.els_p4()[cms2.hyp_trilep_third_index()[i_hyp]];
  }

  mass = -1.;

  for ( int track = 0;
	track < (int)cms2.trks_trk_p4().size();
	++track ) {

    // exclude track from first lepton
    if ( abs(cms2.hyp_trilep_second_type()[i_hyp]) == 1 ) {
      if ( cms2.mus_trkidx()[cms2.hyp_trilep_first_index()[i_hyp]] == track ) continue;
    } else {
      if ( cms2.els_trkidx()[cms2.hyp_trilep_first_index()[i_hyp]] == track ) continue;
    }

    // exclude track from second lepton
    if ( abs(cms2.hyp_trilep_second_type()[i_hyp]) == 1 ) {
      if ( cms2.mus_trkidx()[cms2.hyp_trilep_second_index()[i_hyp]] == track ) continue;
    } else {
      if ( cms2.els_trkidx()[cms2.hyp_trilep_second_index()[i_hyp]] == track ) continue;
    }

    // exclude track from third lepton
    if ( abs(cms2.hyp_trilep_second_type()[i_hyp]) == 1 ) {
      if ( cms2.mus_trkidx()[cms2.hyp_trilep_third_index()[i_hyp]] == track ) continue;
    } else {
      if ( cms2.els_trkidx()[cms2.hyp_trilep_third_index()[i_hyp]] == track ) continue;
    }

    if ( passTrackIsolation(track) && cms2.trks_trk_p4()[track].Pt() >= 20. ) {
      LorentzVector z = lepton + cms2.trks_trk_p4()[track];
      if ( fabs(z.mass() - 91) <= fabs(mass - 91 ) )
	mass = z.mass();
    }
  }

  if ( inZmassWindow(mass) ) return true;

  return false;

}

