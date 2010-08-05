#include <assert.h>
#include <algorithm>
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TDatabasePDG.h"
#include "TSystem.h"
#include "electronSelections.cc"
#include "muonSelections.cc"
#include "metSelections.cc"
//#include "CMS2.cc"
#include "CMS2.h"
#include "electronSelectionsParameters.cc"

using namespace tas;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;
bool isGoodLeptonNoIso(int id, int lepIdx);//, bool used0wrtPV = false);
bool isGoodLeptonwIso(int id, int lepIdx);//, bool used0wrtPV  = false);
bool isGoodHypNoIso(int hypIdx);//, bool used0wrtPV = false);
bool isGoodHypwIso(int hypIdx);//, bool used0wrtPV = false);
bool isGoodDilHypJet(LorentzVector jetp4, unsigned int& hypIdx, double ptCut, double absEtaCut, double dRCut, bool muJetClean);
std::pair<float,float> getMet(string& algo, unsigned int hypIdx, std::string prefix);
bool inZmassWindow (float mass);
bool passTriggersMu9orLisoE15(int dilType);
unsigned int selectHypByHighestSumPt(const vector<unsigned int> &v_goodHyps);
bool isFakeDenominatorElectron_v1(unsigned int lepIdx);
bool isFakeDenominatorElectron_v2(unsigned int lepIdx);
double getd0wrtPV(LorentzVector p4, float d0);
void correctTcMETForHypMus(unsigned int hypIdx, double& met, double& metPhi);
bool passesCaloJetID(LorentzVector jetp4);
bool goodEGTrigger5July2010 (bool mc);
int nHLTObjects(string arg);
LorentzVector p4HLTObject(string arg, int) ;
bool isSpikeElectron(const unsigned int index);

//spike cleaning. Lives here temporarily
bool isSpikeElectron(const unsigned int index) {

    const int scidx = cms2.els_scindex()[index];
    bool isSpike = false;
    if (scidx != -1) {
        //subtract twice max since max is in both 1x3 and 3x1, and we want neither
        const float r4 = (cms2.scs_e1x3()[scidx] + cms2.scs_e3x1()[scidx] - 2*cms2.scs_eMax()[scidx])/cms2.scs_eMax()[scidx];
        if (r4 < 0.05) isSpike = true;
    }

    return isSpike;

}


//--------------------------------------------------
// EG trigger selection from 5 July 2010
// data: Photon10 OR Electron 10 OR Photon 15
// MC  : E10 OR Photon15
// (Note: in data Photon10 is prescaled at some point)
// (Note: in MC Photon15 has a different L1 requirement, so we 
//        use Photon10 and we tighten the pt threshold)
//--------------------------------------------------
bool goodEGTrigger5July2010 (bool mc) {
  if (mc) {
 
    int e10 = nHLTObjects("HLT_Ele10_LW_L1R");
    if (e10 != 0) return true;

    int p10 = nHLTObjects("HLT_Photon10_L1R");
    for (int i=0; i<p10; i++) {
      LorentzVector v = p4HLTObject("HLT_Photon10_L1R", i);
      if (v.pt() > 15.) return true;
    }
 
  } else {  // data now
    int e10 = nHLTObjects("HLT_Ele10_LW_L1R");
    if (e10 != 0) return true;

    int p10 = nHLTObjects("HLT_Photon10_L1R");
    if (p10 != 0) return true;

    p10 = nHLTObjects("HLT_Photon10_Cleaned_L1R");
    if (p10 != 0) return true;

    int p15 = nHLTObjects("HLT_Photon15_L1R");
    if (p15 != 0) return true;

    p15 = nHLTObjects("HLT_Photon15_Cleaned_L1R");
    if (p15 != 0) return true;

    int e15 = nHLTObjects("HLT_Ele15_LW_L1R");
    if(e15 != 0) return true;

 
  }
  return false;
}


//-----------------------------------------------------
// Returns the nth object that passes a given trigger
// (n starts from 0)
///----------------------------------------------------
LorentzVector p4HLTObject( string arg, int objNumber){
 
  TString HLTTrigger( arg );
  int trigIndx = -1;
  vector<TString>::const_iterator begin_it = cms2.hlt_trigNames().begin();
  vector<TString>::const_iterator end_it = cms2.hlt_trigNames().end();
  vector<TString>::const_iterator found_it = find(begin_it, end_it, HLTTrigger );
  if( (found_it != end_it) ){
    trigIndx = found_it - begin_it;
    //cout << "p4HLTObject: Found Trigger: " << arg << endl;
  }
  else {
    cout << "p4HLTObject: Cannot find Trigger: " << arg << endl;
    gSystem->Exit(1);
  }

  int nobj = cms2.hlt_trigObjs_p4().at(trigIndx).size();
  if (nobj == 0 ) {
    cout << "ERROR: nobj == 0" << endl;
    gSystem->Exit(1);
  }

  if (objNumber > (nobj-1)) {
    cout << "ERROR: requested object number " << objNumber << " but we only have " << nobj <<endl;
    gSystem->Exit(1);
  }

  return cms2.hlt_trigObjs_p4().at(trigIndx).at(objNumber);

}
//--------------------------------------------------------
// Returns the number of objects passing a given trigger
// Returns zero if the trigger failed
// Returns -1 if the trigger passed but no onjects were found
//--------------------------------------------------------
int nHLTObjects( string arg ){

  // put the trigger name into a string
  TString HLTTrigger( arg );

  // Did the trigger pass?
  if ( !(cms2.passHLTTrigger(HLTTrigger)) ) return 0;

  // The trigger passed, see how many associated objects there are
  int trigIndx = -1;
  vector<TString>::const_iterator begin_it = cms2.hlt_trigNames().begin();
  vector<TString>::const_iterator end_it = cms2.hlt_trigNames().end();
  vector<TString>::const_iterator found_it = find(begin_it, end_it, HLTTrigger );
  if( (found_it != end_it) ){
    trigIndx = found_it - begin_it;
    //cout << "nHLTObjects: Found Trigger: " << arg << endl;
  }
  else {
    cout << "nHLTObjects: Cannot find Trigger " << arg << endl;
    return 0;
  }

  int nobj = 0;
  for( unsigned int i=0; i < cms2.hlt_trigObjs_p4().at(trigIndx).size(); i++ ){
    nobj++;
    //cout << "\t" << i << ", (pt, eta, phi): " << hlt_trigObjs_p4().at(trigIndx).at(i).pt() << " "
    //              << hlt_trigObjs_p4().at(trigIndx).at(i).eta() << " " << hlt_trigObjs_p4().at(trigIndx).at(i).phi() << endl;
  }

  // cout << " Number of jets = " << njets << endl;

  if (nobj == 0) return -1;
  return nobj;
}

/******************************************************************************************/     
// good lepton (either mu or electron, no isolation cuts)
/******************************************************************************************/
bool isGoodLeptonNoIso(int id, int lepIdx) {//, bool used0wrtPV) {

  if(abs(id) == 11) {
    
    const cuts_t elIDcuts =	(1ll<<ELEIP_400)		|								   
      (1ll<<ELENOMUON_010)		|
      (1ll<<ELENOTCONV_HITPATTERN)	|
      (1ll<<ELENOTCONV_DISTDCOT002)	|
      (1ll<<ELESCET_010)		|
      (1ll<<ELEPT_010)		|
      (1ll<<ELEETA_250)		|
      (1ll<<ELESEED_ECAL);
    
    bool isSpike = isSpikeElectron(lepIdx);
    unsigned int answer_vbtf = electronId_VBTF(lepIdx, VBTF_35X_90);
    bool elsvbtf90_ = ( ( answer_vbtf & (1ll<<ELEID_ID) ) == (1ll<<ELEID_ID) );
    return (pass_electronSelection(lepIdx, elIDcuts) && elsvbtf90_ && !isSpike);
  }

  if(abs(id) == 13) {
    if ( mus_p4()[lepIdx].pt() < 5.) {
      std::cout << "muonID ERROR: requested muon is too low pt,  Abort." << std::endl;
      return false;
    }
    if ( TMath::Abs(cms2.mus_p4()[lepIdx].eta()) > 2.5)  return false; // eta cut
    if (cms2.mus_gfit_chi2().at(lepIdx)/cms2.mus_gfit_ndof().at(lepIdx) >= 10) return false; //glb fit chisq
    if (((cms2.mus_type().at(lepIdx)) & (1<<1)) == 0)    return false; // global muon
    if (((cms2.mus_type().at(lepIdx)) & (1<<2)) == 0)    return false; // tracker muon
    if (cms2.mus_validHits().at(lepIdx) < 11)            return false; // # of tracker hits
    if (cms2.mus_gfit_validSTAHits().at(lepIdx) == 0)    return false; // Glb fit must have hits in mu chambers
    if (TMath::Abs(cms2.mus_d0corr().at(lepIdx)) > 0.02) return false; // d0 from beamspot
  }
  

  return true;
}

/******************************************************************************************/     
// isolated lepton (either mu or electron)
/******************************************************************************************/
bool isGoodLeptonwIso(int id, int lepIdx) { //, bool used0wrtPV) {

  
  if(!isGoodLeptonNoIso(id, lepIdx))
    return false;

  // 11 is a electron
  if(abs(id)== 11) {
       const cuts_t elISOcuts =   (1ll<<ELEISO_REL015);
       if (!pass_electronSelection(lepIdx, elISOcuts))
	    return false;
  }

  // 13 is a muon
  if(abs(id) == 13) 
    if(muonIsoValue(lepIdx) > 0.15)   return false;
  
  return true;
}

/******************************************************************************************/     
// are the leptons in the hypothesis good (all cuts but isolation?)
/******************************************************************************************/
bool isGoodHypNoIso(int hypIdx) {//, bool used0wrtPV) {
  
  
  if(!isGoodLeptonNoIso(hyp_lt_id()[hypIdx], hyp_lt_index()[hypIdx]))//, used0wrtPV)
     return false;
  if(!isGoodLeptonNoIso(hyp_ll_id()[hypIdx], hyp_ll_index()[hypIdx]))//, used0wrtPV)
    return false;

  return true;
}

/******************************************************************************************/     
// are the leptons in the hypothesis isolated?
/******************************************************************************************/     
bool isGoodHypwIso(int hypIdx) {


  if(!isGoodLeptonwIso(hyp_lt_id()[hypIdx], hyp_lt_index()[hypIdx]))
    return false;
  if(!isGoodLeptonwIso(hyp_ll_id()[hypIdx], hyp_ll_index()[hypIdx]))
    return false;


  return true;
}

/******************************************************************************************/     
// is it a good jet?
/******************************************************************************************/     
bool isGoodDilHypJet(LorentzVector jetp4, unsigned int& hypIdx, double ptCut, double absEtaCut, double dRCut, bool muJetClean){
		     
  if(jetp4.Pt() < ptCut)
    return false;  
  if(fabs(jetp4.Eta()) > absEtaCut)
    return false;
  
  double dR_ll = ROOT::Math::VectorUtil::DeltaR(hyp_ll_p4()[hypIdx],jetp4);
  double dR_lt = ROOT::Math::VectorUtil::DeltaR(hyp_lt_p4()[hypIdx],jetp4);
  
  if (abs(hyp_ll_id()[hypIdx]) == 11){
    if (dR_ll < dRCut) return false;
  }
  if (abs(hyp_lt_id()[hypIdx]) == 11){
    if (dR_lt < dRCut) return false;
  }

  if (muJetClean){
    if (abs(hyp_ll_id()[hypIdx]) == 13){
      if (dR_ll < dRCut) return false;
    }
    if (abs(hyp_lt_id()[hypIdx]) == 13){
      if (dR_lt < dRCut) return false;
    }
  }

  return true;
}

/******************************************************************************************/     
// if jpt jet, does it pass jet ID
/******************************************************************************************/     
bool passesCaloJetID(LorentzVector jetp4)
{
     int jet_idx = -1;
     double minDR = 999;

     for (unsigned int i = 0; i < jets_p4().size(); i++)
     {
	  double deltaR = ROOT::Math::VectorUtil::DeltaR(jetp4, jets_p4()[i]);

	  if (deltaR < minDR)
	  {
	       minDR = deltaR;
	       jet_idx = i;
	  }
     }

     if (jet_idx < 0)
	  return false;

     if (jets_emFrac()[jet_idx] < 0.01 || jets_fHPD()[jet_idx] > 0.98 || jets_n90Hits()[jet_idx] < 2)
	  return false;

     return true;
}

/******************************************************************************************/     
//return the MET and the MET phi instead of a bool because the MT2 needs it
/******************************************************************************************/     
std::pair<float,float> getMet(string& algo, unsigned int hypIdx, std::string prefix) {
  
  if(algo != "tcMET" && algo != "muCorMET" && algo != "pfMET") {
    cout << algo << "IS NOT A RECOGNIZED MET ALGORITHM!!!!! PLEASE CHECK YOUR CODE!!!";
    return make_pair(-99999., -99999.);
  }
  if(algo == "tcMET") {
    double tcmet = evt_tcmet();
    double tcmetPhi = evt_tcmetPhi();
    correctTcMETForHypMus(hypIdx, tcmet, tcmetPhi);
    return make_pair(tcmet, tcmetPhi);
  }
  if(algo == "muCorMET")
    return make_pair(evt_metMuonCorr(), evt_metMuonCorrPhi());
  if(algo == "pfMET")
    return make_pair(evt_pfmet(), evt_pfmetPhi());
  
  
  return make_pair(-99999., -99999);
  
}


/******************************************************************************************/     
//trigger requirement
/******************************************************************************************/         
bool passTriggersMu9orLisoE15(int dilType) {
  
  //TString method
  bool hlt_ele15_lw_l1r = cms2.passHLTTrigger("HLT_Ele15_SW_L1R");
  bool hltMu9           = cms2.passHLTTrigger("HLT_Mu9");
  
  if (dilType == 0 && ! (hltMu9) ) return false;
  if ((dilType == 1 || dilType == 2) && ! (hltMu9 || hlt_ele15_lw_l1r)) return false;
  if (dilType == 3 && ! hlt_ele15_lw_l1r) return false;     

  return true;
}

/******************************************************************************************/     
//hypothesis disabmiguation
/******************************************************************************************/     
unsigned int selectHypByHighestSumPt(const vector<unsigned int> &v_goodHyps) {
  
  float maxSumPt = 0;
  unsigned int bestHypIdx = 0;
  for(unsigned int i = 0; i < v_goodHyps.size(); i++) {
    
    unsigned int index = v_goodHyps.at(i);
    float sumPt = hyp_lt_p4()[index].Pt() + hyp_ll_p4()[index].Pt();
    if( sumPt > maxSumPt) {
      maxSumPt = sumPt;
      bestHypIdx = index;
    }
  }

  return bestHypIdx;

}


/******************************************************************************************/     
//electron FO v1
/******************************************************************************************/     
bool isFakeDenominatorElectron_v1(unsigned int lepIdx) {
  if (fabs(els_p4()[lepIdx].Eta()) > 2.5)    return false;
  if (els_p4()[lepIdx].Pt() < 20.)           return false;
  if (!electronId_noMuon(lepIdx))            return false;
  if (isFromConversionPartnerTrack(lepIdx))  return false;
//  if (electronIsolation_relsusy_cand1(lepIdx, true) > 0.40) return false;

  return true;
  
}

/******************************************************************************************/     
//electron FO v1
/******************************************************************************************/     
bool isFakeDenominatorElectron_v2(unsigned int lepIdx) {

  if (fabs(els_p4()[lepIdx].Eta()) > 2.5)    return false;
  if (els_p4()[lepIdx].Pt() < 20.)           return false;
  if (!electronId_noMuon(lepIdx))            return false;
  if (isFromConversionPartnerTrack(lepIdx))  return false;
//  if (electronIsolation_relsusy_cand1(lepIdx, true) > 0.10) return false;

  return true;

}


/******************************************************************************************/     
//electron impact parameter with respect to the highest sumpt pv
/******************************************************************************************/     
double getd0wrtPV(LorentzVector p4, float d0) {

   
  double max_sumpt = -1;
  int i_max = -1;
  assert(cms2.vtxs_sumpt().size() == cms2.vtxs_isFake().size());
  assert(cms2.vtxs_sumpt().size() == cms2.vtxs_position().size());
  assert(cms2.vtxs_sumpt().size() == cms2.vtxs_covMatrix().size());
  for (unsigned int i = 0; i < cms2.vtxs_sumpt().size(); ++i) {
    if (cms2.vtxs_isFake().at(i))
	       continue;
    if (cms2.vtxs_sumpt().at(i) > max_sumpt) {
      max_sumpt = cms2.vtxs_sumpt().at(i);
      i_max = i;
    }
  }

   if (i_max != -1) {
     const double bx = vtxs_position().at(i_max).x();
     const double by = vtxs_position().at(i_max).y();
     double phi = p4.phi();
     double d0vtx = d0 - bx * sin(phi) + by * cos(phi);
     return d0vtx;
   }

   
   cout << "did not find a PV!!!" << endl;
   return 99999;


 }


//*****************************************************************************************
//correct MET for hyp mus that are not used in MET correction
//*****************************************************************************************
void correctTcMETForHypMus(unsigned int hypIdx, double& met, double& metPhi){
  if (cms2.hyp_type()[hypIdx] ==3) return;
  double lmetx = met*cos(metPhi);
  double lmety = met*sin(metPhi);

  unsigned int i_lt = cms2.hyp_lt_index()[hypIdx];
  unsigned int i_ll = cms2.hyp_ll_index()[hypIdx];
  if (abs(cms2.hyp_lt_id()[hypIdx])==13){
    if(cms2.mus_tcmet_flag()[i_lt] == 0){
      lmetx += cms2.mus_met_deltax()[i_lt] - cms2.mus_p4()[i_lt].x();
      lmety += cms2.mus_met_deltay()[i_lt] - cms2.mus_p4()[i_lt].y();
    } else if (cms2.mus_tcmet_flag()[i_lt] == 4){
	 lmetx += -cms2.mus_tcmet_deltax()[i_lt] + cms2.mus_met_deltax()[i_lt] - cms2.mus_p4()[i_lt].x() + cms2.trks_trk_p4()[cms2.mus_trkidx()[i_lt]].x(); 
	 lmety += -cms2.mus_tcmet_deltay()[i_lt] + cms2.mus_met_deltay()[i_lt] - cms2.mus_p4()[i_lt].y() + cms2.trks_trk_p4()[cms2.mus_trkidx()[i_lt]].y(); 
    }
  }
  if (abs(cms2.hyp_ll_id()[hypIdx])==13){
    if(cms2.mus_tcmet_flag()[i_ll] == 0){ 
      lmetx+= cms2.mus_met_deltax()[i_ll] - cms2.mus_p4()[i_ll].x(); 
      lmety+= cms2.mus_met_deltay()[i_ll] - cms2.mus_p4()[i_ll].y(); 
    } else if (cms2.mus_tcmet_flag()[i_ll] == 4){ 
	 lmetx+= - cms2.mus_tcmet_deltax()[i_ll] + cms2.mus_met_deltax()[i_ll] - cms2.mus_p4()[i_ll].x() + cms2.trks_trk_p4()[cms2.mus_trkidx()[i_ll]].x();  
	 lmety+= - cms2.mus_tcmet_deltay()[i_ll] + cms2.mus_met_deltay()[i_ll] - cms2.mus_p4()[i_ll].y() + cms2.trks_trk_p4()[cms2.mus_trkidx()[i_ll]].y();  
    } 
  }
  met = sqrt(lmetx*lmetx+lmety*lmety);
  metPhi = atan2(lmety,lmetx);

  return;
}
