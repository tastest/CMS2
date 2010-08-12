#include <assert.h>
#include <algorithm>
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TDatabasePDG.h"
#include "electronSelections.h"
#include "electronSelectionsParameters.h"
#include "muonSelections.h"
#include "metSelections.h"
#include "triggerUtils.h"
#include "CMS2.h"

#include "ttbarSelections.h"

using namespace tas;

/******************************************************************************************/     
// good lepton (either mu or electron, no isolation cuts)
/******************************************************************************************/
bool isGoodLeptonNoIso(int id, int lepIdx, bool applyAlignmentCorrection, bool removedEtaCutInEndcap) {

  //electrons
  if(abs(id) == 11) {
    
    const cuts_t elIDcuts =     (1ll<<ELEIP_400)                |                                                                  
      (1ll<<ELENOMUON_010)              |
      (1ll<<ELENOTCONV_HITPATTERN)      |
      (1ll<<ELENOTCONV_DISTDCOT002)     |
      (1ll<<ELESCET_010)                |
      (1ll<<ELEPT_020)          |
      (1ll<<ELEETA_250)         |
      (1ll<<ELENOSPIKE_SWISS005) |
      (1ll<<ELESEED_ECAL);
    
    unsigned int answer_vbtf = electronId_VBTF(lepIdx, VBTF_35X_90, applyAlignmentCorrection, removedEtaCutInEndcap);
    bool elsvbtf90_ = ( ( answer_vbtf & (1ll<<ELEID_ID) ) == (1ll<<ELEID_ID) );
    return (pass_electronSelection(lepIdx, elIDcuts, applyAlignmentCorrection, removedEtaCutInEndcap) && elsvbtf90_);
  }

  //muons
  if(abs(id) == 13) {

    if(muonIdNotIsolated(lepIdx, NominalTTbarV2)) 
      return true;

    return false;
      
  }
  

  return false;
}


/******************************************************************************************/     
// isolated lepton (either mu or electron)
/******************************************************************************************/
bool isGoodLeptonwIso(int id, int lepIdx, bool applyAlignmentCorrection, bool removedEtaCutInEndcap) {

  if(!isGoodLeptonNoIso(id, lepIdx, applyAlignmentCorrection, removedEtaCutInEndcap))
    return false;

  // 11 is a electron
  if(abs(id)== 11) {
       const cuts_t elISOcuts =   (1ll<<ELEISO_REL015);
       if (!pass_electronSelection(lepIdx, elISOcuts, applyAlignmentCorrection, removedEtaCutInEndcap))
            return false;
  }

  // 13 is a muon
  if(abs(id) == 13) 
    if(muonIsoValue(lepIdx) > 0.15)  return false;
        
  return true;
}


/******************************************************************************************/     
// are the leptons in the hypothesis good (all cuts but isolation?)
/******************************************************************************************/
bool isGoodHypNoIso(int hypIdx, bool applyAlignmentCorrection, bool removedEtaCutInEndcap) {
 
  if(!isGoodLeptonNoIso(hyp_lt_id()[hypIdx], hyp_lt_index()[hypIdx], applyAlignmentCorrection, removedEtaCutInEndcap))//, used0wrtPV)
     return false;
  if(!isGoodLeptonNoIso(hyp_ll_id()[hypIdx], hyp_ll_index()[hypIdx], applyAlignmentCorrection, removedEtaCutInEndcap))//, used0wrtPV)
    return false;

  return true;
}


/******************************************************************************************/     
// are the leptons in the hypothesis isolated?
/******************************************************************************************/     
bool isGoodHypwIso(int hypIdx, bool applyAlignmentCorrection, bool removedEtaCutInEndcap) {

  if(!isGoodLeptonwIso(hyp_lt_id()[hypIdx], hyp_lt_index()[hypIdx], applyAlignmentCorrection, removedEtaCutInEndcap))
    return false;
  if(!isGoodLeptonwIso(hyp_ll_id()[hypIdx], hyp_ll_index()[hypIdx], applyAlignmentCorrection, removedEtaCutInEndcap))
    return false;


  return true;
}


/******************************************************************************************/     
// is it a good jet?
/******************************************************************************************/     
bool isGoodDilHypJet(LorentzVector jetp4, unsigned int& hypIdx, double ptCut, double absEtaCut, double dRCut, bool muJetClean) {

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
//return the MET and the MET phi, correcting for mus that are not corrected for by default
/******************************************************************************************/     
std::pair<float,float> getMet(const string algo, unsigned int hypIdx) {
  
  if(algo != "tcMET" && algo != "muCorMET" && algo != "pfMET" && algo != "tcMET35X") {
    cout << algo << "IS NOT A RECOGNIZED MET ALGORITHM!!!!! PLEASE CHECK YOUR CODE!!!";
    return make_pair(-99999., -99999.);
  }

  
  if(algo == "tcMET") {

    float tcmetX = evt_tcmet()*cos(evt_tcmetPhi());
    float tcmetY = evt_tcmet()*sin(evt_tcmetPhi());
    
    if(abs(hyp_lt_id()[hypIdx]) == 13)
      fixMetForThisMuon(hyp_lt_index().at(hypIdx), tcmetX, tcmetY, usingTcMet);
    if(abs(hyp_ll_id()[hypIdx]) == 13)
      fixMetForThisMuon(hyp_ll_index().at(hypIdx), tcmetX, tcmetY, usingTcMet);
    
    return make_pair(sqrt(tcmetX*tcmetX + tcmetY*tcmetY), atan2(tcmetY, tcmetX));
  }


  if(algo == "tcMET35X") {

    float tcmetX = evt35X_tcmet()*cos(evt35X_tcmetPhi());
    float tcmetY = evt35X_tcmet()*sin(evt35X_tcmetPhi());
    
    if(abs(hyp_lt_id()[hypIdx]) == 13)
      fixMetForThisMuon(hyp_lt_index().at(hypIdx), tcmetX, tcmetY, usingTcMet35X);
    if(abs(hyp_ll_id()[hypIdx]) == 13)
      fixMetForThisMuon(hyp_ll_index().at(hypIdx), tcmetX, tcmetY, usingTcMet35X);
    
    return make_pair(sqrt(tcmetX*tcmetX + tcmetY*tcmetY), atan2(tcmetY, tcmetX));
  }

  
  if(algo == "muCorMET") {

    float metX = evt_metMuonCorr()*cos(evt_metMuonCorrPhi());
    float metY = evt_metMuonCorr()*sin(evt_metMuonCorrPhi());
    
    if(abs(hyp_lt_id()[hypIdx]) == 13)
      fixMetForThisMuon(hyp_lt_index().at(hypIdx), metX, metY, usingCaloMet);
    if(abs(hyp_ll_id()[hypIdx]) == 13)
      fixMetForThisMuon(hyp_ll_index().at(hypIdx), metX, metY, usingCaloMet);

    return make_pair(sqrt(metX*metX + metY*metY), atan2(metY, metX));
  }
  
  //nothing to do here because they're perfect
  if(algo == "pfMET") 
    return make_pair(evt_pfmet(), evt_pfmetPhi());
  
  
  return make_pair(-99999., -99999);
  
}


/*****************************************************************************************/
//hypothesis disambiguation. Returns the hypothesis that has the highest sum Pt
/*****************************************************************************************/
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



/*****************************************************************************************/
//passes the EGamma triggers
/*****************************************************************************************/
bool passEGTrigger(bool mc) {

  if (mc) {
 
    int e10 = nHLTObjects("HLT_Ele10_LW_L1R");
    for (int i=0; i<e10; i++) {
      LorentzVector p4 = p4HLTObject("HLT_Ele10_LW_L1R", i);
      if(p4.Pt() > 15.) return true;
    }

  } else {  // data now

    
    if(evt_run() < 138000) {
      int e10 = nHLTObjects("HLT_Ele10_LW_L1R");
      for (int i=0; i<e10; i++) {
        LorentzVector p4 = p4HLTObject("HLT_Ele10_LW_L1R", i);
        if(p4.Pt() > 15.) return true;
      }
    }

    if(evt_run() >= 138000 && evt_run() < 141900) {
      int e15 = nHLTObjects("HLT_Ele15_LW_L1R");
      if(e15 != 0) 
        return true;
    }

    if(evt_run() >= 141900) {
      int e15 = nHLTObjects("HLT_Ele15_SW_L1R");
      if(e15 != 0) 
        return true;
    }

 
  }
  return false;
}

/*****************************************************************************************/
//get the impact parameter wrt the PV
/*****************************************************************************************/
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








