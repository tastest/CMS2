
#include <assert.h>
#include <algorithm>
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TDatabasePDG.h"
#include "../CORE/electronSelections.h"
#include "../CORE/muonSelections.h"
#include "../CORE/metSelections.h"
#include "CMS2.h"

#include "ttbarSelections.h"

using namespace tas;

/******************************************************************************************/     
// good lepton (either mu or electron, no isolation cuts)
/******************************************************************************************/
bool isGoodLeptonNoIso(int id, int lepIdx) {

  if(abs(id) == 11) {
    if (!pass_electronSelection(lepIdx, electronSelection_ttbar)) return false;
  }

  if(abs(id) == 13) 
    if(!muonIdNotIsolated(lepIdx, NominalTTbar))    return false;
  

  return true;
}

/******************************************************************************************/     
// isolated lepton (either mu or electron)
/******************************************************************************************/
bool isGoodLeptonwIso(int id, int lepIdx) {

  if(!isGoodLeptonNoIso(id, lepIdx)) return false;

  // 11 is a electron
  if(abs(id)== 11) {
    if (!pass_electronSelection(lepIdx, electronSelection_ttbar_iso)) return false;
  }
  // 13 is a muon
  if(abs(id) == 13) 
    if(!muonId(lepIdx, NominalTTbar))               return false;
  
  return true;
}

/******************************************************************************************/     
// are the leptons in the hypothesis good (all cuts but isolation?)
/******************************************************************************************/
bool isGoodHypNoIso(int hypIdx) {
  
  if(!isGoodLeptonNoIso(hyp_lt_id()[hypIdx], hyp_lt_index()[hypIdx]))
     return false;
  if(!isGoodLeptonNoIso(hyp_ll_id()[hypIdx], hyp_ll_index()[hypIdx]))
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
//return the MET and the MET phi instead of a bool because the MT2 needs it
/******************************************************************************************/     
std::pair<float,float> getMet(string& algo) {
  
  if(algo != "tcMET" && algo != "muCorMET" && algo != "pfMET") {
    cout << algo << "IS NOT A RECOGNIZED MET ALGORITHM!!!!! PLEASE CHECK YOUR CODE!!!";
    return make_pair(-99999., -99999.);
  }
  if(algo == "tcMET") {
    metStruct cortcmet = correctedTCMET();
    return make_pair(cortcmet.met, cortcmet.metphi);
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

    //float isoTk_lt = leptonTrkIsolationTTDil08(id_lt, i_lt);
    //float isoTk_ll = leptonTrkIsolationTTDil08(id_ll, i_ll);
    float isoTk_lt, isoTk_ll;
    float isoCal_lt, isoCal_ll;
    if(abs(id_lt) == 11) {
      isoTk_lt = els_p4()[i_lt].Pt()/(els_tkJuraIso()[i_lt]+els_p4()[i_lt].Pt());
      isoCal_lt = els_p4()[i_lt].Pt()/(els_hcalIso()[i_lt]+els_ecalJuraIso()[i_lt]+els_p4()[i_lt].Pt());
    } else {
      isoTk_lt = mus_p4()[i_lt].Pt()/(mus_iso03_sumPt()[i_lt]+mus_p4()[i_lt].Pt());
      isoCal_lt = mus_p4()[i_lt].Pt()/(mus_iso03_emEt()[i_lt]+mus_iso03_hadEt()[i_lt]+mus_p4()[i_lt].Pt());
    }
    if(abs(id_ll) == 11) {
      isoTk_ll = els_p4()[i_ll].Pt()/(els_tkJuraIso()[i_ll]+els_p4()[i_ll].Pt());
      isoCal_ll = els_p4()[i_ll].Pt()/(els_hcalIso()[i_ll]+els_ecalJuraIso()[i_ll]+els_p4()[i_ll].Pt());
    } else {
      isoTk_ll = mus_p4()[i_ll].Pt()/(mus_iso03_sumPt()[i_ll]+mus_p4()[i_ll].Pt());
      isoCal_ll = mus_p4()[i_ll].Pt()/(mus_iso03_emEt()[i_ll]+mus_iso03_hadEt()[i_ll]+mus_p4()[i_ll].Pt());
    }

  
    
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
    if(pt_lt > 20.)
      hypWeight_lt += (1. - 20./pt_lt*20./pt_lt);
    else
      hypWeight_lt += (1. - 10./pt_lt*10./pt_lt);
    if(pt_ll > 20.)
      hypWeight_ll += (1. - 20./pt_ll*20./pt_ll);
    else
      hypWeight_ll += (1. - 10./pt_ll*10./pt_ll);
    
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

  result = maxWeightIndex;
  return result;
}

//
// return the hyp index of the pair with the highest sum pt
//

unsigned int eventDilIndexByWeightTTDil10(const std::vector<unsigned int> &hyp_index_selected)
{

    float           max_pt = 0.0;
    float           this_pt = 0.0;
    unsigned int    max_index = 0;

    for (unsigned int i = 0; i < hyp_index_selected.size(); ++i) {

        // get the pt of this lepton pair
        this_pt = cms2.hyp_ll_p4()[hyp_index_selected[i]].Pt() 
                    + cms2.hyp_lt_p4()[hyp_index_selected[i]].Pt();

        // if the pt of this lepton pair is 
        // larger than the previous largest
        // then remember its index in the hyp block
        if (this_pt > max_pt) {
            max_pt = this_pt;
            max_index = hyp_index_selected[i];
        }   
    }

    return max_index;

}



/******************************************************************************************/
// MET cut - MET > 20 (OF) MET > 30 (SF)
/******************************************************************************************/
bool passMetAsIs_OF20_SF30(float met, int hyp_type)
{
  if  (hyp_type == 0 || hyp_type == 3) {
    if (met < 30) return false;
  }
  
  if (hyp_type == 1 || hyp_type == 2) {
    if (met < 20) return false;
  }
  return true;
}  

