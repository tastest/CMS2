#include <assert.h>
#include <algorithm>
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TDatabasePDG.h"
#include "electronSelections.cc"
#include "electronSelectionsParameters.cc"
#include "muonSelections.cc"
#include "metSelections.cc"
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
      (1ll<<ELESCET_015)                |
      (1ll<<ELEPT_020)          |
      (1ll<<ELEETA_250)         |      
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
  
  if(algo != "tcMET" && algo != "muCorMET" && algo != "pfMET" && algo != "tcMET35X" && algo != "tcMET_looper") {
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
/*
  if(algo == "tcMET_looper") {

    metStruct tcmetStruct = correctedTCMET();
    float tcmetX = tcmetStruct.metx;
    float tcmetY = tcmetStruct.mety;
    
    if(abs(hyp_lt_id()[hypIdx]) == 13)
      fixMetForThisMuon(hyp_lt_index().at(hypIdx), tcmetX, tcmetY, usingTcMet_looper);
    if(abs(hyp_ll_id()[hypIdx]) == 13)
      fixMetForThisMuon(hyp_ll_index().at(hypIdx), tcmetX, tcmetY, usingTcMet_looper);
    
    return make_pair(sqrt(tcmetX*tcmetX + tcmetY*tcmetY), atan2(tcmetY, tcmetX));
  }
*/

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
      if (p4.Pt() > 15.) return true;
    }

  } else {  // data now

    
    if (cms2.evt_run() < 138000) {
      int e10 = nHLTObjects("HLT_Ele10_LW_L1R");
      for (int i=0; i<e10; i++) {
        LorentzVector p4 = p4HLTObject("HLT_Ele10_LW_L1R", i);
        if (p4.Pt() > 15.) return true;
      }
    }

    if (cms2.evt_run() >= 138000 && cms2.evt_run() < 141900) {
      int e15 = nHLTObjects("HLT_Ele15_LW_L1R");
      if (e15 != 0) 
        return true;
    }

    if (cms2.evt_run() >= 141900 && cms2.evt_run() <= 144000) {
      int e15 = nHLTObjects("HLT_Ele15_SW_L1R");
      if (e15 != 0) 
	return true;
    }   
    
    if (cms2.evt_run() > 144000 && cms2.evt_run() <= 144114) {

      int e15caloId = nHLTObjects("HLT_Ele15_SW_CaloEleId_L1R");
      if (e15caloId != 0)
	return true;
      
      int e20 = nHLTObjects("HLT_Ele20_SW_L1R");
      if (e20 != 0)
	return true;
      
      int ed10 = nHLTObjects("HLT_DoubleEle10_SW_L1R");
      if(ed10 != 0)
	return true;
    }

    if(cms2.evt_run() > 146000 && cms2.evt_run() <= 147120) {

      if(nHLTObjects("HLT_DoubleEle10_SW_L1R") != 0 )
	return true;
      
      if(nHLTObjects("HLT_Ele17_SW_CaloEleId_L1R") != 0)
	return true;
    }

    if(cms2.evt_run() > 147120 && cms2.evt_run() <= 148100 ) {
    
      if(nHLTObjects("HLT_DoubleEle15_SW_L1R_v1") != 0)
	return true;

      if(nHLTObjects("HLT_Ele17_SW_TightCaloEleId_SC8HE_L1R_v1") != 0)
	return true;
      
      if(nHLTObjects("HLT_Ele17_SW_TightEleId_L1R") != 0)
	return true;

    }

    if(cms2.evt_run()  > 148100) {
      if(nHLTObjects("HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v2") != 0)
	return true;

      if(nHLTObjects("HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v1") != 0) 
	return true;
      
      if(nHLTObjects("HLT_Ele22_SW_TighterEleId_L1R_v3") != 0)
	return true;

      if(nHLTObjects("HLT_Ele22_SW_TighterEleId_L1R_v2") != 0)
	return true;
      
      if(nHLTObjects("HLT_Ele17_SW_TighterEleIdIsol_L1R_v3") != 0) 
	return true;

      if(nHLTObjects("HLT_Ele17_SW_TighterEleIdIsol_L1R_v2") != 0) 
	return true;
    }
    
  return false;
  
  }
}

/*****************************************************************************************/
//passes the EGamma triggers
/*****************************************************************************************/
bool passMuTrigger(bool mc) {
  

  if(mc) {
    if(passHLTTrigger("HLT_Mu9"))
      return true;

    
  } else {

    if(cms2.evt_run() <= 145000) {
      if(nHLTObjects("HLT_Mu9") != 0)
	return true;
    }

    if(cms2.evt_run() > 145000 && cms2.evt_run() <= 147120) {
      if(nHLTObjects("HLT_Mu11") != 0)
	return true;
    }
  
    if(cms2.evt_run() > 147120) {
      if(nHLTObjects("HLT_Mu15_v1") != 0)
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

/*****************************************************************************************/
//get the number of jets passing btag discriminator cuts
// takes as arguments a vector of indices, the jet algorithm
// and the btag discriminator
// the working points are hard coded
/*****************************************************************************************/
int getNbtags(const vector<unsigned int> v_jetsIdx, const string jetAlgo, const string bTagDiscriminator) {

  if(jetAlgo != "jptJets" && jetAlgo != "caloJets" && jetAlgo != "pfJets") {
    cout << "Unknown jet algorithm. Returning spurious value" << endl;
    return -9999;
  }

  if(bTagDiscriminator != "trackCountingHighEffBJetTag" &&
     bTagDiscriminator != "simpleSecondaryVertexHighEffBJetTag" && 
     bTagDiscriminator != "simpleSecondaryVertexHighPurBJetTag") {
    cout << "Unknown bTag Discriminator. Returning spurious value" << endl;
    return -9999;
  }

  if(jetAlgo == "jptJets") {
    int ntags = 0;
    if(bTagDiscriminator == "trackCountingHighEffBJetTag") {
      for(unsigned int i = 0; i < v_jetsIdx.size(); i++) {
	if(jets_trackCountingHighEffBJetTag()[v_jetsIdx.at(i)] > 1.7)
	  ntags++;
      }
      return ntags;
    }else if(bTagDiscriminator == "simpleSecondaryVertexHighEffBJetTag") {
      for(unsigned int i = 0; i < v_jetsIdx.size(); i++) {
	if(jets_simpleSecondaryVertexHighEffBJetTag()[v_jetsIdx.at(i)] > 1.74)
	  ntags++;
      }
      return ntags;
    } else if(bTagDiscriminator == "simpleSecondaryVertexHighPurBJetTag") {
      for(unsigned int i = 0; i < v_jetsIdx.size(); i++) {
	if(jets_simpleSecondaryVertexHighPurBJetTags()[v_jetsIdx.at(i)] > 2)
	  ntags++;
      }
      return ntags;
    }
  }


  if(jetAlgo == "caloJets") {		
    int ntags = 0;
    if(bTagDiscriminator == "trackCountingHighEffBJetTag") {
      for(unsigned int i = 0; i < v_jetsIdx.size(); i++) {
	if(jets_trackCountingHighEffBJetTag()[v_jetsIdx.at(i)] > 1.7)
	  ntags++;
      }
      return ntags;
    }else if(bTagDiscriminator == "simpleSecondaryVertexHighEffBJetTag") {
      for(unsigned int i = 0; i < v_jetsIdx.size(); i++) {
	if(jets_simpleSecondaryVertexHighEffBJetTag()[v_jetsIdx.at(i)] > 1.74)
	  ntags++;
      }
      return ntags;
    } else if(bTagDiscriminator == "simpleSecondaryVertexHighPurBJetTag") {
      for(unsigned int i = 0; i < v_jetsIdx.size(); i++) {
	if(jets_simpleSecondaryVertexHighPurBJetTags()[v_jetsIdx.at(i)] > 2)
	  ntags++;
      }
      return ntags;
    }
  }


  
  if(jetAlgo == "pfJets") {		
    int ntags = 0;
    if(bTagDiscriminator == "trackCountingHighEffBJetTag") {
      for(unsigned int i = 0; i < v_jetsIdx.size(); i++) {
	if(pfjets_trackCountingHighEffBJetTag()[v_jetsIdx.at(i)] > 1.7)
	  ntags++;
      }
      return ntags;
    }else if(bTagDiscriminator == "simpleSecondaryVertexHighEffBJetTag") {
      for(unsigned int i = 0; i < v_jetsIdx.size(); i++) {
	if(pfjets_simpleSecondaryVertexHighEffBJetTag()[v_jetsIdx.at(i)] > 1.74)
	  ntags++;
      }
      return ntags;
    } else if(bTagDiscriminator == "simpleSecondaryVertexHighPurBJetTag") {
      for(unsigned int i = 0; i < v_jetsIdx.size(); i++) {
	if(pfjets_simpleSecondaryVertexHighPurBJetTags()[v_jetsIdx.at(i)] > 2)
	  ntags++;
      }
      return ntags;
    }
  }


  
  //should not be here
  std::string name;
  name = system("whoami");
  cout << "You're somewhere you're not supposed to be. Naught " << name << "!" << endl;
  return -9999;

}

      
