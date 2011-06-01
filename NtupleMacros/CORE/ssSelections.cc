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
#include "eventSelections.h"

#include "ssSelections.h"


/****************************************************************
***                                                           ***
***                                                           ***
***                                                           ***
***          2011 Selections                                  ***
***                                                           ***
***                                                           ***
***                                                           ***
****************************************************************/

/******************************************************************************************/     
// 2011 good lepton
/******************************************************************************************/     
bool isGoodLepton(int id, int idx, int vidx)
{
    // electrons
    if (abs(id) == 11)
        return (pass_electronSelection(idx, electronSelection_ssV4_noIso, false, false, vidx));

    // muons
    if (abs(id) == 13)
        return (muonIdNotIsolated(idx, NominalSSv4, vidx));

    return false;
}


/******************************************************************************************/     
// 2011 isolated lepton
/******************************************************************************************/     
bool isIsolatedLepton(int id, int idx, int vidx)
{
    // electrons
    if (abs(id) == 11)
        return (pass_electronSelection(idx, electronSelection_ssV4_iso));

    // muons
    if (abs(id) == 13)
        return (muonIsoValue(idx, true) < 0.15);

    return false;
}


/******************************************************************************************/     
// 2011 numerator lepton
/******************************************************************************************/     
bool isNumeratorLepton(int id, int idx, int vidx)
{
    return (isGoodLepton(id, idx, vidx) && isIsolatedLepton(id, idx, vidx));
}


/******************************************************************************************/     
// 2011 numerator hypothesis
/******************************************************************************************/     
bool isNumeratorHypothesis(int idx, int vidx)
{
    if (!isNumeratorLepton(cms2.hyp_lt_id().at(idx), cms2.hyp_lt_index().at(idx), vidx))
        return false;
    if (!isNumeratorLepton(cms2.hyp_ll_id().at(idx), cms2.hyp_ll_index().at(idx), vidx))
        return false;

    return true;
}


/******************************************************************************************/     
// 2011 denominator lepton
/******************************************************************************************/     
bool isDenominatorLepton(int id, int idx, int vidx)
{
    // electrons
    if (abs(id) == 11)
        return (pass_electronSelection(idx, electronSelectionFOV4_ssVBTF80_v3, false, false, vidx));

    // muons
    if (abs(id) == 13)
        return (muonId(idx, muonSelectionFO_ssV4, vidx));

    return false;
}


/*****************************************************************************************/
// extra Z veto
/*****************************************************************************************/
bool makesExtraZ (int hypIdx, bool applyAlignmentCorrection, bool removedEtaCutInEndcap)
{
 	 LorentzVector ltp4	= cms2.hyp_lt_p4()[hypIdx];
	 int           ltid	= abs(cms2.hyp_lt_id()[hypIdx]);
	 int           ltch = cms2.hyp_lt_charge()[hypIdx];

	 LorentzVector llp4	= cms2.hyp_ll_p4()[hypIdx];
	 int           llid	= abs(cms2.hyp_ll_id()[hypIdx]);
	 int           llch = cms2.hyp_ll_charge()[hypIdx];

	 if (ltid == 13 || llid == 13)
	 {
		  for (unsigned int imu = 0; imu < cms2.mus_p4().size(); imu++)
		  {
			   if (cms2.mus_p4()[imu].pt() < 10.)
					continue;

			   if (!isGoodLeptonwIsoSS(13, imu, applyAlignmentCorrection, removedEtaCutInEndcap) || !isGoodLeptonNoIsoSS(13, imu, applyAlignmentCorrection, removedEtaCutInEndcap))
					continue;

               if (ltid == 13 && imu == cms2.hyp_lt_index()[hypIdx])
                   continue;

               if (llid == 13 && imu == cms2.hyp_ll_index()[hypIdx])
                   continue;

			   if (ltid == 13 && ltch * cms2.mus_charge()[imu] < 0)
			   {
					LorentzVector temp = ltp4 + cms2.mus_p4()[imu];
					if (inZmassWindow(sqrt(temp.M2())))
						 return true;
			   }
		  
			   if (llid == 13 && llch * cms2.mus_charge()[imu] < 0)
			   {
					LorentzVector temp = llp4 + cms2.mus_p4()[imu];
					if (inZmassWindow(sqrt(temp.M2())))
						 return true;
			   }
		  }
	 }

	 if (ltid == 11 || llid == 11)
	 {
		  for (unsigned int iel = 0; iel < cms2.els_p4().size(); iel++)
		  {
			   if (cms2.els_p4()[iel].pt() < 10.)
					continue;

			   if (!isGoodLeptonwIsoSS(11, iel, applyAlignmentCorrection, removedEtaCutInEndcap) || !isGoodLeptonNoIsoSS(11, iel, applyAlignmentCorrection, removedEtaCutInEndcap))
					continue;

               if (ltid == 11 && iel == cms2.hyp_lt_index()[hypIdx])
                   continue;

               if (llid == 11 && iel == cms2.hyp_ll_index()[hypIdx])
                   continue;

			   if (ltid == 11 && ltch * cms2.els_charge()[iel] < 0)
			   {
					LorentzVector temp = ltp4 + cms2.els_p4()[iel];
					if (inZmassWindow(sqrt(temp.M2())))
						 return true;
			   }
		  
			   if (llid == 11 && llch * cms2.els_charge()[iel] < 0)
			   {
					LorentzVector temp = llp4 + cms2.els_p4()[iel];
					if (inZmassWindow(sqrt(temp.M2())))
						 return true;
			   }
		  }
	 }

	 return false;
}

/*****************************************************************************************/
// require electron GSF, CTF and SC charges agree
/*****************************************************************************************/
bool passThreeChargeRequirement(int elIdx)
{
	 int trk_idx = cms2.els_trkidx()[elIdx];

	 if (trk_idx >= 0)
	 {
		  if (cms2.els_sccharge()[elIdx] == cms2.els_trk_charge()[elIdx] && cms2.els_trk_charge()[elIdx] == cms2.trks_charge()[trk_idx])			   
			   return true;
	 }

	 return false;
}

/******************************************************************************************/     
// is it a good jet?
/******************************************************************************************/     
bool isGoodJet(LorentzVector jetp4, int hypIdx, double ptCut, double absEtaCut, double dRCut, bool muJetClean) {

	 if(jetp4.Pt() < ptCut)
		  return false;  
	 if(fabs(jetp4.Eta()) > absEtaCut)
		  return false;
  
	 float deltaR_lt = ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[hypIdx], jetp4);
	 float deltaR_ll = ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[hypIdx], jetp4);
	 
	 if (deltaR_lt < dRCut || deltaR_ll < dRCut)
		  return false;

	 for (unsigned int elidx = 0; elidx < cms2.els_p4().size(); elidx++)
	 {
		  if (cms2.els_p4()[elidx].pt() < 10.)
			   continue;
		  if (fabs(cms2.els_p4()[elidx].eta()) > 2.4)
			   continue;

		  if (!isGoodLeptonNoIsoSS(11, elidx))
			   continue;
		  if (!isGoodLeptonwIsoSS(11, elidx))
			   continue;

		  float deltaR = ROOT::Math::VectorUtil::DeltaR(cms2.els_p4()[elidx], jetp4);
		  if (deltaR < dRCut)
			   return false;
	 }

	 if (muJetClean)
	 {
		  for (unsigned int muidx = 0; muidx < cms2.mus_p4().size(); muidx++)
		  {
			   if (cms2.mus_p4()[muidx].pt() < 5.)
					continue;
			   if (fabs(cms2.mus_p4()[muidx].eta()) > 2.4)
					continue;

			   if (!isGoodLeptonNoIsoSS(13, muidx))
					continue;
			   if (!isGoodLeptonwIsoSS(13, muidx))
					continue;

			   float deltaR = ROOT::Math::VectorUtil::DeltaR(cms2.mus_p4()[muidx], jetp4);
			   if (deltaR < dRCut)
					return false;
		  }
	 }

	 return true;
}









/****************************************************************
***                                                           ***
***                                                           ***
***                                                           ***
***          2010 Selections                                  ***
***                                                           ***
***                                                           ***
***                                                           ***
****************************************************************/

/******************************************************************************************/     
// good lepton (either mu or electron, no isolation cuts)
/******************************************************************************************/
bool isGoodLeptonNoIsoSS (int id, int lepIdx, bool applyAlignmentCorrection, bool removedEtaCutInEndcap)
{
	 // electrons
	 if (abs(id) == 11)
		  return (pass_electronSelection(lepIdx, electronSelection_ssV2_NoIso, applyAlignmentCorrection, removedEtaCutInEndcap));

	 // muons
	 if (abs(id) == 13)
		  return (muonIdNotIsolated(lepIdx, NominalSS)) ;

	 return false;
}


/******************************************************************************************/     
// isolated lepton (either mu or electron), no other cuts applied
/******************************************************************************************/
bool isGoodLeptonwIsoSS (int id, int lepIdx, bool applyAlignmentCorrection, bool removedEtaCutInEndcap)
{
	 // electrons
	 if (abs(id)== 11)
		  return (pass_electronSelection(lepIdx, electronSelection_ss_Iso, applyAlignmentCorrection, removedEtaCutInEndcap));

	 // muons
	 if (abs(id) == 13) 
		  return (muonIsoValue(lepIdx) < 0.1);

	 return false;
}


/******************************************************************************************/     
// are the leptons in the hypothesis good (all cuts but isolation?)
/******************************************************************************************/
bool isGoodHypNoIsoSS (int hypIdx, bool applyAlignmentCorrection, bool removedEtaCutInEndcap)
{
	 if(!isGoodLeptonNoIsoSS(cms2.hyp_lt_id()[hypIdx], cms2.hyp_lt_index()[hypIdx], applyAlignmentCorrection, removedEtaCutInEndcap))
		  return false;
	 if(!isGoodLeptonNoIsoSS(cms2.hyp_ll_id()[hypIdx], cms2.hyp_ll_index()[hypIdx], applyAlignmentCorrection, removedEtaCutInEndcap))
		  return false;

	 return true;
}


/******************************************************************************************/     
// are the leptons in the hypothesis isolated? (no other cuts applied)
/******************************************************************************************/     
bool isGoodHypwIsoSS (int hypIdx, bool applyAlignmentCorrection, bool removedEtaCutInEndcap)
{
	 if(!isGoodLeptonwIsoSS(cms2.hyp_lt_id()[hypIdx], cms2.hyp_lt_index()[hypIdx], applyAlignmentCorrection, removedEtaCutInEndcap))
		  return false;
	 if(!isGoodLeptonwIsoSS(cms2.hyp_ll_id()[hypIdx], cms2.hyp_ll_index()[hypIdx], applyAlignmentCorrection, removedEtaCutInEndcap))
		  return false;

	 return true;
}

/*****************************************************************************************/
//passes the EGamma triggers
/*****************************************************************************************/
bool passEGTriggerSS (unsigned int hypIdx, bool mc)
{
	 unsigned int hypType = cms2.hyp_type()[hypIdx];	 

	 if (mc)
	 {
		  if (cms2.passHLTTrigger("HLT_Ele10_LW_L1R"))
			   return true;

		  if (cms2.passHLTTrigger("HLT_Ele10_SW_L1R"))
			   return true;

		  if (cms2.passHLTTrigger("HLT_Ele17_SW_L1R"))
			   return true;

		  if (cms2.passHLTTrigger("HLT_Ele10_LW_EleId_L1R"))
			   return true;

		  if (cms2.passHLTTrigger("HLT_Ele15_LW_L1R"))
			   return true;

		  if (cms2.passHLTTrigger("HLT_DoubleEle5_SW_L1R"))
			   return true;
	 } 
	 else // data now
	 {
		  if(cms2.evt_run() < 138000)
		  {
			   if (cms2.passHLTTrigger("HLT_Ele10_LW_L1R"))
					return true;

			   if (cms2.passHLTTrigger("HLT_Ele15_LW_L1R"))
					return true;

			   if (hypType == 3)
					if (cms2.passHLTTrigger("HLT_DoubleEle5_SW_L1R"))
						 return true;
		  }

		  if(cms2.evt_run() >= 138000 && cms2.evt_run() < 141900)
		  {
			   if (cms2.passHLTTrigger("HLT_Ele10_LW_EleId_L1R"))
					return true;

			   if (cms2.passHLTTrigger("HLT_Ele15_LW_L1R"))
					return true;

			   if (hypType == 3)
					if (cms2.passHLTTrigger("HLT_DoubleEle5_SW_L1R"))
						 return true;
		  }

		  if (cms2.evt_run() >= 141900)
		  {
			   if (cms2.passHLTTrigger("HLT_Ele10_SW_EleId_L1R"))
					return true;

			   if (cms2.passHLTTrigger("HLT_Ele15_SW_CaloEleId_L1R"))
					return true;

			   if (cms2.passHLTTrigger("HLT_Ele15_SW_EleId_L1R"))
					return true;

			   if (cms2.passHLTTrigger("HLT_Ele17_SW_LooseEleId_L1R"))
					return true;

			   if (cms2.passHLTTrigger("HLT_Ele17_SW_CaloEleId_L1R"))
					return true;

			   if (cms2.passHLTTrigger("HLT_Ele17_SW_EleId_L1R"))
					return true;

			   if (cms2.passHLTTrigger("HLT_Ele17_SW_TightCaloEleId_SC8HE_L1R_v1"))
					return true;

			   if (cms2.passHLTTrigger("HLT_Ele17_SW_TightEleIdIsol_L1R_v1"))
					return true;

			   if (hypType == 3)
			   {
					if (cms2.passHLTTrigger("HLT_DoubleEle10_SW_L1R"))
						 return true;

					if (cms2.passHLTTrigger("HLT_DoubleEle15_SW_L1R_v1"))
						 return true;
			   }
		  }
	 }

	 return false;
}

/*****************************************************************************************/
//passes the muon triggers
/*****************************************************************************************/
bool passMuTriggerSS(unsigned int hypIdx, bool mc)
{
	 if (mc)
	 {
		  if (cms2.passHLTTrigger("HLT_Mu9"))
			   return true;

		  if (cms2.passHLTTrigger("HLT_Mu11"))
			   return true;

		  if (cms2.hyp_type()[hypIdx] == 0)
			   if (cms2.passHLTTrigger("HLT_DoubleMu3"))
					return true;
	 }
	 else {
		  if (cms2.passHLTTrigger("HLT_Mu9"))
			   return true;

		  if (cms2.passHLTTrigger("HLT_Mu11"))
			   return true;

		  if (cms2.passHLTTrigger("HLT_Mu15_v1"))
			   return true;

		  if (cms2.hyp_type()[hypIdx] == 0)
		  {
			   if (cms2.passHLTTrigger("HLT_DoubleMu3"))
					return true;
			   if (cms2.passHLTTrigger("HLT_DoubleMu5_v1"))
					return true;
		  }
	 
		  if (cms2.hyp_type()[hypIdx] == 1 || cms2.hyp_type()[hypIdx] == 2)
		  {
			   if (cms2.passHLTTrigger("HLT_Mu5_Ele9_v1"))
					return true;
			   if (cms2.passHLTTrigger("HLT_Mu8_Ele8_v1"))
					return true;
		  }
	 }

	 return false;
}

/******************************************************************************************/     
// good lepton (either mu or electron, no isolation cuts), d0 corrected for PV
/******************************************************************************************/
bool isGoodLeptonNoIsoSSd0PV (int id, int lepIdx, bool applyAlignmentCorrection, bool removedEtaCutInEndcap)
{
	 // electrons
	 if (abs(id) == 11)
		  return (pass_electronSelection(lepIdx, electronSelection_ssV2d0PV_NoIso, applyAlignmentCorrection, removedEtaCutInEndcap));

	 // muons
	 if (abs(id) == 13)
		  return (muonIdNotIsolated(lepIdx, NominalSSd0PV)) ;

	 return false;
}


/******************************************************************************************/     
// are the leptons in the hypothesis good (all cuts but isolation?), d0 corrected for PV
/******************************************************************************************/
bool isGoodHypNoIsoSSd0PV (int hypIdx, bool applyAlignmentCorrection, bool removedEtaCutInEndcap)
{
	 if(!isGoodLeptonNoIsoSSd0PV(cms2.hyp_lt_id()[hypIdx], cms2.hyp_lt_index()[hypIdx], applyAlignmentCorrection, removedEtaCutInEndcap))
		  return false;
	 if(!isGoodLeptonNoIsoSSd0PV(cms2.hyp_ll_id()[hypIdx], cms2.hyp_ll_index()[hypIdx], applyAlignmentCorrection, removedEtaCutInEndcap))
		  return false;

	 return true;
}




/******************************************************************************************/     
// good lepton (either mu or electron, no isolation cuts), no d0 cut
/******************************************************************************************/
bool isGoodLeptonNoIsoSSnod0 (int id, int lepIdx, bool applyAlignmentCorrection, bool removedEtaCutInEndcap)
{
	 // electrons
	 if (abs(id) == 11)
		  return (pass_electronSelection(lepIdx, electronSelection_ssV2nod0_NoIso, applyAlignmentCorrection, removedEtaCutInEndcap));

	 // muons
	 if (abs(id) == 13)
		  return (muonIdNotIsolated(lepIdx, NominalSSnod0)) ;

	 return false;
}


/******************************************************************************************/     
// are the leptons in the hypothesis good (all cuts but isolation?), no d0 cut
/******************************************************************************************/
bool isGoodHypNoIsoSSnod0 (int hypIdx, bool applyAlignmentCorrection, bool removedEtaCutInEndcap)
{
	 if(!isGoodLeptonNoIsoSSnod0(cms2.hyp_lt_id()[hypIdx], cms2.hyp_lt_index()[hypIdx], applyAlignmentCorrection, removedEtaCutInEndcap))
		  return false;
	 if(!isGoodLeptonNoIsoSSnod0(cms2.hyp_ll_id()[hypIdx], cms2.hyp_ll_index()[hypIdx], applyAlignmentCorrection, removedEtaCutInEndcap))
		  return false;

	 return true;
}




/******************************************************************************************/     
// good lepton (either mu or electron, no isolation cuts), no triple charge consistency requirement
/******************************************************************************************/
bool isGoodLeptonNoIsoSSnoTripleChargeReq (int id, int lepIdx, bool applyAlignmentCorrection, bool removedEtaCutInEndcap)
{
	 // electrons
	 if (abs(id) == 11)
		  return (pass_electronSelection(lepIdx, electronSelection_ssV2noTripleChargeReq_NoIso, applyAlignmentCorrection, removedEtaCutInEndcap));

	 // muons
	 if (abs(id) == 13)
		  return (muonIdNotIsolated(lepIdx, NominalSS)) ;

	 return false;
}


/******************************************************************************************/     
// are the leptons in the hypothesis good (all cuts but isolation?), no triple charge consistency requirement
/******************************************************************************************/
bool isGoodHypNoIsoSSnoTripleChargeReq (int hypIdx, bool applyAlignmentCorrection, bool removedEtaCutInEndcap)
{
	 if(!isGoodLeptonNoIsoSSnoTripleChargeReq(cms2.hyp_lt_id()[hypIdx], cms2.hyp_lt_index()[hypIdx], applyAlignmentCorrection, removedEtaCutInEndcap))
		  return false;
	 if(!isGoodLeptonNoIsoSSnoTripleChargeReq(cms2.hyp_ll_id()[hypIdx], cms2.hyp_ll_index()[hypIdx], applyAlignmentCorrection, removedEtaCutInEndcap))
		  return false;

	 return true;
}


/******************************************************************************************/     
// is muon a fakeable object
/******************************************************************************************/
bool isFakeableMuon(int index, bool useTightIsolation)
{
	 if (useTightIsolation)
		  return muonId(index, muonSelectionFO_mu_ss);
	 else
		  return muonId(index, muonSelectionFO_mu_ss_iso10);
}
