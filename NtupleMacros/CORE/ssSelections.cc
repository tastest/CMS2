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
#include "CMS2.h"

#include "ssSelections.h"

using namespace tas;

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
// isolated lepton (either mu or electron)
/******************************************************************************************/
bool isGoodLeptonwIsoSS (int id, int lepIdx, bool applyAlignmentCorrection, bool removedEtaCutInEndcap)
{
	 // electrons
	 if (abs(id)== 11)
		  return (pass_electronSelection(lepIdx, electronSelection_ssV2, applyAlignmentCorrection, removedEtaCutInEndcap));

	 // muons
	 if (abs(id) == 13) 
		  return (muonId(lepIdx, NominalSS));

	 return false;
}


/******************************************************************************************/     
// are the leptons in the hypothesis good (all cuts but isolation?)
/******************************************************************************************/
bool isGoodHypNoIsoSS (int hypIdx, bool applyAlignmentCorrection, bool removedEtaCutInEndcap)
{
	 if(!isGoodLeptonNoIsoSS(hyp_lt_id()[hypIdx], hyp_lt_index()[hypIdx], applyAlignmentCorrection, removedEtaCutInEndcap))
		  return false;
	 if(!isGoodLeptonNoIsoSS(hyp_ll_id()[hypIdx], hyp_ll_index()[hypIdx], applyAlignmentCorrection, removedEtaCutInEndcap))
		  return false;

	 return true;
}


/******************************************************************************************/     
// are the leptons in the hypothesis isolated?
/******************************************************************************************/     
bool isGoodHypwIsoSS (int hypIdx, bool applyAlignmentCorrection, bool removedEtaCutInEndcap)
{
	 if(!isGoodLeptonwIsoSS(hyp_lt_id()[hypIdx], hyp_lt_index()[hypIdx], applyAlignmentCorrection, removedEtaCutInEndcap))
		  return false;
	 if(!isGoodLeptonwIsoSS(hyp_ll_id()[hypIdx], hyp_ll_index()[hypIdx], applyAlignmentCorrection, removedEtaCutInEndcap))
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
		  if (passHLTTrigger("HLT_Ele10_LW_L1R"))
			   return true;

		  if (passHLTTrigger("HLT_Ele10_LW_EleId_L1R"))
			   return true;

		  if (passHLTTrigger("HLT_Ele15_LW_L1R"))
			   return true;

		  if (passHLTTrigger("HLT_DoubleEle5_SW_L1R"))
			   return true;
	 } 
	 else // data now
	 {
		  if(evt_run() < 138000)
		  {
			   if (passHLTTrigger("HLT_Ele10_LW_L1R"))
					return true;

			   if (passHLTTrigger("HLT_Ele15_LW_L1R"))
					return true;

			   if (hypType == 3)
					if (passHLTTrigger("HLT_DoubleEle5_SW_L1R"))
						 return true;
		  }

		  if(evt_run() >= 138000 && evt_run() < 141900)
		  {
			   if (passHLTTrigger("HLT_Ele10_LW_EleId_L1R"))
					return true;

			   if (passHLTTrigger("HLT_Ele15_LW_L1R"))
					return true;

			   if (hypType == 3)
					if (passHLTTrigger("HLT_DoubleEle5_SW_L1R"))
						 return true;
		  }

		  if (cms2.evt_run() >= 141900)
		  {
			   if (passHLTTrigger("HLT_Ele10_SW_EleId_L1R"))
					return true;

			   if (passHLTTrigger("HLT_Ele15_SW_CaloEleId_L1R"))
					return true;

			   if (passHLTTrigger("HLT_Ele15_SW_EleId_L1R"))
					return true;

			   if (hypType == 3)
					if (passHLTTrigger("HLT_DoubleEle10_SW_L1R"))
						 return true;
		  }
	 }

	 return false;
}

/*****************************************************************************************/
//passes the muon triggers
/*****************************************************************************************/
bool passMuTriggerSS(unsigned int hypIdx)
{
	 if (passHLTTrigger("HLT_Mu9"))
		  return true;

	 if (cms2.hyp_type()[hypIdx] == 0)
		  if (passHLTTrigger("HLT_DoubleMu3"))
			   return true;

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

			   if (!isGoodLeptonwIsoSS(13, imu, applyAlignmentCorrection, removedEtaCutInEndcap))
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

			   if (!isGoodLeptonwIsoSS(11, iel, applyAlignmentCorrection, removedEtaCutInEndcap))
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
bool isGoodJet(LorentzVector jetp4, double ptCut, double absEtaCut, double dRCut, bool muJetClean) {

	 if(jetp4.Pt() < ptCut)
		  return false;  
	 if(fabs(jetp4.Eta()) > absEtaCut)
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
			   if (cms2.mus_p4()[muidx].pt() < 10.)
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
