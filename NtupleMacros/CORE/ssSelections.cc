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
#include "utilities.h"

#include "ssSelections.h"

using namespace samesign;

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
bool samesign::isGoodLepton(int id, int idx, int vidx)
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
bool samesign::isIsolatedLepton(int id, int idx, int vidx)
{
    // electrons
    if (abs(id) == 11)
        return (pass_electronSelection(idx, electronSelection_ssV4_iso));

    // muons
    if (abs(id) == 13)
        return (muonIsoValue(idx, false) < 0.15);

    return false;
}

/******************************************************************************************/     
// 2011 numerator lepton
/******************************************************************************************/     
bool samesign::isNumeratorLepton(int id, int idx, int vidx)
{
  return (isGoodLepton(id, idx, vidx) && isIsolatedLepton(id, idx, vidx));
}


/******************************************************************************************/     
// 2011 numerator hypothesis
/******************************************************************************************/     
bool samesign::isNumeratorHypothesis(int idx, int vidx)
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
bool samesign::isDenominatorLepton(int id, int idx, int vidx)
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
bool samesign::overlapsOtherNNHypInZ(int idx){
  int vidx = -1; // use default now // hypsFromSameVtx2011(idx,1.0,true, false);
  bool result = false;
  int nHyps = cms2.hyp_lt_p4().size();
  for (int iH = 0; iH< nHyps; ++iH){
    if (iH == idx || !hypsOverlap(idx,iH)
	|| abs(cms2.hyp_lt_id()[iH])!= abs(cms2.hyp_ll_id()[iH]) || cms2.hyp_lt_id()[iH]*cms2.hyp_ll_id()[iH] > 0) 
      continue;
    if (! isNumeratorHypothesis(iH,vidx)) continue;
    if (cms2.hyp_p4()[iH].mass2()>0 && fabs(cms2.hyp_p4()[iH].mass() - 91)< 15){ 
      result = true; break;
    }
  }
  return result;
}

/*****************************************************************************************/
// require electron GSF, CTF and SC charges agree
/*****************************************************************************************/
bool samesign::passThreeChargeRequirement(int elIdx)
{
	 int trk_idx = cms2.els_trkidx()[elIdx];

	 if (trk_idx >= 0)
	 {
		  if (cms2.els_sccharge()[elIdx] == cms2.els_trk_charge()[elIdx] && cms2.els_trk_charge()[elIdx] == cms2.trks_charge()[trk_idx])			   
			   return true;
	 }

	 return false;
}

/*****************************************************************************************/
// get sumpt, skip jets overlapping with numerator e/mu with pt>x (defaults are 10/5 GeV)
/*****************************************************************************************/
float samesign::sumJetPt(int idx, enum JetType type, double deltaR, double min_pt, double max_eta, double mu_minpt, double ele_minpt){
  vector<LorentzVector> selJets = getJets(idx, true, type, JETS_CLEAN_NONE, deltaR, min_pt, max_eta);
  int nJets = selJets.size();
  if (nJets == 0 ) return 0.;

  int ivtx = -1;// hypsFromSameVtx2011(idx,1.0, true, false);

  vector<int> jetFlags(nJets, 0);

  int nMus = cms2.mus_type().size();
  for (int iMu = 0; iMu < nMus; ++iMu){    
    if (cms2.mus_p4()[iMu].pt() > mu_minpt && isNumeratorLepton(13,iMu,ivtx)){
      for (int iJ = 0; iJ < nJets; ++iJ){
	float dR = ROOT::Math::VectorUtil::DeltaR(cms2.mus_p4()[iMu], selJets[iJ]);
	if (dR < deltaR) jetFlags[iJ]++;
      }
    }
  }

  int nEls = cms2.els_type().size();
  for (int iEl = 0; iEl < nEls; ++iEl){
    if (cms2.els_p4()[iEl].pt() > ele_minpt && isNumeratorLepton(11,iEl,ivtx)){
      for (int iJ = 0; iJ < nJets; ++iJ){
	float dR = ROOT::Math::VectorUtil::DeltaR(cms2.els_p4()[iEl], selJets[iJ]);
	if (dR < deltaR) jetFlags[iJ]++;
      }
    }
  }

  float result = 0.0;

  for (int iJ = 0; iJ < nJets; ++iJ){
    if (jetFlags[iJ] == 0) result += selJets[iJ].pt();
  }

  return result;
}
