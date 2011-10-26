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
#include "ssSelections.h"
#include "jetSelections.h"
#include "triggerUtils.h"
#include "eventSelections.h"
#include "utilities.h"

struct jet_pt_gt {
     bool operator () (const LorentzVector &v1, const LorentzVector &v2) 
	  {
	       return v1.pt() > v2.pt();
	  }
};

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
        return (pass_electronSelection(idx, electronSelection_ssV6_noIso, false, false, vidx));

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
        return (pass_electronSelection(idx, electronSelection_ssV6_iso));

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
        return (pass_electronSelection(idx, electronSelectionFOV6_ssVBTF80_v3, false, false, vidx) && electronIsolation_rel_v1(idx, true) < 0.60);

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
// get jets and perform overlap removal with numerator e/mu with pt > x (defaults are 10/5 GeV)
/*****************************************************************************************/
std::vector<LorentzVector> samesign::getJets(int idx, enum JetType type, double deltaR, double min_pt, double max_eta, double mu_minpt, double ele_minpt) {

    std::vector<LorentzVector> tmp_jets = getJets(idx, true, type, JETS_CLEAN_HYP_E_MU, deltaR, 0., max_eta);

    // ok, now perform the rest of the lepton overlap removal
    // and the impose the pt requirement after applying the
    // extra corrections
    std::vector<LorentzVector> final_jets;
    for (unsigned int jidx = 0; jidx < tmp_jets.size(); jidx++) {

        bool jetIsLep = false;

        LorentzVector vjet = tmp_jets.at(jidx);
        if (vjet.pt() < min_pt)
            continue;

        for (unsigned int eidx = 0; eidx < cms2.els_p4().size(); eidx++) {
            if (cms2.els_p4().at(eidx).pt() < ele_minpt)
                continue;
            if (!isNumeratorLepton(11, eidx))
                continue;

            if (ROOT::Math::VectorUtil::DeltaR(vjet, cms2.els_p4().at(eidx)) > deltaR)
                continue;

            jetIsLep = true;
            break;
        }

        if (jetIsLep) continue;

        for (unsigned int midx = 0; midx < cms2.mus_p4().size(); midx++) {
            if (cms2.mus_p4().at(midx).pt() < mu_minpt)
                continue;
            if (!isNumeratorLepton(13, midx))
                continue;

            if (ROOT::Math::VectorUtil::DeltaR(vjet, cms2.mus_p4().at(midx)) > deltaR)
                continue;                

            jetIsLep = true;
            break;
        }            

        if (jetIsLep) continue;

        final_jets.push_back(vjet);
    }

    sort(final_jets.begin(), final_jets.end(), jet_pt_gt());
    return final_jets;    
}

/*****************************************************************************************/
// get jets and apply an on-the-fly JEC
/*****************************************************************************************/
std::vector<LorentzVector> samesign::getJets(int idx, FactorizedJetCorrector* jet_corrector, enum JetType type, double deltaR, double min_pt, double max_eta, double mu_minpt, double ele_minpt) {

    std::vector<LorentzVector> tmp_jets = samesign::getJets(idx, type, deltaR, 0., max_eta, mu_minpt, ele_minpt);

    // now impose the pt requirement after applying the extra corrections
    std::vector<LorentzVector> final_jets;
    for (unsigned int jidx = 0; jidx < tmp_jets.size(); jidx++) {

        LorentzVector vjet = tmp_jets.at(jidx);
        float jet_cor = jetCorrection(vjet, jet_corrector);
        vjet *= jet_cor;
        if (vjet.pt() < min_pt)
            continue;

        final_jets.push_back(vjet);
    }

    sort(final_jets.begin(), final_jets.end(), jet_pt_gt());
    return final_jets;
}


/*****************************************************************************************/
// get sumpt, skip jets overlapping with numerator e/mu with pt>x (defaults are 10/5 GeV)
/*****************************************************************************************/
float samesign::sumJetPt(int idx, enum JetType type, double deltaR, double min_pt, double max_eta, double mu_minpt, double ele_minpt){
    std::vector<LorentzVector> good_jets = samesign::getJets(idx, type, deltaR, min_pt, max_eta, mu_minpt, ele_minpt);
    unsigned int nJets = good_jets.size();
    if (nJets == 0)
        return 0.;

    float sumpt = 0.;
    for (unsigned int idx = 0; idx < nJets; idx++) {
        sumpt += good_jets.at(idx).pt();
    }

    return sumpt;
}

/*****************************************************************************************/
// same as above, but allowing use of on-the-fly JEC corrections
/*****************************************************************************************/
float samesign::sumJetPt(int idx, FactorizedJetCorrector* jet_corrector, enum JetType type, double deltaR, double min_pt, double max_eta, double mu_minpt, double ele_minpt) {

    std::vector<LorentzVector> good_jets = samesign::getJets(idx, jet_corrector, type, deltaR, min_pt, max_eta, mu_minpt, ele_minpt);
    unsigned int nJets = good_jets.size();
    if (nJets == 0)
        return 0.;

    float sumpt = 0.;
    for (unsigned int idx = 0; idx < nJets; idx++) {
        sumpt += good_jets.at(idx).pt();
    }

    return sumpt;
}

/*****************************************************************************************/
// get sumpt, skip jets overlapping with numerator e/mu with pt>x (defaults are 10/5 GeV)
/*****************************************************************************************/
int samesign::nJets(int idx, enum JetType type, double deltaR, double min_pt, double max_eta, double mu_minpt, double ele_minpt) {
    
    std::vector<LorentzVector> good_jets = samesign::getJets(idx, type, deltaR, min_pt, max_eta, mu_minpt, ele_minpt);
    return good_jets.size();
}

/*****************************************************************************************/
// same as above, but allowing use of on-the-fly JEC corrections
/*****************************************************************************************/
int samesign::nJets(int idx, FactorizedJetCorrector* jet_corrector, enum JetType type, double deltaR, double min_pt, double max_eta, double mu_minpt, double ele_minpt) {

    std::vector<LorentzVector> good_jets = samesign::getJets(idx, jet_corrector, type, deltaR, min_pt, max_eta, mu_minpt, ele_minpt);
    return good_jets.size();
}


/*****************************************************************************************/
// number of good vertices in the event
/*****************************************************************************************/
int samesign::numberOfGoodVertices() {
    
    int ngv = 0;
    for (unsigned int vidx = 0; vidx < cms2.vtxs_position().size(); vidx++) {
        if (isGoodVertex(vidx))
            ++ngv;
    }

    return ngv;
}
