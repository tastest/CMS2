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
#include "susySelections.h"
#include "jetcorr/FactorizedJetCorrector.h"

using namespace wp2012;

struct jet_pt_gt {
    bool operator () (const LorentzVector &v1, const LorentzVector &v2) 
    {
        return v1.pt() > v2.pt();
    }
};

/////////////////////////////////////////////////////////////////
///                                                           ///
///                                                           ///
///                                                           ///
///          2012 Selections                                  ///
///                                                           ///
///                                                           ///
///                                                           ///
/////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////     
// 2012 good lepton
////////////////////////////////////////////////////////////////////////////////////////////     
bool samesign::isGoodLepton(int id, int idx)
{
    // electrons
    if (abs(id) == 11) {
        if (!pass_electronSelection(idx, electronSelection_ssV7_noIso, false, false))
            return false;
        if ((cms2.els_fiduciality().at(idx) & (1<<ISEB)) == (1<<ISEB))
            return (cms2.els_hOverE().at(idx) < 0.10);
        else
            return (cms2.els_hOverE().at(idx) < 0.075);
    }

    // muons
    if (abs(id) == 13)
        return (muonIdNotIsolated(idx, NominalSSv5));

    return false;
}


////////////////////////////////////////////////////////////////////////////////////////////     
// 2012 isolated lepton
////////////////////////////////////////////////////////////////////////////////////////////     
bool samesign::isIsolatedLepton(int id, int idx)
{
    // electrons
    if (abs(id) == 11)
        return (samesign::electronIsolationPF2012(idx) < 0.09);

    // muons
    if (abs(id) == 13)
        return (muonIsoValuePF2012_deltaBeta(idx) < 0.10);

    return false;
}


////////////////////////////////////////////////////////////////////////////////////////////     
// 2012 numerator lepton
////////////////////////////////////////////////////////////////////////////////////////////     
bool samesign::isNumeratorLepton(int id, int idx)
{
    return (samesign::isGoodLepton(id, idx) && samesign::isIsolatedLepton(id, idx));
}


////////////////////////////////////////////////////////////////////////////////////////////     
// 2012 numerator hypothesis
////////////////////////////////////////////////////////////////////////////////////////////     
bool samesign::isNumeratorHypothesis(int idx)
{
    if (!samesign::isNumeratorLepton(cms2.hyp_lt_id().at(idx), cms2.hyp_lt_index().at(idx)))
        return false;
    if (!samesign::isNumeratorLepton(cms2.hyp_ll_id().at(idx), cms2.hyp_ll_index().at(idx)))
        return false;

    return true;
}

////////////////////////////////////////////////////////////////////////////////////////////     
// 2012 denominator lepton
////////////////////////////////////////////////////////////////////////////////////////////     
bool samesign::isDenominatorLepton(int id, int idx)
{
    // electrons
    if (abs(id) == 11)
        return (pass_electronSelection(idx, electronSelectionFOV7_v3) && samesign::electronIsolationPF2012(idx) < 0.60);

    // muons
    if (abs(id) == 13)
        return (muonId(idx, muonSelectionFO_ssV5));

    return false;
}


///////////////////////////////////////////////////////////////////////////////////////////
// require electron GSF, CTF and SC charges agree
///////////////////////////////////////////////////////////////////////////////////////////
bool samesign::passThreeChargeRequirement(int elIdx)
{
    int trk_idx = cms2.els_trkidx().at(elIdx);

    if (trk_idx >= 0)
    {
        if (cms2.els_sccharge().at(elIdx) == cms2.els_trk_charge().at(elIdx) && cms2.els_trk_charge().at(elIdx) == cms2.trks_charge().at(trk_idx))			   
            return true;
    }

    return false;
}

///////////////////////////////////////////////////////////////////////////////////////////
// calculate PF-based isolation for electrons with rho*Aeff correction
///////////////////////////////////////////////////////////////////////////////////////////
float samesign::electronIsolationPF2012(int idx)
{
    float etaAbs = fabs(cms2.els_etaSC()[idx]);
    float pt     = cms2.els_p4()[idx].pt();

    // get effective area
    float AEff = 0.;
    if (etaAbs <= 1.0) AEff = 0.10;
    else if (etaAbs > 1.0 && etaAbs <= 1.479) AEff = 0.12;
    else if (etaAbs > 1.479 && etaAbs <= 2.0) AEff = 0.085;
    else if (etaAbs > 2.0 && etaAbs <= 2.2) AEff = 0.11;
    else if (etaAbs > 2.2 && etaAbs <= 2.3) AEff = 0.12;
    else if (etaAbs > 2.3 && etaAbs <= 2.4) AEff = 0.12;
    else if (etaAbs > 2.4) AEff = 0.13;

    // pf iso
    // calculate from the ntuple for now...
    //float pfiso_ch = cms2.els_iso03_pf2012_ch().at(idx);
    //float pfiso_em = cms2.els_iso03_pf2012_em().at(idx);
    //float pfiso_nh = cms2.els_iso03_pf2012_nh().at(idx);
    float pfiso_ch = cms2.els_iso03_pf2012ext_ch().at(idx);
    float pfiso_em = cms2.els_iso03_pf2012ext_em().at(idx);
    float pfiso_nh = cms2.els_iso03_pf2012ext_nh().at(idx);

    // rho
    float rhoPrime = std::max(cms2.evt_kt6pf_foregiso_rho(), float(0.0));
    float pfiso_n = std::max(pfiso_em + pfiso_nh - rhoPrime * AEff, float(0.0));
    float pfiso = (pfiso_ch + pfiso_n) / pt;

    return pfiso;
}


///////////////////////////////////////////////////////////////////////////////////////////
// passes dilepton trigger
///////////////////////////////////////////////////////////////////////////////////////////
bool samesign::passesTrigger(int hyp_type)
{
    //----------------------------------------
    // no trigger requirements applied to MC
    //----------------------------------------
  
    if (!cms2.evt_isRealData())
        return true; 
  
    //---------------------------------
    // triggers for dilepton datasets
    //---------------------------------
    //mm
    if (hyp_type == 0) {
        if( passUnprescaledHLTTriggerPattern("HLT_Mu17_Mu8_v" ) )   return true;
    }

    //em
    else if ((hyp_type == 1 || hyp_type == 2)) {
        if( passUnprescaledHLTTriggerPattern("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v") )   return true;
        if( passUnprescaledHLTTriggerPattern("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v") )   return true;
    }
    
    //ee
    else if (hyp_type == 3) {
        if( passUnprescaledHLTTriggerPattern("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v") ) return true;
    }

    return false;
}


///////////////////////////////////////////////////////////////////////////////////////////
// extra Z veto for b-tagged same sign analysis
///////////////////////////////////////////////////////////////////////////////////////////
bool samesign::makesExtraZ(int idx, bool apply_id_iso) {

    std::vector<unsigned int> ele_idx;
    std::vector<unsigned int> mu_idx;

    int lt_id           = cms2.hyp_lt_id().at(idx);
    int ll_id           = cms2.hyp_ll_id().at(idx);
    unsigned int lt_idx = cms2.hyp_lt_index().at(idx);
    unsigned int ll_idx = cms2.hyp_ll_index().at(idx);

    (abs(lt_id) == 11) ? ele_idx.push_back(lt_idx) : mu_idx.push_back(lt_idx);
    (abs(ll_id) == 11) ? ele_idx.push_back(ll_idx) : mu_idx.push_back(ll_idx);

    if (ele_idx.size() + mu_idx.size() != 2) {
        std::cout << "ERROR: don't have 2 leptons in hypothesis!!!  Exiting" << std::endl;
        return false;
    }
        
    if (ele_idx.size() > 0) {
        for (unsigned int eidx = 0; eidx < cms2.els_p4().size(); eidx++) {

            bool is_hyp_lep = false;
            for (unsigned int vidx = 0; vidx < ele_idx.size(); vidx++) {
                if (eidx == ele_idx.at(vidx))
                    is_hyp_lep = true;                
            }
            if (is_hyp_lep)
                continue;

            if (fabs(cms2.els_p4().at(eidx).eta()) > 2.4)
            {
                continue;
            }

            if (cms2.els_p4().at(eidx).pt() < 10.)
            {
                continue;
            }


            if (apply_id_iso) {
                float iso_val = samesign::electronIsolationPF2012(eidx);
                if (iso_val > 0.2)
                {
                    continue;
                }
                
                electronIdComponent_t passAllVetoCuts = DETAIN | DPHIIN | SIGMAIETAIETA | HOE | D0VTX | DZVTX;
                electronIdComponent_t vetoid = electronId_WP2012(eidx, VETO);
                if ((passAllVetoCuts & vetoid) != passAllVetoCuts)
                {
                    continue;
                }
            }

            for (unsigned int vidx = 0; vidx < ele_idx.size(); vidx++) {

                if (cms2.els_charge().at(eidx) * cms2.els_charge().at(ele_idx.at(vidx)) > 0)
                {
                    continue;
                }

                LorentzVector zp4 = cms2.els_p4().at(eidx) + cms2.els_p4().at(ele_idx.at(vidx));
                float zcandmass = sqrt(fabs(zp4.mass2()));
                if (fabs(zcandmass-91.) < 15.)
                {
                    return true;
                }
            }
        }        
    }

    if (mu_idx.size() > 0) {
        for (unsigned int midx = 0; midx < cms2.mus_p4().size(); midx++) {

            bool is_hyp_lep = false;
            for (unsigned int vidx = 0; vidx < mu_idx.size(); vidx++) {
                if (midx == mu_idx.at(vidx))
                    is_hyp_lep = true;                
            }
            if (is_hyp_lep)
                continue;

            if (fabs(cms2.mus_p4().at(midx).eta()) > 2.4)
                continue;

            if (cms2.mus_p4().at(midx).pt() < 10.)
                continue;

            if (apply_id_iso) {
                float iso_val = muonIsoValuePF2012_deltaBeta(midx);
                if (iso_val > 0.2)
                    continue;
                
                if (!cms2.mus_pid_PFMuon().at(midx))
                    continue;

                bool is_global  = ((cms2.mus_type().at(midx) & (1<<1)) == (1<<1));
                bool is_tracker = ((cms2.mus_type().at(midx) & (1<<2)) == (1<<2));
                if (!is_global && !is_tracker)
                    continue;
            }

            for (unsigned int vidx = 0; vidx < mu_idx.size(); vidx++) {

                if (cms2.mus_charge().at(midx) * cms2.mus_charge().at(mu_idx.at(vidx)) > 0)
                    continue;

                LorentzVector zp4 = cms2.mus_p4().at(midx) + cms2.mus_p4().at(mu_idx.at(vidx));
                float zcandmass = sqrt(fabs(zp4.mass2()));
                if (fabs(zcandmass-91.) < 15.)
                    return true;
            }
        }
    }

    return false;
}


////////////////////////////////////////////////////////////////////////////////////////////////
// get jets and perform overlap removal with numerator e/mu with pt > x (defaults are 20/20 GeV)
////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<LorentzVector> samesign::getJets(int idx, enum JetType type, float deltaR, float min_pt, float max_eta, float mu_minpt, float ele_minpt, float rescale, int systFlag) {

  std::vector<LorentzVector> tmp_jets = getJets(idx, true, type, JETS_CLEAN_HYP_E_MU, deltaR, 0., max_eta, (double) rescale, systFlag);

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
            if (!samesign::isNumeratorLepton(11, eidx))
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
            if (!samesign::isNumeratorLepton(13, midx))
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


///////////////////////////////////////////////////////////////////////////////////////////
// get jets and apply an on-the-fly JEC
///////////////////////////////////////////////////////////////////////////////////////////
std::vector<LorentzVector> samesign::getJets(int idx, FactorizedJetCorrector* jet_corrector, enum JetType type, float deltaR, float min_pt, float max_eta, float mu_minpt, float ele_minpt, float rescale, int systFlag) {

  std::vector<bool> tmp_jet_flags = samesign::getJetFlags(idx, type, deltaR, 0., max_eta, mu_minpt, ele_minpt, rescale, systFlag);

    // now impose the pt requirement after applying the extra corrections
    std::vector<LorentzVector> final_jets;
    for (unsigned int jidx = 0; jidx < tmp_jet_flags.size(); jidx++) {

        if (!tmp_jet_flags.at(jidx))
            continue;

        LorentzVector vjet = cms2.pfjets_p4().at(jidx);
        jet_corrector->setRho(cms2.evt_ww_rho_vor());
        jet_corrector->setJetA(cms2.pfjets_area().at(jidx));
        jet_corrector->setJetPt(cms2.pfjets_p4().at(jidx).pt());
        jet_corrector->setJetEta(cms2.pfjets_p4().at(jidx).eta());        
        float jet_cor = jet_corrector->getCorrection();
        vjet *= jet_cor * rescale;
	if (systFlag != 0) {
	  float c = getJetMetSyst(systFlag, vjet.pt(), vjet.eta());
	  vjet *= c;
	}
        if (vjet.pt() < min_pt)
            continue;

        final_jets.push_back(vjet);
    }

    sort(final_jets.begin(), final_jets.end(), jet_pt_gt());
    return final_jets;
}


///////////////////////////////////////////////////////////////////////////////////////////
// get jets and perform overlap removal with numerator e/mu with pt > x (defaults are 20/20 GeV)
///////////////////////////////////////////////////////////////////////////////////////////
std::vector<bool> samesign::getJetFlags(int idx, enum JetType type, float deltaR, float min_pt, float max_eta, float mu_minpt, float ele_minpt, float rescale, int systFlag) {

  std::vector<bool> tmp_jet_flags = getJetFlags((unsigned int)idx, type, JETS_CLEAN_HYP_E_MU, (double)deltaR, 0., (double)max_eta, (double)rescale, systFlag);

    // ok, now perform the rest of the lepton overlap removal
    // and the impose the pt requirement after applying the
    // extra corrections
    std::vector<bool> final_jets;
    for (unsigned int jidx = 0; jidx < tmp_jet_flags.size(); jidx++) {

        if (!tmp_jet_flags.at(jidx)) {
            final_jets.push_back(tmp_jet_flags.at(jidx));
            continue;
        }

        LorentzVector vjet;
        if (cms2.evt_isRealData()) 
	  vjet = cms2.pfjets_p4().at(jidx) * cms2.pfjets_corL1FastL2L3residual().at(jidx) * rescale;
        else
            vjet = cms2.pfjets_p4().at(jidx) * cms2.pfjets_corL1FastL2L3().at(jidx) * rescale;

	if (systFlag != 0) {
	  float c = getJetMetSyst(systFlag, vjet.pt(), vjet.eta());
	  vjet *= c;
	}
        if (vjet.pt() < min_pt) {
            final_jets.push_back(false);
            continue;
        }

        bool jetIsLep = false;

        for (unsigned int eidx = 0; eidx < cms2.els_p4().size(); eidx++) {
            if (cms2.els_p4().at(eidx).pt() < ele_minpt)
                continue;
            if (!samesign::isNumeratorLepton(11, eidx))
                continue;

            if (ROOT::Math::VectorUtil::DeltaR(vjet, cms2.els_p4().at(eidx)) > deltaR)
                continue;

            jetIsLep = true;
            break;
        }

        if (jetIsLep) {
            final_jets.push_back(false);
            continue;
        }

        for (unsigned int midx = 0; midx < cms2.mus_p4().size(); midx++) {
            if (cms2.mus_p4().at(midx).pt() < mu_minpt)
                continue;
            if (!samesign::isNumeratorLepton(13, midx))
                continue;

            if (ROOT::Math::VectorUtil::DeltaR(vjet, cms2.mus_p4().at(midx)) > deltaR)
                continue;                

            jetIsLep = true;
            break;
        }            

        if (jetIsLep) {
            final_jets.push_back(false);
            continue;   
        }

        final_jets.push_back(true);
    }

    return final_jets;    
}


///////////////////////////////////////////////////////////////////////////////////////////
// get jets and apply an on-the-fly JEC and perform overlap removal with numerator
// e/mu with pt > x (defaults are 20/20 GeV)
///////////////////////////////////////////////////////////////////////////////////////////
std::vector<bool> samesign::getJetFlags(int idx, FactorizedJetCorrector* jet_corrector, enum JetType type, float deltaR, float min_pt, float max_eta, float mu_minpt, float ele_minpt, float rescale, int systFlag) {

  std::vector<bool> tmp_jet_flags = samesign::getJetFlags(idx, type, deltaR, 0., max_eta, mu_minpt, ele_minpt, rescale, systFlag);

    assert(tmp_jet_flags.size() == cms2.pfjets_p4().size());

    // now impose the pt requirement after applying the extra corrections
    std::vector<bool> final_jets;
    for (unsigned int jidx = 0; jidx < tmp_jet_flags.size(); jidx++) {

        if (!tmp_jet_flags.at(jidx)) {
            final_jets.push_back(tmp_jet_flags.at(jidx));
            continue;
        }

        LorentzVector vjet = cms2.pfjets_p4().at(jidx);
        jet_corrector->setRho(cms2.evt_ww_rho_vor());
        jet_corrector->setJetA(cms2.pfjets_area().at(jidx));
        jet_corrector->setJetPt(cms2.pfjets_p4().at(jidx).pt());
        jet_corrector->setJetEta(cms2.pfjets_p4().at(jidx).eta());        
        float jet_cor = jet_corrector->getCorrection();
        vjet *= jet_cor * rescale;
	if (systFlag != 0) {
	  float c = getJetMetSyst(systFlag, vjet.pt(), vjet.eta());
	  vjet *= c;
	}
        if (vjet.pt() < min_pt) {
            final_jets.push_back(false);
            continue;   
        }

        final_jets.push_back(true);
    }

    return final_jets;    
}


///////////////////////////////////////////////////////////////////////////////////////////
// get sumpt, skip jets overlapping with numerator e/mu with pt>x (defaults are 20/20 GeV)
///////////////////////////////////////////////////////////////////////////////////////////
float samesign::sumJetPt(int idx_arg, enum JetType type, float deltaR, float min_pt, float max_eta, float mu_minpt, float ele_minpt, float rescale, int systFlag) {
  std::vector<LorentzVector> good_jets = samesign::getJets(idx_arg, type, deltaR, min_pt, max_eta, mu_minpt, ele_minpt, rescale, systFlag);
    unsigned int nJets = good_jets.size();
    if (nJets == 0) return 0.0;
    float sumpt = 0.0;
    for (unsigned int idx = 0; idx < nJets; idx++) {
        sumpt += good_jets.at(idx).pt();
    }
    return sumpt;
}


///////////////////////////////////////////////////////////////////////////////////////////
// same as above, but allowing use of on-the-fly JEC corrections
///////////////////////////////////////////////////////////////////////////////////////////
float samesign::sumJetPt(int idx_arg, FactorizedJetCorrector* jet_corrector, enum JetType type, float deltaR, float min_pt, float max_eta, float mu_minpt, float ele_minpt, float rescale, int systFlag) {
  std::vector<LorentzVector> good_jets = samesign::getJets(idx_arg, jet_corrector, type, deltaR, min_pt, max_eta, mu_minpt, ele_minpt, rescale, systFlag);
    unsigned int nJets = good_jets.size();
    if (nJets == 0) return 0.0;
    float sumpt = 0.0;
    for (unsigned int idx = 0; idx < nJets; idx++) {
        sumpt += good_jets.at(idx).pt();
    }
    return sumpt;
}


///////////////////////////////////////////////////////////////////////////////////////////
// get sumpt, skip jets overlapping with numerator e/mu with pt>x (defaults are 20/20 GeV)
///////////////////////////////////////////////////////////////////////////////////////////
int samesign::nJets(int idx, enum JetType type, float deltaR, float min_pt, float max_eta, float mu_minpt, float ele_minpt, float rescale, int systFlag) {
    
  std::vector<LorentzVector> good_jets = samesign::getJets(idx, type, deltaR, min_pt, max_eta, mu_minpt, ele_minpt, rescale, systFlag);
    return good_jets.size();
}


///////////////////////////////////////////////////////////////////////////////////////////
// same as above, but allowing use of on-the-fly JEC corrections
///////////////////////////////////////////////////////////////////////////////////////////
int samesign::nJets(int idx, FactorizedJetCorrector* jet_corrector, enum JetType type, float deltaR, float min_pt, float max_eta, float mu_minpt, float ele_minpt, float rescale, int systFlag) {

  std::vector<LorentzVector> good_jets = samesign::getJets(idx, jet_corrector, type, deltaR, min_pt, max_eta, mu_minpt, ele_minpt, rescale, systFlag);
    return good_jets.size();
}


///////////////////////////////////////////////////////////////////////////////////////////
// get b-tagged jets and perform overlap removal with numerator e/mu with pt > x (defaults are 20/20 GeV)
///////////////////////////////////////////////////////////////////////////////////////////
std::vector<bool> samesign::getBtaggedJetFlags(int idx, enum JetType type, enum BtagType btag_type, float deltaR, float min_pt, float max_eta, float mu_minpt, float ele_minpt, float rescale, int systFlag)
{
  std::vector<bool> tmp_jet_flags = getBtaggedJetFlags(idx, type, JETS_CLEAN_HYP_E_MU, btag_type, deltaR, 0., max_eta, (double) rescale, systFlag);

    // ok, now perform the rest of the lepton overlap removal
    // and the impose the pt requirement after applying the
    // extra corrections
    std::vector<bool> final_jets;
    for (unsigned int jidx = 0; jidx < tmp_jet_flags.size(); jidx++) {

        if (!tmp_jet_flags.at(jidx)) {
            final_jets.push_back(tmp_jet_flags.at(jidx));
            continue;
        }

        LorentzVector vjet;
        if (cms2.evt_isRealData()) 
            vjet = cms2.pfjets_p4().at(jidx) * cms2.pfjets_corL1FastL2L3residual().at(jidx) * rescale;
        else
            vjet = cms2.pfjets_p4().at(jidx) * cms2.pfjets_corL1FastL2L3().at(jidx) * rescale;
	if (systFlag != 0) {
	  float c = getJetMetSyst(systFlag, vjet.pt(), vjet.eta());
	  vjet *= c;
	}
        if (vjet.pt() < min_pt) {
            final_jets.push_back(false);
            continue;
        }

        bool jetIsLep = false;

        for (unsigned int eidx = 0; eidx < cms2.els_p4().size(); eidx++) {
            if (cms2.els_p4().at(eidx).pt() < ele_minpt)
                continue;
            if (!samesign::isNumeratorLepton(11, eidx))
                continue;

            if (ROOT::Math::VectorUtil::DeltaR(vjet, cms2.els_p4().at(eidx)) > deltaR)
                continue;

            jetIsLep = true;
            break;
        }

        if (jetIsLep) {
            final_jets.push_back(false);
            continue;
        }

        for (unsigned int midx = 0; midx < cms2.mus_p4().size(); midx++) {
            if (cms2.mus_p4().at(midx).pt() < mu_minpt)
                continue;
            if (!samesign::isNumeratorLepton(13, midx))
                continue;

            if (ROOT::Math::VectorUtil::DeltaR(vjet, cms2.mus_p4().at(midx)) > deltaR)
                continue;                

            jetIsLep = true;
            break;
        }            

        if (jetIsLep) {
            final_jets.push_back(false);
            continue;   
        }

        final_jets.push_back(true);
    }

    return final_jets;    
}

///////////////////////////////////////////////////////////////////////////////////////////
// get b=tagged jets and apply an on-the-fly JEC and perform overlap removal with numerator
// e/mu with pt > x (defaults are 20/20 GeV)
///////////////////////////////////////////////////////////////////////////////////////////
std::vector<bool> samesign::getBtaggedJetFlags(int idx, FactorizedJetCorrector* jet_corrector, enum JetType type, enum BtagType btag_type, float deltaR, float min_pt, float max_eta, float mu_minpt, float ele_minpt, float rescale, int systFlag)
{
  std::vector<bool> tmp_jet_flags = samesign::getBtaggedJetFlags(idx, type, btag_type, deltaR, 0., max_eta, mu_minpt, ele_minpt, rescale, systFlag);

    assert(tmp_jet_flags.size() == cms2.pfjets_p4().size());

    // now impose the pt requirement after applying the extra corrections
    std::vector<bool> final_jets;
    for (unsigned int jidx = 0; jidx < tmp_jet_flags.size(); jidx++) {

        if (!tmp_jet_flags.at(jidx)) {
            final_jets.push_back(tmp_jet_flags.at(jidx));
            continue;
        }

        LorentzVector vjet = cms2.pfjets_p4().at(jidx);
        jet_corrector->setRho(cms2.evt_ww_rho_vor());
        jet_corrector->setJetA(cms2.pfjets_area().at(jidx));
        jet_corrector->setJetPt(cms2.pfjets_p4().at(jidx).pt());
        jet_corrector->setJetEta(cms2.pfjets_p4().at(jidx).eta());        
        float jet_cor = jet_corrector->getCorrection();
        vjet *= jet_cor * rescale;
	if (systFlag != 0) {
	  float c = getJetMetSyst(systFlag, vjet.pt(), vjet.eta());
	  vjet *= c;
	}
        if (vjet.pt() < min_pt) {
            final_jets.push_back(false);
            continue;   
        }

        final_jets.push_back(true);
    }

    return final_jets;
}


///////////////////////////////////////////////////////////////////////////////////////////
// get jets and perform overlap removal with numerator e/mu with pt > x (defaults are 10/5 GeV)
///////////////////////////////////////////////////////////////////////////////////////////
std::vector<LorentzVector> samesign::getBtaggedJets(int idx, enum JetType type, enum BtagType btag_type, float deltaR, float min_pt, float max_eta, float mu_minpt, float ele_minpt, float rescale, int systFlag) {

  std::vector<LorentzVector> tmp_jets = getBtaggedJets(idx, true, type, JETS_CLEAN_HYP_E_MU, btag_type, deltaR, 0., max_eta, (double) rescale, systFlag);

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
            if (!samesign::isNumeratorLepton(11, eidx))
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
            if (!samesign::isNumeratorLepton(13, midx))
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


///////////////////////////////////////////////////////////////////////////////////////////
// get jets and apply an on-the-fly JEC and perform overlap removal with numerator
// e/mu with pt > x (defaults are 10/5 GeV)
///////////////////////////////////////////////////////////////////////////////////////////
std::vector<LorentzVector> samesign::getBtaggedJets(int idx, FactorizedJetCorrector* jet_corrector, enum JetType type, enum BtagType btag_type, float deltaR, float min_pt, float max_eta, float mu_minpt, float ele_minpt, float rescale, int systFlag) {

  std::vector<bool> tmp_jet_flags = samesign::getBtaggedJetFlags(idx, type, btag_type, deltaR, 0., max_eta, mu_minpt, ele_minpt, rescale, systFlag);

    // now impose the pt requirement after applying the extra corrections
    std::vector<LorentzVector> final_jets;
    for (unsigned int jidx = 0; jidx < tmp_jet_flags.size(); jidx++) {

        LorentzVector vjet = cms2.pfjets_p4().at(jidx);
        jet_corrector->setRho(cms2.evt_ww_rho_vor());
        jet_corrector->setJetA(cms2.pfjets_area().at(jidx));
        jet_corrector->setJetPt(cms2.pfjets_p4().at(jidx).pt());
        jet_corrector->setJetEta(cms2.pfjets_p4().at(jidx).eta());        
        float jet_cor = jet_corrector->getCorrection();
        vjet *= jet_cor * rescale;
	if (systFlag != 0) {
	  float c = getJetMetSyst(systFlag, vjet.pt(), vjet.eta());
	  vjet *= c;
	}
        if (vjet.pt() < min_pt)
            continue;

        final_jets.push_back(vjet);
    }

    sort(final_jets.begin(), final_jets.end(), jet_pt_gt());
    return final_jets;
    
}


///////////////////////////////////////////////////////////////////////////////////////////
// get sumpt, skip jets overlapping with numerator e/mu with pt>x (defaults are 10/5 GeV)
///////////////////////////////////////////////////////////////////////////////////////////
int samesign::nBtaggedJets(int idx, enum JetType type, enum BtagType btag_type, float deltaR, float min_pt, float max_eta, float mu_minpt, float ele_minpt, float rescale, int systFlag) {

  std::vector<LorentzVector> good_jets = samesign::getBtaggedJets(idx, type, btag_type, deltaR, min_pt, max_eta, mu_minpt, ele_minpt, rescale, systFlag);
    return good_jets.size();
}


///////////////////////////////////////////////////////////////////////////////////////////
// same as above, but allowing use of on-the-fly JEC corrections
///////////////////////////////////////////////////////////////////////////////////////////
int samesign::nBtaggedJets(int idx, FactorizedJetCorrector* jet_corrector, enum JetType type, enum BtagType btag_type, float deltaR, float min_pt, float max_eta, float mu_minpt, float ele_minpt, float rescale, int systFlag) {

  std::vector<LorentzVector> good_jets = samesign::getBtaggedJets(idx, jet_corrector, type, btag_type, deltaR, min_pt, max_eta, mu_minpt, ele_minpt, rescale, systFlag);
    return good_jets.size();    
}


///////////////////////////////////////////////////////////////////////////////////////////
// get sumpt, skip jets overlapping with numerator e/mu with pt>x (defaults are 20/20 GeV)
///////////////////////////////////////////////////////////////////////////////////////////
int samesign::nJets(int idx, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type, enum JetType type, float deltaR, float min_pt, float max_eta, float mu_minpt, float ele_minpt)
{
    std::vector<LorentzVector> good_jets = samesign::getJets(idx, jet_unc, scale_type, type, deltaR, min_pt, max_eta, mu_minpt, ele_minpt);
    return good_jets.size();
}

///////////////////////////////////////////////////////////////////////////////////////////
// same as above, but allowing use of on-the-fly JEC corrections
///////////////////////////////////////////////////////////////////////////////////////////
int samesign::nJets(int idx, FactorizedJetCorrector* jet_corrector, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type, enum JetType type, float deltaR, float min_pt, float max_eta, float mu_minpt, float ele_minpt)
{
    std::vector<LorentzVector> good_jets = samesign::getJets(idx, jet_corrector, jet_unc, scale_type, type, deltaR, min_pt, max_eta, mu_minpt, ele_minpt);
    return good_jets.size();
}

 ///////////////////////////////////////////////////////////////////////////////////////////	 
 // get jets and perform overlap removal with numerator e/mu with pt > x (defaults are 20/20 GeV)	 
 ///////////////////////////////////////////////////////////////////////////////////////////	 
 std::vector<LorentzVector> samesign::getJets(int idx, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type, enum JetType type, float deltaR, float min_pt, float max_eta, float mu_minpt, float ele_minpt)	 
 {	 
     std::vector<LorentzVector> tmp_jets = samesign::getJets(idx, type, deltaR, 0., max_eta, mu_minpt, ele_minpt);	 
 	 
     std::vector<LorentzVector> final_jets;	 
     for (unsigned int jidx = 0; jidx < tmp_jets.size(); jidx++) {	 
 	 
         LorentzVector vjet = tmp_jets.at(jidx);	 
         jet_unc->setJetPt(vjet.pt());	 
         jet_unc->setJetEta(vjet.eta());	 
         float jet_cor = jet_unc->getUncertainty(true);	 
         vjet *= (1. + jet_cor * scale_type);	 
         if (vjet.pt() < min_pt)	 
             continue;	 
 	 
         final_jets.push_back(vjet);	 
     }	 
 	 
     sort(final_jets.begin(), final_jets.end(), jet_pt_gt());	 
     return final_jets;	 
 }	 
 	 
 	 
 ///////////////////////////////////////////////////////////////////////////////////////////	 
 // get jets and apply an on-the-fly JEC and perform overlap removal with numerator	 
 // e/mu with pt > x (defaults are 10/5 GeV)	 
 ///////////////////////////////////////////////////////////////////////////////////////////	 
 std::vector<LorentzVector> samesign::getJets(int idx, FactorizedJetCorrector* jet_corrector, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type,  enum JetType type, float deltaR, float min_pt, float max_eta, float mu_minpt, float ele_minpt)	 
 {	 
     std::vector<LorentzVector> tmp_jets = samesign::getJets(idx, jet_corrector, type, deltaR, 0., max_eta, mu_minpt, ele_minpt);	 
 	 
     std::vector<LorentzVector> final_jets;	 
     for (unsigned int jidx = 0; jidx < tmp_jets.size(); jidx++) {	 
 	 
         LorentzVector vjet = tmp_jets.at(jidx);	 
         jet_unc->setJetPt(vjet.pt());	 
         jet_unc->setJetEta(vjet.eta());	 
         float jet_cor = jet_unc->getUncertainty(true);	 
         vjet *= (1. + jet_cor * scale_type);	 
         if (vjet.pt() < min_pt)	 
             continue;	 
 	 
         final_jets.push_back(vjet);	 
     }	 
 	 
     sort(final_jets.begin(), final_jets.end(), jet_pt_gt());	 
     return final_jets;	 
 }	 
 	 
 	 
 ///////////////////////////////////////////////////////////////////////////////////////////	 
 // get sumpt, skip jets overlapping with numerator e/mu with pt>x (defaults are 20/20 GeV)	 
 ///////////////////////////////////////////////////////////////////////////////////////////	 
 float samesign::sumJetPt(int idx_arg, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type, enum JetType type, float deltaR, float min_pt, float max_eta, float mu_minpt, float ele_minpt)	 
 {	 
     std::vector<LorentzVector> good_jets = samesign::getJets(idx_arg, jet_unc, scale_type, type, deltaR, min_pt, max_eta, mu_minpt, ele_minpt);	 
     unsigned int nJets = good_jets.size();	 
     if (nJets == 0) return 0.0;	 
     float sumpt = 0.0;	 
     for (unsigned int idx = 0; idx < nJets; idx++) {	 
         sumpt += good_jets.at(idx).pt();	 
     }	 
     return sumpt;	 
 }


 ///////////////////////////////////////////////////////////////////////////////////////////	 
 // get jets and perform overlap removal with numerator e/mu with pt > x (defaults are 20/20 GeV)	 
 ///////////////////////////////////////////////////////////////////////////////////////////	 
 std::vector<LorentzVector> samesign::getBtaggedJets(int idx, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type, enum JetType type, enum BtagType btag_type, float deltaR, float min_pt, float max_eta, float mu_minpt, float ele_minpt)	 
 {	 
     std::vector<LorentzVector> tmp_jets = samesign::getBtaggedJets(idx, type, btag_type, deltaR, 0., max_eta, mu_minpt, ele_minpt);	 
 	 
     std::vector<LorentzVector> final_jets;	 
     for (unsigned int jidx = 0; jidx < tmp_jets.size(); jidx++) {	 
 	 
         LorentzVector vjet = tmp_jets.at(jidx);	 
         jet_unc->setJetPt(vjet.pt());	 
         jet_unc->setJetEta(vjet.eta());	 
         float jet_cor = jet_unc->getUncertainty(true);	 
         vjet *= (1. + jet_cor * scale_type);	 
         if (vjet.pt() < min_pt)	 
             continue;	 
 	 
         final_jets.push_back(vjet);	 
     }	 
 	 
     sort(final_jets.begin(), final_jets.end(), jet_pt_gt());	 
     return final_jets;	 
 }	 
 	 
 	 
 ///////////////////////////////////////////////////////////////////////////////////////////	 
 // get jets and apply an on-the-fly JEC and perform overlap removal with numerator	 
 // e/mu with pt > x (defaults are 20/20 GeV)	 
 ///////////////////////////////////////////////////////////////////////////////////////////	 
 std::vector<LorentzVector> samesign::getBtaggedJets(int idx, FactorizedJetCorrector* jet_corrector, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type, enum JetType type, enum BtagType btag_type, float deltaR, float min_pt, float max_eta, float mu_minpt, float ele_minpt)	 
 {	 
     std::vector<LorentzVector> tmp_jets = samesign::getBtaggedJets(idx, jet_corrector, type, btag_type, deltaR, 0., max_eta, mu_minpt, ele_minpt);	 
 	 
     std::vector<LorentzVector> final_jets;	 
     for (unsigned int jidx = 0; jidx < tmp_jets.size(); jidx++) {	 
 	 
         LorentzVector vjet = tmp_jets.at(jidx);	 
         jet_unc->setJetPt(vjet.pt());	 
         jet_unc->setJetEta(vjet.eta());	 
         float jet_cor = jet_unc->getUncertainty(true);	 
         vjet *= (1. + jet_cor * scale_type);	 
         if (vjet.pt() < min_pt)	 
             continue;	 
 	 
         final_jets.push_back(vjet);	 
     }	 
 	 
     sort(final_jets.begin(), final_jets.end(), jet_pt_gt());	 
     return final_jets;	 
 }	 
 	 
 	 
 ///////////////////////////////////////////////////////////////////////////////////////////	 
 // get sumpt, skip jets overlapping with numerator e/mu with pt>x (defaults are 10/5 GeV)	 
 ///////////////////////////////////////////////////////////////////////////////////////////	 
 int samesign::nBtaggedJets(int idx, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type, enum JetType type, enum BtagType btag_type, float deltaR, float min_pt, float max_eta, float mu_minpt, float ele_minpt)	 
 {	 
     std::vector<LorentzVector> good_btags = samesign::getBtaggedJets(idx, jet_unc, scale_type, type, btag_type, deltaR, min_pt, max_eta, mu_minpt, ele_minpt);	 
     return good_btags.size();	 
 }	 

 	 
 ///////////////////////////////////////////////////////////////////////////////////////////	 
 // same as above, but allowing use of on-the-fly JEC corrections	 
 ///////////////////////////////////////////////////////////////////////////////////////////	 
 int samesign::nBtaggedJets(int idx, FactorizedJetCorrector* jet_corrector, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type, enum JetType type, enum BtagType btag_type, float deltaR, float min_pt, float max_eta, float mu_minpt, float ele_minpt)	 
 {	 
     std::vector<LorentzVector> good_btags = samesign::getBtaggedJets(idx, jet_corrector, jet_unc, scale_type, type, btag_type, deltaR, min_pt, max_eta, mu_minpt, ele_minpt);	 
     return good_btags.size();	 
 }	 

 	 
 ///////////////////////////////////////////////////////////////////////////////////////////	 
 // get vector of good els p4s	 
 ///////////////////////////////////////////////////////////////////////////////////////////	 
 std::vector<LorentzVector> samesign::getGoodElectrons()	 
 {	 
     std::vector<LorentzVector> good_els;	 
     for (unsigned int idx = 0; idx < cms2.els_p4().size(); idx++) {	 
 	 
         if (cms2.els_p4().at(idx).pt() < 20.) continue;	 
         if (fabs(cms2.els_p4().at(idx).eta()) > 2.4) continue;	 
 	 
         if (!isNumeratorLepton(11, idx)) continue;	 
	 	 
         good_els.push_back(cms2.els_p4().at(idx));	 
     }	 
	 	 
     sort(good_els.begin(), good_els.end(), jet_pt_gt());	 
     return good_els;	 
 }	 

	 	 
 ///////////////////////////////////////////////////////////////////////////////////////////	 
 // get vector of good mus p4s	 
 ///////////////////////////////////////////////////////////////////////////////////////////	 
 std::vector<LorentzVector> samesign::getGoodMuons()	 
 {	 
     std::vector<LorentzVector> good_mus;	 
     for (unsigned int idx = 0; idx < cms2.mus_p4().size(); idx++) {	 
	 	 
         if (cms2.mus_p4().at(idx).pt() < 20.) continue;	 
         if (fabs(cms2.mus_p4().at(idx).eta()) > 2.4) continue;	 
	 	 
         if (!isNumeratorLepton(13, idx)) continue;	 
	 	 
         good_mus.push_back(cms2.mus_p4().at(idx));	 
     }	 
	 	 
     sort(good_mus.begin(), good_mus.end(), jet_pt_gt());	 
     return good_mus;	 
 }





/////////////////////////////////////////////////////////////////
///                                                           ///
///                                                           ///
///                                                           ///
///          2011 Selections                                  ///
///                                                           ///
///                                                           ///
///                                                           ///
/////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////     
// 2011 good lepton
////////////////////////////////////////////////////////////////////////////////////////////     
bool samesign2011::isGoodLepton(int id, int idx)
{
    // electrons
    if (abs(id) == 11)
        return (pass_electronSelection(idx, electronSelection_ssV6_noIso, false, false));

    // muons
    if (abs(id) == 13)
        return (muonIdNotIsolated(idx, NominalSSv4));

    return false;
}


////////////////////////////////////////////////////////////////////////////////////////////     
// 2011 isolated lepton
////////////////////////////////////////////////////////////////////////////////////////////     
bool samesign2011::isIsolatedLepton(int id, int idx, enum IsolationType iso_type)
{
    // electrons
    if (abs(id) == 11) {
        if (iso_type == DET_ISO)
            return (pass_electronSelection(idx, electronSelection_ssV6_iso));
        else if (iso_type == COR_DET_ISO)
            return (electronIsolation_cor_rel_v1(idx, true) < 0.10);
        else if (iso_type == TIGHT_DET_ISO)
            return (electronIsolation_rel_v1(idx, true) < 0.10);
    }

    // muons
    if (abs(id) == 13) {
        if (iso_type == DET_ISO)
            return (muonIsoValue(idx, false) < 0.15);
        else if (iso_type == COR_DET_ISO)
            return (muonCorIsoValue(idx, false) < 0.10);
        else if (iso_type == TIGHT_DET_ISO)
            return (muonIsoValue(idx, false) < 0.10);
    }

    return false;
}

////////////////////////////////////////////////////////////////////////////////////////////     
// 2011 numerator lepton
////////////////////////////////////////////////////////////////////////////////////////////     
bool samesign2011::isNumeratorLepton(int id, int idx, enum IsolationType iso_type)
{
    return (samesign2011::isGoodLepton(id, idx) && samesign2011::isIsolatedLepton(id, idx, iso_type));
}


////////////////////////////////////////////////////////////////////////////////////////////     
// 2011 numerator hypothesis
////////////////////////////////////////////////////////////////////////////////////////////     
bool samesign2011::isNumeratorHypothesis(int idx, enum IsolationType iso_type)
{
    if (!samesign2011::isNumeratorLepton(cms2.hyp_lt_id().at(idx), cms2.hyp_lt_index().at(idx), iso_type))
        return false;
    if (!samesign2011::isNumeratorLepton(cms2.hyp_ll_id().at(idx), cms2.hyp_ll_index().at(idx), iso_type))
        return false;

    return true;
}


////////////////////////////////////////////////////////////////////////////////////////////     
// 2011 denominator lepton
////////////////////////////////////////////////////////////////////////////////////////////     
bool samesign2011::isDenominatorLepton(int id, int idx, enum IsolationType iso_type)
{
    // electrons
    if (abs(id) == 11) {
        if (iso_type == DET_ISO || iso_type == TIGHT_DET_ISO)
            return (pass_electronSelection(idx, electronSelectionFOV6_ssVBTF80_v3, false, false) && electronIsolation_rel_v1(idx, true) < 0.60);
        else if (iso_type == COR_DET_ISO)
            return (pass_electronSelection(idx, electronSelectionFOV6_ssVBTF80_v3, false, false) && electronIsolation_cor_rel_v1(idx, true) < 0.60);            
    }

    // muons
    if (abs(id) == 13) {
        if (iso_type == DET_ISO || iso_type == TIGHT_DET_ISO)
            return (muonId(idx, muonSelectionFO_ssV4));
        else if (iso_type == COR_DET_ISO)
            return (muonIdNotIsolated(idx, muonSelectionFO_ssV4) && muonCorIsoValue(idx, false) < 0.40);
    }

    return false;
}


///////////////////////////////////////////////////////////////////////////////////////////
// extra Z veto
///////////////////////////////////////////////////////////////////////////////////////////
bool samesign2011::overlapsOtherNNHypInZ(int idx, enum IsolationType iso_type){
    bool result = false;
    int nHyps = cms2.hyp_lt_p4().size();
    for (int iH = 0; iH< nHyps; ++iH){
        if (iH == idx || !hypsOverlap(idx,iH)
            || abs(cms2.hyp_lt_id()[iH])!= abs(cms2.hyp_ll_id()[iH]) || cms2.hyp_lt_id()[iH]*cms2.hyp_ll_id()[iH] > 0) 
            continue;
        if (! samesign2011::isNumeratorHypothesis(iH,iso_type)) continue;
        if (cms2.hyp_p4()[iH].mass2()>0 && fabs(cms2.hyp_p4()[iH].mass() - 91)< 15){ 
            result = true; break;
        }
    }
    return result;
}

///////////////////////////////////////////////////////////////////////////////////////////
// require electron GSF, CTF and SC charges agree
///////////////////////////////////////////////////////////////////////////////////////////
bool samesign2011::passThreeChargeRequirement(int elIdx)
{
    int trk_idx = cms2.els_trkidx()[elIdx];

    if (trk_idx >= 0)
    {
        if (cms2.els_sccharge()[elIdx] == cms2.els_trk_charge()[elIdx] && cms2.els_trk_charge()[elIdx] == cms2.trks_charge()[trk_idx])			   
            return true;
    }

    return false;
}

///////////////////////////////////////////////////////////////////////////////////////////
// get jets and perform overlap removal with numerator e/mu with pt > x (defaults are 10/5 GeV)
///////////////////////////////////////////////////////////////////////////////////////////
std::vector<LorentzVector> samesign2011::getJets(int idx, enum JetType type, double deltaR, double min_pt, double max_eta, double mu_minpt, double ele_minpt, enum IsolationType iso_type, double rescale) {

    std::vector<LorentzVector> tmp_jets = getJets(idx, true, type, JETS_CLEAN_HYP_E_MU, deltaR, 0., max_eta, rescale);

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
            if (!samesign2011::isNumeratorLepton(11, eidx, iso_type))
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
            if (!samesign2011::isNumeratorLepton(13, midx, iso_type))
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

///////////////////////////////////////////////////////////////////////////////////////////
// get jets and apply an on-the-fly JEC
///////////////////////////////////////////////////////////////////////////////////////////
std::vector<LorentzVector> samesign2011::getJets(int idx, FactorizedJetCorrector* jet_corrector, enum JetType type, double deltaR, double min_pt, double max_eta, double mu_minpt, double ele_minpt, enum IsolationType iso_type, double rescale) {

    std::vector<LorentzVector> tmp_jets = samesign2011::getJets(idx, type, deltaR, 0., max_eta, mu_minpt, ele_minpt, iso_type);

    // now impose the pt requirement after applying the extra corrections
    std::vector<LorentzVector> final_jets;
    for (unsigned int jidx = 0; jidx < tmp_jets.size(); jidx++) {

        LorentzVector vjet = tmp_jets.at(jidx);
        float jet_cor = jetCorrection(vjet, jet_corrector);
        vjet *= jet_cor * rescale;
        if (vjet.pt() < min_pt)
            continue;

        final_jets.push_back(vjet);
    }

    sort(final_jets.begin(), final_jets.end(), jet_pt_gt());
    return final_jets;
}


///////////////////////////////////////////////////////////////////////////////////////////
// get jets and perform overlap removal with numerator e/mu with pt > x (defaults are 10/5 GeV)
///////////////////////////////////////////////////////////////////////////////////////////
std::vector<bool> samesign2011::getJetFlags(int idx, enum JetType type, double deltaR, double min_pt, double max_eta, double mu_minpt, double ele_minpt, enum IsolationType iso_type, double rescale) {

    std::vector<bool> tmp_jet_flags = getJetFlags((unsigned int)idx, type, JETS_CLEAN_HYP_E_MU, deltaR, 0., max_eta, rescale);

    // ok, now perform the rest of the lepton overlap removal
    // and the impose the pt requirement after applying the
    // extra corrections
    std::vector<bool> final_jets;
    for (unsigned int jidx = 0; jidx < tmp_jet_flags.size(); jidx++) {

        if (!tmp_jet_flags.at(jidx)) {
            final_jets.push_back(tmp_jet_flags.at(jidx));
            continue;
        }

        LorentzVector vjet = cms2.pfjets_p4().at(jidx) * cms2.pfjets_corL1FastL2L3().at(jidx) * rescale;
        if (vjet.pt() < min_pt) {
            final_jets.push_back(false);
            continue;
        }

        bool jetIsLep = false;

        for (unsigned int eidx = 0; eidx < cms2.els_p4().size(); eidx++) {
            if (cms2.els_p4().at(eidx).pt() < ele_minpt)
                continue;
            if (!samesign2011::isNumeratorLepton(11, eidx, iso_type))
                continue;

            if (ROOT::Math::VectorUtil::DeltaR(vjet, cms2.els_p4().at(eidx)) > deltaR)
                continue;

            jetIsLep = true;
            break;
        }

        if (jetIsLep) {
            final_jets.push_back(false);
            continue;
        }

        for (unsigned int midx = 0; midx < cms2.mus_p4().size(); midx++) {
            if (cms2.mus_p4().at(midx).pt() < mu_minpt)
                continue;
            if (!samesign2011::isNumeratorLepton(13, midx, iso_type))
                continue;

            if (ROOT::Math::VectorUtil::DeltaR(vjet, cms2.mus_p4().at(midx)) > deltaR)
                continue;                

            jetIsLep = true;
            break;
        }            

        if (jetIsLep) {
            final_jets.push_back(false);
            continue;   
        }

        final_jets.push_back(true);
    }

    return final_jets;    
}

///////////////////////////////////////////////////////////////////////////////////////////
// get jets and apply an on-the-fly JEC and perform overlap removal with numerator
// e/mu with pt > x (defaults are 10/5 GeV)
///////////////////////////////////////////////////////////////////////////////////////////
std::vector<bool> samesign2011::getJetFlags(int idx, FactorizedJetCorrector* jet_corrector, enum JetType type, double deltaR, double min_pt, double max_eta, double mu_minpt, double ele_minpt, enum IsolationType iso_type, double rescale) {

    std::vector<bool> tmp_jet_flags     = samesign2011::getJetFlags(idx, type, deltaR, 0., max_eta, mu_minpt, ele_minpt, iso_type);
    std::vector<LorentzVector> tmp_jets = samesign2011::getJets(idx, type, deltaR, 0., max_eta, mu_minpt, ele_minpt, iso_type);

    // now impose the pt requirement after applying the extra corrections
    std::vector<bool> final_jets;
    for (unsigned int jidx = 0; jidx < tmp_jet_flags.size(); jidx++) {

        if (!tmp_jet_flags.at(jidx)) {
            final_jets.push_back(tmp_jet_flags.at(jidx));
            continue;
        }

        LorentzVector vjet = tmp_jets.at(jidx);
        float jet_cor = jetCorrection(vjet, jet_corrector);
        vjet *= jet_cor * rescale;
        if (vjet.pt() < min_pt) {
            final_jets.push_back(false);
            continue;   
        }

        final_jets.push_back(true);
    }

    return final_jets;    
}


///////////////////////////////////////////////////////////////////////////////////////////
// get sumpt, skip jets overlapping with numerator e/mu with pt>x (defaults are 10/5 GeV)
///////////////////////////////////////////////////////////////////////////////////////////
float samesign2011::sumJetPt(int idx_arg, enum JetType type, double deltaR, double min_pt, double max_eta, double mu_minpt, double ele_minpt, enum IsolationType iso_type, double rescale) {
    std::vector<LorentzVector> good_jets = samesign2011::getJets(idx_arg, type, deltaR, min_pt, max_eta, mu_minpt, ele_minpt, iso_type, rescale);
    unsigned int nJets = good_jets.size();
    if (nJets == 0) return 0.0;
    float sumpt = 0.0;
    for (unsigned int idx = 0; idx < nJets; idx++) {
        sumpt += good_jets.at(idx).pt();
    }
    return sumpt;
}

///////////////////////////////////////////////////////////////////////////////////////////
// same as above, but allowing use of on-the-fly JEC corrections
///////////////////////////////////////////////////////////////////////////////////////////
float samesign2011::sumJetPt(int idx_arg, FactorizedJetCorrector* jet_corrector, enum JetType type, double deltaR, double min_pt, double max_eta, double mu_minpt, double ele_minpt, enum IsolationType iso_type, double rescale) {
    std::vector<LorentzVector> good_jets = samesign2011::getJets(idx_arg, jet_corrector, type, deltaR, min_pt, max_eta, mu_minpt, ele_minpt, iso_type, rescale);
    unsigned int nJets = good_jets.size();
    if (nJets == 0) return 0.0;
    float sumpt = 0.0;
    for (unsigned int idx = 0; idx < nJets; idx++) {
        sumpt += good_jets.at(idx).pt();
    }
    return sumpt;
}

///////////////////////////////////////////////////////////////////////////////////////////
// get sumpt, skip jets overlapping with numerator e/mu with pt>x (defaults are 10/5 GeV)
///////////////////////////////////////////////////////////////////////////////////////////
int samesign2011::nJets(int idx, enum JetType type, double deltaR, double min_pt, double max_eta, double mu_minpt, double ele_minpt, enum IsolationType iso_type, double rescale) {
    
    std::vector<LorentzVector> good_jets = samesign2011::getJets(idx, type, deltaR, min_pt, max_eta, mu_minpt, ele_minpt, iso_type, rescale);
    return good_jets.size();
}

///////////////////////////////////////////////////////////////////////////////////////////
// same as above, but allowing use of on-the-fly JEC corrections
///////////////////////////////////////////////////////////////////////////////////////////
int samesign2011::nJets(int idx, FactorizedJetCorrector* jet_corrector, enum JetType type, double deltaR, double min_pt, double max_eta, double mu_minpt, double ele_minpt, enum IsolationType iso_type, double rescale) {

    std::vector<LorentzVector> good_jets = samesign2011::getJets(idx, jet_corrector, type, deltaR, min_pt, max_eta, mu_minpt, ele_minpt, iso_type, rescale);
    return good_jets.size();
}


///////////////////////////////////////////////////////////////////////////////////////////
// passes dilepton trigger
///////////////////////////////////////////////////////////////////////////////////////////
bool samesign2011::passesTrigger(bool is_data, int hyp_type, bool is_high_pt) {
    
    //----------------------------------------
    // no trigger requirements applied to MC
    //----------------------------------------
  
    if (!is_data)
        return true; 
  
    //---------------------------------
    // triggers for lepton-HT datasets
    //---------------------------------
  
    if (!is_high_pt) {
  
        //mm
        if (hyp_type == 0) {
            if( passUnprescaledHLTTriggerPattern("HLT_DoubleMu3_HT150_v") )   return true;
            if( passUnprescaledHLTTriggerPattern("HLT_DoubleMu3_HT160_v") )   return true;
            if( passUnprescaledHLTTriggerPattern("HLT_DoubleMu3_HT200_v") )   return true;
            if( passUnprescaledHLTTriggerPattern("HLT_DoubleMu5_HT150_v") )   return true;
            if( passUnprescaledHLTTriggerPattern("HLT_DoubleMu5_Mass4_HT150_v") )   return true;
            if( passUnprescaledHLTTriggerPattern("HLT_DoubleMu5_Mass8_HT150_v") )   return true;
            if( passUnprescaledHLTTriggerPattern("HLT_DoubleMu8_Mass8_HT150_v") )   return true;
            if( passUnprescaledHLTTriggerPattern("HLT_DoubleMu8_Mass8_HT200_v") )   return true;
            if( passUnprescaledHLTTriggerPattern("HLT_DoubleTkIso10Mu5_Mass8_HT150_v") )   return true;
            if( passUnprescaledHLTTriggerPattern("HLT_DoubleTkIso10Mu5_Mass8_HT200_v") )   return true;
        }
    
        //em
        else if (hyp_type == 1 || hyp_type == 2) {
            if( passUnprescaledHLTTriggerPattern("HLT_Mu3_Ele8_CaloIdL_TrkIdVL_HT150_v") )     return true; 
            if( passUnprescaledHLTTriggerPattern("HLT_Mu3_Ele8_CaloIdT_TrkIdVL_HT150_v") )     return true;
            if( passUnprescaledHLTTriggerPattern("HLT_Mu3_Ele8_CaloIdL_TrkIdVL_HT160_v") )     return true; 
            if( passUnprescaledHLTTriggerPattern("HLT_Mu3_Ele8_CaloIdT_TrkIdVL_HT160_v") )     return true;
            if( passUnprescaledHLTTriggerPattern("HLT_Mu5_Ele8_CaloIdT_TrkIdVL_Mass4_HT150_v") )     return true;
            if( passUnprescaledHLTTriggerPattern("HLT_Mu5_Ele8_CaloIdT_TrkIdVL_Mass8_HT150_v") )     return true;
            if( passUnprescaledHLTTriggerPattern("HLT_Mu8_Ele8_CaloIdT_TrkIdVL_Mass8_HT150_v") )     return true;
            if( passUnprescaledHLTTriggerPattern("HLT_Mu8_Ele8_CaloIdT_TrkIdVL_Mass8_HT200_v") )     return true;
            if( passUnprescaledHLTTriggerPattern("HLT_TkIso10Mu5_Ele8_CaloIdT_CaloIsoVVL_TrkIdVL_Mass8_HT150_v") )     return true;
            if( passUnprescaledHLTTriggerPattern("HLT_TkIso10Mu5_Ele8_CaloIdT_CaloIsoVVL_TrkIdVL_Mass8_HT200_v") )     return true;
        }
    
        //ee
        else if (hyp_type == 3) {
            if( passUnprescaledHLTTriggerPattern("HLT_DoubleEle8_CaloIdL_TrkIdVL_HT150_v") )   return true;
            if( passUnprescaledHLTTriggerPattern("HLT_DoubleEle8_CaloIdT_TrkIdVL_HT150_v") )   return true;
            if( passUnprescaledHLTTriggerPattern("HLT_DoubleEle8_CaloIdL_TrkIdVL_HT160_v") )   return true;
            if( passUnprescaledHLTTriggerPattern("HLT_DoubleEle8_CaloIdT_TrkIdVL_HT160_v") )   return true;
            if( passUnprescaledHLTTriggerPattern("HLT_DoubleEle8_CaloIdT_TrkIdVL_Mass4_HT150_v") )   return true;
            if( passUnprescaledHLTTriggerPattern("HLT_DoubleEle8_CaloIdT_TrkIdVL_Mass8_HT150_v") )   return true;
            if( passUnprescaledHLTTriggerPattern("HLT_DoubleEle8_CaloIdT_TrkIdVL_Mass8_HT200_v") )   return true;
        }
    }
  
    //---------------------------------
    // triggers for dilepton datasets
    //---------------------------------
  
    else {
  
        //mm
        if (hyp_type == 0) {
            if( passUnprescaledHLTTriggerPattern("HLT_DoubleMu7_v") )   return true;
            if( passUnprescaledHLTTriggerPattern("HLT_Mu13_Mu7_v" ) )   return true;
            if( passUnprescaledHLTTriggerPattern("HLT_Mu13_Mu8_v" ) )   return true;
            if( passUnprescaledHLTTriggerPattern("HLT_Mu17_Mu8_v" ) )   return true;
        }
    
        //em
        else if (hyp_type == 1 || hyp_type == 2) {
            if( passUnprescaledHLTTriggerPattern("HLT_Mu17_Ele8_CaloIdL_v") )   return true;
            if( passUnprescaledHLTTriggerPattern("HLT_Mu8_Ele17_CaloIdL_v") )   return true;
            if( passUnprescaledHLTTriggerPattern("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v") )   return true;
            if( passUnprescaledHLTTriggerPattern("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v") )   return true;

        }
    
        //ee
        else if (hyp_type == 3) {
            if( passUnprescaledHLTTriggerPattern("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v") )                                   return true;
            if( passUnprescaledHLTTriggerPattern("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v") ) return true;
            if( passUnprescaledHLTTriggerPattern("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v") ) return true;
        }                                     
    }        
  
    return false;
}

///////////////////////////////////////////////////////////////////////////////////////////
// get jets and perform overlap removal with numerator e/mu with pt > x (defaults are 10/5 GeV)
///////////////////////////////////////////////////////////////////////////////////////////
std::vector<LorentzVector> samesign2011::getBtaggedJets(int idx, enum JetType type, enum BtagType btag_type, double deltaR, double min_pt, double max_eta, double mu_minpt, double ele_minpt, enum IsolationType iso_type, double rescale) {

    std::vector<LorentzVector> tmp_jets = getBtaggedJets(idx, true, type, JETS_CLEAN_HYP_E_MU, btag_type, deltaR, 0., max_eta, rescale);

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
            if (!samesign2011::isNumeratorLepton(11, eidx, iso_type))
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
            if (!samesign2011::isNumeratorLepton(13, midx, iso_type))
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

///////////////////////////////////////////////////////////////////////////////////////////
// get jets and apply an on-the-fly JEC and perform overlap removal with numerator
// e/mu with pt > x (defaults are 10/5 GeV)
///////////////////////////////////////////////////////////////////////////////////////////
std::vector<LorentzVector> samesign2011::getBtaggedJets(int idx, FactorizedJetCorrector* jet_corrector, enum JetType type, enum BtagType btag_type, double deltaR, double min_pt, double max_eta, double mu_minpt, double ele_minpt, enum IsolationType iso_type, double rescale) {

    std::vector<LorentzVector> tmp_jets = samesign2011::getBtaggedJets(idx, type, btag_type, deltaR, 0., max_eta, mu_minpt, ele_minpt, iso_type);

    // now impose the pt requirement after applying the extra corrections
    std::vector<LorentzVector> final_jets;
    for (unsigned int jidx = 0; jidx < tmp_jets.size(); jidx++) {

        LorentzVector vjet = tmp_jets.at(jidx);
        float jet_cor = jetCorrection(vjet, jet_corrector);
        vjet *= jet_cor * rescale;
        if (vjet.pt() < min_pt)
            continue;

        final_jets.push_back(vjet);
    }

    sort(final_jets.begin(), final_jets.end(), jet_pt_gt());
    return final_jets;
    
}


///////////////////////////////////////////////////////////////////////////////////////////
// get sumpt, skip jets overlapping with numerator e/mu with pt>x (defaults are 10/5 GeV)
///////////////////////////////////////////////////////////////////////////////////////////
int samesign2011::nBtaggedJets(int idx, enum JetType type, enum BtagType btag_type, double deltaR, double min_pt, double max_eta, double mu_minpt, double ele_minpt, enum IsolationType iso_type, double rescale) {

    std::vector<LorentzVector> good_jets = samesign2011::getBtaggedJets(idx, type, btag_type, deltaR, min_pt, max_eta, mu_minpt, ele_minpt, iso_type, rescale);
    return good_jets.size();
}

///////////////////////////////////////////////////////////////////////////////////////////
// same as above, but allowing use of on-the-fly JEC corrections
///////////////////////////////////////////////////////////////////////////////////////////
int samesign2011::nBtaggedJets(int idx, FactorizedJetCorrector* jet_corrector, enum JetType type, enum BtagType btag_type, double deltaR, double min_pt, double max_eta, double mu_minpt, double ele_minpt, enum IsolationType iso_type, double rescale) {

    std::vector<LorentzVector> good_jets = samesign2011::getBtaggedJets(idx, jet_corrector, type, btag_type, deltaR, min_pt, max_eta, mu_minpt, ele_minpt, iso_type, rescale);
    return good_jets.size();    
}

///////////////////////////////////////////////////////////////////////////////////////////
// extra Z veto for b-tagged same sign analysis
///////////////////////////////////////////////////////////////////////////////////////////
bool samesign2011::makesExtraZ(int idx, enum IsolationType iso_type, bool apply_id_iso) {

    std::vector<unsigned int> ele_idx;
    std::vector<unsigned int> mu_idx;


    int lt_id   = cms2.hyp_lt_id().at(idx);
    int ll_id   = cms2.hyp_ll_id().at(idx);
    unsigned int lt_idx  = cms2.hyp_lt_index().at(idx);
    unsigned int ll_idx  = cms2.hyp_ll_index().at(idx);

    (abs(lt_id) == 11) ? ele_idx.push_back(lt_idx) : mu_idx.push_back(lt_idx);
    (abs(ll_id) == 11) ? ele_idx.push_back(ll_idx) : mu_idx.push_back(ll_idx);

    if (ele_idx.size() + mu_idx.size() != 2) {
        std::cout << "ERROR: don't have 2 leptons in hypothesis!!!  Exiting" << std::endl;
        return false;
    }
        
    if (ele_idx.size() > 0) {
        for (unsigned int eidx = 0; eidx < cms2.els_p4().size(); eidx++) {

            bool is_hyp_lep = false;
            for (unsigned int vidx = 0; vidx < ele_idx.size(); vidx++) {
                if (eidx == ele_idx.at(vidx))
                    is_hyp_lep = true;                
            }
            if (is_hyp_lep)
                continue;

            if (fabs(cms2.els_p4().at(eidx).eta()) > 2.5)
                continue;

            if (cms2.els_p4().at(eidx).pt() < 10.)
                continue;

            if (apply_id_iso) {
                float iso_val = (iso_type == DET_ISO || iso_type == TIGHT_DET_ISO) ? electronIsolation_rel_v1(eidx, true) : electronIsolation_cor_rel_v1(eidx, true);
                if (iso_val > 0.2)
                    continue;
                
                if (!electronId_VBTF(eidx, VBTF_95_NOHOEEND))
                    continue;                
            }

            for (unsigned int vidx = 0; vidx < ele_idx.size(); vidx++) {

                if (cms2.els_charge().at(eidx) * cms2.els_charge().at(ele_idx.at(vidx)) > 0)
                    continue;

                LorentzVector zp4 = cms2.els_p4().at(eidx) + cms2.els_p4().at(ele_idx.at(vidx));
                float zcandmass = sqrt(fabs(zp4.mass2()));
                if (fabs(zcandmass-91.) < 15.)
                    return true;
            }
        }        
    }

    if (mu_idx.size() > 0) {
        for (unsigned int midx = 0; midx < cms2.mus_p4().size(); midx++) {

            bool is_hyp_lep = false;
            for (unsigned int vidx = 0; vidx < mu_idx.size(); vidx++) {
                if (midx == mu_idx.at(vidx))
                    is_hyp_lep = true;                
            }
            if (is_hyp_lep)
                continue;

            if (fabs(cms2.mus_p4().at(midx).eta()) > 2.5)
                continue;

            if (cms2.mus_p4().at(midx).pt() < 10.)
                continue;

            if (apply_id_iso) {
                float iso_val = (iso_type == DET_ISO || iso_type == TIGHT_DET_ISO) ? muonIsoValue(midx, false) : muonCorIsoValue(midx, false);
                if (iso_val > 0.2)
                    continue;
                
                if (!muonIdNotIsolated(midx, OSGeneric_v4))
                    continue;                
            }

            for (unsigned int vidx = 0; vidx < mu_idx.size(); vidx++) {

                if (cms2.mus_charge().at(midx) * cms2.mus_charge().at(mu_idx.at(vidx)) > 0)
                    continue;

                LorentzVector zp4 = cms2.mus_p4().at(midx) + cms2.mus_p4().at(mu_idx.at(vidx));
                float zcandmass = sqrt(fabs(zp4.mass2()));
                if (fabs(zcandmass-91.) < 15.)
                    return true;
            }
        }
    }

    return false;
}
