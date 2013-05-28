#include <assert.h>
#include <iostream>
#include <algorithm>
#include "Math/Point2Dfwd.h"
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TDatabasePDG.h"
#include "electronSelections.h"
#include "electronSelectionsParameters.h"
#include "muonSelections.h"
#include "metSelections.h"
#include "osSelections.h"
#include "jetSelections.h"
#include "trackSelections.h"
#include "MITConversionUtilities.h"
#include "triggerUtils.h"
#include "eventSelections.h"
#include "utilities.h"
#include "susySelections.h"
#include "jetcorr/FactorizedJetCorrector.h"
#include "jetcorr/JetCorrectionUncertainty.h"
#include "Math/Vector2D.h"
#include "TRandom3.h"

using namespace wp2012;

struct SortByPt 
{
    bool operator () (const LorentzVector& lhs, const LorentzVector& rhs) 
    {
        return lhs.pt() > rhs.pt();
    }

    bool operator () (const std::pair<LorentzVector, unsigned int>& lhs, const std::pair<LorentzVector, unsigned int>& rhs) 
    {
        return lhs.first.pt() > rhs.first.pt();
    }

    bool operator () (const std::pair<LorentzVector, float>& lhs, const std::pair<LorentzVector, float>& rhs) 
    {
        return lhs.first.pt() > rhs.first.pt();
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
// 2012 good lepton (passes ID)
////////////////////////////////////////////////////////////////////////////////////////////     
bool os2012::isGoodLepton(int id, int idx)
{
    // electrons
    if (abs(id) == 11)
    {
        if (!pass_electronSelection(idx, electronSelection_osV1_noIso))
        {
            return false;
        }
    }

    // muons
    if (abs(id) == 13)
    {
        return (muonIdNotIsolated(idx, NominalOSv1));
    }

    return false;
}


////////////////////////////////////////////////////////////////////////////////////////////     
// 2012 isolated lepton
////////////////////////////////////////////////////////////////////////////////////////////     
bool os2012::isIsolatedLepton(int id, int idx)
{
    // electrons
    if (abs(id) == 11)
        return (os2012::electronIsolationPF2012(idx) < 0.15);

    // muons
    if (abs(id) == 13)
        return (muonIsoValuePF2012_deltaBeta(idx) < 0.15);

    return false;
}


////////////////////////////////////////////////////////////////////////////////////////////     
// 2012 lepton isolation value
////////////////////////////////////////////////////////////////////////////////////////////     
double os2012::leptonIsolation(int id, int idx)
{
    // electrons
    if (abs(id) == 11)
    {
        return os2012::electronIsolationPF2012(idx);
    }

    // muons
    if (abs(id) == 13)
    {
        return muonIsoValuePF2012_deltaBeta(idx);
    }

    return -999999.0;
}


////////////////////////////////////////////////////////////////////////////////////////////     
// 2012 effective area 
////////////////////////////////////////////////////////////////////////////////////////////     

float os2012::EffectiveArea03(int id, int idx)
{
    if (abs(id)!=11)
        return -999990.0;

    // use SC eta
    float eta = fabs(cms2.els_etaSC().at(idx));

    // get effective area from electronSelections.h
    return fastJetEffArea03_v2(eta);    // 2013
}

float os2012::EffectiveArea04(int id, int idx)
{
    if (abs(id)!=11)
        return -999990.0;

    // use SC eta
    float eta = fabs(cms2.els_etaSC().at(idx));

    // get effective area from electronSelections.h
    return fastJetEffArea04_v2(eta);    // 2013
}


////////////////////////////////////////////////////////////////////////////////////////////     
// 2012 numerator lepton (passes ID and isolation)
////////////////////////////////////////////////////////////////////////////////////////////     
bool os2012::isNumeratorLepton(int id, int idx)
{
    return (os2012::isGoodLepton(id, idx) && os2012::isIsolatedLepton(id, idx));
}


////////////////////////////////////////////////////////////////////////////////////////////     
// 2012 numerator hypothesis (passes ID and isolation)
////////////////////////////////////////////////////////////////////////////////////////////     
bool os2012::isNumeratorHypothesis(int idx)
{
    if (!os2012::isNumeratorLepton(cms2.hyp_lt_id().at(idx), cms2.hyp_lt_index().at(idx)))
    {
        return false;
    }
    if (!os2012::isNumeratorLepton(cms2.hyp_ll_id().at(idx), cms2.hyp_ll_index().at(idx)))
    {
        return false;
    }

    return true;
}

///////////////////////////////////////////////////////////////////////////////////////////
// calculate PF-based isolation for electrons with rho*Aeff correction
///////////////////////////////////////////////////////////////////////////////////////////
float os2012::electronIsolationPF2012(int idx)
{
    return os2012::electronIsolationPF2012_cone03(idx);
}

float os2012::electronIsolationPF2012_cone03(int idx)
{
    // electron pT
    const float pt = cms2.els_p4().at(idx).pt();

    // get effective area
    const float AEff = EffectiveArea03(11, idx);

    // pf iso
    const float pfiso_ch = cms2.els_iso03_pf2012ext_ch().at(idx);
    const float pfiso_em = cms2.els_iso03_pf2012ext_em().at(idx);
    const float pfiso_nh = cms2.els_iso03_pf2012ext_nh().at(idx);

    // rho
    const float rhoPrime = std::max(cms2.evt_kt6pf_foregiso_rho(), 0.0f);
    const float pfiso_n = std::max(pfiso_em + pfiso_nh - rhoPrime * AEff, 0.0f);
    const float pfiso = (pfiso_ch + pfiso_n) / pt;

    return pfiso;
}

float os2012::electronIsolationPF2012_cone04(int idx)
{
    // electron pT
    const float pt = cms2.els_p4().at(idx).pt();

    // get effective area
    const float AEff = EffectiveArea04(11, idx);

    // pf iso
    // calculate from the ntuple for now...
    const float pfiso_ch = cms2.els_iso04_pf2012ext_ch().at(idx);
    const float pfiso_em = cms2.els_iso04_pf2012ext_em().at(idx);
    const float pfiso_nh = cms2.els_iso04_pf2012ext_nh().at(idx);

    // rho
    const float rhoPrime = std::max(cms2.evt_kt6pf_foregiso_rho(), 0.0f);
    const float pfiso_n = std::max(pfiso_em + pfiso_nh - rhoPrime * AEff, 0.0f);
    const float pfiso = (pfiso_ch + pfiso_n) / pt;

    return pfiso;
}


///////////////////////////////////////////////////////////////////////////////////////////
// passes dilepton trigger
///////////////////////////////////////////////////////////////////////////////////////////

// analysis type:
//   0 --> use high pT analysis triggers
//   anything else will return false

bool os2012::passesTrigger(int hyp_type, int analysis_type)
{
    //----------------------------------------
    // no trigger requirements applied to MC
    //----------------------------------------

    if (!cms2.evt_isRealData())
        return true; 

    switch(analysis_type)
    {
        case 0: return passesTriggerHighPt(hyp_type); break;
        default: return false;
    }

    return false;
}

bool os2012::passesTriggerHighPt(int hyp_type)
{
    //----------------------------------------
    // no trigger requirements applied to MC
    //----------------------------------------

    if (!cms2.evt_isRealData())
        return true; 

    //---------------------------------
    // triggers for dilepton datasets
    //---------------------------------

    // mm
    if (hyp_type == 0) {
        if (passUnprescaledHLTTriggerPattern("HLT_Mu17_Mu8_v")) {return true;}
    }

    // em
    else if ((hyp_type == 1 || hyp_type == 2)) {
        if (passUnprescaledHLTTriggerPattern("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v")) {return true;}
        if (passUnprescaledHLTTriggerPattern("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v")) {return true;}
    }

    // ee
    else if (hyp_type == 3) {
        if (passUnprescaledHLTTriggerPattern("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v")) {return true;}
    }

    return false;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// 2012 get jets and perform overlap removal with numerator e/mu with pt > x (defaults are 20/20 GeV)
////////////////////////////////////////////////////////////////////////////////////////////////
    
// JEC taken from ntuple
std::vector<LorentzVector> os2012::getJets(int idx, enum JetType type, float deltaR, float min_pt, float max_eta, float min_lep_pt, int systFlag, bool sort_by_pt) {

    std::vector<LorentzVector> tmp_jets = getJets(idx, /*sort=*/true, type, JETS_CLEAN_HYP_E_MU, deltaR, 0.0, max_eta, (double) /*rescale=*/1., systFlag);

    // ok, now perform the rest of the lepton overlap removal
    // and the impose the pt requirement after applying the
    // extra corrections
    std::vector<LorentzVector> final_jets;
    for (unsigned int jidx = 0; jidx < tmp_jets.size(); jidx++) {

        bool jetIsLep = false;

        LorentzVector vjet = tmp_jets.at(jidx);
        if (vjet.pt() < min_pt)
            continue;

        std::vector<LorentzVector> ele_p4s = os2012::getGoodElectrons(min_lep_pt);
        for (unsigned int eidx = 0; eidx < ele_p4s.size(); eidx++) {
            if (ROOT::Math::VectorUtil::DeltaR(vjet, ele_p4s.at(eidx)) > deltaR)
                continue;

            jetIsLep = true;
            break;
        }

        if (jetIsLep) continue;

        std::vector<LorentzVector> mu_p4s = os2012::getGoodMuons(min_lep_pt);
        for (unsigned int midx = 0; midx < mu_p4s.size(); midx++) {
            if (ROOT::Math::VectorUtil::DeltaR(vjet, mu_p4s.at(midx)) > deltaR)
                continue;                

            jetIsLep = true;
            break;
        }            

        if (jetIsLep) continue;

        final_jets.push_back(vjet);
    }

    if (sort_by_pt)
        sort(final_jets.begin(), final_jets.end(), SortByPt());
    return final_jets;    
}


// JEC applied otf
std::vector<LorentzVector> os2012::getJets(int idx, FactorizedJetCorrector* jet_corrector, enum JetType type, float deltaR, float min_pt, float max_eta, float min_lep_pt, int systFlag, bool sort_by_pt)
{
    std::vector<bool> tmp_jet_flags = os2012::getJetFlags(idx, type, deltaR, /*min_pt=*/0.0, max_eta, min_lep_pt, systFlag);

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
        vjet *= jet_cor;
        if (systFlag != 0) {
            float c = getJetMetSyst(systFlag, vjet.pt(), vjet.eta());
            vjet *= c;
        }
        if (vjet.pt() < min_pt)
        {
            continue;
        }

        final_jets.push_back(vjet);
    }

    if (sort_by_pt)
        sort(final_jets.begin(), final_jets.end(), SortByPt());
    return final_jets;
}

// JEC uncertainty applied otf
std::vector<LorentzVector> os2012::getJets(int idx, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type, enum JetType type, float deltaR, float min_pt, float max_eta, float min_lep_pt, bool sort_by_pt)
{    
    std::vector<bool> tmp_jet_flags = os2012::getJetFlags(idx, type, deltaR, /*min_pt=*/0.0, max_eta, min_lep_pt);

    // now impose the pt requirement after applying the extra corrections
    std::vector<LorentzVector> final_jets;
    for (unsigned int jidx = 0; jidx < tmp_jet_flags.size(); jidx++) {

        if (!tmp_jet_flags.at(jidx))
            continue;

        const float jet_cor = (cms2.evt_isRealData() ? cms2.pfjets_corL1FastL2L3residual().at(jidx) :  cms2.pfjets_corL1FastL2L3().at(jidx));
        LorentzVector vjet  = cms2.pfjets_p4().at(jidx) * jet_cor; 
        jet_unc->setJetPt(vjet.pt());    
        jet_unc->setJetEta(vjet.eta());  
        const float jet_cor_unc = jet_unc->getUncertainty(true);     
        vjet *= (1.0 + jet_cor_unc * scale_type);    
        if (vjet.pt() < min_pt)
        {
            continue;
        }

        final_jets.push_back(vjet);
    }

    if (sort_by_pt)
        sort(final_jets.begin(), final_jets.end(), SortByPt());
    return final_jets;
}

// JEC AND JEC uncertainty applied otf
std::vector<LorentzVector> os2012::getJets(int idx, FactorizedJetCorrector* jet_corrector, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type,  enum JetType type, float deltaR, float min_pt, float max_eta, float min_lep_pt, bool sort_by_pt)
{    
    std::vector<bool> tmp_jet_flags = os2012::getJetFlags(idx, type, deltaR, /*min_pt=*/0.0, max_eta, min_lep_pt);

    // now impose the pt requirement after applying the extra corrections
    std::vector<LorentzVector> final_jets;
    for (unsigned int jidx = 0; jidx < tmp_jet_flags.size(); jidx++) {

        if (!tmp_jet_flags.at(jidx))
            continue;

        jet_corrector->setRho(cms2.evt_ww_rho_vor());
        jet_corrector->setJetA(cms2.pfjets_area().at(jidx));
        jet_corrector->setJetPt(cms2.pfjets_p4().at(jidx).pt());
        jet_corrector->setJetEta(cms2.pfjets_p4().at(jidx).eta());        
        const float jet_cor = jet_corrector->getCorrection();
        LorentzVector vjet = cms2.pfjets_p4().at(jidx) * jet_cor;
        jet_unc->setJetPt(vjet.pt());    
        jet_unc->setJetEta(vjet.eta());  
        const float jet_cor_unc = jet_unc->getUncertainty(true);     
        vjet *= (1.0 + jet_cor_unc * scale_type);    
        if (vjet.pt() < min_pt)
        {
            continue;
        }

        final_jets.push_back(vjet);
    }
     
    if (sort_by_pt)
        sort(final_jets.begin(), final_jets.end(), SortByPt());  
    return final_jets;   
}    
     

///////////////////////////////////////////////////////////////////////////////////////////////////////////
// 2012 get jet flags and perform overlap removal with numerator e/mu with pt > x (defaults are 20/20 GeV)
///////////////////////////////////////////////////////////////////////////////////////////////////////////

// JEC taken from ntuple
std::vector<bool> os2012::getJetFlags(int idx, enum JetType type, float deltaR, float min_pt, float max_eta, float min_lep_pt, int systFlag)
{
    std::vector<bool> tmp_jet_flags = getJetFlags((unsigned int)idx, type, JETS_CLEAN_HYP_E_MU, (double)deltaR, /*min_pt=*/0.0, (double)max_eta, /*rescale=*/1., systFlag);

    // ok, now perform the rest of the lepton overlap removal
    // and the impose the pt requirement after applying the
    // extra corrections
    std::vector<bool> final_jets;
    for (unsigned int jidx = 0; jidx < tmp_jet_flags.size(); jidx++) {

        if (!tmp_jet_flags.at(jidx)) {
            final_jets.push_back(tmp_jet_flags.at(jidx));
            continue;
        }

        const float jet_cor = (cms2.evt_isRealData() ? cms2.pfjets_corL1FastL2L3residual().at(jidx) :  cms2.pfjets_corL1FastL2L3().at(jidx));
        LorentzVector vjet  = cms2.pfjets_p4().at(jidx) * jet_cor; 
        if (systFlag != 0) {
            float c = getJetMetSyst(systFlag, vjet.pt(), vjet.eta());
            vjet *= c;
        }
        if (vjet.pt() < min_pt) {
            final_jets.push_back(false);
            continue;
        }

        bool jetIsLep = false;

        std::vector<LorentzVector> ele_p4s = os2012::getGoodElectrons(min_lep_pt);
        for (unsigned int eidx = 0; eidx < ele_p4s.size(); eidx++) {
            if (ROOT::Math::VectorUtil::DeltaR(vjet, ele_p4s.at(eidx)) > deltaR)
                continue;

            jetIsLep = true;
            break;
        }

        if (jetIsLep) {
            final_jets.push_back(false);
            continue;
        }

        std::vector<LorentzVector> mu_p4s = os2012::getGoodMuons(min_lep_pt);
        for (unsigned int midx = 0; midx < mu_p4s.size(); midx++) {
            if (ROOT::Math::VectorUtil::DeltaR(vjet, mu_p4s.at(midx)) > deltaR)
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


// JEC applied otf
std::vector<bool> os2012::getJetFlags(int idx, FactorizedJetCorrector* jet_corrector, enum JetType type, float deltaR, float min_pt, float max_eta, float min_lep_pt, int systFlag)
{
    std::vector<bool> tmp_jet_flags = os2012::getJetFlags(idx, type, deltaR, /*min_pt=*/0.0, max_eta, min_lep_pt, systFlag);

    assert(tmp_jet_flags.size() == cms2.pfjets_p4().size());

    // now impose the pt requirement after applying the extra corrections
    std::vector<bool> final_jets;
    for (unsigned int jidx = 0; jidx < tmp_jet_flags.size(); jidx++) {

        if (!tmp_jet_flags.at(jidx)) {
            final_jets.push_back(tmp_jet_flags.at(jidx));
            continue;
        }

        jet_corrector->setRho(cms2.evt_ww_rho_vor());
        jet_corrector->setJetA(cms2.pfjets_area().at(jidx));
        jet_corrector->setJetPt(cms2.pfjets_p4().at(jidx).pt());
        jet_corrector->setJetEta(cms2.pfjets_p4().at(jidx).eta());        
        float jet_cor = jet_corrector->getCorrection();
        LorentzVector vjet = cms2.pfjets_p4().at(jidx) * jet_cor;
        if (systFlag != 0) 
        {
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


// JEC uncertainty applied otf
std::vector<bool> os2012::getJetFlags(int idx, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type, enum JetType type, float deltaR, float min_pt, float max_eta, float min_lep_pt)     
{
    std::vector<bool> tmp_jet_flags = os2012::getJetFlags(idx, type, deltaR, /*min_pt=*/0.0, max_eta, min_lep_pt);

    assert(tmp_jet_flags.size() == cms2.pfjets_p4().size());

    // now impose the pt requirement after applying the extra corrections
    std::vector<bool> final_jets;
    for (unsigned int jidx = 0; jidx < tmp_jet_flags.size(); jidx++) {

        if (!tmp_jet_flags.at(jidx)) {
            final_jets.push_back(tmp_jet_flags.at(jidx));
            continue;
        }
        const float jet_cor = (cms2.evt_isRealData() ? cms2.pfjets_corL1FastL2L3residual().at(jidx) :  cms2.pfjets_corL1FastL2L3().at(jidx));
        LorentzVector vjet  = cms2.pfjets_p4().at(jidx) * jet_cor; 
        jet_unc->setJetPt(vjet.pt());    
        jet_unc->setJetEta(vjet.eta());  
        const float jet_cor_unc = jet_unc->getUncertainty(true);     
        vjet *= (1.0 + jet_cor_unc * scale_type);    
        if (vjet.pt() < min_pt) {
            final_jets.push_back(false);
            continue;   
        }

        final_jets.push_back(true);
    }

    return final_jets;    
}


// JEC AND JEC uncertainty applied otf
std::vector<bool> os2012::getJetFlags(int idx, FactorizedJetCorrector* jet_corrector, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type,  enum JetType type, float deltaR, float min_pt, float max_eta, float min_lep_pt)
{
    std::vector<bool> tmp_jet_flags = os2012::getJetFlags(idx, type, deltaR, /*min_pt=*/0.0, max_eta, min_lep_pt);

    assert(tmp_jet_flags.size() == cms2.pfjets_p4().size());

    // now impose the pt requirement after applying the extra corrections
    std::vector<bool> final_jets;
    for (unsigned int jidx = 0; jidx < tmp_jet_flags.size(); jidx++) {

        if (!tmp_jet_flags.at(jidx)) {
            final_jets.push_back(tmp_jet_flags.at(jidx));
            continue;
        }
        jet_corrector->setRho(cms2.evt_ww_rho_vor());
        jet_corrector->setJetA(cms2.pfjets_area().at(jidx));
        jet_corrector->setJetPt(cms2.pfjets_p4().at(jidx).pt());
        jet_corrector->setJetEta(cms2.pfjets_p4().at(jidx).eta());        
        float jet_cor = jet_corrector->getCorrection();
        LorentzVector vjet = cms2.pfjets_p4().at(jidx) * jet_cor;
        jet_unc->setJetPt(vjet.pt());    
        jet_unc->setJetEta(vjet.eta());  
        const float jet_cor_unc = jet_unc->getUncertainty(true);     
        vjet *= (1.0 + jet_cor_unc * scale_type);    
        if (vjet.pt() < min_pt) {
            final_jets.push_back(false);
            continue;   
        }

        final_jets.push_back(true);
    }

    return final_jets;    
}


///////////////////////////////////////////////////////////////////////////////////////////
// 2012 get sumpt, skip jets overlapping with numerator e/mu with pt>x (defaults are 20/20 GeV)
///////////////////////////////////////////////////////////////////////////////////////////

// JEC taken from ntuple
float os2012::sumJetPt(int idx_arg, enum JetType type, float deltaR, float min_pt, float max_eta, float min_lep_pt, int systFlag) {
    std::vector<LorentzVector> good_jets = os2012::getJets(idx_arg, type, deltaR, min_pt, max_eta, min_lep_pt, systFlag);
    unsigned int nJets = good_jets.size();
    if (nJets == 0) return 0.0;
    float sumpt = 0.0;
    for (unsigned int idx = 0; idx < nJets; idx++) {
        sumpt += good_jets.at(idx).pt();
    }
    return sumpt;
}


// JEC applied otf
float os2012::sumJetPt(int idx_arg, FactorizedJetCorrector* jet_corrector, enum JetType type, float deltaR, float min_pt, float max_eta, float min_lep_pt, int systFlag) {
    std::vector<LorentzVector> good_jets = os2012::getJets(idx_arg, jet_corrector, type, deltaR, min_pt, max_eta, min_lep_pt, systFlag);
    unsigned int nJets = good_jets.size();
    if (nJets == 0) return 0.0;
    float sumpt = 0.0;
    for (unsigned int idx = 0; idx < nJets; idx++) {
        sumpt += good_jets.at(idx).pt();
    }
    return sumpt;
}


// JEC uncertainty applied otf
float os2012::sumJetPt(int idx_arg, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type, enum JetType type, float deltaR, float min_pt, float max_eta, float min_lep_pt)    
{    
    std::vector<LorentzVector> good_jets = os2012::getJets(idx_arg, jet_unc, scale_type, type, deltaR, min_pt, max_eta, min_lep_pt);  
    unsigned int nJets = good_jets.size();   
    if (nJets == 0) return 0.0;  
    float sumpt = 0.0;   
    for (unsigned int idx = 0; idx < nJets; idx++) {     
        sumpt += good_jets.at(idx).pt();     
    }    
    return sumpt;    
}


// JEC AND JEC uncertainty applied otf
float os2012::sumJetPt(int idx_arg, FactorizedJetCorrector* jet_corrector, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type, enum JetType type, float deltaR, float min_pt, float max_eta, float min_lep_pt)
{    
    std::vector<LorentzVector> good_jets = os2012::getJets(idx_arg, jet_corrector, jet_unc, scale_type, type, deltaR, min_pt, max_eta, min_lep_pt);
    unsigned int nJets = good_jets.size();
    if (nJets == 0) return 0.0;  
    float sumpt = 0.0;   
    for (unsigned int idx = 0; idx < nJets; idx++) {     
        sumpt += good_jets.at(idx).pt();     
    }    
    return sumpt;    
}


///////////////////////////////////////////////////////////////////////////////////////////
// 2012 get njets, skip jets overlapping with numerator e/mu with pt>x (defaults are 20/20 GeV)
///////////////////////////////////////////////////////////////////////////////////////////

// JEC taken from ntuple
int os2012::nJets(int idx, enum JetType type, float deltaR, float min_pt, float max_eta, float min_lep_pt, int systFlag)
{
    std::vector<LorentzVector> good_jets = os2012::getJets(idx, type, deltaR, min_pt, max_eta, min_lep_pt, systFlag);
    return good_jets.size();
}


// JEC applied otf
int os2012::nJets(int idx, FactorizedJetCorrector* jet_corrector, enum JetType type, float deltaR, float min_pt, float max_eta, float min_lep_pt, int systFlag)
{
    std::vector<LorentzVector> good_jets = os2012::getJets(idx, jet_corrector, type, deltaR, min_pt, max_eta, min_lep_pt, systFlag);
    return good_jets.size();
}


// JEC uncertainty applied otf
int os2012::nJets(int idx, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type, enum JetType type, float deltaR, float min_pt, float max_eta, float min_lep_pt)
{
    std::vector<LorentzVector> good_jets = os2012::getJets(idx, jet_unc, scale_type, type, deltaR, min_pt, max_eta, min_lep_pt);
    return good_jets.size();
}


// JEC AND JEC uncertainty applied otf
int os2012::nJets(int idx, FactorizedJetCorrector* jet_corrector, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type, enum JetType type, float deltaR, float min_pt, float max_eta, float min_lep_pt)
{
    std::vector<LorentzVector> good_jets = os2012::getJets(idx, jet_corrector, jet_unc, scale_type, type, deltaR, min_pt, max_eta, min_lep_pt);
    return good_jets.size();
}


///////////////////////////////////////////////////////////////////////////////////////////
// 2012 get b-tagged jets and perform overlap removal with numerator e/mu with pt > x (defaults are 20/20 GeV)
///////////////////////////////////////////////////////////////////////////////////////////

// JEC taken from ntuple
std::vector<LorentzVector> os2012::getBtaggedJets(int idx, enum JetType type, enum BtagType btag_type, float deltaR, float min_pt, float max_eta, float min_lep_pt, int systFlag, bool sort_by_pt)
{
    std::vector<LorentzVector> tmp_jets = getBtaggedJets(idx, true, type, JETS_CLEAN_HYP_E_MU, btag_type, deltaR, /*min_pt=*/0.0, max_eta, /*rescale=*/1., systFlag);

    // ok, now perform the rest of the lepton overlap removal
    // and the impose the pt requirement after applying the
    // extra corrections
    std::vector<LorentzVector> final_jets;
    for (unsigned int jidx = 0; jidx < tmp_jets.size(); jidx++) {

        bool jetIsLep = false;

        LorentzVector vjet = tmp_jets.at(jidx);
        if (vjet.pt() < min_pt)
            continue;

        std::vector<LorentzVector> ele_p4s = os2012::getGoodElectrons(min_lep_pt);
        for (unsigned int eidx = 0; eidx < ele_p4s.size(); eidx++) {
            if (ROOT::Math::VectorUtil::DeltaR(vjet, ele_p4s.at(eidx)) > deltaR)
                continue;

            jetIsLep = true;
            break;
        }

        if (jetIsLep) continue;

        std::vector<LorentzVector> mu_p4s = os2012::getGoodMuons(min_lep_pt);
        for (unsigned int midx = 0; midx < mu_p4s.size(); midx++) {
            if (ROOT::Math::VectorUtil::DeltaR(vjet, mu_p4s.at(midx)) > deltaR)
                continue;

            jetIsLep = true;
            break;
        }            

        if (jetIsLep) continue;

        final_jets.push_back(vjet);
    }


    if (sort_by_pt)
        sort(final_jets.begin(), final_jets.end(), SortByPt());
    return final_jets;        
}


// JEC applied otf
std::vector<LorentzVector> os2012::getBtaggedJets(int idx, FactorizedJetCorrector* jet_corrector, enum JetType type, enum BtagType btag_type, float deltaR, float min_pt, float max_eta, float min_lep_pt, int systFlag, bool sort_by_pt)
{
    std::vector<bool> tmp_jet_flags = os2012::getBtaggedJetFlags(idx, type, btag_type, deltaR, /*min_pt=*/0.0, max_eta, min_lep_pt, systFlag);

    // now impose the pt requirement after applying the extra corrections
    std::vector<LorentzVector> final_jets;
    for (unsigned int jidx = 0; jidx < tmp_jet_flags.size(); jidx++) {

        if (!tmp_jet_flags.at(jidx)) continue;

        jet_corrector->setRho(cms2.evt_ww_rho_vor());
        jet_corrector->setJetA(cms2.pfjets_area().at(jidx));
        jet_corrector->setJetPt(cms2.pfjets_p4().at(jidx).pt());
        jet_corrector->setJetEta(cms2.pfjets_p4().at(jidx).eta());        
        float jet_cor = jet_corrector->getCorrection();
        LorentzVector vjet = cms2.pfjets_p4().at(jidx) * jet_cor;
        if (systFlag != 0) {
            float c = getJetMetSyst(systFlag, vjet.pt(), vjet.eta());
            vjet *= c;
        }
        if (vjet.pt() < min_pt)
        {
            continue;
        }
        final_jets.push_back(vjet);
    }

    if (sort_by_pt)
        sort(final_jets.begin(), final_jets.end(), SortByPt());
    return final_jets;
}


// JEC uncertainty applied otf
std::vector<LorentzVector> os2012::getBtaggedJets(int idx, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type, enum JetType type, enum BtagType btag_type, float deltaR, float min_pt, float max_eta, float min_lep_pt, bool sort_by_pt)
{    
    std::vector<bool> tmp_jet_flags = os2012::getBtaggedJetFlags(idx, type, btag_type, deltaR, /*min_pt=*/0.0, max_eta, min_lep_pt);

    // now impose the pt requirement after applying the extra corrections
    std::vector<LorentzVector> final_jets;
    for (unsigned int jidx = 0; jidx < tmp_jet_flags.size(); jidx++) {

        if (!tmp_jet_flags.at(jidx)) continue;

        const float jet_cor = (cms2.evt_isRealData() ? cms2.pfjets_corL1FastL2L3residual().at(jidx) :  cms2.pfjets_corL1FastL2L3().at(jidx));
        LorentzVector vjet  = cms2.pfjets_p4().at(jidx) * jet_cor; 
        jet_unc->setJetPt(vjet.pt());    
        jet_unc->setJetEta(vjet.eta());  
        const float jet_cor_unc = jet_unc->getUncertainty(true);     
        vjet *= (1.0 + jet_cor_unc * scale_type);    
        if (vjet.pt() < min_pt)
        {
            continue;
        }
        final_jets.push_back(vjet);
    }

    if (sort_by_pt)
        sort(final_jets.begin(), final_jets.end(), SortByPt());
    return final_jets;
}    
     
 
// JEC AND JEC uncertainty applied otf
std::vector<LorentzVector> os2012::getBtaggedJets(int idx, FactorizedJetCorrector* jet_corrector, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type, enum JetType type, enum BtagType btag_type, float deltaR, float min_pt, float max_eta, float min_lep_pt, bool sort_by_pt)
{    
    std::vector<bool> tmp_jet_flags = os2012::getBtaggedJetFlags(idx, type, btag_type, deltaR, /*min_pt=*/0.0, max_eta, min_lep_pt);

    // now impose the pt requirement after applying the extra corrections
    std::vector<LorentzVector> final_jets;
    for (unsigned int jidx = 0; jidx < tmp_jet_flags.size(); jidx++) {

        if (!tmp_jet_flags.at(jidx)) continue;

        jet_corrector->setRho(cms2.evt_ww_rho_vor());
        jet_corrector->setJetA(cms2.pfjets_area().at(jidx));
        jet_corrector->setJetPt(cms2.pfjets_p4().at(jidx).pt());
        jet_corrector->setJetEta(cms2.pfjets_p4().at(jidx).eta());        
        float jet_cor = jet_corrector->getCorrection();
        LorentzVector vjet = cms2.pfjets_p4().at(jidx) * jet_cor;
        jet_unc->setJetPt(vjet.pt());    
        jet_unc->setJetEta(vjet.eta());  
        const float jet_cor_unc = jet_unc->getUncertainty(true);     
        vjet *= (1.0 + jet_cor_unc * scale_type);    
        if (vjet.pt() < min_pt)
        {
            continue;
        }
        final_jets.push_back(vjet);
    }

    if (sort_by_pt)
        sort(final_jets.begin(), final_jets.end(), SortByPt());
    return final_jets;
}    
     
     
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 2012 get b-tagged jet flags and perform overlap removal with numerator e/mu with pt > x (defaults are 20/20 GeV)
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// JEC taken from ntuple
std::vector<bool> os2012::getBtaggedJetFlags(int idx, enum JetType type, enum BtagType btag_type, float deltaR, float min_pt, float max_eta, float min_lep_pt, int systFlag)
{
    std::vector<bool> tmp_jet_flags       = getBtaggedJetFlags(idx, type, JETS_CLEAN_HYP_E_MU, btag_type, deltaR, /*min_pt=*/0.0, max_eta, /*rescale=*/1., systFlag);

    // ok, now perform the rest of the lepton overlap removal
    // and the impose the pt requirement after applying the
    // extra corrections
    std::vector<bool> final_jets;
    for (unsigned int jidx = 0; jidx < tmp_jet_flags.size(); jidx++) {

        if (!tmp_jet_flags.at(jidx)) {
            final_jets.push_back(tmp_jet_flags.at(jidx));
            continue;
        }

        const float jet_cor = (cms2.evt_isRealData() ? cms2.pfjets_corL1FastL2L3residual().at(jidx) :  cms2.pfjets_corL1FastL2L3().at(jidx));
        LorentzVector vjet  = cms2.pfjets_p4().at(jidx) * jet_cor; 
        if (systFlag != 0) {
            float c = getJetMetSyst(systFlag, vjet.pt(), vjet.eta());
            vjet *= c;
        }
        if (vjet.pt() < min_pt) {
            final_jets.push_back(false);
            continue;
        }

        bool jetIsLep = false;

        std::vector<LorentzVector> ele_p4s = os2012::getGoodElectrons(min_lep_pt);
        for (unsigned int eidx = 0; eidx < ele_p4s.size(); eidx++) {
            if (ROOT::Math::VectorUtil::DeltaR(vjet, ele_p4s.at(eidx)) > deltaR)
                continue;

            jetIsLep = true;
            break;
        }

        if (jetIsLep) {
            final_jets.push_back(false);
            continue;
        }

        std::vector<LorentzVector> mu_p4s = os2012::getGoodMuons(min_lep_pt);
        for (unsigned int midx = 0; midx < mu_p4s.size(); midx++) {
            if (ROOT::Math::VectorUtil::DeltaR(vjet, mu_p4s.at(midx)) > deltaR)
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


// JEC applied otf
std::vector<bool> os2012::getBtaggedJetFlags(int idx, FactorizedJetCorrector* jet_corrector, enum JetType type, enum BtagType btag_type, float deltaR, float min_pt, float max_eta, float min_lep_pt, int systFlag)
{
    std::vector<bool> tmp_jet_flags = os2012::getBtaggedJetFlags(idx, type, btag_type, deltaR, 0.0, max_eta, min_lep_pt, systFlag);

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
        vjet *= jet_cor;
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


// JEC uncertainty applied otf
std::vector<bool> os2012::getBtaggedJetFlags(int idx, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type, enum JetType type, enum BtagType btag_type, float deltaR, float min_pt, float max_eta, float min_lep_pt)
{
    std::vector<bool> tmp_jet_flags = os2012::getBtaggedJetFlags(idx, type, btag_type, deltaR, 0.0, max_eta, min_lep_pt);

    assert(tmp_jet_flags.size() == cms2.pfjets_p4().size());

    // now impose the pt requirement after applying the extra corrections
    std::vector<bool> final_jets;
    for (unsigned int jidx = 0; jidx < tmp_jet_flags.size(); jidx++) {

        if (!tmp_jet_flags.at(jidx)) {
            final_jets.push_back(tmp_jet_flags.at(jidx));
            continue;
        }
        const float jet_cor = (cms2.evt_isRealData() ? cms2.pfjets_corL1FastL2L3residual().at(jidx) :  cms2.pfjets_corL1FastL2L3().at(jidx));
        LorentzVector vjet  = cms2.pfjets_p4().at(jidx) * jet_cor; 
        jet_unc->setJetPt(vjet.pt());    
        jet_unc->setJetEta(vjet.eta());  
        const float jet_cor_unc = jet_unc->getUncertainty(true);     
        vjet *= (1.0 + jet_cor_unc * scale_type);    
        if (vjet.pt() < min_pt) {
            final_jets.push_back(false);
            continue;   
        }
        final_jets.push_back(true);
    }

    return final_jets;
}


// JEC AND JEC uncertainty applied otf
std::vector<bool> os2012::getBtaggedJetFlags(int idx, FactorizedJetCorrector* jet_corrector, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type, enum JetType type, enum BtagType btag_type, float deltaR, float min_pt, float max_eta, float min_lep_pt)
{
    std::vector<bool> tmp_jet_flags = os2012::getBtaggedJetFlags(idx, type, btag_type, deltaR, 0.0, max_eta, min_lep_pt);

    assert(tmp_jet_flags.size() == cms2.pfjets_p4().size());

    // now impose the pt requirement after applying the extra corrections
    std::vector<bool> final_jets;
    for (unsigned int jidx = 0; jidx < tmp_jet_flags.size(); jidx++) {

        if (!tmp_jet_flags.at(jidx)) {
            final_jets.push_back(tmp_jet_flags.at(jidx));
            continue;
        }
        jet_corrector->setRho(cms2.evt_ww_rho_vor());
        jet_corrector->setJetA(cms2.pfjets_area().at(jidx));
        jet_corrector->setJetPt(cms2.pfjets_p4().at(jidx).pt());
        jet_corrector->setJetEta(cms2.pfjets_p4().at(jidx).eta());        
        const float jet_cor = jet_corrector->getCorrection();
        LorentzVector vjet = cms2.pfjets_p4().at(jidx) * jet_cor;
        jet_unc->setJetPt(vjet.pt());    
        jet_unc->setJetEta(vjet.eta());  
        const float jet_cor_unc = jet_unc->getUncertainty(true);     
        vjet *= (1.0 + jet_cor_unc * scale_type);    
        if (vjet.pt() < min_pt) {
            final_jets.push_back(false);
            continue;   
        }
        final_jets.push_back(true);
    }

    return final_jets;
}


///////////////////////////////////////////////////////////////////////////////////////////
// 2012 get sumpt, skip jets overlapping with numerator e/mu with pt>x (defaults are 20/20 GeV)
///////////////////////////////////////////////////////////////////////////////////////////

// JEC taken from ntuple
int os2012::nBtaggedJets(int idx, enum JetType type, enum BtagType btag_type, float deltaR, float min_pt, float max_eta, float min_lep_pt, int systFlag)
{
    std::vector<LorentzVector> good_jets = os2012::getBtaggedJets(idx, type, btag_type, deltaR, min_pt, max_eta, min_lep_pt, systFlag);
    return good_jets.size();
}


// JEC applied otf
int os2012::nBtaggedJets(int idx, FactorizedJetCorrector* jet_corrector, enum JetType type, enum BtagType btag_type, float deltaR, float min_pt, float max_eta, float min_lep_pt, int systFlag)
{
    std::vector<LorentzVector> good_jets = os2012::getBtaggedJets(idx, jet_corrector, type, btag_type, deltaR, min_pt, max_eta, min_lep_pt, systFlag);
    return good_jets.size();    
}


// JEC uncertainty applied otf
int os2012::nBtaggedJets(int idx, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type, enum JetType type, enum BtagType btag_type, float deltaR, float min_pt, float max_eta, float min_lep_pt)     
{    
    std::vector<LorentzVector> good_btags = os2012::getBtaggedJets(idx, jet_unc, scale_type, type, btag_type, deltaR, min_pt, max_eta, min_lep_pt);
    return good_btags.size();    
}    

 
// JEC AND JEC uncertainty applied otf
int os2012::nBtaggedJets(int idx, FactorizedJetCorrector* jet_corrector, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type, enum JetType type, enum BtagType btag_type, float deltaR, float min_pt, float max_eta, float min_lep_pt)  
{    
    std::vector<LorentzVector> good_btags = os2012::getBtaggedJets(idx, jet_corrector, jet_unc, scale_type, type, btag_type, deltaR, min_pt, max_eta, min_lep_pt);
    return good_btags.size();    
}    


///////////////////////////////////////////////////////////////////////////////////////////
// 2012 get jet b-tag discriminators, skip jets overlapping with numerator e/mu with pt>x (defaults are 20/20 GeV)
///////////////////////////////////////////////////////////////////////////////////////////

// JEC taken from ntuple
std::vector<float> os2012::getJetBtagDiscriminators(int idx, enum JetType type, enum BtagType btag_type, float deltaR, float min_pt, float max_eta, float min_lep_pt, int systFlag, bool sort_by_pt)
{
    std::vector<LorentzVector> tmp_jet_p4s = os2012::getJets(idx, type, /*deltaR*/0.0, /*min_pt=*/-9999., /*max_eta=*/9999., /*min_lep_pt=*/9999., systFlag, /*sort_by_pt=*/false);
    std::vector<bool> tmp_jet_flags        = os2012::getJetFlags(idx, type, deltaR, min_pt, max_eta, min_lep_pt, systFlag);

    assert(tmp_jet_p4s.size() == tmp_jet_flags.size());

    std::vector<float> tmp_btag_disc = getJetBtagDiscriminators(btag_type, type);

    assert(tmp_jet_flags.size() == tmp_btag_disc.size());

    std::vector<float> ret;
    if (sort_by_pt)
    {
        std::vector<std::pair<LorentzVector, float> > tmp_jet_p4_disc;
        for (unsigned int jidx = 0; jidx < tmp_jet_flags.size(); jidx++)
        {
            if (!tmp_jet_flags.at(jidx)) continue;
            tmp_jet_p4_disc.push_back(std::make_pair(tmp_jet_p4s.at(jidx), tmp_btag_disc.at(jidx)));
        }
        
        sort(tmp_jet_p4_disc.begin(), tmp_jet_p4_disc.end(), SortByPt());
        
        for (unsigned int bidx = 0; bidx < tmp_jet_p4_disc.size(); bidx++)
        {
            ret.push_back(tmp_jet_p4_disc.at(bidx).second);
        }
    }
    else
    {
        for (unsigned int jidx = 0; jidx < tmp_jet_flags.size(); jidx++)
        {
            if (!tmp_jet_flags.at(jidx)) continue;
            ret.push_back(tmp_btag_disc.at(jidx));
        }
    }

    return ret;
}

// JEC applied otf
std::vector<float> os2012::getJetBtagDiscriminators(int idx, FactorizedJetCorrector* jet_corrector, enum JetType type, enum BtagType btag_type, float deltaR, float min_pt, float max_eta, float min_lep_pt, int systFlag, bool sort_by_pt)
{
    std::vector<LorentzVector> tmp_jet_p4s = os2012::getJets(idx, jet_corrector, type, /*deltaR*/0.0, /*min_pt=*/-9999., /*max_eta=*/9999., /*min_lep_pt=*/9999., systFlag, /*sort_by_pt=*/false);
    std::vector<bool> tmp_jet_flags        = os2012::getJetFlags(idx, jet_corrector, type, deltaR, min_pt, max_eta, min_lep_pt, systFlag);

    assert(tmp_jet_p4s.size() == tmp_jet_flags.size());

    std::vector<float> tmp_btag_disc = getJetBtagDiscriminators(btag_type, type);

    assert(tmp_jet_flags.size() == tmp_btag_disc.size());

    std::vector<float> ret;
    if (sort_by_pt)
    {
        std::vector<std::pair<LorentzVector, float> > tmp_jet_p4_disc;
        for (unsigned int jidx = 0; jidx < tmp_jet_flags.size(); jidx++)
        {
            if (!tmp_jet_flags.at(jidx)) continue;
            tmp_jet_p4_disc.push_back(std::make_pair(tmp_jet_p4s.at(jidx), tmp_btag_disc.at(jidx)));
        }
        
        sort(tmp_jet_p4_disc.begin(), tmp_jet_p4_disc.end(), SortByPt());
        
        for (unsigned int bidx = 0; bidx < tmp_jet_p4_disc.size(); bidx++)
        {
            ret.push_back(tmp_jet_p4_disc.at(bidx).second);
        }
    }
    else
    {
        for (unsigned int jidx = 0; jidx < tmp_jet_flags.size(); jidx++)
        {
            if (!tmp_jet_flags.at(jidx)) continue;
            ret.push_back(tmp_btag_disc.at(jidx));
        }
    }

    return ret;    
}

// JEC uncertainty applied otf
std::vector<float> os2012::getJetBtagDiscriminators(int idx, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type, enum JetType type, enum BtagType btag_type, float deltaR, float min_pt, float max_eta, float min_lep_pt, bool sort_by_pt)
{
    std::vector<LorentzVector> tmp_jet_p4s = os2012::getJets(idx, jet_unc, scale_type, type, /*deltaR*/0.0, /*min_pt=*/-9999., /*max_eta=*/9999., /*min_lep_pt=*/9999., /*sort_by_pt=*/false);
    std::vector<bool> tmp_jet_flags        = os2012::getJetFlags(idx, jet_unc, scale_type, type, deltaR, min_pt, max_eta, min_lep_pt);

    assert(tmp_jet_p4s.size() == tmp_jet_flags.size());

    std::vector<float> tmp_btag_disc = getJetBtagDiscriminators(btag_type, type);

    assert(tmp_jet_flags.size() == tmp_btag_disc.size());

    std::vector<float> ret;
    if (sort_by_pt)
    {
        std::vector<std::pair<LorentzVector, float> > tmp_jet_p4_disc;
        for (unsigned int jidx = 0; jidx < tmp_jet_flags.size(); jidx++)
        {
            if (!tmp_jet_flags.at(jidx)) continue;
            tmp_jet_p4_disc.push_back(std::make_pair(tmp_jet_p4s.at(jidx), tmp_btag_disc.at(jidx)));
        }
        
        sort(tmp_jet_p4_disc.begin(), tmp_jet_p4_disc.end(), SortByPt());
        
        for (unsigned int bidx = 0; bidx < tmp_jet_p4_disc.size(); bidx++)
        {
            ret.push_back(tmp_jet_p4_disc.at(bidx).second);
        }
    }
    else
    {
        for (unsigned int jidx = 0; jidx < tmp_jet_flags.size(); jidx++)
        {
            if (!tmp_jet_flags.at(jidx)) continue;
            ret.push_back(tmp_btag_disc.at(jidx));
        }
    }

    return ret;    
}

// JEC AND JEC uncertainty applied otf
std::vector<float> os2012::getJetBtagDiscriminators(int idx, FactorizedJetCorrector* jet_corrector, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type, enum JetType type, enum BtagType btag_type, float deltaR, float min_pt, float max_eta, float min_lep_pt, bool sort_by_pt)
{
    std::vector<LorentzVector> tmp_jet_p4s = os2012::getJets(idx, jet_corrector, jet_unc, scale_type, type, /*deltaR*/0.0, /*min_pt=*/-9999., /*max_eta=*/9999., /*min_lep_pt=*/9999., /*sort_by_pt=*/false);
    std::vector<bool> tmp_jet_flags        = os2012::getJetFlags(idx, jet_corrector, jet_unc, scale_type, type, deltaR, min_pt, max_eta, min_lep_pt);

    assert(tmp_jet_p4s.size() == tmp_jet_flags.size());

    std::vector<float> tmp_btag_disc = getJetBtagDiscriminators(btag_type, type);

    assert(tmp_jet_flags.size() == tmp_btag_disc.size());

    std::vector<float> ret;
    if (sort_by_pt)
    {
        std::vector<std::pair<LorentzVector, float> > tmp_jet_p4_disc;
        for (unsigned int jidx = 0; jidx < tmp_jet_flags.size(); jidx++)
        {
            if (!tmp_jet_flags.at(jidx)) continue;
            tmp_jet_p4_disc.push_back(std::make_pair(tmp_jet_p4s.at(jidx), tmp_btag_disc.at(jidx)));
        }
        
        sort(tmp_jet_p4_disc.begin(), tmp_jet_p4_disc.end(), SortByPt());
        
        for (unsigned int bidx = 0; bidx < tmp_jet_p4_disc.size(); bidx++)
        {
            ret.push_back(tmp_jet_p4_disc.at(bidx).second);
        }
    }
    else
    {
        for (unsigned int jidx = 0; jidx < tmp_jet_flags.size(); jidx++)
        {
            if (!tmp_jet_flags.at(jidx)) continue;
            ret.push_back(tmp_btag_disc.at(jidx));
        }
    }

    return ret;    
}



///////////////////////////////////////////////////////////////////////////////////////////
// 2012 rescale the jet energy resolution (JER) 
// function that returns the sigma(pT)*pT of the MC jets. those numbers are old,
// but at least they are in an understandable format.
///////////////////////////////////////////////////////////////////////////////////////////
float os2012::getErrPt(const float pt, const float eta)
{
    float InvPerr2;
    float N = 0.0;
    float S = 0.0;
    float C = 0.0;
    float m = 0.0;
    if(fabs(eta) < 0.5)
    {
        N = 3.96859;
        S = 0.18348;
        C = 0.;
        m = 0.62627;
    } 
    else if(fabs(eta) < 1.0)
    {
        N = 3.55226;
        S = 0.24026;
        C = 0.;
        m = 0.52571;
    } 
    else if(fabs(eta) < 1.5) 
    {
        N = 4.54826;
        S = 0.22652;
        C = 0.;
        m = 0.58963;
    }
    else if(fabs(eta) < 2.0)
    {
        N = 4.62622;
        S = 0.23664;
        C = 0.;
        m = 0.48738;
    }
    else if(fabs(eta) < 2.5)
    {
        N = 2.53324;
        S = 0.34306;
        C = 0.;
        m = 0.28662;
    }
    else if(fabs(eta) < 3.0)
    {
        N = -3.33814;
        S = 0.73360;
        C = 0.;
        m = 0.08264;
    }
    else if(fabs(eta) < 5.0)
    {
        N = 2.95397;
        S = 0.11619;
        C = 0.;
        m = 0.96086;
    }
    
    // this is the absolute resolution (squared), not sigma(pt)/pt
    // so have to multiply by pt^2, thats why m+1 instead of m-1
    InvPerr2 = (N * fabs(N) ) + (S * S) * pow(pt, m+1) + (C * C) * pt * pt;

    return sqrt(InvPerr2);
}


///////////////////////////////////////////////////////////////////////////////////////////
// function to get the jer scale factors. values taken from the twiki: 
// https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
///////////////////////////////////////////////////////////////////////////////////////////
float os2012::getJERScale(const float jet_eta)
{
    const float aeta = fabs(jet_eta);
    if      (aeta < 0.5) {return 1.052;}
    else if (aeta < 1.1) {return 1.057;}
    else if (aeta < 1.7) {return 1.096;}
    else if (aeta < 2.3) {return 1.134;}
    else                 {return 1.288;}

}

// rescaled the jet p4s, met, met_phi and ht scaling up the JER
void os2012::smearJETScaleJetsMetHt(std::vector<LorentzVector>& vjets_p4, float& met, float& met_phi, float& ht, const unsigned int seed)
{
    static TRandom3 random(seed);
    float new_ht = 0;

    // rescale the jets/met/ht
    ROOT::Math::XYVector cmet(met*cos(met_phi), met*sin(met_phi));
    std::vector<LorentzVector> new_vjets_p4;
    for (size_t jidx = 0; jidx != vjets_p4.size(); jidx++)
    {
        random.SetSeed(seed*(jidx+1));
    
        // rescale the jet pt
        const LorentzVector& jet_p4 = vjets_p4.at(jidx);
        const float jer_scale       = getJERScale(jet_p4.eta());
        const float sigma_mc        = getErrPt(jet_p4.pt(), jet_p4.eta())/jet_p4.pt();
        const float jet_rescaled    = random.Gaus(1.0, sqrt(jer_scale*jer_scale-1.0)*sigma_mc);
        LorentzVector new_jet_p4    = (jet_p4 * jet_rescaled);

        // propogate to the met
        ROOT::Math::XYVector old_jet(jet_p4.px(), jet_p4.py());
        ROOT::Math::XYVector new_jet(new_jet_p4.px(), new_jet_p4.py());
        cmet = cmet - new_jet + old_jet;

        // check that the new jets pass the min pt cut
        if (new_jet_p4.pt() < 40.0)
        {
            continue;
        }

        // add to the new ht
        new_ht += new_jet_p4.pt();

        // return the new jet
        new_vjets_p4.push_back(LorentzVector(new_jet_p4));
    }
    
    // set the new met
    met     = cmet.r();
    met_phi = cmet.phi();

    // set the new pt
    ht = new_ht;

    // set the new jets
    vjets_p4 = new_vjets_p4;

    // done
    return;
}       
            
void os2012::smearJETScaleJetsMetHt
(
    std::vector<LorentzVector>& vjets_p4, 
    float& met,
    float& met_phi,
    float& ht,
    int idx,
    enum JetType type,
    const unsigned int seed,
    float deltaR,
    float min_pt,
    float max_eta,
    float lep_minpt
)
{
    static TRandom3 random(seed);
    float new_ht = 0;

    // rescale the jets/met/ht
    ROOT::Math::XYVector cmet(met*cos(met_phi), met*sin(met_phi));
    std::vector<LorentzVector> new_vjets_p4;
    std::vector<LorentzVector> tmp_vjets_p4 = os2012::getJets(idx, type, deltaR, /*min_pt=*/10, /*max_eta=*/2.4, lep_minpt);
    for (size_t jidx = 0; jidx != tmp_vjets_p4.size(); jidx++)
    {
        random.SetSeed(seed*(jidx+1));

        // rescale the jet pt
        const LorentzVector& jet_p4 = tmp_vjets_p4.at(jidx);
        const float jer_scale       = getJERScale(jet_p4.eta());
        const float sigma_mc        = getErrPt(jet_p4.pt(), jet_p4.eta())/jet_p4.pt();
        const float arg             = sqrt(jer_scale*jer_scale-1.0)*sigma_mc;
        const float jet_rescaled    = random.Gaus(1.0, arg);
        LorentzVector new_jet_p4    = (jet_p4 * jet_rescaled);

        // propogate to the met
        ROOT::Math::XYVector old_jet(jet_p4.px(), jet_p4.py());
        ROOT::Math::XYVector new_jet(new_jet_p4.px(), new_jet_p4.py());
        cmet = cmet - new_jet + old_jet;

        // check that the new jets pass the min pt cut
        if (not (new_jet_p4.pt() > min_pt and fabs(new_jet_p4.eta()) < max_eta))
        {
            continue;
        }

        // add to the new ht
        new_ht += new_jet_p4.pt();

        // return the new jet
        new_vjets_p4.push_back(LorentzVector(new_jet_p4));
    }
    
    // set the new met
    met     = cmet.r();
    met_phi = cmet.phi();

    // set the new pt
    ht = new_ht;

    // set the new jets
    vjets_p4 = new_vjets_p4;

    // done
    return;
}


// semar JER for jets
void os2012::smearJETScaleJets(std::vector<LorentzVector>& vjets_p4, const unsigned int seed)
{
    static TRandom3 random(seed);

    // rescale the b-tagged jets
    std::vector<LorentzVector> new_vjets_p4;
    for (size_t jidx = 0; jidx != vjets_p4.size(); jidx++)
    {
        random.SetSeed(seed*(jidx+1));

        // rescale the jet pt
        const LorentzVector& jet_p4 = vjets_p4.at(jidx);
        const float jer_scale       = getJERScale(jet_p4.eta());
        const float sigma_mc        = getErrPt(jet_p4.pt(), jet_p4.eta())/jet_p4.pt();
        const float jet_rescaled    = random.Gaus(1.0, sqrt(jer_scale*jer_scale-1.0)*sigma_mc);
        LorentzVector new_jet_p4    = (jet_p4 * jet_rescaled);

        // check that the new jets pass the min pt cut
        if (new_jet_p4.pt() < 40.0)
        {
            continue;
        }

        // return the new jet
        new_vjets_p4.push_back(new_jet_p4);
    }
    vjets_p4 = new_vjets_p4;

    // done
    return;
}       


///////////////////////////////////////////////////////////////////////////////////////////
// 2012 rescale the MET by scaling up/down the unclustered erngy 
///////////////////////////////////////////////////////////////////////////////////////////
float os2012::scaleMET
(
    const float met,
    const float met_phi,
    int idx,
    enum JetType type,
    float deltaR,
    float min_pt,
    float max_eta,
    float lep_minpt,
    const int scale_type,
    const float scale
)
{
    using namespace tas;
    typedef ROOT::Math::Polar2DVectorF Polar2D;

    // 2D vector to keep the sums
    Polar2D jets;
    Polar2D leps;

    // sum up the jets 
    std::vector<LorentzVector> vjets_p4 = os2012::getJets(idx, type, deltaR, min_pt, max_eta, lep_minpt);
    for (size_t jidx = 0; jidx != vjets_p4.size(); jidx++)
    {
        jets += Polar2D(vjets_p4.at(jidx).pt(), vjets_p4.at(jidx).phi());
    }
    // sum up the electrons
    for (size_t eidx = 0; eidx < els_p4().size(); eidx++)
    {
        if (els_p4().at(eidx).pt() < lep_minpt)     {continue;}
        if (!os2012::isNumeratorLepton(11, eidx)) {continue;}
        leps += Polar2D(els_p4().at(eidx).pt(), els_p4().at(eidx).phi());
    }
    // sum up the muons
    for (size_t midx = 0; midx < mus_p4().size(); midx++)
    {
        if (mus_p4().at(midx).pt() < lep_minpt)      {continue;}
        if (!os2012::isNumeratorLepton(13, midx)) {continue;}
        leps += Polar2D(mus_p4().at(midx).pt(), mus_p4().at(midx).phi());
    }
    
    // subtract the leptons and jet contributions
    Polar2D umet(met, met_phi);
    umet = umet + leps + jets;

    // scale the unclustered energy by 10%
    umet.SetR(umet.r() * (1.0 + scale_type * scale));

    // resum the met
    Polar2D new_met = umet - jets - leps;
    return new_met.r();
}

     
///////////////////////////////////////////////////////////////////////////////////////////  
// 2012 get vector of good els p4s   
///////////////////////////////////////////////////////////////////////////////////////////  
std::vector<LorentzVector> os2012::getGoodElectrons(const float ptcut)     
{    
    std::vector<LorentzVector> good_els;     
    for (unsigned int idx = 0; idx < cms2.els_p4().size(); idx++) {  
     
        if (cms2.els_p4().at(idx).pt() < ptcut) continue;    
        if (fabs(cms2.els_p4().at(idx).eta()) > 2.4) continue;
        if (!os2012::isNumeratorLepton(11, idx)) continue;   
         
        good_els.push_back(cms2.els_p4().at(idx));   
    }    
         
    sort(good_els.begin(), good_els.end(), SortByPt());  
    return good_els;     
}    


std::vector<std::pair<LorentzVector, unsigned int> > os2012::getNumeratorElectrons(const float ptcut)
{    
    std::vector<std::pair<LorentzVector, unsigned int> > good_els;   
    for (unsigned int idx = 0; idx < cms2.els_p4().size(); idx++)
    {    
        if (cms2.els_p4().at(idx).pt() < ptcut)      {continue;}
        if (fabs(cms2.els_p4().at(idx).eta()) > 2.4) {continue;}
        if (!os2012::isNumeratorLepton(11, idx))   {continue;}
        good_els.push_back(std::make_pair(cms2.els_p4().at(idx), idx));  
    }    
         
    sort(good_els.begin(), good_els.end(), SortByPt());  
    return good_els;     
}    

         
         
///////////////////////////////////////////////////////////////////////////////////////////  
// 2012 get vector of good mus p4s   
///////////////////////////////////////////////////////////////////////////////////////////  
std::vector<LorentzVector> os2012::getGoodMuons(const float ptcut)     
{    
    std::vector<LorentzVector> good_mus;     
    for (unsigned int idx = 0; idx < cms2.mus_p4().size(); idx++) {  
         
        if (cms2.mus_p4().at(idx).pt() < ptcut) continue;
        if (!os2012::isNumeratorLepton(13, idx)) continue;
        good_mus.push_back(cms2.mus_p4().at(idx));   
    }    
         
    sort(good_mus.begin(), good_mus.end(), SortByPt());  
    return good_mus;     
}

std::vector<std::pair<LorentzVector, unsigned int> > os2012::getNumeratorMuons(const float ptcut)
{    
    std::vector<std::pair<LorentzVector, unsigned int> > good_mus;   
    for (unsigned int idx = 0; idx < cms2.mus_p4().size(); idx++)
    {    
        if (cms2.mus_p4().at(idx).pt() < ptcut)      {continue;}
        if (!os2012::isNumeratorLepton(13, idx))   {continue;}
        good_mus.push_back(std::make_pair(cms2.mus_p4().at(idx), idx));  
    }    
    sort(good_mus.begin(), good_mus.end(), SortByPt());  
    return good_mus;     
}    
