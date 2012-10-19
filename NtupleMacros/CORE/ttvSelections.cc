#include <assert.h>
#include <iostream>
#include <algorithm>

//core headers
#include "ttvSelections.h"
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

//Root headers
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TDatabasePDG.h"

using namespace ttv;

struct jet_pt_gt {
    bool operator () (const LorentzVector &v1, const LorentzVector &v2) 
    {
        return v1.pt() > v2.pt();
    }
};

bool ttv::passesTrigger (int hyp_type)
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

bool ttv::isGoodLepton (int id, int idx, LeptonType::value_type lep_type)
{
    if (abs(id) == 11)
    {
        if (lep_type == ttv::LeptonType::LOOSE)
        {
		  if( fabs(cms2.els_etaSC()[idx]) > 1.4442 && fabs(cms2.els_etaSC()[idx]) < 1.566 )    return false;
		  if( fabs(cms2.els_p4()[idx].eta()) > 2.4 )                                          return false;
		  if( overlapMuon(idx, ttv::LeptonType::LOOSE, 10.0) )                                                         return false; // is 10 GeV, eta 2.4, and deltaR 0.4 good?
		  electronIdComponent_t answer_loose_2012 = electronId_WP2012_v2(idx, LOOSE);
		  if ((answer_loose_2012 & wp2012::PassWP2012CutsNoIso) == wp2012::PassWP2012CutsNoIso) return true; 
		  return false;
        }
        else if (lep_type == ttv::LeptonType::TIGHT)
        {
		  return pass_electronSelection(idx, electronSelection_ssV7_noIso);
        }
		else if (lep_type == ttv::LeptonType::LOOSEMVA)
		  {
			std::cout<<"Warning: in ttvSelections bool ttv::isGoodLepton, enum LOOSEMVA not supported for electrons yet."<<std::endl;
		  }
		else if (lep_type == ttv::LeptonType::TIGHTMVA)
		  {
			std::cout<<"Warning: in ttvSelections bool ttv::isGoodLepton, enum TIGHTMVA not supported for electrons yet."<<std::endl;
		  }
		else
		  {
			std::cout<<"Warning: in ttvSelections bool ttv::isGoodLepton, unsupported enum "<<lep_type<<" for electrons."<<std::endl;
		  }
    }
    else if (abs(id) == 13)
    {
        if (lep_type == ttv::LeptonType::LOOSE)
        {
		  return (muonIdNotIsolated(idx, NominalTTZ_loose_v1));
        }
        else if (lep_type == ttv::LeptonType::TIGHT)
        {
		  return (muonIdNotIsolated(idx, NominalTTZ_tight_v1));
        }
		else
		  {
			std::cout<<"Warning: in ttvSelections bool ttv::isGoodLepton, unsupported enum "<<lep_type<<" for muons."<<endl;
		  }
    }

    return false;
}


bool ttv::isIsolatedLepton (int id, int idx, LeptonType::value_type lep_type)
{
    if (abs(id) == 11)
    {
        if (lep_type == ttv::LeptonType::LOOSE)
        {
		  return electronIsoValuePF2012_FastJetEffArea_v2( idx ) < 0.15;
        }
        else if (lep_type == ttv::LeptonType::TIGHT)
        {
		  return samesign::electronIsolationPF2012(idx) < 0.09;
        }
    }
    else if (abs(id) == 13)
    {
        if (lep_type == ttv::LeptonType::LOOSE)
        {
            return (muonIsoValuePF2012_deltaBeta(idx) < 0.15);
        }
        else if (lep_type == ttv::LeptonType::TIGHT)
        {
            return (muonIsoValuePF2012_deltaBeta(idx) < 0.10);
        }
    }

    return false;
}

bool ttv::isNumeratorLepton (int id, int idx, LeptonType::value_type lep_type)
{
    return (ttv::isGoodLepton(id, idx, lep_type) && ttv::isIsolatedLepton(id, idx, lep_type));
}

bool ttv::isDenominatorLepton (int id, int idx, LeptonType::value_type lep_type)
{
    if (abs(id) == 11)
    {
        if (lep_type == ttv::LeptonType::LOOSE)
        {
		  // still need to do this one
        }
        else if (lep_type == ttv::LeptonType::TIGHT)
        {
		  return pass_electronSelection(idx, electronSelectionFOV7_v1);
        }
    }
    else if (abs(id) == 13)
    {
        if (lep_type == ttv::LeptonType::LOOSE)
        {
            return (muonId(idx, NominalTTZ_looseFO_v1));
        }
        else if (lep_type == ttv::LeptonType::TIGHT)
        {
            return (muonId(idx, NominalTTZ_tightFO_v1));
        }
    }

    return false;
}

bool overlapMuon(int idx, LeptonType::value_type lep_type, float pt, float eta, float deltaR){
  for( unsigned int muidx = 0; muidx < cms2.mus_p4().size(); muidx++){
	if (cms2.mus_p4()[muidx].pt() > pt &&
		abs(cms2.mus_p4()[muidx].eta()) < eta &&
		ROOT::Math::VectorUtil::DeltaR(cms2.mus_p4()[muidx], cms2.els_p4()[idx]) < deltaR){ 
	  return true; 
	}
  }
  return false;
}

std::vector<LorentzVector> ttv::getJets(std::vector<LorentzVector>& leps, enum JetType type, float deltaR, float min_pt, float max_eta, float rescale, int systFlag)
{

    std::vector<LorentzVector> tmp_jets = getJets(999999, true, type, JETS_CLEAN_NONE, deltaR, min_pt, max_eta, (double) rescale, systFlag);

    // ok, now perform the rest of the lepton overlap removal
    // and the impose the pt requirement after applying the
    // extra corrections
    std::vector<LorentzVector> final_jets;

    for (unsigned int jidx = 0; jidx < tmp_jets.size(); jidx++)
    {
        LorentzVector vjet = tmp_jets.at(jidx);
        bool jetIsLep = false;
        for (unsigned int lidx = 0; lidx < leps.size(); lidx++)
        {
            if (ROOT::Math::VectorUtil::DeltaR(vjet, leps.at(lidx)) < deltaR)
            {
                jetIsLep = true;
                break;
            }
        }

        if (jetIsLep) continue;        

        final_jets.push_back(vjet);
    }

    sort(final_jets.begin(), final_jets.end(), jet_pt_gt());
    return final_jets;    
}

std::vector<LorentzVector> ttv::getJets(std::vector<LorentzVector>& leps, FactorizedJetCorrector* jet_corrector, enum JetType type, float deltaR, float min_pt, float max_eta, float rescale, int systFlag)
{
    std::vector<bool> tmp_jet_flags = ttv::getJetFlags(leps, type, deltaR, 0., max_eta, 1., 0);

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
        double jet_cor = jet_corrector->getCorrection();

        double unc = (systFlag == 0) ? 1. : getJetMetSyst(systFlag, vjet.pt(), vjet.eta());        
        vjet *= jet_cor * rescale * unc;
        if (vjet.pt() < min_pt) continue;

        final_jets.push_back(vjet);
    }

    sort(final_jets.begin(), final_jets.end(), jet_pt_gt());
    return final_jets;
}

std::vector<bool> ttv::getJetFlags(std::vector<LorentzVector>& leps, enum JetType type, float deltaR, float min_pt, float max_eta, float rescale, int systFlag)
{
    std::vector<bool> tmp_jet_flags = getJetFlags(999999, type, JETS_CLEAN_NONE, (double)deltaR, min_pt, (double)max_eta, (double)rescale, systFlag);
    std::vector<LorentzVector> tmp_jets = getJets(999999, false, type, JETS_CLEAN_NONE, deltaR, min_pt, max_eta, (double) rescale, systFlag);

    // ok, now perform the rest of the lepton overlap removal
    // and the impose the pt requirement after applying the
    // extra corrections
    std::vector<bool> final_jets;
    for (unsigned int jidx = 0; jidx < tmp_jet_flags.size(); jidx++) {

        if (!tmp_jet_flags.at(jidx)) {
            final_jets.push_back(tmp_jet_flags.at(jidx));
            continue;
        }
        
        LorentzVector vjet = tmp_jets.at(jidx);
        bool jetIsLep = false;
        for (unsigned int lidx = 0; lidx < leps.size(); lidx++)
        {
            if (ROOT::Math::VectorUtil::DeltaR(vjet, leps.at(lidx)) < deltaR)
            {
                jetIsLep = true;
                break;
            }
        }
            
        final_jets.push_back(!jetIsLep);
    }

    return final_jets;
}

std::vector<bool> ttv::getJetFlags(std::vector<LorentzVector>& leps, FactorizedJetCorrector* jet_corrector, enum JetType type, float deltaR, float min_pt, float max_eta, float rescale, int systFlag)
{
    std::vector<bool> tmp_jet_flags = ttv::getJetFlags(leps, type, deltaR, 0., max_eta, 1., 0);
    std::vector<LorentzVector> *tmp_jets = NULL;
    if (type <= JETS_TYPE_PF_UNCORR)
        tmp_jets = &cms2.pfjets_p4();
    else if (type <= JETS_TYPE_CALO_UNCORR)
        tmp_jets = &cms2.jets_p4();
    else
        tmp_jets = &cms2.genjets_p4();

    assert(tmp_jets->size() == tmp_jet_flags.size());

    // ok, now perform the rest of the lepton overlap removal
    // and the impose the pt requirement after applying the
    // extra corrections
    std::vector<bool> final_jets;
    for (unsigned int jidx = 0; jidx < tmp_jet_flags.size(); jidx++) {

        if (!tmp_jet_flags.at(jidx)) {
            final_jets.push_back(tmp_jet_flags.at(jidx));
            continue;
        }
        
        LorentzVector vjet = tmp_jets->at(jidx);
        jet_corrector->setRho(cms2.evt_ww_rho_vor());
        jet_corrector->setJetA(cms2.pfjets_area().at(jidx));
        jet_corrector->setJetPt(cms2.pfjets_p4().at(jidx).pt());
        jet_corrector->setJetEta(cms2.pfjets_p4().at(jidx).eta());        
        double jet_cor = jet_corrector->getCorrection();

        double unc = (systFlag == 0) ? 1. : getJetMetSyst(systFlag, vjet.pt(), vjet.eta());        
        vjet *= jet_cor * rescale * unc;
        bool is_good_jet = (vjet.pt() >= min_pt);

        final_jets.push_back(is_good_jet);
    }

    return final_jets;
}
    
float ttv::sumJetPt(std::vector<LorentzVector>& leps, enum JetType type, float deltaR, float min_pt, float max_eta, float rescale, int systFlag)
{
    std::vector<LorentzVector> tmp_jets = ttv::getJets(leps, type, deltaR, min_pt, max_eta, rescale, systFlag);
    float ht = 0.;
    for (unsigned int jidx = 0; jidx < tmp_jets.size(); jidx++)
        ht += tmp_jets.at(jidx).pt();

    return ht;
}


float ttv::sumJetPt(std::vector<LorentzVector>& leps, FactorizedJetCorrector* jet_corrector, enum JetType type, float deltaR, float min_pt, float max_eta, float rescale, int systFlag)
{
    std::vector<LorentzVector> tmp_jets = ttv::getJets(leps, jet_corrector, type, deltaR, min_pt, max_eta, rescale, systFlag);
    float ht = 0.;
    for (unsigned int jidx = 0; jidx < tmp_jets.size(); jidx++)
        ht += tmp_jets.at(jidx).pt();

    return ht;
}

int ttv::nJets(std::vector<LorentzVector>& leps, enum JetType type, float deltaR, float min_pt, float max_eta, float rescale, int systFlag)
{
    std::vector<LorentzVector> tmp_jets = ttv::getJets(leps, type, deltaR, min_pt, max_eta, rescale, systFlag);
    return tmp_jets.size();
}

int ttv::nJets(std::vector<LorentzVector>& leps, FactorizedJetCorrector* jet_corrector, enum JetType type, float deltaR, float min_pt, float max_eta, float rescale, int systFlag)
{
    std::vector<LorentzVector> tmp_jets = ttv::getJets(leps, jet_corrector, type, deltaR, min_pt, max_eta, rescale, systFlag);
    return tmp_jets.size();
}

std::vector<LorentzVector> ttv::getBtaggedJets(std::vector<LorentzVector>& leps, enum JetType type, enum BtagType btag_type, float deltaR, float min_pt, float max_eta, float rescale, int systFlag)
{
    std::vector<LorentzVector> tmp_jets = getBtaggedJets(999999, true, type, JETS_CLEAN_NONE, btag_type, deltaR, min_pt, max_eta, (double) rescale, systFlag);

    // ok, now perform the rest of the lepton overlap removal
    // and the impose the pt requirement after applying the
    // extra corrections
    std::vector<LorentzVector> final_jets;
    for (unsigned int jidx = 0; jidx < tmp_jets.size(); jidx++)
    {
        LorentzVector vjet = tmp_jets.at(jidx);
        bool jetIsLep = false;
        for (unsigned int lidx = 0; lidx < leps.size(); lidx++)
        {
            if (ROOT::Math::VectorUtil::DeltaR(vjet, leps.at(lidx)) < deltaR)
            {
                jetIsLep = true;
                break;
            }
        }

        if (jetIsLep) continue;        

        final_jets.push_back(vjet);
    }

    sort(final_jets.begin(), final_jets.end(), jet_pt_gt());
    return final_jets;    
}

std::vector<LorentzVector> ttv::getBtaggedJets(std::vector<LorentzVector>& leps, FactorizedJetCorrector* jet_corrector, enum JetType type, enum BtagType btag_type, float deltaR, float min_pt, float max_eta, float rescale, int systFlag)
{
    std::vector<bool> tmp_jet_flags = ttv::getBtaggedJetFlags(leps, type, btag_type, deltaR, 0., max_eta, 1., 0);

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
        double jet_cor = jet_corrector->getCorrection();

        double unc = (systFlag == 0) ? 1. : getJetMetSyst(systFlag, vjet.pt(), vjet.eta());        
        vjet *= jet_cor * rescale * unc;
        if (vjet.pt() < min_pt) continue;

        final_jets.push_back(vjet);
    }

    sort(final_jets.begin(), final_jets.end(), jet_pt_gt());
    return final_jets;
}

std::vector<bool> ttv::getBtaggedJetFlags(std::vector<LorentzVector>& leps, enum JetType type, enum BtagType btag_type, float deltaR, float min_pt, float max_eta, float rescale, int systFlag)
{
  std::vector<bool> tmp_jet_flags = getBtaggedJetFlags(999999, type, JETS_CLEAN_NONE, btag_type, (double)deltaR, (double)min_pt, (double)max_eta, (double)rescale, systFlag);
    std::vector<LorentzVector> tmp_jets = getBtaggedJets(999999, false, type, JETS_CLEAN_NONE, btag_type, deltaR, min_pt, max_eta, (double) rescale, systFlag);

    // ok, now perform the rest of the lepton overlap removal
    // and the impose the pt requirement after applying the
    // extra corrections
    std::vector<bool> final_jets;
    for (unsigned int jidx = 0; jidx < tmp_jet_flags.size(); jidx++) {

        if (!tmp_jet_flags.at(jidx)) {
            final_jets.push_back(tmp_jet_flags.at(jidx));
            continue;
        }
        
        LorentzVector vjet = tmp_jets.at(jidx);
        bool jetIsLep = false;
        for (unsigned int lidx = 0; lidx < leps.size(); lidx++)
        {
            if (ROOT::Math::VectorUtil::DeltaR(vjet, leps.at(lidx)) < deltaR)
            {
                jetIsLep = true;
                break;
            }
        }
            
        final_jets.push_back(!jetIsLep);
    }

    return final_jets;
}

std::vector<bool> ttv::getBtaggedJetFlags(std::vector<LorentzVector>& leps, FactorizedJetCorrector* jet_corrector, enum JetType type, enum BtagType btag_type, float deltaR, float min_pt, float max_eta, float rescale, int systFlag)
{
  std::vector<bool> tmp_jet_flags = ttv::getBtaggedJetFlags(leps, type, btag_type, deltaR, 0., max_eta, 1., 0);
    std::vector<LorentzVector> *tmp_jets = NULL;
    if (type <= JETS_TYPE_PF_UNCORR)
        tmp_jets = &cms2.pfjets_p4();
    else if (type <= JETS_TYPE_CALO_UNCORR)
        tmp_jets = &cms2.jets_p4();
    else
        tmp_jets = &cms2.genjets_p4();

    assert(tmp_jets->size() == tmp_jet_flags.size());

    // ok, now perform the rest of the lepton overlap removal
    // and the impose the pt requirement after applying the
    // extra corrections
    std::vector<bool> final_jets;
    for (unsigned int jidx = 0; jidx < tmp_jet_flags.size(); jidx++) {

        if (!tmp_jet_flags.at(jidx)) {
            final_jets.push_back(tmp_jet_flags.at(jidx));
            continue;
        }
        
        LorentzVector vjet = tmp_jets->at(jidx);
        jet_corrector->setRho(cms2.evt_ww_rho_vor());
        jet_corrector->setJetA(cms2.pfjets_area().at(jidx));
        jet_corrector->setJetPt(cms2.pfjets_p4().at(jidx).pt());
        jet_corrector->setJetEta(cms2.pfjets_p4().at(jidx).eta());        
        double jet_cor = jet_corrector->getCorrection();

        double unc = (systFlag == 0) ? 1. : getJetMetSyst(systFlag, vjet.pt(), vjet.eta());        
        vjet *= jet_cor * rescale * unc;
        bool is_good_jet = (vjet.pt() >= min_pt);

        final_jets.push_back(is_good_jet);
    }

    return final_jets;
}

int ttv::nBtaggedJets(std::vector<LorentzVector>& leps, enum JetType type, enum BtagType btag_type, float deltaR, float min_pt, float max_eta, float rescale, int systFlag)
{
    std::vector<LorentzVector> tmp_jets = ttv::getBtaggedJets(leps, type, btag_type, deltaR, min_pt, max_eta, rescale, systFlag);
    return tmp_jets.size();
}

int ttv::nBtaggedJets(std::vector<LorentzVector>& leps, FactorizedJetCorrector* jet_corrector, enum JetType type, enum BtagType btag_type, float deltaR, float min_pt, float max_eta, float rescale, int systFlag)
{
    std::vector<LorentzVector> tmp_jets = ttv::getBtaggedJets(leps, type, btag_type, deltaR, min_pt, max_eta, rescale, systFlag);
    return tmp_jets.size();
}

