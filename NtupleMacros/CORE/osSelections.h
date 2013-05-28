#ifndef OSSELECTIONS_H
#define OSSELECTIONS_H

#include <vector>
#include <string>
#include "Math/LorentzVector.h"
#include "CMS2.h"
#include "jetSelections.h"
#include "electronSelections.h"
#include "jetcorr/JetCorrectionUncertainty.h"

/////////////////////////////////////////////////////////////////
///                                                           ///
///                                                           ///
///                                                           ///
///          2012 Selections                                  ///
///                                                           ///
///                                                           ///
///                                                           ///
/////////////////////////////////////////////////////////////////

namespace os2012
{
    ////////////////////////////////////////////////////////////////////////////////////////////     
    // 2012 good lepton (passes ID)
    ////////////////////////////////////////////////////////////////////////////////////////////     
    bool isGoodLepton(int id, int idx);


    ////////////////////////////////////////////////////////////////////////////////////////////     
    // 2012 isolated lepton
    ////////////////////////////////////////////////////////////////////////////////////////////     
    bool isIsolatedLepton(int id, int idx);


    ////////////////////////////////////////////////////////////////////////////////////////////     
    // 2012 lepton isolation value
    ////////////////////////////////////////////////////////////////////////////////////////////     
    double leptonIsolation(int id, int idx);


    ////////////////////////////////////////////////////////////////////////////////////////////     
    // 2012 effective area 
    ////////////////////////////////////////////////////////////////////////////////////////////     
    float EffectiveArea03(int id, int idx);
    float EffectiveArea04(int id, int idx);


    ////////////////////////////////////////////////////////////////////////////////////////////     
    // 2012 numerator lepton (passes ID and isolation)
    ////////////////////////////////////////////////////////////////////////////////////////////     
    bool isNumeratorLepton(int id, int idx);


    ////////////////////////////////////////////////////////////////////////////////////////////     
    // 2012 numerator hypothesis (passes ID and isolation)
    ////////////////////////////////////////////////////////////////////////////////////////////     
    bool isNumeratorHypothesis(int idx);


    ///////////////////////////////////////////////////////////////////////////////////////////
    // 2012 calculate PF-based isolation for electrons with rho*Aeff correction
    ///////////////////////////////////////////////////////////////////////////////////////////
    float electronIsolationPF2012_cone03(int idx);  // uses ∆R < 0.3
    float electronIsolationPF2012_cone04(int idx);  // uses ∆R < 0.4
    float electronIsolationPF2012(int idx);         // wrapper ∆R < 0.3 version which is used in the analysis


    ///////////////////////////////////////////////////////////////////////////////////////////
    // 2012 passes dilepton trigger
    ///////////////////////////////////////////////////////////////////////////////////////////

    // analysis type:
    //   0 --> use high pT analysis triggers
    //   anything else will return false

    bool passesTrigger(int hyp_type, int analysis_type);
    bool passesTriggerHighPt(int hyp_type);

    ///////////////////////////////////////////////////////////////////////////////////////////
    // JET CONSTANTS - DEFINE HERE TO MAKE EASY TO CHANGE    
    ///////////////////////////////////////////////////////////////////////////////////////////
    // const float deltaR_cut     = 0.4;
    // const float min_pt_cut     = 30.;
    // const float min_lep_pt_cut = 20.;
    // const int systFlag_val     = 0;
    // const bool sort_by_pt_val  = true;


    ///////////////////////////////////////////////////////////////////////////////////////////
    // 2012 get jets and perform overlap removal with numerator e/mu with pt > x (defaults are 20/20 GeV)
    ///////////////////////////////////////////////////////////////////////////////////////////

    // JEC taken from ntuple
    std::vector<LorentzVector> getJets(int idx, enum JetType type, float deltaR = 0.4, float min_pt = 30.0, float max_eta = 2.4, float min_lep_pt = 20.0, int systFlag = 0, bool sort_by_pt = true);

    // JEC applied otf
    std::vector<LorentzVector> getJets(int idx, FactorizedJetCorrector* jet_corrector, enum JetType type, float deltaR = 0.4, float min_pt = 30.0, float max_eta = 2.4, float min_lep_pt = 20.0, int systFlag = 0, bool sort_by_pt = true);

    // JEC uncertainty applied otf
    std::vector<LorentzVector> getJets(int idx, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type, enum JetType type, float deltaR = 0.4, float min_pt = 30.0, float max_eta = 2.4, float min_lep_pt = 20.0, bool sort_by_pt = true);	 

    // JEC AND JEC uncertainty applied otf
    std::vector<LorentzVector> getJets(int idx, FactorizedJetCorrector* jet_corrector, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type,  enum JetType type, float deltaR = 0.4, float min_pt = 30.0, float max_eta = 2.4, float min_lep_pt = 20.0, bool sort_by_pt = true);	 


    ///////////////////////////////////////////////////////////////////////////////////////////
    // 2012 get jet flags and perform overlap removal with numerator e/mu with pt > x (defaults are 20/20 GeV)
    ///////////////////////////////////////////////////////////////////////////////////////////

    // JEC taken from ntuple
    std::vector<bool> getJetFlags(int idx, enum JetType type, float deltaR = 0.4, float min_pt = 30.0, float max_eta = 2.4, float min_lep_pt = 20.0, int systFlag = 0);

    // JEC applied otf
    std::vector<bool> getJetFlags(int idx, FactorizedJetCorrector* jet_corrector, enum JetType type, float deltaR = 0.4, float min_pt = 30.0, float max_eta = 2.4, float min_lep_pt = 20.0, int systFlag = 0);

    // JEC uncertainty applied otf
    std::vector<bool> getJetFlags(int idx, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type, enum JetType type, float deltaR = 0.4, float min_pt = 30.0, float max_eta = 2.4, float min_lep_pt = 20.0);

    // JEC AND JEC uncertainty applied otf
    std::vector<bool> getJetFlags(int idx, FactorizedJetCorrector* jet_corrector, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type,  enum JetType type, float deltaR = 0.4, float min_pt = 30.0, float max_eta = 2.4, float min_lep_pt = 20.0);


    ///////////////////////////////////////////////////////////////////////////////////////////
    // 2012 get sumpt, skip jets overlapping with numerator e/mu with pt>x (defaults are 20/20 GeV)
    ///////////////////////////////////////////////////////////////////////////////////////////

    // JEC taken from ntuple
    float sumJetPt(int idx, enum JetType type, float deltaR = 0.4, float min_pt = 30.0, float max_eta = 2.4, float min_lep_pt = 20.0, int systFlag = 0);

    // JEC applied otf
    float sumJetPt(int idx, FactorizedJetCorrector* jet_corrector, enum JetType type, float deltaR = 0.4, float min_pt = 30.0, float max_eta = 2.4, float min_lep_pt = 20.0, int systFlag = 0);

    // JEC uncertainty applied otf
    float sumJetPt(int idx, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type, enum JetType type, float deltaR = 0.4, float min_pt = 30.0, float max_eta = 2.4, float min_lep_pt = 20.0);	 

    // JEC AND JEC uncertainty applied otf
    float sumJetPt(int idx, FactorizedJetCorrector* jet_corrector, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type, enum JetType type, float deltaR = 0.4, float min_pt = 30.0, float max_eta = 2.4, float min_lep_pt = 20.0);	 


    ///////////////////////////////////////////////////////////////////////////////////////////
    // 2012 get njets, skip jets overlapping with numerator e/mu with pt>x (defaults are 20/20 GeV)
    ///////////////////////////////////////////////////////////////////////////////////////////

    // JEC taken from ntuple
    int nJets(int idx, enum JetType type, float deltaR = 0.4, float min_pt = 30.0, float max_eta = 2.4, float min_lep_pt = 20.0, int systFlag = 0);

    // JEC applied otf
    int nJets(int idx, FactorizedJetCorrector* jet_corrector, enum JetType type, float deltaR = 0.4, float min_pt = 30.0, float max_eta = 2.4, float min_lep_pt = 20.0, int systFlag = 0);

    // JEC uncertainty applied otf
    int nJets(int idx, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type, enum JetType type, float deltaR = 0.4, float min_pt = 30.0, float max_eta = 2.4, float min_lep_pt = 20.0);	 

    // JEC AND JEC uncertainty applied otf
    int nJets(int idx, FactorizedJetCorrector* jet_corrector, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type, enum JetType type, float deltaR = 0.4, float min_pt = 30.0, float max_eta = 2.4, float min_lep_pt = 20.0);	 


    ///////////////////////////////////////////////////////////////////////////////////////////
    // 2012 get b-tagged jets and perform overlap removal with numerator e/mu with pt > x (defaults are 20/20 GeV)
    ///////////////////////////////////////////////////////////////////////////////////////////

    // JEC taken from ntuple
    std::vector<LorentzVector> getBtaggedJets(int idx, enum JetType type, enum BtagType btag_type, float deltaR = 0.4, float min_pt = 30.0, float max_eta = 2.4, float min_lep_pt = 20.0, int systFlag = 0, bool sort_by_pt = true);

    // JEC applied otf
    std::vector<LorentzVector> getBtaggedJets(int idx, FactorizedJetCorrector* jet_corrector, enum JetType type, enum BtagType btag_type, float deltaR = 0.4, float min_pt = 30.0, float max_eta = 2.4, float min_lep_pt = 20.0, int systFlag = 0, bool sort_by_pt = true);

    // JEC uncertainty applied otf
    std::vector<LorentzVector> getBtaggedJets(int idx, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type, enum JetType type, enum BtagType btag_type, float deltaR = 0.4, float min_pt = 30.0, float max_eta = 2.4, float min_lep_pt = 20.0, bool sort_by_pt = true);

    // JEC AND JEC uncertainty applied otf
    std::vector<LorentzVector> getBtaggedJets(int idx, FactorizedJetCorrector* jet_corrector, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type, enum JetType type, enum BtagType btag_type, float deltaR = 0.4, float min_pt = 30.0, float max_eta = 2.4, float min_lep_pt = 20.0, bool sort_by_pt = true);


    ///////////////////////////////////////////////////////////////////////////////////////////
    // 2012 get b-tagged jet flags and perform overlap removal with numerator e/mu with pt > x (defaults are 20/20 GeV)
    ///////////////////////////////////////////////////////////////////////////////////////////

    // JEC taken from ntuple
    std::vector<bool> getBtaggedJetFlags(int idx, enum JetType type, enum BtagType btag_type, float deltaR = 0.4, float min_pt = 30.0, float max_eta = 2.4, float min_lep_pt = 20.0, int systFlag = 0);

    // JEC applied otf
    std::vector<bool> getBtaggedJetFlags(int idx, FactorizedJetCorrector* jet_corrector, enum JetType type, enum BtagType btag_type, float deltaR = 0.4, float min_pt = 30.0, float max_eta = 2.4, float min_lep_pt = 20.0, int systFlag = 0);

    // JEC uncertainty applied otf
    std::vector<bool> getBtaggedJetFlags(int idx, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type, enum JetType type, enum BtagType btag_type, float deltaR = 0.4, float min_pt = 30.0, float max_eta = 2.4, float min_lep_pt = 20.0);

    // JEC AND JEC uncertainty applied otf
    std::vector<bool> getBtaggedJetFlags(int idx, FactorizedJetCorrector* jet_corrector, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type, enum JetType type, enum BtagType btag_type, float deltaR = 0.4, float min_pt = 30.0, float max_eta = 2.4, float min_lep_pt = 20.0);


    ///////////////////////////////////////////////////////////////////////////////////////////
    // 2012 get sumpt, skip jets overlapping with numerator e/mu with pt>x (defaults are 20/20 GeV)
    ///////////////////////////////////////////////////////////////////////////////////////////

    // JEC taken from ntuple
    int nBtaggedJets(int idx, enum JetType type, enum BtagType btag_type, float deltaR = 0.4, float min_pt = 30.0, float max_eta = 2.4, float min_lep_pt = 20.0, int systFlag = 0);

    // JEC applied otf
    int nBtaggedJets(int idx, FactorizedJetCorrector* jet_corrector, enum JetType type, enum BtagType btag_type, float deltaR = 0.4, float min_pt = 30.0, float max_eta = 2.4, float min_lep_pt = 20.0, int systFlag = 0);

    // JEC uncertainty applied otf
    int nBtaggedJets(int idx, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type, enum JetType type, enum BtagType btag_type, float deltaR = 0.4, float min_pt = 30.0, float max_eta = 2.4, float min_lep_pt = 20.0);	 

    // JEC AND JEC uncertainty applied otf
    int nBtaggedJets(int idx, FactorizedJetCorrector* jet_corrector, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type, enum JetType type, enum BtagType btag_type, float deltaR = 0.4, float min_pt = 30.0, float max_eta = 2.4, float min_lep_pt = 20.0);	 


    ///////////////////////////////////////////////////////////////////////////////////////////
    // 2012 get jet b-tag discriminators, skip jets overlapping with numerator e/mu with pt>x (defaults are 20/20 GeV)
    ///////////////////////////////////////////////////////////////////////////////////////////

    // JEC taken from ntuple
    std::vector<float> getJetBtagDiscriminators(int idx, enum JetType type, enum BtagType btag_type, float deltaR = 0.4, float min_pt = 30.0, float max_eta = 2.4, float min_lep_pt = 20.0, int systFlag = 0, bool sort_by_pt = true);

    // JEC applied otf
    std::vector<float> getJetBtagDiscriminators(int idx, FactorizedJetCorrector* jet_corrector, enum JetType type, enum BtagType btag_type, float deltaR = 0.4, float min_pt = 30.0, float max_eta = 2.4, float min_lep_pt = 20.0, int systFlag = 0, bool sort_by_pt = true);

    // JEC uncertainty applied otf
    std::vector<float> getJetBtagDiscriminators(int idx, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type, enum JetType type, enum BtagType btag_type, float deltaR = 0.4, float min_pt = 30.0, float max_eta = 2.4, float min_lep_pt = 20.0, bool sort_by_pt = true);	 

    // JEC AND JEC uncertainty applied otf
    std::vector<float> getJetBtagDiscriminators(int idx, FactorizedJetCorrector* jet_corrector, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type, enum JetType type, enum BtagType btag_type, float deltaR = 0.4, float min_pt = 30.0, float max_eta = 2.4, float min_lep_pt = 20.0, bool sort_by_pt = true);	 


    ///////////////////////////////////////////////////////////////////////////////////////////
    // 2012 rescale the jet energy resolution (JER) 
    ///////////////////////////////////////////////////////////////////////////////////////////
    void smearJETScaleJetsMetHt(std::vector<LorentzVector>& vjets_p4, float& met, float& met_phi, float& ht, const unsigned int seed);
    void smearJETScaleJets(std::vector<LorentzVector>& vjets_p4, const unsigned int seed);
    void smearJETScaleJetsMetHt
    (
        std::vector<LorentzVector>& vjets_p4, 
        float& met,
        float& met_phi,
        float& ht, 
        int idx,
        enum JetType type,
        const unsigned int seed,
        float deltaR = 0.4,
        float min_pt = 30.0,
        float max_eta = 2.4,
        float lep_minpt = 20.0
    );

    ///////////////////////////////////////////////////////////////////////////////////////////
    // 2012 rescale the MET by scaling up/down the unclustered erngy 
    ///////////////////////////////////////////////////////////////////////////////////////////
    float scaleMET
    (
        const float met,
        const float met_phi,
        int idx,
        enum JetType type,
        float deltaR = 0.4,
        float min_pt = 30.0,
        float max_eta = 2.4,
        float lep_minpt = 20.0,
        const int scale_type = 0,
        const float scale = 0.1
    );

    ///////////////////////////////////////////////////////////////////////////////////////////	 
    // 2012 get vector of good els p4s	 
    ///////////////////////////////////////////////////////////////////////////////////////////	 
    std::vector<LorentzVector> getGoodElectrons(const float ptcut = 20.0f);	 
    std::vector<std::pair<LorentzVector, unsigned int> > getNumeratorElectrons(const float ptcut = 20.0f);	 


    ///////////////////////////////////////////////////////////////////////////////////////////	 
    // 2012 get vector of good mus p4s	 
    ///////////////////////////////////////////////////////////////////////////////////////////	 
    std::vector<LorentzVector> getGoodMuons(const float ptcut = 20.0f);
    std::vector<std::pair<LorentzVector, unsigned int> > getNumeratorMuons(const float ptcut = 20.0f);	 
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    // 2012 rescale the jet energy resolution (JER) 
    // function that returns the sigma(pT)*pT of the MC jets. those numbers are old,
    // but at least they are in an understandable format.
    ///////////////////////////////////////////////////////////////////////////////////////////
    float getErrPt(const float pt, const float eta);


    ///////////////////////////////////////////////////////////////////////////////////////////
    // function to get the jer scale factors. values taken from the twiki: 
    // https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
    ///////////////////////////////////////////////////////////////////////////////////////////
    float getJERScale(const float jet_eta);

} // namespace os2012

#endif //OSSELECTIONS_H
