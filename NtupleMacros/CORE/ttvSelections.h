#ifndef TTVSELECTIONS_H
#define TTVSELECTIONS_H
#include <vector>
#include "jetSelections.h"

class FactorizedJetCorrector;

namespace ttv
{
    struct LeptonType
    {
        enum value_type
        {
            LOOSE, // eg. lepton from Z
            TIGHT, // eg. lepton from W
			LOOSEMVA,
			TIGHTMVA,
            static_size
        };
    };

    struct LeptonInfo
    {
        int id;
        int idx;        
    };

    struct ZleptonInfo
    {
        LeptonInfo lep1;
        LeptonInfo lep2;
    };

    struct TTZleptonInfo
    {
        ZleptonInfo z;
        LeptonInfo w;
    };
    
    bool passesTrigger (int hyp_type);

    bool isGoodLepton        (int id, int idx, LeptonType::value_type lep_type);
    bool isIsolatedLepton    (int id, int idx, LeptonType::value_type lep_type);
    bool isNumeratorLepton   (int id, int idx, LeptonType::value_type lep_type);
    bool isDenominatorLepton (int id, int idx, LeptonType::value_type lep_type);
  bool overlapMuon(int idx, LeptonType::value_type lep_type, float pt=10., float eta=2.4, float deltaR=0.4);

    std::vector<LorentzVector> getJets(std::vector<LorentzVector>& leps, enum JetType type, float deltaR = 0.4, float min_pt = 30., float max_eta = 2.4, float rescale = 1.0, int systFlag = 0);
    std::vector<LorentzVector> getJets(std::vector<LorentzVector>& leps, FactorizedJetCorrector* jet_corrector, enum JetType type, float deltaR = 0.4, float min_pt = 30., float max_eta = 2.4, float rescale = 1.0, int systFlag = 0);

    std::vector<bool> getJetFlags(std::vector<LorentzVector>& leps, enum JetType type, float deltaR = 0.4, float min_pt = 30., float max_eta = 2.4, float rescale = 1.0, int systFlag = 0);
    std::vector<bool> getJetFlags(std::vector<LorentzVector>& leps, FactorizedJetCorrector* jet_corrector, enum JetType type, float deltaR = 0.4, float min_pt = 30., float max_eta = 2.4, float rescale = 1.0, int systFlag = 0);
    
    float sumJetPt(std::vector<LorentzVector>& leps, enum JetType type, float deltaR = 0.4, float min_pt = 30., float max_eta = 2.4, float rescale = 1.0, int systFlag = 0);
    float sumJetPt(std::vector<LorentzVector>& leps, FactorizedJetCorrector* jet_corrector, enum JetType type, float deltaR = 0.4, float min_pt = 30., float max_eta = 2.4, float rescale = 1.0, int systFlag = 0);

    int nJets(std::vector<LorentzVector>& leps, enum JetType type, float deltaR = 0.4, float min_pt = 30., float max_eta = 2.4, float rescale = 1.0, int systFlag = 0);
    int nJets(std::vector<LorentzVector>& leps, FactorizedJetCorrector* jet_corrector, enum JetType type, float deltaR = 0.4, float min_pt = 30., float max_eta = 2.4, float rescale = 1.0, int systFlag = 0);

    std::vector<LorentzVector> getBtaggedJets(std::vector<LorentzVector>& leps, enum JetType type, enum BtagType btag_type, float deltaR = 0.4, float min_pt = 30., float max_eta = 2.4, float rescale = 1.0, int systFlag = 0);
    std::vector<LorentzVector> getBtaggedJets(std::vector<LorentzVector>& leps, FactorizedJetCorrector* jet_corrector, enum JetType type, enum BtagType btag_type, float deltaR = 0.4, float min_pt = 30., float max_eta = 2.4, float rescale = 1.0, int systFlag = 0);

    std::vector<bool> getBtaggedJetFlags(std::vector<LorentzVector>& leps, enum JetType type, enum BtagType btag_type, float deltaR = 0.4, float min_pt = 30., float max_eta = 2.4, float rescale = 1.0, int systFlag = 0);
    std::vector<bool> getBtaggedJetFlags(std::vector<LorentzVector>& leps, FactorizedJetCorrector* jet_corrector, enum JetType type, enum BtagType btag_type, float deltaR = 0.4, float min_pt = 30., float max_eta = 2.4, float rescale = 1.0, int systFlag = 0);

    int nBtaggedJets(std::vector<LorentzVector>& leps, enum JetType type, enum BtagType btag_type, float deltaR = 0.4, float min_pt = 30., float max_eta = 2.4, float rescale = 1.0, int systFlag = 0);
    int nBtaggedJets(std::vector<LorentzVector>& leps, FactorizedJetCorrector* jet_corrector, enum JetType type, enum BtagType btag_type, float deltaR = 0.4, float min_pt = 30., float max_eta = 2.4, float rescale = 1.0, int systFlag = 0);


}

#endif
