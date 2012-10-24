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
        LeptonInfo (int id_, int idx_) : id(id_), idx(idx_) {};
        LeptonInfo () : id(0), idx(-1) {};
        int id;
        int idx;
        LorentzVector p4();
    };

    struct DileptonInfo
    {
        DileptonInfo (const LeptonInfo& lep1_, const LeptonInfo& lep2_) : lep1(lep1_), lep2(lep2_) {};
        DileptonInfo (int id1_, int idx1_, int id2_, int idx2_);
        DileptonInfo (int hyp_index_);
        DileptonInfo () : lep1(), lep2() {};
        LeptonInfo lep1;
        LeptonInfo lep2;
        int hyp_index;
        float sumpt ();
        LorentzVector p4 ();
        float mass ();        
    };

    struct TrileptonInfo
    {
        TrileptonInfo (const DileptonInfo& z_, const LeptonInfo& w_) : z(z_), w(w_) {};
        TrileptonInfo (int id1_, int idx1_, int id2_, int idx2_, int id3_, int idx3_);
        TrileptonInfo (int hyp_index_, int id3_, int idx3_);
        TrileptonInfo () : z(), w() {};
        DileptonInfo z;
        LeptonInfo w;
        float sumpt ();
    };
    
    bool passesTrigger (int hyp_type);
    bool passesTrigger (TrileptonInfo& trilep_info);

    bool isGoodLepton        (int id, int idx, LeptonType::value_type lep_type);
    bool isIsolatedLepton    (int id, int idx, LeptonType::value_type lep_type);
    bool isNumeratorLepton   (int id, int idx, LeptonType::value_type lep_type);
    bool isDenominatorLepton (int id, int idx, LeptonType::value_type lep_type);
    bool overlapMuon         (int idx, LeptonType::value_type lep_type, float pt = 20., float eta = 2.4, float deltaR = 0.1);

    std::vector<LorentzVector> getGoodLeptons(LeptonType::value_type lep_type, float pt = 20., float mu_eta = 2.4, float ele_eta = 2.5);
    float getTrigMVAThreshold(int idx, ttv::LeptonType::value_type lep_type);
    float getNonTrigMVAThreshold(int idx, ttv::LeptonType::value_type lep_type);

    std::vector<LorentzVector> getJets(std::vector<LorentzVector>& leps, enum JetType type, float deltaR = 0.5, float min_pt = 15., float max_eta = 2.4, float rescale = 1.0, int systFlag = 0);
    std::vector<LorentzVector> getJets(std::vector<LorentzVector>& leps, FactorizedJetCorrector* jet_corrector, enum JetType type, float deltaR = 0.5, float min_pt = 15., float max_eta = 2.4, float rescale = 1.0, int systFlag = 0);

    std::vector<bool> getJetFlags(std::vector<LorentzVector>& leps, enum JetType type, float deltaR = 0.5, float min_pt = 15., float max_eta = 2.4, float rescale = 1.0, int systFlag = 0);
    std::vector<bool> getJetFlags(std::vector<LorentzVector>& leps, FactorizedJetCorrector* jet_corrector, enum JetType type, float deltaR = 0.5, float min_pt = 15., float max_eta = 2.4, float rescale = 1.0, int systFlag = 0);
    
    float sumJetPt(std::vector<LorentzVector>& leps, enum JetType type, float deltaR = 0.5, float min_pt = 15., float max_eta = 2.4, float rescale = 1.0, int systFlag = 0);
    float sumJetPt(std::vector<LorentzVector>& leps, FactorizedJetCorrector* jet_corrector, enum JetType type, float deltaR = 0.5, float min_pt = 15., float max_eta = 2.4, float rescale = 1.0, int systFlag = 0);

    int nJets(std::vector<LorentzVector>& leps, enum JetType type, float deltaR = 0.5, float min_pt = 15., float max_eta = 2.4, float rescale = 1.0, int systFlag = 0);
    int nJets(std::vector<LorentzVector>& leps, FactorizedJetCorrector* jet_corrector, enum JetType type, float deltaR = 0.5, float min_pt = 15., float max_eta = 2.4, float rescale = 1.0, int systFlag = 0);

    std::vector<LorentzVector> getBtaggedJets(std::vector<LorentzVector>& leps, enum JetType type, enum BtagType btag_type, float deltaR = 0.5, float min_pt = 15., float max_eta = 2.4, float rescale = 1.0, int systFlag = 0);
    std::vector<LorentzVector> getBtaggedJets(std::vector<LorentzVector>& leps, FactorizedJetCorrector* jet_corrector, enum JetType type, enum BtagType btag_type, float deltaR = 0.5, float min_pt = 15., float max_eta = 2.4, float rescale = 1.0, int systFlag = 0);

    std::vector<bool> getBtaggedJetFlags(std::vector<LorentzVector>& leps, enum JetType type, enum BtagType btag_type, float deltaR = 0.5, float min_pt = 15., float max_eta = 2.4, float rescale = 1.0, int systFlag = 0);
    std::vector<bool> getBtaggedJetFlags(std::vector<LorentzVector>& leps, FactorizedJetCorrector* jet_corrector, enum JetType type, enum BtagType btag_type, float deltaR = 0.5, float min_pt = 15., float max_eta = 2.4, float rescale = 1.0, int systFlag = 0);

    int nBtaggedJets(std::vector<LorentzVector>& leps, enum JetType type, enum BtagType btag_type, float deltaR = 0.5, float min_pt = 15., float max_eta = 2.4, float rescale = 1.0, int systFlag = 0);
    int nBtaggedJets(std::vector<LorentzVector>& leps, FactorizedJetCorrector* jet_corrector, enum JetType type, enum BtagType btag_type, float deltaR = 0.5, float min_pt = 15., float max_eta = 2.4, float rescale = 1.0, int systFlag = 0);


}

#endif
