#ifndef dilepbabymaker_h
#define dilepbabymaker_h

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

#include "Math/LorentzVector.h"
#include "Math/PxPyPzE4D.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

class dilepbabymaker
{
    public:
        dilepbabymaker() {};
        ~dilepbabymaker() {
            delete babyFile_;
            delete babyTree_;
        };
        void MakeBabyNtuple (const char *);
        void InitBabyNtuple ();

        bool PassTriggerGroup(const std::vector<std::string> &triggers, const LorentzVector &obj);
        bool PassTriggerGroup(const std::vector<std::string> &triggers);

        void SetEventLevelInfo();
        void NewRun();

        void FillBabyNtuple ();
        void CloseBabyNtuple ();
        void ScanChain (const char *, const char *, int nEvents = -1);

    private:
        //
        // BABY NTUPLE VARIABLES
        //
        TFile *babyFile_;
        TTree *babyTree_;

        std::vector<std::string> triggers_e_;
        std::vector<std::string> triggers_m_;
        std::vector<std::string> triggers_ee_;
        std::vector<std::string> triggers_mm_;
        std::vector<std::string> triggers_em_;

        // event stuff
        char    dataset_[200];
        Int_t   run_;
        Int_t   ls_;
        Int_t   evt_;
        Int_t   evt_clean082010_;
        Int_t   evt_clean102010_;
        Int_t   isdata_;
        Int_t   nvtx_;
        Float_t scale1fb_;
        Float_t pthat_;
        Int_t   hyp_type_;
        Float_t pfmet_;
        Float_t tcmet_;
        Float_t calotcmet_;
        Float_t proj_pfmet_;
        Float_t proj_tcmet_;
        Int_t   ntrks_;
        Int_t   njets_;
        Int_t   njets25_;
        Int_t   njetsSS_;
        Float_t jet1pt_;
        Float_t jet2pt_;
        Float_t jet3pt_;
        Float_t sumjetpt_;
        Float_t sumjetpt25_;
        Float_t sumjetptSS_;
        Float_t jet1eta_;
        Float_t jet2eta_;
        Float_t jet3eta_;
        Float_t jet1phi_;
        Float_t jet2phi_;
        Float_t jet3phi_;
        Float_t jetmass_;
        Float_t pfmth_;
        Float_t tcmth_;
        Int_t  jet1isBtag_;
        Int_t  jet2isBtag_;
        Int_t  jet3isBtag_;
        Float_t dphipfmetjet_;
        Float_t dphitcmetjet_;
        Float_t deltaphi_;	
        Int_t ntchelbtags_;
        Int_t nssvhembtags_;
        Int_t nssvhetbtags_;
        Int_t nssvhptbtags_;
        Float_t pfmeff_;
        Float_t tcmeff_;
        Float_t intLumiPerLS_;

        // lepton stuff
        Int_t   ngoodlep_;
        Int_t   ngoodmus_;
        Int_t   ngoodels_;
        Int_t   ngenels_;
        Int_t   ngenmus_;
        Int_t   ngentaus_;
        Float_t dilpt_;
        Float_t mass_;
        Int_t   eormu1_;
        Int_t   type1_;
        Float_t pt1_;
        Float_t eta1_;
        Float_t phi1_;
        Float_t iso1_;
        Float_t d0corr1_;
        Float_t d0vtx1_;
        Float_t dphipfmet1_;
        Float_t dphitcmet1_;
        Float_t drjet1_;
        Int_t   mcid1_;
        Int_t   mcmotherid1_;
        Int_t   eormu2_;
        Int_t   type2_;
        Float_t pt2_;
        Float_t eta2_;
        Float_t phi2_;
        Float_t iso2_;
        Float_t d0corr2_;
        Float_t d0vtx2_;
        Float_t dphipfmet2_;
        Float_t dphitcmet2_;
        Float_t drjet2_;
        Int_t   mcid2_;
        Int_t   mcmotherid2_;
        Float_t mt2_;
        Float_t mt2j_;
        Bool_t  extraZveto_;
        Float_t trkIso1_;
        Float_t ecalIso1_;
        Float_t hcalIso1_;
        Float_t trkIso2_;
        Float_t ecalIso2_;
        Float_t hcalIso2_;

        // muon stuff
        Bool_t  mu1_muonidfull_;
        Bool_t  mu1_muonid_;
        Bool_t  mu1_muonidfullV1_;
        Bool_t  mu1_muonidV1_;
        Int_t   mu1_goodmask_;
        Float_t mu1_gfitchi2_;
        Int_t   mu1_siHits_;
        Int_t   mu1_saHits_;
        Float_t mu1_emVetoDep_;
        Float_t mu1_hadVetoDep_;
        Bool_t  mu1_isPFmuon_;

        Bool_t  mu1_cosmic_;
        Bool_t  mu2_cosmic_;
        Bool_t  mu2_muonidfull_;
        Bool_t  mu2_muonid_;
        Bool_t  mu2_muonidfullV1_;
        Bool_t  mu2_muonidV1_;
        Int_t   mu2_goodmask_;
        Float_t mu2_gfitchi2_;
        Int_t   mu2_siHits_;
        Int_t   mu2_saHits_;
        Float_t mu2_emVetoDep_;
        Float_t mu2_hadVetoDep_;
        Bool_t  mu2_isPFmuon_;

        // electron stuff
        Bool_t  e1_cand01full_;
        Bool_t  e1_cand01_;
        Bool_t  e1_vbtf90full_;
        Bool_t  e1_vbtf90_;
        Bool_t  e1_vbtf85_;
        Bool_t  e1_vbtf80_;
        Bool_t  e1_vbtf70_;
        Float_t e1_scet_;
        Float_t e1_eopin_;
        Float_t e1_hoe_;
        Float_t e1_dphiin_;
        Float_t e1_detain_;
        Float_t e1_e25Me55_;
        Float_t e1_sigieie_; // sigmaietaieta
        Float_t e1_eMe55_;
        Int_t   e1_nmHits_;
        Float_t e1_dcot_;
        Float_t e1_dist_;
        Float_t e1_drmu_;
        Bool_t  e1_isspike_;
        Int_t   e1_ctfCharge_;
        Int_t   e1_scCharge_;
        Int_t   e1_gsfCharge_;
        Float_t e1_fbrem_;

        Bool_t  e2_cand01full_;
        Bool_t  e2_cand01_;
        Bool_t  e2_vbtf90full_;
        Bool_t  e2_vbtf90_;
        Bool_t  e2_vbtf85_;
        Bool_t  e2_vbtf80_;
        Bool_t  e2_vbtf70_;
        Float_t e2_scet_;
        Float_t e2_eopin_;
        Float_t e2_hoe_;
        Float_t e2_dphiin_;
        Float_t e2_detain_;
        Float_t e2_e25Me55_;
        Float_t e2_sigieie_; // sigmaietaieta
        Float_t e2_eMe55_;
        Int_t   e2_nmHits_;
        Float_t e2_dcot_;
        Float_t e2_dist_;
        Float_t e2_drmu_;
        Bool_t  e2_isspike_;
        Int_t   e2_ctfCharge_;
        Int_t   e2_scCharge_;
        Int_t   e2_gsfCharge_;
        Float_t e2_fbrem_;

        // triggers
        Int_t trg_single_mu_;
        Int_t trg_single_e_;
        Int_t trg_double_mu_;
        Int_t trg_double_e_;
        Int_t trg_cross_emu_;

};

#endif
