#ifndef babymaker_h
#define babymaker_h

#include "TFile.h"
#include "TTree.h"

class TChain;

class babymaker
{
    public:
        babymaker() {};
        ~babymaker() {
            delete babyFile_;
            delete babyTree_;
        };
        void MakeBabyNtuple (const char *);
        void InitBabyNtuple ();
        void FillBabyNtuple ();
        void CloseBabyNtuple ();
        void ScanChain (const char *, const char *, int nEvents = -1);

    private:
        //
        // BABY NTUPLE VARIABLES
        //
        TFile *babyFile_;
        TTree *babyTree_;

        // event stuff
        Int_t   run_;
        Int_t   ls_;
        Int_t   evt_;
        Float_t pfmet_;
        Float_t tcmet_;
        Int_t   njets_;
        Float_t jet1pt_;
        Float_t dphipfmetjet_;
        Float_t dphitcmetjet_;
        Int_t   ntrks_;

        // lepton stuff
        Int_t   eormu_;
        Int_t   type_;
        Float_t pt_;
        Float_t eta_;
        Float_t phi_;
        Float_t iso_;
        Float_t d0corr_;
        Float_t d0vtx_;
        Float_t dphipfmet_;
        Float_t dphitcmet_;
        Float_t drjet_;
        Float_t mt_;

        // muon stuff
        Bool_t  mu_muonid_;
        Int_t   mu_goodmask_;
        Float_t mu_gfitchi2_;
        Bool_t  mu_cosmic_;

        // electron stuff
        Bool_t  e_cand01_;
        Float_t e_eopin_;
        Float_t e_hoe_;
        Float_t e_dphiin_;
        Float_t e_detain_;
        Float_t e_eMe55_;
        Int_t   e_nmHits_;
        Float_t e_dcot_;
        Float_t e_dist_;
        Float_t e_drmu_;
};

#endif
