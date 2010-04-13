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
        Int_t   njets_;
        Float_t jet1pt_;
        Float_t dphimetjet_;
        Int_t   eormu_;

        // lepton stuff
        Int_t   type_;
        Float_t pt_;
        Float_t iso_;
        Float_t d0corr_;
        Float_t dphimet_;
        Float_t drjet_;

        // muon stuff
        Bool_t  mu_muonid_;
        Int_t   mu_goodmask_;
        Float_t mu_gfitchi2_;

        // electron stuff
        Bool_t  e_cand01_;
        Float_t e_eopin_;
        Float_t e_hoe_;
        Float_t e_dphiin_;
        Float_t e_detain_;
        Float_t e_eMe55_;
        Int_t   e_nmHits_;
};

#endif
