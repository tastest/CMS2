#ifndef twinmaker_h
#define twinmaker_h

#include "TFile.h"
#include "TTree.h"

class TChain;

class twinmaker
{
    public:
        twinmaker() {};
        ~twinmaker() {
            delete twinFile_;
            delete twinTree_;
        };
        void MakeTwinNtuple (const char *);
        void InitTwinNtuple ();
        void FillTwinNtuple ();
        void CloseTwinNtuple ();
        void ScanChain (const char *, const char *, int nEvents = -1);

    private:
        //
        // TWIN NTUPLE VARIABLES
        //
        TFile *twinFile_;
        TTree *twinTree_;

        // event stuff
        Int_t   run_;
        Int_t   ls_;
        Int_t   evt_;
        Int_t   njets_;
	Int_t   hyp_type_;
        Float_t pfmet_;
	Float_t tcmet_;
        Float_t jet1pt_;
        Float_t dphipfmetjet_;
        Float_t dphitcmetjet_;
	Float_t dilpt_;
	Float_t deltaphi_;	

        // lepton stuff
        Int_t   eormu1_;
	Int_t   eormu2_;
        Int_t   type1_;
        Int_t   type2_;
        Float_t pt1_;
        Float_t pt2_;
	Float_t eta1_;
	Float_t eta2_;
	Float_t phi1_;
	Float_t phi2_;
        Float_t iso1_;
        Float_t iso2_;
        Float_t d0corr1_;
        Float_t d0corr2_;
        Float_t d0vtx1_;
        Float_t d0vtx2_;
        Float_t dphipfmet1_;
        Float_t dphipfmet2_;
        Float_t dphitcmet1_;
        Float_t dphitcmet2_;
        Float_t drjet1_;
        Float_t drjet2_;
	Float_t mass_;

        // muon stuff
        Bool_t  mu1_muonid_;
        Bool_t  mu2_muonid_;
        Int_t   mu1_goodmask_;
        Int_t   mu2_goodmask_;
        Float_t mu1_gfitchi2_;
        Float_t mu2_gfitchi2_;

        // electron stuff
        Bool_t  e1_cand01_;
        Bool_t  e2_cand01_;
        Int_t   e1_nmHits_;
        Int_t   e2_nmHits_;
        Float_t e1_eopin_;
        Float_t e2_eopin_;
        Float_t e1_hoe_;
        Float_t e2_hoe_;
        Float_t e1_dphiin_;
        Float_t e2_dphiin_;
        Float_t e1_detain_;
        Float_t e2_detain_;
        Float_t e1_eMe55_;
        Float_t e2_eMe55_;
	Float_t e1_drmu_;
	Float_t e2_drmu_;
};

#endif
