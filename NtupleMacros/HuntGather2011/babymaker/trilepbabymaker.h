#ifndef trilepbabymaker_h
#define trilepbabymaker_h

#include "TFile.h"
#include "TTree.h"

class TChain;

class trilepbabymaker
{
    public:
        trilepbabymaker() {};
        ~trilepbabymaker() {
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
        char    dataset_[200];
        Int_t   run_;
        Int_t   ls_;
        Int_t   evt_;
        Int_t   isdata_;
		Int_t   nvtx_;
		Float_t scale1fb_;
		Float_t pthat_;
        Int_t   hyp_type_;
        Float_t pfmet_;
        Float_t tcmet_;
		Float_t calotcmet_;
        Int_t   ntrks_;
        Int_t   njets_;
        Float_t jet1pt_;
        Float_t jet2pt_;
        Float_t jet3pt_;
        Float_t sumjetpt_;
        Float_t jet1eta_;
        Float_t jet2eta_;
        Float_t jet3eta_;
        Float_t jet1phi_;
        Float_t jet2phi_;
        Float_t jet3phi_;
        Float_t jetmass_;
        Bool_t  jet1isBtag_;
        Bool_t  jet2isBtag_;
        Bool_t  jet3isBtag_;
        Float_t dphipfmetjet_;
        Float_t dphitcmetjet_;
        Int_t   neffbtags_;
        Int_t   npurbtags_;
        Int_t   ntceffbtags_;
        Int_t   ntcpurbtags_;
		Float_t pfmeff_;
		Float_t tcmeff_;

        // lepton stuff
		Int_t   ngoodlep_;
		Int_t   ngoodmus_;
		Int_t   ngoodels_;
        Int_t   ngenels_;
        Int_t   ngenmus_;
        Int_t   ngentaus_;
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

        Int_t   eormu3_;
        Int_t   type3_;
        Float_t pt3_;
        Float_t eta3_;
        Float_t phi3_;
        Float_t iso3_;
        Float_t d0corr3_;
        Float_t d0vtx3_;
        Float_t dphipfmet3_;
        Float_t dphitcmet3_;
        Float_t drjet3_;

        // muon stuff
        Bool_t  mu1_muonidfull_;
        Bool_t  mu1_muonid_;
        Bool_t  mu1_muonidfullV1_;
        Bool_t  mu1_muonidV1_;
        Int_t   mu1_goodmask_;
        Float_t mu1_gfitchi2_;
        Bool_t  mu1_cosmic_;
		Int_t   mu1_siHits_;
		Int_t   mu1_saHits_;
		Float_t mu1_emVetoDep_;
		Float_t mu1_hadVetoDep_;
		Bool_t  mu1_isPFmuon_;

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

        Bool_t  mu3_cosmic_;
        Bool_t  mu3_muonidfull_;
        Bool_t  mu3_muonid_;
        Bool_t  mu3_muonidfullV1_;
        Bool_t  mu3_muonidV1_;
        Int_t   mu3_goodmask_;
        Float_t mu3_gfitchi2_;
		Int_t   mu3_siHits_;
		Int_t   mu3_saHits_;
		Float_t mu3_emVetoDep_;
		Float_t mu3_hadVetoDep_;
		Bool_t  mu3_isPFmuon_;

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
	
        Bool_t  e3_cand01full_;
        Bool_t  e3_cand01_;
        Bool_t  e3_vbtf90full_;
        Bool_t  e3_vbtf90_;
        Bool_t  e3_vbtf85_;
        Bool_t  e3_vbtf80_;
        Bool_t  e3_vbtf70_;
        Float_t e3_scet_;
        Float_t e3_eopin_;
        Float_t e3_hoe_;
        Float_t e3_dphiin_;
        Float_t e3_detain_;
        Float_t e3_e25Me55_;
        Float_t e3_sigieie_; // sigmaietaieta
        Float_t e3_eMe55_;
        Int_t   e3_nmHits_;
        Float_t e3_dcot_;
        Float_t e3_dist_;
        Float_t e3_drmu_;
        Bool_t  e3_isspike_;
        Int_t   e3_ctfCharge_;
        Int_t   e3_scCharge_;
        Int_t   e3_gsfCharge_;
}
;

#endif
