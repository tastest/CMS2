#ifndef emubabymaker_h
#define emubabymaker_h

#include "TFile.h"
#include "TTree.h"

class TChain;

class emubabymaker
{
public:
	 emubabymaker() {};
	 ~emubabymaker() {
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
	 Float_t pfmet_;
	 Float_t tcmet_;
	 Int_t   ntrks_;
	 Int_t   njets_;
	 Int_t   njetsClean_;
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
	 Bool_t  jet1passesID_;
	 Bool_t  jet2passesID_;
	 Bool_t  jet3passesID_;
	 Bool_t  jet1isBtag_;
	 Bool_t  jet2isBtag_;
	 Bool_t  jet3isBtag_;
	 Float_t dphipfmetjet_;
	 Float_t dphitcmetjet_;
	 Int_t   neffbtags_;
	 Int_t   npurbtags_;

	 // lepton stuff
	 Int_t   ngoodlep_;
	 Int_t   ngoodmus_;
     Int_t   ngoodels_;
	 Int_t   eormu_;
	 Int_t   type_;
     Int_t   ngenels_;
     Int_t   ngenmus_;
     Int_t   ngentaus_;
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
	 Float_t tcmt_;
	 Float_t sf_mass_; // hyp_mass if found in ee or mumu

	 // muon stuff
	 Bool_t  mu_muonidfull_;
	 Bool_t  mu_muonid_;
	 Bool_t  mu_muonidfullV1_;
	 Bool_t  mu_muonidV1_;
	 Int_t   mu_goodmask_;
	 Float_t mu_gfitchi2_;
	 Bool_t  mu_cosmic_;
	 Int_t   mu_siHits_;
	 Int_t   mu_saHits_;
	 Float_t mu_emVetoDep_;
	 Float_t mu_hadVetoDep_;

	 // electron stuff
	 Bool_t  e_cand01full_;
	 Bool_t  e_cand01_;
	 Bool_t  e_vbtf90fullAlign_;
	 Bool_t  e_vbtf90full_;
	 Bool_t  e_vbtf90_;
	 Bool_t  e_vbtf85_;
	 Bool_t  e_vbtf80_;
	 Bool_t  e_vbtf70_;
	 Float_t e_scet_;
	 Float_t e_eopin_;
	 Float_t e_hoe_;
	 Float_t e_dphiin_;
	 Float_t e_detain_;
	 Float_t e_e25Me55_;
	 Float_t e_sigieie_; // sigmaietaieta
	 Float_t e_eMe55_;
	 Int_t   e_nmHits_;
	 Float_t e_dcot_;
	 Float_t e_dist_;
	 Float_t e_drmu_;
	 Bool_t  e_isspike_;
	 Int_t   e_ctfCharge_;
	 Int_t   e_gsfCharge_;
	 Int_t   e_scCharge_;
};

#endif
