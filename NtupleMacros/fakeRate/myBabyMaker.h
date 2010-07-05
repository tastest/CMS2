#ifndef myBabyMaker_h
#define myBabyMaker_h

#include "TFile.h"
#include "TTree.h"

//class TChain;

class myBabyMaker
{
    public:
        myBabyMaker() {};
        ~myBabyMaker() {
            delete babyFile_;
            delete babyTree_;
        };
        void MakeBabyNtuple (const char *);
        void InitBabyNtuple ();
        void FillBabyNtuple ();
        void CloseBabyNtuple ();
        void ScanChain (TChain *, const char *, int);

    private:
        //
        // BABY NTUPLE VARIABLES
        //
        TFile *babyFile_;
        TTree *babyTree_;

        // event identification
        Int_t   run_;
        Int_t   ls_;
        Int_t   evt_;

	// Lepton pt and eta and phi
	Float_t pt_;
	Float_t eta_;
	Float_t phi_;
	Float_t scet_;
	Int_t id_; // \pm 11 or \pm 13

	// tcmet
	Float_t tcmet_;

	// did it pass the jet trigger and is this lepton unbiased
	// 0=fail 1="pass but biased" 2="pass and unbiased"  -1="pass but cant find jet trg obj"
	Int_t  hlt15u_; // HLT_Jet15U
	Int_t  hlt30u_; // HLT_Jet30U
	Int_t  hlt50u_; // HLT_Jet50U
	Int_t  l16u_;   // HLT_L1Jet6U
	Int_t  l110u_;   // HLT_L1Jet10U

	// What type of electron selection did it pass
	Bool_t v1_;
	Bool_t v2_;
	Bool_t v3_;
	Bool_t num_;

	// What electron trigger did it pass and is the electron matched to a egamma trg object
	// 0=fail 1="pass but no match" 2="pass and matched" -1="pass but egamm trg obj missing"
	Int_t ph10_;  // HLT_Photon10_L1R or HLT_Photon10_Cleaned_L1R
	Int_t ph15_;  // HLT_Photon15_L1R or HLT_Photon15_Cleaned_L1R
	Int_t el10_;  // HLT_Ele10_LW_L1R
	Int_t eg5_;   // HLT_L1SingleEG5 
	Int_t eg8_;   // HLT_L1SingleEG


	// What muon trigger did it pass
	// 0=fail 1="pass but no match" 2="pass and matched" -1="pass but muon trg obj missing"
	Int_t mu9_; // HLT_Mu9
	Int_t mu5_; // HLT_Mu5
	Int_t mu3_; // HLT_Mu3

	//  Minimm dR to the closest eg object
	Float_t drph10_;
	Float_t drph15_;
	Float_t drel10_;
	Float_t dreg5_;
	Float_t dreg8_;

	//  Minimm dR to the closest HLT mu object
	Float_t drmu9_;
	Float_t drmu5_;
	Float_t drmu3_;
	
	// Btag information
	Int_t nbjet_; // number of btagged jet pt>15
	Float_t dRbNear_; // dR between lepton and closest such jet
	Float_t dRbFar_; // dR between lepton and farthest such jet


};

#endif
