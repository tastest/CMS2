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
        void ScanChain (TChain *, const char *, bool, int);

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
	Float_t tcmetphi_;

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
	Bool_t num_;   // {electronSelection_ttbarV1 && (!isSpikeElectron()}  NominalTTbarV2   

  Bool_t v1SS_;
  Bool_t v2SS_;
  Bool_t v3SS_;
	Bool_t numSS_; // electronSelection_ss Nominal

  Bool_t v1SSAug9_;
  Bool_t v2SSAug9_;
  Bool_t v3SSAug9_;
	Bool_t numSSAug9_; // electronSelection_ss Nominal

	Bool_t numv1_; // electronSelection_ttbarV1 NominalTTbar

	Bool_t numAug9_; 
  Bool_t v1Aug9_;
  Bool_t v2Aug9_;
  Bool_t v3Aug9_;

  Bool_t fo_04_;
  Bool_t fo_10_;

	// What electron trigger did it pass and is the electron matched to a egamma trg object
	// 0=fail 1="pass but no match" 2="pass and matched" -1="pass but egamm trg obj missing"
	Int_t ph10_;    // HLT_Photon10_L1R or HLT_Photon10_Cleaned_L1R
	Int_t ph15_;    // HLT_Photon15_L1R or HLT_Photon15_Cleaned_L1R
	Int_t ph20_;    // HLT_Photon20_Cleaned_L1R

	Int_t el10_lw_;    // HLT_Ele10_LW_L1R
	Int_t el10_sw_;    // HLT_Ele10_SW_L1R

	Int_t el10_lw_id_; // HLT_Ele10_LW_EleId_L1R
	Int_t el10_sw_id_; // HLT_Ele10_SW_EleId_L1R

	Int_t el15_lw_;    // HLT_Ele15_LW_L1R
	Int_t el15_sw_;    // HLT_Ele15_SW_L1R

	Int_t el15_lw_id_; // HLT_Ele15_LW_EleId_L1R
	Int_t el15_sw_id_; // HLT_Ele15_SW_EleId_L1R

	Int_t el15_sw_cid_; // HLT_Ele15_SW_CaloEleId_L1R

	// What muon trigger did it pass
	// 0=fail 1="pass but no match" 2="pass and matched" -1="pass but muon trg obj missing"
	Int_t mu9_; // HLT_Mu9
	Int_t mu7_; // HLT_Mu7
	Int_t mu5_; // HLT_Mu5

	//  Minimm dR to the closest eg object
  Float_t drph10_;    // HLT_Photon10_L1R or HLT_Photon10_Cleaned_L1R
  Float_t drph15_;    // HLT_Photon15_L1R or HLT_Photon15_Cleaned_L1R
  Float_t drph20_;    // HLT_Photon20_Cleaned_L1R

  Float_t drel10_lw_;    // HLT_Ele10_LW_L1R
  Float_t drel10_sw_;    // HLT_Ele10_SW_L1R

  Float_t drel10_lw_id_; // HLT_Ele10_LW_EleId_L1R
  Float_t drel10_sw_id_; // HLT_Ele10_SW_EleId_L1R

  Float_t drel15_lw_;    // HLT_Ele15_LW_L1R
  Float_t drel15_sw_;    // HLT_Ele15_SW_L1R

  Float_t drel15_lw_id_; // HLT_Ele15_LW_EleId_L1R
  Float_t drel15_sw_id_; // HLT_Ele15_SW_EleId_L1R

  Float_t drel15_sw_cid_; // HLT_Ele15_SW_CaloEleId_L1R

	//  Minimm dR to the closest HLT mu object
	Float_t drmu9_;
	Float_t drmu7_;
	Float_t drmu5_;
	
	// Btag information
	Int_t nbjet_; // number of btagged jet pt>15
	Float_t dRbNear_; // dR between lepton and closest such jet
	Float_t dRbFar_; // dR between lepton and farthest such jet


	// Information to do offline jet trigger selection
	Float_t ptj1_;        // highest pt jet well separated from the lepton
	Float_t ptj1_b2b_;    // highest pt jet away frmo lepton by dR >= 1.0 and dPhi > 2.5
	Float_t dphij1_b2b_;  // dphi between lepton and jet for jets away from lepton by dR >= 1.0
	Int_t   nj1_;         // number of jets above 10 GeV and away from lepton by dR >= 1.0

	// transverse W mass
	Float_t mt_;

};

#endif
