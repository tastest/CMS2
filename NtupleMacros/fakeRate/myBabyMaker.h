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

  // pfmet
  Float_t pfmet_;
  Float_t pfmetphi_;

  // did it pass the jet trigger and is this lepton unbiased
  // 0=fail 1="pass but biased" 2="pass and unbiased"  -1="pass but cant find jet trg obj"
  Int_t  hlt15u_; // HLT_Jet15U
  Int_t  hlt30u_; // HLT_Jet30U
  Int_t  hlt50u_; // HLT_Jet50U
  Int_t  l16u_;   // HLT_L1Jet6U
  Int_t  l110u_;   // HLT_L1Jet10U

  // What type of selection did it pass
  // v1xx_  v2xx_ v3xx_ are electron FO
  // fo_xx_ are muon fakeable objects
  // numxx_ are numerator selection (both muons and electrons)
  Bool_t v1_;  // electronSelectionFO_el_ttbarV1_v1
  Bool_t v2_;  // electronSelectionFO_el_ttbarV1_v2
  Bool_t v3_;  // electronSelectionFO_el_ttbarV1_v3
  Bool_t num_; // {electronSelection_ttbarV1 && (!isSpikeElectron()} (ele)  
               // NominalTTbarV2 (muons)
  Bool_t numv1_; // electronSelection_ttbarV1 (ele) NominalTTbar (muons)

  Bool_t v1SS_;  // electronSelectionFO_ssVBTF80_v1
  Bool_t v2SS_;  // electronSelectionFO_ssVBTF80_v2
  Bool_t v3SS_;  // electronSelectionFO_ssVBTF80_v3
  Bool_t numSS_; // electronSelection_ss (eletrons) Nominal (muons)
  Bool_t numNomSS_; // Nominal with SS cuts (muons)

  Bool_t v1SSAug9_;  // electronSelectionFO_ssVBTF80_v1, isData, true
  Bool_t v2SSAug9_;  // electronSelectionFO_ssVBTF80_v2, isData, true
  Bool_t v3SSAug9_;  // electronSelectionFO_ssVBTF80_v3, isData, true
  Bool_t numSSAug9_; // electronSelection_ss (electrons) Not filled for muons

  Bool_t v1Aug9_;  // identical to v1_
  Bool_t v2Aug9_;  // identical to v2_
  Bool_t v3Aug9_;  // electronSelectionFO_el_ttbarV1_v3, isData, true
  Bool_t numAug9_; // (electronSelection_ttbarV1, isData, true) && (!isSpikeElectron(iEl)) (ele)
                   // NominalTTbarV2 (muons)

  Bool_t numOct6_; // not filled for muons; electronSelection_ttbarV1_pass5 (for ele)
  Bool_t v1Oct6_;  // identical to v1_
  Bool_t v2Oct6_;  // identical to v2_
  Bool_t v3Oct6_;  // electronSelectionFO_el_ttbarV1_v3

  Bool_t numOSOct18_; // not filled for muons; electronSelection_el_OSV1 (for electrons)
  Bool_t v1OSOct18_;  // electronSelectionFO_el_OSV1_v1
  Bool_t v2OSOct18_;  // electronSelectionFO_el_OSV1_v2
  Bool_t v3OSOct18_;  // electronSelectionFO_el_OSV1_v2

  Bool_t numSSOct18_; // not filled for muons; electronSelection_ss, false, false
  Bool_t v1SSOct18_;  // electronSelectionFO_ssVBTF80_v1, false, false
  Bool_t v2SSOct18_;  // electronSelectionFO_ssVBTF80_v2, false, false
  Bool_t v3SSOct18_;  // electronSelectionFO_ssVBTF80_v3, false, false

  Bool_t numSSV2_; // not filled for muons; electronSelection_ss, false, false
  Bool_t v1SSV2_;  // electronSelectionFO_ssVBTF80_v1, false, false
  Bool_t v2SSV2_;  // electronSelectionFO_ssVBTF80_v2, false, false
  Bool_t v3SSV2_;  // electronSelectionFO_ssVBTF80_v3, false, false

  Bool_t v1_wwV0_;  // electronSelectionFO_el_wwV0_v1
  Bool_t v2_wwV0_;  // electronSelectionFO_el_wwV0_v2
  Bool_t v3_wwV0_;  // electronSelectionFO_el_wwV0_v3
  Bool_t v4_wwV0_;  // electronSelectionFO_el_wwV0_v4
  Bool_t num_wwV0_; // electronSelection_wwV0 | NominalWWV0 (muons)

  Bool_t v1_wwV0b_;  // electronSelectionFO_el_wwV0b_v1
  Bool_t v2_wwV0b_;  // electronSelectionFO_el_wwV0b_v2
  Bool_t v3_wwV0b_;  // electronSelectionFO_el_wwV0b_v3
  Bool_t v4_wwV0b_;  // electronSelectionFO_el_wwV0b_v4
  Bool_t num_wwV0b_; // electronSelection_wwV0b | NominalWWV0 (muons)


  Bool_t fo_04_;  // muonSelectionFO_mu_ttbar
  Bool_t fo_10_;  // muonSelectionFO_mu_ttbar_iso10

  Bool_t fo_muss04_;  // muonSelectionFO_mu_ss
  Bool_t fo_muss10_;  // muonSelectionFO_mu_ss_iso10

  Bool_t fo_wwV0_04_;  // muonSelectionFO_mu_ww
  Bool_t fo_wwV0_10_;  // muonSelectionFO_mu_ww_iso10

  // What electron trigger did it pass and is the electron matched to a egamma trg object
  // 0=fail 1="pass but no match" 2="pass and matched" -1="pass but egamm trg obj missing"
  Int_t ph10_;    // HLT_Photon10_L1R or HLT_Photon10_Cleaned_L1R
  Int_t ph15_;    // HLT_Photon15_L1R or HLT_Photon15_Cleaned_L1R
  Int_t ph20_;    // HLT_Photon20_Cleaned_L1R

  Int_t el10_lw_;     // HLT_Ele10_LW_L1R
  Int_t el10_sw_;     // HLT_Ele10_SW_L1R

  Int_t el10_lw_id_;  // HLT_Ele10_LW_EleId_L1R
  Int_t el10_sw_id_;  // HLT_Ele10_SW_EleId_L1R

  Int_t el15_lw_;     // HLT_Ele15_LW_L1R
  Int_t el15_sw_;     // HLT_Ele15_SW_L1R

  Int_t el15_lw_id_;  // HLT_Ele15_LW_EleId_L1R
  Int_t el15_sw_id_;  // HLT_Ele15_SW_EleId_L1R

  Int_t el15_sw_cid_; // HLT_Ele15_SW_CaloEleId_L1R

  Int_t el20_sw_;     // HLT_Ele20_SW_L1R
  Int_t el25_sw_;     // HLT_Ele25_SW_L1R

  Int_t el17_sw_;     // HLT_Ele17_SW_L1R
  Int_t el17_iso_;    // HLT_Ele17_Isol_L1R
  Int_t el17_loose_;  // HLT_Ele17_SW_LooseEleId_L1R
  Int_t el17_sw_cid_; // HLT_Ele17_SW_CaloEleId_L1R
  Int_t el17_sw_id_;  // HLT_Ele17_SW_EleId_L1R
  Int_t el17_tiso_;   // HLT_Ele17_SW_TightEleIdIsol_L1R_v1

  Int_t Del10_sw_;    // HLT_DoubleEle10_SW_L1R

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

  Float_t drel20_sw_;     // HLT_Ele20_SW_L1R
  Float_t drel25_sw_;     // HLT_Ele25_SW_L1R

  Float_t drel17_sw_;     // HLT_Ele17_SW_L1R
  Float_t drel17_iso_;    // HLT_Ele17_Isol_L1R
  Float_t drel17_loose_;  // HLT_Ele17_SW_LooseEleId_L1R
  Float_t drel17_sw_cid_; // HLT_Ele17_SW_CaloEleId_L1R
  Float_t drel17_sw_id_;  // HLT_Ele17_SW_EleId_L1R
  Float_t drel17_tiso_;   // HLT_Ele17_SW_TightEleIdIsol_L1R_v1

  Float_t drDel10_sw_;    // HLT_DoubleEle10_SW_L1R

  // What muon trigger did it pass
  // 0=fail 1="pass but no match" 2="pass and matched" -1="pass but muon trg obj missing"
  Int_t mu15_; // HLT_Mu15_v1
  Int_t mu11_; // HLT_Mu11
  Int_t mu9_;  // HLT_Mu9
  Int_t mu7_;  // HLT_Mu7
  Int_t mu5_;  // HLT_Mu5

  //  Minimm dR to the closest HLT mu object
  Float_t drmu11_;
  Float_t drmu15_;
  Float_t drmu9_;
  Float_t drmu7_;
  Float_t drmu5_;
  
  // Btag information
  Int_t nbjet_; // number of btagged jet pt>15
  Float_t dRbNear_; // dR between lepton and closest such jet
  Float_t dRbFar_; // dR between lepton and farthest such jet

  // Btag PF Corrected information
  Int_t nbpfcjet_; // number of btagged jet pt>15
  Float_t dRbpfcNear_; // dR between lepton and closest such jet
  Float_t dRbpfcFar_; // dR between lepton and farthest such jet


  // Information to do offline jet trigger selection
  Float_t ptj1_;        // highest pt jet well separated from the lepton
  Float_t ptj1_b2b_;    // highest pt jet away frmo lepton by dR >= 1.0 and dPhi > 2.5
  Float_t dphij1_b2b_;  // dphi between lepton and jet for jets away from lepton by dR >= 1.0
  Int_t   nj1_;         // number of jets above 10 GeV and away from lepton by dR >= 1.0
  Float_t ptpfj1_;        // highest pt pfjet well separated from the lepton
  Float_t ptpfj1_b2b_;    // highest pt pfjet away frmo lepton by dR >= 1.0 and dPhi > 2.5
  Float_t dphipfj1_b2b_;  // dphi between lepton and pfjet for pfjets away from lepton by dR >= 1.0
  Int_t   npfj1_;         // number of pfjets above 10 GeV and away from lepton by dR >= 1.0

  // Same for PF Corrected jets

  Float_t ptpfcj1_; // highest pt jet well separated from the lepton
  Float_t ptpfcj1_b2b_;    // highest pt jet away frmo lepton by dR >= 1.0 and dPhi > 2.5
  Float_t dphipfcj1_b2b_;  // dphi between lepton and jet for jets away from lepton by dR >= 1.0
  Int_t   npfcj1_;         // number of jets above 10 GeV and away from lepton by dR >= 1.0
  Bool_t  btagpfc_; 

  // transverse W mass
  Float_t mt_;
  Float_t pfmt_;

  // do the 3 electron charges agree?
  Bool_t q3_;

 // Missing hit info
  Int_t els_exp_innerlayers_;
  Int_t els_exp_innerlayers39X_;

  //Some MC informatio added 16 Sep 2010
  Int_t mcid_;        // els_mc_id or mus_mc_id
  Int_t mcmotherid_;  // els_mc_motherid or mus_mc_motherid

};

#endif
