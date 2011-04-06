#ifndef myBabyMaker_h
#define myBabyMaker_h

#include "TFile.h"
#include "TTree.h"

class myBabyMaker {

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
      
  // BABY NTUPLE VARIABLES
  TFile *babyFile_;
  TTree *babyTree_;
     
  /////////////////////////// 
  // Event Information     //
  ///////////////////////////

    // Basic Event Information
    Int_t   run_;
    Int_t   ls_;
    Int_t   evt_;
    Float_t weight_;

    // Pileup - PUSummaryInfoMaker
    Int_t         pu_nPUvertices_;
    vector<float> pu_zpositions_;
    vector<float> pu_sumptlowpt_;
    vector<float> pu_sumpthighpt_;
    vector<float> pu_instLumi_;
    vector<int>   pu_ntrkslowpt_;
    vector<int>   pu_ntrkshighpt_;
  
    // Pileup - VertexMaker
    Int_t                   evt_nvtxs_;
    vector<float>           vtxs_xError_;
    vector<float>           vtxs_yError_;
    vector<float>           vtxs_zError_;
    vector<float>           vtxs_chi2_;
    vector<float>           vtxs_ndof_;
    vector<float>           vtxs_sumpt_;
    vector<int>             vtxs_isFake_;
    vector<int>             vtxs_isValid_;
    vector<int>             vtxs_tracksSize_;
    vector< vector<float> > vtxs_covMatrix_;
    vector<LorentzVector>   vtxs_position_;
  
    // Pileup - VertexMaker
    Int_t                   evt_ndavtxs_;
    vector<float>           davtxs_xError_;
    vector<float>           davtxs_yError_;
    vector<float>           davtxs_zError_;
    vector<float>           davtxs_chi2_;
    vector<float>           davtxs_ndof_;
    vector<float>           davtxs_sumpt_;
    vector<int>             davtxs_isFake_;
    vector<int>             davtxs_isValid_;
    vector<int>             davtxs_tracksSize_;
    vector< vector<float> > davtxs_covMatrix_;
    vector<LorentzVector>   davtxs_position_;

  /////////////////////////// 
  // End Event Information //
  ///////////////////////////



  //////////////////////////
  // Lepton Variables     //
  //////////////////////////

    // Lepton pt and eta and phi
    Float_t pt_;
    Float_t eta_;
    Float_t phi_;
    Float_t scet_;
    Int_t   id_;  // \pm 11 or \pm 13
  
    // tcmet
    Float_t tcmet_;
    Float_t tcmetphi_;
  
    // pfmet
    Float_t pfmet_;
    Float_t pfmetphi_;
  
    // isolation
    Float_t iso_;

  //////////////////////////
  // End Lepton Variables //
  //////////////////////////



  //////////////////////////////////////////////////////
  // Fake Rate Numerator & Denominator Selections     //
  //////////////////////////////////////////////////////

    //////////
    // 2011 //
    //////////
  
    //////////
    // 2010 //
    //////////
  
    // ttbar
  
      // electrons ( not filled for muons )
      Bool_t numOct6_;      // electronSelection_ttbarV1_pass5
      Bool_t v1Oct6_;       // electronSelectionFO_el_ttbarV1_v1_pass5
      Bool_t v2Oct6_;       // electronSelectionFO_el_ttbarV1_v1_pass5
      Bool_t v3Oct6_;       // electronSelectionFO_el_ttbarV1_v1_pass5
  
      // muons     ( not filled for electrons )
      Bool_t num_;          // NominalTTbarV2
      Bool_t fo_04_;        // muonSelectionFO_mu_ttbar
      Bool_t fo_10_;        // muonSelectionFO_mu_ttbar_iso10
  
    // Same Sign Susy

      // electrons ( not filled for muons )
      Bool_t numSSV2_;      // electronSelection_ssV2, false, false && !isSpikeElectron
      Bool_t v1SSV2_;       // electronSelectionFOV2_ssVBTF80_v1, false, false
      Bool_t v2SSV2_;       // electronSelectionFOV2_ssVBTF80_v2, false, false
      Bool_t v3SSV2_;       // electronSelectionFOV2_ssVBTF80_v3, false, false

      // muons
      Bool_t numNomSSv2_;   // NominalSSv2
      Bool_t fo_mussV2_04_; // muonSelectionFO_mu_ssV2
      Bool_t fo_mussV2_10_; // muonSelectionFO_mu_ssV2_iso10
  
    // Opposite Sign Susy
    Bool_t num_OSGv1_;      // OSGeneric_v1 (muons) | nothing for ele
    Bool_t num_OSZv1_;      // OSZ_v1 (muons)       | nothing for ele
    Bool_t numOSOct18_;     // not filled for muons; electronSelection_el_OSV1 (for electrons)
    Bool_t v1OSOct18_;      // electronSelectionFO_el_OSV1_v1
    Bool_t v2OSOct18_;      // electronSelectionFO_el_OSV1_v2
    Bool_t v3OSOct18_;      // electronSelectionFO_el_OSV1_v2
  
    // WW
    Bool_t num_wwV1_;       // electronSelection_wwV1 | NominalWWV1 (muons)
    Bool_t v1_wwV1_;        // electronSelectionFO_el_wwV1_v1
    Bool_t v2_wwV1_;        // electronSelectionFO_el_wwV1_v2
    Bool_t v3_wwV1_;        // electronSelectionFO_el_wwV1_v3
    Bool_t v4_wwV1_;        // electronSelectionFO_el_wwV1_v4
    Bool_t fo_wwV1_04_;     // muonSelectionFO_mu_wwV1
    Bool_t fo_wwV1_10_;     // muonSelectionFO_mu_wwV1_iso10
    Bool_t fo_wwV1_10_d0_;  // muonSelectionFO_mu_wwV1_iso10_d0

  //////////////////////////////////////////////////////
  // End Fake Rate Numerator & Denominator Selections //
  //////////////////////////////////////////////////////



  ///////////////////////  
  // 2011 Triggers     //
  ///////////////////////

    // Triggers & HLT matching
    // 0  = fail 
    // 1  = "pass but no match" 
    // 2  = "pass and matched"
    // -1 = "pass but egamm trg obj missing"

    // 2011 Trigger Documenation Rules 
    //
    //  1. The trigger variable name = the trigger name with:
    //      "HLT_Ele" -> "ele"  (electrons)
    //      "HLT_Mu"  -> "mu"   (muons)
    //  2. Each trigger variable name should be commented with the trigger name
    //  3. Delta R to the closest trigger object should be stored for each trigger by prepending "dr_" to the trigger variable name

    // Electrons
    Int_t ele8_v2_;                                               // HLT_Ele8_v2
    Int_t ele8_CaloIdL_TrkIdVL_v2_;                               // HLT_Ele8_CaloIdL_TrkIdVL_v2
    Int_t ele8_CaloIdL_CaloIsoVL_Jet40_v2_;                       // HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v2
    Int_t ele8_CaloIdL_CaloIsoVL_v2_;                             // HLT_Ele8_CaloIdL_CaloIsoVL_v2
    Int_t ele17_CaloIdL_CaloIsoVL_v2_;                            // HLT_Ele17_CaloIdL_CaloIsoVL_v2  
    Int_t photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v2_;      // HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v2

    Float_t dr_ele8_v2_;                                          // HLT_Ele8_v2
    Float_t dr_ele8_CaloIdL_TrkIdVL_v2_;                          // HLT_Ele8_CaloIdL_TrkIdVL_v2
    Float_t dr_ele8_CaloIdL_CaloIsoVL_Jet40_v2_;                  // HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v2
    Float_t dr_ele8_CaloIdL_CaloIsoVL_v2_;                        // HLT_Ele8_CaloIdL_CaloIsoVL_v2
    Float_t dr_ele17_CaloIdL_CaloIsoVL_v2_;                       // HLT_Ele17_CaloIdL_CaloIsoVL_v2  
    Float_t dr_photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v2_; // HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v2

    // Muons
    Int_t mu3_v3_;                                                // HLT_Mu3_v3
    Int_t mu5_v3_;                                                // HLT_Mu5_v3
    Int_t mu8_v1_;                                                // HLT_Mu8_v1
    Int_t mu12_v1_;                                               // HLT_Mu12_v1
    Int_t mu15_v2_;                                               // HLT_Mu15_v2
    Int_t mu20_v1_;                                               // HLT_Mu20_v1
    Int_t mu24_v1_;                                               // HLT_Mu24_v1
    Int_t mu30_v1_;                                               // HLT_Mu30_v1
    Int_t mu8_Jet40_v3_;                                          // HLT_Mu8_Jet40_v3

    Float_t dr_mu3_v3_;                                           // HLT_Mu3_v3
    Float_t dr_mu5_v3_;                                           // HLT_Mu5_v3
    Float_t dr_mu8_v1_;                                           // HLT_Mu8_v1
    Float_t dr_mu12_v1_;                                          // HLT_Mu12_v1
    Float_t dr_mu15_v2_;                                          // HLT_Mu15_v2
    Float_t dr_mu20_v1_;                                          // HLT_Mu20_v1
    Float_t dr_mu24_v1_;                                          // HLT_Mu24_v1
    Float_t dr_mu30_v1_;                                          // HLT_Mu30_v1
    Float_t dr_mu8_Jet40_v3_;                                     // HLT_Mu8_Jet40_v3

  ///////////////////////  
  // End 2011 Triggers //
  ///////////////////////



  ///////////////////////  
  // 2010 Triggers     //
  ///////////////////////

    // Jet Triggers
  
    // did it pass the jet trigger and is this lepton unbiased
    // 0=fail 1="pass but biased" 2="pass and unbiased"  -1="pass but cant find jet trg obj"
    Int_t  hlt15u_; // HLT_Jet15U
    Int_t  hlt30u_; // HLT_Jet30U
    Int_t  hlt50u_; // HLT_Jet50U
    Int_t  l16u_;   // HLT_L1Jet6U
    Int_t  l110u_;   // HLT_L1Jet10U
  
    // Electron Triggers
  
    // What electron trigger did it pass and is the electron matched to a egamma trg object
    // 0=fail 1="pass but no match" 2="pass and matched" -1="pass but egamm trg obj missing"
    Int_t ph10_;    // HLT_Photon10_L1R or HLT_Photon10_Cleaned_L1R
    Int_t ph15_;    // HLT_Photon15_L1R or HLT_Photon15_Cleaned_L1R
    Int_t ph20_;    // HLT_Photon20_Cleaned_L1R
    Int_t el10_lw_;     // HLT_Ele10_LW_L1R
    Int_t el10_sw_;     // HLT_Ele10_SW_L1R
    Int_t el10_sw_v2_;  // HLT_Ele10_SW_L1R_v2
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
    Int_t el17_sw_v2_;  // HLT_Ele17_SW_L1R_v2
    Int_t el17_iso_;    // HLT_Ele17_Isol_L1R
    Int_t el17_loose_;  // HLT_Ele17_SW_LooseEleId_L1R
    Int_t el17_sw_cid_; // HLT_Ele17_SW_CaloEleId_L1R
    Int_t el17_sw_id_;  // HLT_Ele17_SW_EleId_L1R
    Int_t el17_tiso_;   // HLT_Ele17_SW_TightEleIdIsol_L1R
    Int_t el17_tiso_v1_;// HLT_Ele17_SW_TightEleIdIsol_L1R_v1
    Int_t Del10_sw_;    // HLT_DoubleEle10_SW_L1R
  
    //  Minimm dR to the closest eg object
    Float_t drph10_;    // HLT_Photon10_L1R or HLT_Photon10_Cleaned_L1R
    Float_t drph15_;    // HLT_Photon15_L1R or HLT_Photon15_Cleaned_L1R
    Float_t drph20_;    // HLT_Photon20_Cleaned_L1R
    Float_t drel10_lw_;    // HLT_Ele10_LW_L1R
    Float_t drel10_sw_;    // HLT_Ele10_SW_L1R
    Float_t drel10_sw_v2_; // HLT_Ele10_SW_L1R_v2
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
    Float_t drel17_sw_v2_;  // HLT_Ele17_SW_L1R_v2
    Float_t drel17_iso_;    // HLT_Ele17_Isol_L1R
    Float_t drel17_loose_;  // HLT_Ele17_SW_LooseEleId_L1R
    Float_t drel17_sw_cid_; // HLT_Ele17_SW_CaloEleId_L1R
    Float_t drel17_sw_id_;  // HLT_Ele17_SW_EleId_L1R
    Float_t drel17_tiso_;   // HLT_Ele17_SW_TightEleIdIsol_L1R
    Float_t drel17_tiso_v1_;// HLT_Ele17_SW_TightEleIdIsol_L1R_v1
    Float_t drDel10_sw_;    // HLT_DoubleEle10_SW_L1R
  
    // Muon Triggers
  
    // What muon trigger did it pass
    // 0=fail 1="pass but no match" 2="pass and matched" -1="pass but muon trg obj missing"
    Int_t mu17_;      // HLT_Mu17_v1
    Int_t mu15_;      // HLT_Mu15_v1
    Int_t mu13_;      // HLT_Mu13_v1
    Int_t mu11_;      // HLT_Mu11
    Int_t mu9_;       // HLT_Mu9
    Int_t mu7_;       // HLT_Mu7
    Int_t mu5_;       // HLT_Mu5
  
    //  Minimm dR to the closest HLT mu object
    Float_t drmu17_;  // HLT_Mu17_v1
    Float_t drmu15_;  // HLT_Mu15_v1
    Float_t drmu13_;  // HLT_Mu13_v1
    Float_t drmu11_;  // HLT_Mu11
    Float_t drmu9_;   // HLT_Mu9
    Float_t drmu7_;   // HLT_Mu7
    Float_t drmu5_;   // HLT_Mu5
  
  ///////////////////////  
  // End 2010 Triggers //
  ///////////////////////



  //////////////
  // Jets     //
  //////////////

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

  //////////////
  // End Jets //
  //////////////



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
