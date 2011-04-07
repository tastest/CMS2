#ifndef myBabyMaker_h
#define myBabyMaker_h

// C++ Includes
#include <vector>

// ROOT Includes
#include "TFile.h"
#include "TTree.h"

// TAS Includes
#include "CMS2.cc"


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

    // Id
    Bool_t  el_id_vbtf80_;
    Bool_t  el_id_vbtf90_;

    // Conversion Rejection
    Bool_t convHitPattern_;   // isFromConversionHitPattern(iEl)
    Bool_t convPartnerTrack_; // isFromConversionPartnerTrack(iEl)
    Bool_t convMIT_;          // isFromConversionMIT(iEl)

    // HT
    float ht_calo_;          
    float ht_calo_L2L3_;     
    float ht_jpt_L2L3_;      
    float ht_pf_;            
    float ht_pf_L2L3_;       
    float ht_pf_L1FastL2L3_;  

  //////////////////////////
  // End Lepton Variables //
  //////////////////////////



  //////////////////////////////////////////////////////
  // Fake Rate Numerator & Denominator Selections     //
  //////////////////////////////////////////////////////

    //////////
    // 2011 //
    //////////
  
      // SS
      Bool_t num_el_ssV3_;
      Bool_t v1_el_ssV3_;
      Bool_t v2_el_ssV3_;
      Bool_t v3_el_ssV3_;

      // WW
      Bool_t num_el_smurfV3_;
      Bool_t v1_el_smurfV1_;
      Bool_t v3_el_smurfV1_;
      Bool_t v4_el_smurfV1_;

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
    Float_t ptj1_;          // highest pt jet well separated from the lepton
    Float_t ptj1_b2b_;      // highest pt jet away frmo lepton by dR >= 1.0 and dPhi > 2.5
    Float_t dphij1_b2b_;    // dphi between lepton and jet for jets away from lepton by dR >= 1.0
    Int_t   nj1_;           // number of jets above 10 GeV and away from lepton by dR >= 1.0
    Float_t ptpfj1_;        // highest pt pfjet well separated from the lepton
    Float_t ptpfj1_b2b_;    // highest pt pfjet away frmo lepton by dR >= 1.0 and dPhi > 2.5
    Float_t dphipfj1_b2b_;  // dphi between lepton and pfjet for pfjets away from lepton by dR >= 1.0
    Int_t   npfj1_;         // number of pfjets above 10 GeV and away from lepton by dR >= 1.0
  
    // Same for PF Corrected jets
    Float_t ptpfcj1_;       // highest pt jet well separated from the lepton
    Float_t ptpfcj1_b2b_;   // highest pt jet away frmo lepton by dR >= 1.0 and dPhi > 2.5
    Float_t dphipfcj1_b2b_; // dphi between lepton and jet for jets away from lepton by dR >= 1.0
    Int_t   npfcj1_;        // number of jets above 10 GeV and away from lepton by dR >= 1.0
    Bool_t  btagpfc_; 

    // Same for PF Corrected jets
    Float_t ptpfcL1Fj1_;       // highest pt jet well separated from the lepton
    Float_t ptpfcL1Fj1_b2b_;   // highest pt jet away frmo lepton by dR >= 1.0 and dPhi > 2.5
    Float_t dphipfcL1Fj1_b2b_; // dphi between lepton and jet for jets away from lepton by dR >= 1.0
    Int_t   npfcL1Fj1_;        // number of jets above 10 GeV and away from lepton by dR >= 1.0
    Bool_t  btagpfcL1F_;


    // Same for PF Corrected jets
    Float_t ptjptcj1_;       // highest pt jet well separated from the lepton
    Float_t ptjptcj1_b2b_;   // highest pt jet away frmo lepton by dR >= 1.0 and dPhi > 2.5
    Float_t dphijptcj1_b2b_; // dphi between lepton and jet for jets away from lepton by dR >= 1.0
    Int_t   njptcj1_;        // number of jets above 10 GeV and away from lepton by dR >= 1.0
    Bool_t  btagjptc_; 

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

  //Some MC informatio added 16 Sep 2010
  Int_t mcid_;        // els_mc_id or mus_mc_id
  Int_t mcmotherid_;  // els_mc_motherid or mus_mc_motherid

};
#endif

//------------------------------------------
// Initialize baby ntuple variables
//------------------------------------------
void myBabyMaker::InitBabyNtuple () {

  /////////////////////////// 
  // Event Information     //
  ///////////////////////////
    
    // 
    run_ = -1;
    ls_  = -1;
    evt_ = -1;
    weight_ = 1.;
  
    // Pileup - PUSummaryInfoMaker
    pu_nPUvertices_ = -1;
    pu_sumptlowpt_.clear();
    pu_sumpthighpt_.clear();
    pu_instLumi_.clear();
    pu_ntrkslowpt_.clear();
    pu_ntrkshighpt_.clear();
  
    // Pileup - VertexMaker
    evt_nvtxs_ = -1;
    vtxs_xError_.clear();
    vtxs_yError_.clear();
    vtxs_zError_.clear();
    vtxs_chi2_.clear();
    vtxs_ndof_.clear();
    vtxs_sumpt_.clear();
    vtxs_isFake_.clear();
    vtxs_isValid_.clear();
    vtxs_tracksSize_.clear();
    vtxs_covMatrix_.clear();
    vtxs_position_.clear();
  
    // Pileup - VertexMaker
    evt_ndavtxs_ = -1;
    davtxs_xError_.clear();
    davtxs_yError_.clear();
    davtxs_zError_.clear();
    davtxs_chi2_.clear();
    davtxs_ndof_.clear();
    davtxs_sumpt_.clear();
    davtxs_isFake_.clear();
    davtxs_isValid_.clear();
    davtxs_tracksSize_.clear();
    davtxs_covMatrix_.clear();
    davtxs_position_.clear();

  /////////////////////////// 
  // End Event Information //
  ///////////////////////////



  //////////////////////////// 
  // Lepton Information     //
  ////////////////////////////

    id_               = -1;
    pt_               = -999.;
    eta_              = -999.;
    phi_              = -999.;
    scet_             = -999.;
    tcmet_            = -999.;
    tcmetphi_         = -999.;
    pfmet_            = -999.;
    pfmetphi_         = -999.;
    iso_              = -999.;
    el_id_vbtf80_     = false;
    el_id_vbtf90_     = false;
    convHitPattern_   = false;
    convPartnerTrack_ = false;
    convMIT_          = false;

    mt_                   = -999;
    pfmt_                 = -999;
    q3_                   = false;
    els_exp_innerlayers_  = 999;
    mcid_                 = 0;
    mcmotherid_           = 0;
      
    // HT
    ht_calo_          = -999;           
    ht_calo_L2L3_     = -999;      
    ht_jpt_L2L3_      = -999;       
    ht_pf_            = -999;            
    ht_pf_L2L3_       = -999;        
    ht_pf_L1FastL2L3_ = -999;  

  //////////////////////////// 
  // End Lepton Information //
  ////////////////////////////



  //////////////////////////////////////////////////////
  // Fake Rate Numerator & Denominator Selections     //
  //////////////////////////////////////////////////////

    //////////
    // 2011 //
    //////////

      // SS
      num_el_ssV3_    = false;
      v1_el_ssV3_     = false;
      v2_el_ssV3_     = false;
      v3_el_ssV3_     = false;

      // WW
      num_el_smurfV3_ = false;
      v1_el_smurfV1_  = false;
      v3_el_smurfV1_  = false;
      v4_el_smurfV1_  = false;

    //////////
    // 2010 //
    //////////

    // ttbar
    numOct6_ = false;
    v1Oct6_  = false;
    v2Oct6_  = false;
    v3Oct6_  = false;
    num_     = false;
    fo_04_   = false;
    fo_10_   = false;

    // SS
    numSSV2_      = false;
    v1SSV2_       = false;
    v2SSV2_       = false;
    v3SSV2_       = false;
    numNomSSv2_   = false;
    fo_mussV2_04_ = false;
    fo_mussV2_10_ = false;

    // OS
    num_OSGv1_  = false;
    num_OSZv1_  = false;
    numOSOct18_ = false;
    v1OSOct18_  = false;
    v2OSOct18_  = false;
    v3OSOct18_  = false;

    // WW
    num_wwV1_ = false;
    v1_wwV1_  = false;
    v2_wwV1_  = false;
    v3_wwV1_  = false;
    v4_wwV1_  = false;
  
    fo_wwV1_04_    = false;
    fo_wwV1_10_    = false;
    fo_wwV1_10_d0_ = false;

  //////////////////////////////////////////////////////
  // End Fake Rate Numerator & Denominator Selections //
  //////////////////////////////////////////////////////

   
  
  ///////////////////////  
  // 2011 Triggers     //
  ///////////////////////

    // Electrons
    ele8_v2_                                             = 0;
    ele8_CaloIdL_TrkIdVL_v2_                             = 0;
    ele8_CaloIdL_CaloIsoVL_Jet40_v2_                     = 0;
    ele8_CaloIdL_CaloIsoVL_v2_                           = 0;
    ele17_CaloIdL_CaloIsoVL_v2_                          = 0;
    photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v2_    = 0;

    dr_ele8_v2_                                          = 99.0; 
    dr_ele8_CaloIdL_TrkIdVL_v2_                          = 99.0; 
    dr_ele8_CaloIdL_CaloIsoVL_Jet40_v2_                  = 99.0; 
    dr_ele8_CaloIdL_CaloIsoVL_v2_                        = 99.0; 
    dr_ele17_CaloIdL_CaloIsoVL_v2_                       = 99.0;    
    dr_photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v2_ = 99.0; 

    // Muons
    mu3_v3_          = 0;  
    mu5_v3_          = 0;  
    mu8_v1_          = 0;  
    mu12_v1_         = 0;  
    mu15_v2_         = 0;  
    mu20_v1_         = 0;  
    mu24_v1_         = 0;  
    mu30_v1_         = 0;  
    mu8_Jet40_v3_    = 0;    

    dr_mu3_v3_       = 99.0;
    dr_mu5_v3_       = 99.0;
    dr_mu8_v1_       = 99.0;
    dr_mu12_v1_      = 99.0; 
    dr_mu15_v2_      = 99.0;
    dr_mu20_v1_      = 99.0;
    dr_mu24_v1_      = 99.0;
    dr_mu30_v1_      = 99.0;
    dr_mu8_Jet40_v3_ = 99.0;

  ///////////////////////  
  // End 2011 Triggers //
  ///////////////////////



  ///////////////////////  
  // 2010 Triggers     //
  ///////////////////////

    // Jets
    hlt15u_ = 0;
    hlt30u_ = 0;
    hlt50u_ = 0;
    l16u_   = 0;
    l110u_  = 0;

    // Electrons
    ph10_           = 0;
    ph15_           = 0;
    ph20_           = 0;
    el10_lw_        = 0;
    el10_sw_        = 0;
    el10_sw_v2_     = 0;
    el10_lw_id_     = 0;
    el10_sw_id_     = 0;
    el15_lw_        = 0;
    el15_sw_        = 0;
    el15_lw_id_     = 0;
    el15_sw_id_     = 0;
    el15_sw_cid_    = 0;
    el20_sw_        = 0;
    el25_sw_        = 0;
    Del10_sw_       = 0;
    el17_sw_        = 0;
    el17_sw_v2_     = 0;
    el17_iso_       = 0;
    el17_loose_     = 0;
    el17_sw_cid_    = 0;
    el17_sw_id_     = 0;
    el17_tiso_      = 0;
    el17_tiso_v1_   = 0;
  
    drph10_         = 99.0;
    drph15_         = 99.0;
    drph20_         = 99.0;
    drel10_lw_      = 99.0;
    drel10_sw_      = 99.0;
    drel10_sw_v2_   = 99.0;
    drel10_lw_id_   = 99.0;
    drel10_sw_id_   = 99.0;
    drel15_lw_      = 99.0;
    drel15_sw_      = 99.0;
    drel15_lw_id_   = 99.0;
    drel15_sw_id_   = 99.0;
    drel15_sw_cid_  = 99.0;
    drel20_sw_      = 99.0;
    drel25_sw_      = 99.0;
    drDel10_sw_     = 99.0;
    drel17_sw_      = 99.0;
    drel17_sw_v2_   = 99.0;
    drel17_iso_     = 99.0;
    drel17_loose_   = 99.0;
    drel17_sw_cid_  = 99.0;
    drel17_sw_id_   = 99.0;
    drel17_tiso_    = 99.0;
    drel17_tiso_v1_ = 99.0;
  
    // Muons
    mu5_    = 0;
    mu7_    = 0;
    mu9_    = 0;
    mu11_   = 0;
    mu13_   = 0;
    mu15_   = 0;
    mu17_   = 0;
  
    drmu5_  = 99.0;
    drmu7_  = 99.0;
    drmu9_  = 99.0;
    drmu11_ = 99.0;
    drmu13_ = 99.0;
    drmu15_ = 99.0;
    drmu17_ = 99.0;
   
  ///////////////////////  
  // End 2010 Triggers //
  ///////////////////////



  //////////////
  // Jets     //
  //////////////

    ptj1_   = 0.;
    nj1_    = 0;
    ptj1_b2b_ = -999.;
    dphij1_b2b_ = -999.;
    ptpfj1_   = 0.;
    npfj1_    = 0;
    ptpfj1_b2b_ = -999.;
    dphipfj1_b2b_ = -999.;
  
    ptpfcj1_   = 0.;
    npfcj1_    = 0;
    ptpfcj1_b2b_ = -999.;
    dphipfcj1_b2b_ = -999.;
    btagpfc_ = false;
  
    //
    nbjet_  = 0;
    dRbNear_ = 99.;
    dRbFar_ = -99.;
    nbpfcjet_  = 0;
    dRbpfcNear_ = 99.;
    dRbpfcFar_ = -99.;

  //////////////
  // End Jets //
  //////////////

}

// Book the baby ntuple
void myBabyMaker::MakeBabyNtuple(const char *babyFilename)
{
    TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
    rootdir->cd();
    babyFile_ = new TFile(Form("%s", babyFilename), "RECREATE");
    babyFile_->cd();
    babyTree_ = new TTree("tree", "A Baby Ntuple");

    /////////////////////////// 
    // Event Information     //
    ///////////////////////////

      babyTree_->Branch("run",          &run_,         "run/I"         );
      babyTree_->Branch("ls",           &ls_,          "ls/I"          );
      babyTree_->Branch("evt",          &evt_,         "evt/I"         );
      babyTree_->Branch("weight",       &weight_,      "weight/F"      );
  
      // Pileup - PUSummaryInfoMaker
      babyTree_->Branch("pu_nPUvertices", &pu_nPUvertices_ );
      babyTree_->Branch("pu_zpositions" , &pu_zpositions_  );
      babyTree_->Branch("pu_sumptlowpt" , &pu_sumptlowpt_  );
      babyTree_->Branch("pu_sumpthighpt", &pu_sumpthighpt_ );
      babyTree_->Branch("pu_instLumi"   , &pu_instLumi_    );
      babyTree_->Branch("pu_ntrkslowpt" , &pu_ntrkslowpt_  );
      babyTree_->Branch("pu_ntrkshighpt", &pu_ntrkshighpt_ );
  
      // Pileup - VertexMaker
      babyTree_->Branch("evt_nvtxs"      , &evt_nvtxs_      );
      babyTree_->Branch("vtxs_xError"    , &vtxs_xError_    );
      babyTree_->Branch("vtxs_yError"    , &vtxs_yError_    );
      babyTree_->Branch("vtxs_zError"    , &vtxs_zError_    );
      babyTree_->Branch("vtxs_chi2"      , &vtxs_chi2_      );
      babyTree_->Branch("vtxs_ndof"      , &vtxs_ndof_      );
      babyTree_->Branch("vtxs_sumpt"     , &vtxs_sumpt_     );
      babyTree_->Branch("vtxs_isFake"    , &vtxs_isFake_    );
      babyTree_->Branch("vtxs_isValid"   , &vtxs_isValid_   );
      babyTree_->Branch("vtxs_tracksSize", &vtxs_tracksSize_);
      babyTree_->Branch("vtxs_covMatrix" , &vtxs_covMatrix_ );
      babyTree_->Branch("vtxs_position"  , &vtxs_position_  );
  
      // Pileup - VertexMaker
      babyTree_->Branch("evt_ndavtxs"      , &evt_ndavtxs_      );
      babyTree_->Branch("davtxs_xError"    , &davtxs_xError_    );
      babyTree_->Branch("davtxs_yError"    , &davtxs_yError_    );
      babyTree_->Branch("davtxs_zError"    , &davtxs_zError_    );
      babyTree_->Branch("davtxs_chi2"      , &davtxs_chi2_      );
      babyTree_->Branch("davtxs_ndof"      , &davtxs_ndof_      );
      babyTree_->Branch("davtxs_sumpt"     , &davtxs_sumpt_     );
      babyTree_->Branch("davtxs_isFake"    , &davtxs_isFake_    );
      babyTree_->Branch("davtxs_isValid"   , &davtxs_isValid_   );
      babyTree_->Branch("davtxs_tracksSize", &davtxs_tracksSize_);
      babyTree_->Branch("davtxs_covMatrix" , &davtxs_covMatrix_ );
      babyTree_->Branch("davtxs_position"  , &davtxs_position_  );

    /////////////////////////// 
    // End Event Information //
    ///////////////////////////
        
        
        
    //////////////////////////// 
    // Lepton Information     //
    ////////////////////////////

      babyTree_->Branch("pt"              , &pt_      , "pt/F"      );
      babyTree_->Branch("eta"             , &eta_     , "eta/F"     );
      babyTree_->Branch("phi"             , &phi_     , "phi/F"     );
      babyTree_->Branch("scet"            , &scet_    , "scet/F"    );
      babyTree_->Branch("tcmet"           , &tcmet_   , "tcmet/F"   );
      babyTree_->Branch("tcmetphi"        , &tcmetphi_, "tcmetphi/F");
      babyTree_->Branch("pfmet"           , &pfmet_   , "pfmet/F"   );
      babyTree_->Branch("pfmetphi"        , &pfmetphi_, "pfmetphi/F");
      babyTree_->Branch("iso"             , &iso_     , "iso/F"     );
      babyTree_->Branch("id"              , &id_      , "id/I"      );
      babyTree_->Branch("el_id_vbtf80"    , &el_id_vbtf80_     );
      babyTree_->Branch("el_id_vbtf90"    , &el_id_vbtf90_     );
      babyTree_->Branch("convHitPattern"  , &convHitPattern_   );
      babyTree_->Branch("convPartnerTrack", &convPartnerTrack_ );
      babyTree_->Branch("convMIT"         , &convMIT_          );

      babyTree_->Branch("mt",          &mt_,         "mt/F"         );
      babyTree_->Branch("pfmt",          &pfmt_,         "pfmt/F"         );
      babyTree_->Branch("q3",          &q3_,         "q3/O"         );
      babyTree_->Branch("els_exp_innerlayers", &els_exp_innerlayers_, "els_exp_innerlayers/I" );
      babyTree_->Branch("mcid",       &mcid_,       "mcid/I"      );
      babyTree_->Branch("mcmotherid", &mcmotherid_, "mcmotherid/I"      );

      // HT
      babyTree_->Branch("ht_calo"         , &ht_calo_          );
      babyTree_->Branch("ht_calo_L2L3"    , &ht_calo_L2L3_     );
      babyTree_->Branch("ht_jpt_L2L3"     , &ht_jpt_L2L3_      );
      babyTree_->Branch("ht_pf"           , &ht_pf_            );
      babyTree_->Branch("ht_pf_L2L3"      , &ht_pf_L2L3_       );
      babyTree_->Branch("ht_pf_L1FastL2L3", &ht_pf_L1FastL2L3_ );

    //////////////////////////// 
    // End Lepton Information //
    ////////////////////////////



    //////////////////////////////////////////////////////
    // Fake Rate Numerator & Denominator Selections     //
    //////////////////////////////////////////////////////

      //////////
      // 2011 //
      //////////

        // SS
        babyTree_->Branch("num_el_ssV3"   , &num_el_ssV3_    );
        babyTree_->Branch("v1_el_ssV3"    , &v1_el_ssV3_     );
        babyTree_->Branch("v2_el_ssV3"    , &v2_el_ssV3_     );
        babyTree_->Branch("v3_el_ssV3"    , &v3_el_ssV3_     );
 
        // OS
 
        // WW
        babyTree_->Branch("num_el_smurfV3", &num_el_smurfV3_ );
        babyTree_->Branch("v1_el_smurfV1" , &v1_el_smurfV1_  );
        babyTree_->Branch("v3_el_smurfV1" , &v3_el_smurfV1_  );
        babyTree_->Branch("v4_el_smurfV1" , &v4_el_smurfV1_  );

      //////////
      // 2010 //
      //////////

      // ttbar
      babyTree_->Branch("numOct6",         &numOct6_,        "numOct6/O"      );
      babyTree_->Branch("v1Oct6",         &v1Oct6_,        "v1Oct6/O"      );
      babyTree_->Branch("v2Oct6",         &v2Oct6_,        "v2Oct6/O"      );
      babyTree_->Branch("v3Oct6",         &v3Oct6_,        "v3Oct6/O"      );
      babyTree_->Branch("num",         &num_,        "num/O"      );
      babyTree_->Branch("fo_04",         &fo_04_,        "fo_04/O"      );
      babyTree_->Branch("fo_10",         &fo_10_,        "fo_10/O"      );

      // SS
      babyTree_->Branch("v1SSV2",         &v1SSV2_,        "v1SSV2/O"      );
      babyTree_->Branch("v2SSV2",         &v2SSV2_,        "v2SSV2/O"      );
      babyTree_->Branch("v3SSV2",         &v3SSV2_,        "v3SSV2/O"      );
      babyTree_->Branch("numSSV2",         &numSSV2_,        "numSSV2/O"      );
      babyTree_->Branch("numNomSSv2",         &numNomSSv2_,        "numNomSSv2/O"      );
      babyTree_->Branch("fo_mussV2_04",         &fo_mussV2_04_,        "fo_mussV2_04/O"      );
      babyTree_->Branch("fo_mussV2_10",         &fo_mussV2_10_,        "fo_mussV2_10/O"      );

      // OS
      babyTree_->Branch("num_OSGv1",         &num_OSGv1_,        "num_OSGv1/O"      );
      babyTree_->Branch("num_OSZv1",         &num_OSZv1_,        "num_OSZv1/O"      );
      babyTree_->Branch("numOSOct18",         &numOSOct18_,        "numOSOct18/O"      );
      babyTree_->Branch("v1OSOct18",         &v1OSOct18_,        "v1OSOct18/O"      );
      babyTree_->Branch("v2OSOct18",         &v2OSOct18_,        "v2OSOct18/O"      );
      babyTree_->Branch("v3OSOct18",         &v3OSOct18_,        "v3OSOct18/O"      );

      // WW
      babyTree_->Branch("num_wwV1",         &num_wwV1_,        "num_wwV1/O"      );
      babyTree_->Branch("v1_wwV1",         &v1_wwV1_,        "v1_wwV1/O"      );
      babyTree_->Branch("v2_wwV1",         &v2_wwV1_,        "v2_wwV1/O"      );
      babyTree_->Branch("v3_wwV1",         &v3_wwV1_,        "v3_wwV1/O"      );
      babyTree_->Branch("v4_wwV1",         &v4_wwV1_,        "v4_wwV1/O"      );
      babyTree_->Branch("fo_wwV1_04",         &fo_wwV1_04_,        "fo_wwV1_04/O"      );
      babyTree_->Branch("fo_wwV1_10",         &fo_wwV1_10_,        "fo_wwV1_10/O"      );
      babyTree_->Branch("fo_wwV1_10_d0",         &fo_wwV1_10_d0_,        "fo_wwV1_10_d0/O"      );

    //////////////////////////////////////////////////////
    // End Fake Rate Numerator & Denominator Selections //
    //////////////////////////////////////////////////////

      // Electrons
      babyTree_->Branch("ele8_v2_"                                            , &ele8_v2_                                             );
      babyTree_->Branch("ele8_CaloIdL_TrkIdVL_v2_"                            , &ele8_CaloIdL_TrkIdVL_v2_                             );
      babyTree_->Branch("ele8_CaloIdL_CaloIsoVL_Jet40_v2_"                    , &ele8_CaloIdL_CaloIsoVL_Jet40_v2_                     );
      babyTree_->Branch("ele8_CaloIdL_CaloIsoVL_v2_"                          , &ele8_CaloIdL_CaloIsoVL_v2_                           );
      babyTree_->Branch("ele17_CaloIdL_CaloIsoVL_v2_"                         , &ele17_CaloIdL_CaloIsoVL_v2_                          );
      babyTree_->Branch("photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v2_"   , &photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v2_    );

      babyTree_->Branch("dr_ele8_v2_"                                         , &dr_ele8_v2_                                          );
      babyTree_->Branch("dr_ele8_CaloIdL_TrkIdVL_v2_"                         , &dr_ele8_CaloIdL_TrkIdVL_v2_                          );
      babyTree_->Branch("dr_ele8_CaloIdL_CaloIsoVL_Jet40_v2_"                 , &dr_ele8_CaloIdL_CaloIsoVL_Jet40_v2_                  );
      babyTree_->Branch("dr_ele8_CaloIdL_CaloIsoVL_v2_"                       , &dr_ele8_CaloIdL_CaloIsoVL_v2_                        );
      babyTree_->Branch("dr_ele17_CaloIdL_CaloIsoVL_v2_"                      , &dr_ele17_CaloIdL_CaloIsoVL_v2_                       );
      babyTree_->Branch("dr_photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v2_", &dr_photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v2_ );

      babyTree_->Branch("mu3_v3"          , &mu3_v3_          );
      babyTree_->Branch("mu5_v3"          , &mu5_v3_          );
      babyTree_->Branch("mu8_v1"          , &mu8_v1_          );
      babyTree_->Branch("mu12_v1"         , &mu12_v1_         );
      babyTree_->Branch("mu15_v2"         , &mu15_v2_         );
      babyTree_->Branch("mu20_v1"         , &mu20_v1_         );
      babyTree_->Branch("mu24_v1"         , &mu24_v1_         );
      babyTree_->Branch("mu30_v1"         , &mu30_v1_         );
      babyTree_->Branch("mu8_Jet40_v3"    , &mu8_Jet40_v3_    );

      babyTree_->Branch("dr_mu3_v3"       , &dr_mu3_v3_       );
      babyTree_->Branch("dr_mu5_v3"       , &dr_mu5_v3_       );
      babyTree_->Branch("dr_mu8_v1"       , &dr_mu8_v1_       );
      babyTree_->Branch("dr_mu12_v1"      , &dr_mu12_v1_      );
      babyTree_->Branch("dr_mu15_v2"      , &dr_mu15_v2_      );
      babyTree_->Branch("dr_mu20_v1"      , &dr_mu20_v1_      );
      babyTree_->Branch("dr_mu24_v1"      , &dr_mu24_v1_      );
      babyTree_->Branch("dr_mu30_v1"      , &dr_mu30_v1_      );
      babyTree_->Branch("dr_mu8_Jet40_v3" , &dr_mu8_Jet40_v3_ );
  
    ///////////////////////  
    // End 2010 Triggers //
    ///////////////////////

      // Jets
      babyTree_->Branch("hlt15u",       &hlt15u_,       "hlt15u/I"      );
      babyTree_->Branch("hlt30u",       &hlt30u_,       "hlt30u/I"      );
      babyTree_->Branch("hlt50u",       &hlt50u_,       "hlt50u/I"      );
      babyTree_->Branch("l16u",         &l16u_,         "l16uu/I"      );
      babyTree_->Branch("l110",         &l110u_,        "l110u/I"      );
  
      // Electrons
      babyTree_->Branch("ph10",       &ph10_,       "ph10/I"      );
      babyTree_->Branch("ph15",       &ph15_,       "ph15/I"      );
      babyTree_->Branch("ph20",       &ph20_,       "ph20/I"      );
      babyTree_->Branch("el10_lw",         &el10_lw_,         "el10_lw/I"      );
      babyTree_->Branch("el10_sw",         &el10_sw_,         "el10_sw/I"      );
      babyTree_->Branch("el10_sw_v2",         &el10_sw_v2_,         "el10_sw_v2/I"      );
      babyTree_->Branch("el10_lw_id",         &el10_lw_id_,         "el10_lw_id/I"      );
      babyTree_->Branch("el10_sw_id",         &el10_sw_id_,         "el10_sw_id/I"      );
      babyTree_->Branch("el15_lw",         &el15_lw_,         "el15_lw/I"      );
      babyTree_->Branch("el15_sw",         &el15_sw_,         "el15_sw/I"      );
      babyTree_->Branch("el15_lw_id",         &el15_lw_id_,         "el15_lw_id/I"      );
      babyTree_->Branch("el15_sw_id",         &el15_sw_id_,         "el15_sw_id/I"      );
      babyTree_->Branch("el15_sw_cid",         &el15_sw_cid_,         "el15_sw_cid/I"      );
      babyTree_->Branch("el20_sw",         &el20_sw_,         "el20_sw/I"      );
      babyTree_->Branch("el25_sw",         &el25_sw_,         "el25_sw/I"      );
      babyTree_->Branch("Del10_sw",         &Del10_sw_,         "Del10_sw/I"      );
      babyTree_->Branch("el17_sw",         &el17_sw_,         "el17_sw/I"      );
      babyTree_->Branch("el17_sw_v2",         &el17_sw_v2_,         "el17_sw_v2/I"      );
      babyTree_->Branch("el17_iso",         &el17_iso_,         "el17_iso/I"      );
      babyTree_->Branch("el17_loose",         &el17_loose_,         "el17_loose/I"      );
      babyTree_->Branch("el17_sw_cid",         &el17_sw_cid_,         "el17_sw_cid/I"      );
      babyTree_->Branch("el17_sw_id",         &el17_sw_id_,         "el17_sw_id/I"      );
      babyTree_->Branch("el17_tiso",         &el17_tiso_,         "el17_tiso/I"      );
      babyTree_->Branch("el17_tiso_v1",         &el17_tiso_v1_,         "el17_tiso_v1/I"      );
  
      babyTree_->Branch("drph10",       &drph10_,       "drph10/F"      );
      babyTree_->Branch("drph15",       &drph15_,       "drph15/F"      );
      babyTree_->Branch("drph20",       &drph20_,       "drph20/F"      );
      babyTree_->Branch("drel10_lw",         &drel10_lw_,         "drel10_lw/F"      );
      babyTree_->Branch("drel10_sw",         &drel10_sw_,         "drel10_sw/F"      );
      babyTree_->Branch("drel10_sw_v2",         &drel10_sw_v2_,         "drel10_sw_v2/F"      );
      babyTree_->Branch("drel10_lw_id",         &drel10_lw_id_,         "drel10_lw_id/F"      );
      babyTree_->Branch("drel10_sw_id",         &drel10_sw_id_,         "drel10_sw_id/F"      );
      babyTree_->Branch("drel15_lw",         &drel15_lw_,         "drel15_lw/F"      );
      babyTree_->Branch("drel15_sw",         &drel15_sw_,         "drel15_sw/F"      );
      babyTree_->Branch("drel15_lw_id",         &drel15_lw_id_,         "drel15_lw_id/F"      );
      babyTree_->Branch("drel15_sw_id",         &drel15_sw_id_,         "drel15_sw_id/F"      );
      babyTree_->Branch("drel15_sw_cid",         &drel15_sw_cid_,         "drel15_sw_cid/F"      );
      babyTree_->Branch("drel20_sw",         &drel20_sw_,         "drel20_sw/F"      );
      babyTree_->Branch("drel25_sw",         &drel25_sw_,         "drel25_sw/F"      );
      babyTree_->Branch("drDel10_sw",         &drDel10_sw_,         "drDel10_sw/F"      );
      babyTree_->Branch("drel17_sw",         &drel17_sw_,         "drel17_sw/F"      );
      babyTree_->Branch("drel17_sw_v2",         &drel17_sw_v2_,         "drel17_sw_v2/F"      );
      babyTree_->Branch("drel17_iso",         &drel17_iso_,         "drel17_iso/F"      );
      babyTree_->Branch("drel17_loose",         &drel17_loose_,         "drel17_loose/F"      );
      babyTree_->Branch("drel17_sw_cid",         &drel17_sw_cid_,         "drel17_sw_cid/F"      );
      babyTree_->Branch("drel17_sw_id",         &drel17_sw_id_,         "drel17_sw_id/F"      );
      babyTree_->Branch("drel17_tiso",         &drel17_tiso_,         "drel17_tiso/F"      );
      babyTree_->Branch("drel17_tiso_v1",         &drel17_tiso_v1_,         "drel17_tiso_v1/F"      );
  
      // Muons
      babyTree_->Branch("mu17",       &mu17_,       "mu17/I"      );
      babyTree_->Branch("mu15",       &mu15_,       "mu15/I"      );
      babyTree_->Branch("mu13",       &mu13_,       "mu13/I"      );
      babyTree_->Branch("mu11",       &mu11_,       "mu11/I"      );
      babyTree_->Branch("mu9",       &mu9_,       "mu9/I"      );
      babyTree_->Branch("mu7",       &mu7_,       "mu7/I"      );
      babyTree_->Branch("mu5",       &mu5_,       "mu5/I"      );
  
      babyTree_->Branch("drmu17",       &drmu17_,       "drmu17/F"      );
      babyTree_->Branch("drmu15",       &drmu15_,       "drmu15/F"      );
      babyTree_->Branch("drmu13",       &drmu13_,       "drmu13/F"      );
      babyTree_->Branch("drmu11",       &drmu11_,       "drmu11/F"      );
      babyTree_->Branch("drmu9",       &drmu9_,       "drmu9/F"      );
      babyTree_->Branch("drmu7",       &drmu7_,       "drmu7/F"      );
      babyTree_->Branch("drmu5",       &drmu5_,       "drmu5/F"      );

    ///////////////////////  
    // End 2010 Triggers //
    ///////////////////////


        
    //////////////
    // Jets     //
    //////////////

      babyTree_->Branch("ptj1"         , &ptj1_         , "ptj1/F"          );
      babyTree_->Branch("nj1"          , &nj1_          , "nj1/I"           );
      babyTree_->Branch("ptj1_b2b"     , &ptj1_b2b_     , "ptj1_b2b/F"      );
      babyTree_->Branch("dphij1_b2b"   , &dphij1_b2b_   , "dphij1_b2b/F"    );
      babyTree_->Branch("ptpfj1"       , &ptpfj1_       , "ptpfj1/F"        );
      babyTree_->Branch("npfj1"        , &npfj1_        , "npfj1/I"         );
      babyTree_->Branch("ptpfj1_b2b"   , &ptpfj1_b2b_   , "ptpfj1_b2b/F"    );
      babyTree_->Branch("dphipfj1_b2b" , &dphipfj1_b2b_ , "dphipfj1_b2b/F"  );
  
      babyTree_->Branch("ptpfcj1"      , &ptpfcj1_      , "ptpfcj1/F"       );
      babyTree_->Branch("npfcj1"       , &npfcj1_       , "npfcj1/I"        );
      babyTree_->Branch("ptpfcj1_b2b"  , &ptpfcj1_b2b_  , "ptpfcj1_b2b/F"   );
      babyTree_->Branch("dphipfcj1_b2b", &dphipfcj1_b2b_, "dphipfcj1_b2b/F" );
      babyTree_->Branch("btagpfc"      , &btagpfc_      , "btagpfc/O"       );

      babyTree_->Branch("nbjet"        , &nbjet_        , "nbjet/I"         );
      babyTree_->Branch("dRNear"       , &dRbNear_      , "dRbNear/F"       );
      babyTree_->Branch("dRFar"        , &dRbFar_       , "dRbFar/F"        );
  
      babyTree_->Branch("nbpfcjet"     , &nbpfcjet_     , "nbpfcjet/I"      );
      babyTree_->Branch("dRpfcNear"    , &dRbpfcNear_   , "dRbpfcNear/F"    );   
      babyTree_->Branch("dRpfcFar"     , &dRbpfcFar_    , "dRbpfcFar/F"     );

    //////////////
    // End Jets //
    //////////////
}

// Fill the baby
void myBabyMaker::FillBabyNtuple() { 
  babyTree_->Fill(); 
}

// Close the baby
void myBabyMaker::CloseBabyNtuple() {
    babyFile_->cd();
    babyTree_->Write();
    babyFile_->Close();
}
