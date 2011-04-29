#ifndef myBabyMaker_h
#define myBabyMaker_h

// C++ Includes

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
    Int_t pu_nPUvertices_;
  
    // Pileup - VertexMaker
    Int_t evt_nvtxs_;
  
    // Pileup - VertexMaker
    Int_t evt_ndavtxs_;

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
    Float_t iso_;             // Isolation ( truncated )
    Float_t iso_nps_;         // Isolation ( truncated with 1 GeV pedestal subtraction in ecal barrel )
    Float_t nt_iso_;          // Isolation ( not truncated )
    Float_t nt_iso_nps_;      // Isolation ( not truncated with 1 GeV pedestal subtraction in ecal barrel )
    Float_t trck_iso_;        // TRK Isolation (truncated )
    Float_t trck_nt_iso_;     // TRK Isolation ( not truncated )
    Float_t ecal_iso_;        // ECAL Isolation ( truncated )
    Float_t ecal_iso_nps_;    // ECAL Isolation ( truncated with 1 GeV pedestal subtraction in ecal barrel )
    Float_t ecal_nt_iso_;     // ECAL Isolation ( not truncated )
    Float_t ecal_nt_iso_nps_; // ECAL Isolation ( not truncated with 1 GeV pedestal subtraction in ecal barrel )
    Float_t hcal_iso_;        // HCAL Isolation ( not truncated )
    Float_t hcal_nt_iso_;     // HCAL Isolation ( truncated )

    // PV
    Float_t d0PV_wwV1_;       // electron_d0PV_wwV1(iEl)
    Float_t dzPV_wwV1_;       // electron_dzPV_wwV1(iEl)


    // Id
    Bool_t closestMuon_;  // true if els_closestMuon().at(index) == -1
    Bool_t el_id_smurfV3_;
    Bool_t el_id_vbtf80_;
    Bool_t el_id_vbtf90_;

    // Conversion Rejection
    Bool_t convHitPattern_;   // isFromConversionHitPattern(iEl)
    Bool_t convPartnerTrack_; // isFromConversionPartnerTrack(iEl)
    Bool_t convMIT_;          // isFromConversionMIT(iEl)
    Bool_t conv0MissHits_;    // true if els_exp_innerlayers().at(index) == 0

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
      
    // Electrons
    Bool_t num_el_ssV3_;
    Bool_t v1_el_ssV3_;
    Bool_t v2_el_ssV3_;
    Bool_t v3_el_ssV3_;
    
    // Muons
    Bool_t numNomSSv3_;   // NominalSSv3
    Bool_t fo_mussV3_04_; // muonSelectionFO_mu_ssV3

    // WW, HWW

    // Electrons
    Bool_t num_el_smurfV3_;
    Bool_t v1_el_smurfV1_;
    Bool_t v2_el_smurfV1_;
    Bool_t v3_el_smurfV1_;
    Bool_t v4_el_smurfV1_;

    // Muons
    Bool_t num_mu_smurfV3_;
    Bool_t fo_mu_smurf_04_;
    Bool_t fo_mu_smurf_10_;


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
    Int_t ele8_vstar_;                                               // HLT_Ele8_v*
    Int_t ele8_CaloIdL_TrkIdVL_vstar_;                               // HLT_Ele8_CaloIdL_TrkIdVL_v*
    Int_t ele8_CaloIdL_CaloIsoVL_Jet40_vstar_;                       // HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v*
    Int_t ele8_CaloIdL_CaloIsoVL_vstar_;                             // HLT_Ele8_CaloIdL_CaloIsoVL_v*
    Int_t ele17_CaloIdL_CaloIsoVL_vstar_;                            // HLT_Ele17_CaloIdL_CaloIsoVL_v*  
    Int_t photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_vstar_;      // HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v*

    Int_t ele8_version_;                                               // HLT_Ele8_v*
    Int_t ele8_CaloIdL_TrkIdVL_version_;                               // HLT_Ele8_CaloIdL_TrkIdVL_v*
    Int_t ele8_CaloIdL_CaloIsoVL_Jet40_version_;                       // HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v*
    Int_t ele8_CaloIdL_CaloIsoVL_version_;                             // HLT_Ele8_CaloIdL_CaloIsoVL_v*
    Int_t ele17_CaloIdL_CaloIsoVL_version_;                            // HLT_Ele17_CaloIdL_CaloIsoVL_v*  
    Int_t photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_version_;      // HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v*

    Float_t dr_ele8_vstar_;                                          // HLT_Ele8_v*
    Float_t dr_ele8_CaloIdL_TrkIdVL_vstar_;                          // HLT_Ele8_CaloIdL_TrkIdVL_v*
    Float_t dr_ele8_CaloIdL_CaloIsoVL_Jet40_vstar_;                  // HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v*
    Float_t dr_ele8_CaloIdL_CaloIsoVL_vstar_;                        // HLT_Ele8_CaloIdL_CaloIsoVL_v*
    Float_t dr_ele17_CaloIdL_CaloIsoVL_vstar_;                       // HLT_Ele17_CaloIdL_CaloIsoVL_v*  
    Float_t dr_photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_vstar_; // HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v*

    // Muons
    Int_t mu3_vstar_;                                                // HLT_Mu3_v*
    Int_t mu5_vstar_;                                                // HLT_Mu5_v*
    Int_t mu8_vstar_;                                                // HLT_Mu8_v*
    Int_t mu12_vstar_;                                               // HLT_Mu12_v*
    Int_t mu15_vstar_;                                               // HLT_Mu15_v*
    Int_t mu20_vstar_;                                               // HLT_Mu20_v*
    Int_t mu24_vstar_;                                               // HLT_Mu24_v*
    Int_t mu30_vstar_;                                               // HLT_Mu30_v*
    Int_t mu8_Jet40_vstar_;                                          // HLT_Mu8_Jet40_v*

    Int_t mu3_version_;                                                // HLT_Mu3_v*
    Int_t mu5_version_;                                                // HLT_Mu5_v*
    Int_t mu8_version_;                                                // HLT_Mu8_v*
    Int_t mu12_version_;                                               // HLT_Mu12_v*
    Int_t mu15_version_;                                               // HLT_Mu15_v*
    Int_t mu20_version_;                                               // HLT_Mu20_v*
    Int_t mu24_version_;                                               // HLT_Mu24_v*
    Int_t mu30_version_;                                               // HLT_Mu30_v*
    Int_t mu8_Jet40_version_;                                          // HLT_Mu8_Jet40_v*

    Float_t dr_mu3_vstar_;                                           // HLT_Mu5_v*
    Float_t dr_mu5_vstar_;                                           // HLT_Mu5_v*
    Float_t dr_mu8_vstar_;                                           // HLT_Mu8_v*
    Float_t dr_mu12_vstar_;                                          // HLT_Mu12_v*
    Float_t dr_mu15_vstar_;                                          // HLT_Mu15_v*
    Float_t dr_mu20_vstar_;                                          // HLT_Mu20_v*
    Float_t dr_mu24_vstar_;                                          // HLT_Mu24_v*
    Float_t dr_mu30_vstar_;                                          // HLT_Mu30_v*
    Float_t dr_mu8_Jet40_vstar_;                                     // HLT_Mu8_Jet40_v*

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
  
    // Pileup - VertexMaker
    evt_nvtxs_ = -1;
  
    // Pileup - VertexMaker
    evt_ndavtxs_ = -1;

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
    iso_nps_          = -999.;
    nt_iso_           = -999.;
    nt_iso_nps_       = -999.;
    trck_iso_         = -999.;
    trck_nt_iso_      = -999.;
    ecal_iso_         = -999.;
    ecal_iso_nps_     = -999.;
    ecal_nt_iso_      = -999.;
    ecal_nt_iso_nps_  = -999.;
    hcal_iso_         = -999.;
    hcal_nt_iso_      = -999.;

    closestMuon_      = false;
    el_id_smurfV3_    = false;
    el_id_vbtf80_     = false;
    el_id_vbtf90_     = false;
    convHitPattern_   = false;
    convPartnerTrack_ = false;
    convMIT_          = false;

    d0PV_wwV1_        = -999.;
    dzPV_wwV1_        = -999.;

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

    // Electrons
    num_el_ssV3_    = false;
    v1_el_ssV3_     = false;
    v2_el_ssV3_     = false;
    v3_el_ssV3_     = false;

    // Muons
    numNomSSv3_     = false;
    fo_mussV3_04_   = false;

    // WW, HWW

    // Electrons
    num_el_smurfV3_ = false;
    v1_el_smurfV1_  = false;
    v2_el_smurfV1_  = false;
    v3_el_smurfV1_  = false;
    v4_el_smurfV1_  = false;

    // Muons
    num_mu_smurfV3_ = false;
    fo_mu_smurf_04_ = false;
    fo_mu_smurf_10_ = false;

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
    ele8_vstar_                                             = 0;
    ele8_CaloIdL_TrkIdVL_vstar_                             = 0;
    ele8_CaloIdL_CaloIsoVL_Jet40_vstar_                     = 0;
    ele8_CaloIdL_CaloIsoVL_vstar_                           = 0;
    ele17_CaloIdL_CaloIsoVL_vstar_                          = 0;
    photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_vstar_    = 0;

    ele8_version_                                           = -1;
    ele8_CaloIdL_TrkIdVL_version_                           = -1;
    ele8_CaloIdL_CaloIsoVL_Jet40_version_                   = -1;
    ele8_CaloIdL_CaloIsoVL_version_                         = -1;
    ele17_CaloIdL_CaloIsoVL_version_                        = -1;
    photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_version_  = -1;

    dr_ele8_vstar_                                          = 99.0; 
    dr_ele8_CaloIdL_TrkIdVL_vstar_                          = 99.0; 
    dr_ele8_CaloIdL_CaloIsoVL_Jet40_vstar_                  = 99.0; 
    dr_ele8_CaloIdL_CaloIsoVL_vstar_                        = 99.0; 
    dr_ele17_CaloIdL_CaloIsoVL_vstar_                       = 99.0;    
    dr_photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_vstar_ = 99.0; 

    // Muons
    mu3_vstar_          = 0;
    mu5_vstar_          = 0;
    mu8_vstar_          = 0;  
    mu12_vstar_         = 0;  
    mu15_vstar_         = 0;  
    mu20_vstar_         = 0;  
    mu24_vstar_         = 0;  
    mu30_vstar_         = 0;  
    mu8_Jet40_vstar_    = 0;    

    mu3_version_        = -1;
    mu5_version_        = -1;
    mu8_version_        = -1;  
    mu12_version_       = -1;  
    mu15_version_       = -1;  
    mu20_version_       = -1;  
    mu24_version_       = -1;  
    mu30_version_       = -1;  
    mu8_Jet40_version_  = -1;    

    dr_mu3_vstar_       = 99.0;
    dr_mu5_vstar_       = 99.0;
    dr_mu8_vstar_       = 99.0;
    dr_mu12_vstar_      = 99.0; 
    dr_mu15_vstar_      = 99.0; 
    dr_mu20_vstar_      = 99.0; 
    dr_mu24_vstar_      = 99.0; 
    dr_mu30_vstar_      = 99.0;
    dr_mu8_Jet40_vstar_ = 99.0;

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

    // Calo Jets
    ptj1_       = 0.;
    ptj1_b2b_   = -999.;
    dphij1_b2b_ = -999.;
    nj1_        = 0;

    // PF Jets
    ptpfj1_       = 0.;
    ptpfj1_b2b_   = -999.;
    dphipfj1_b2b_ = -999.;
    npfj1_        = 0;

    // PF L2L3 Corrected jets
    ptpfcj1_        = 0.;
    ptpfcj1_b2b_    = -999.;
    dphipfcj1_b2b_  = -999.;
    npfcj1_         = 0;
    btagpfc_        = false;

    // PF L1FastL2L3 Corrected jets
    ptpfcL1Fj1_        = 0.;
    ptpfcL1Fj1_b2b_    = -999.;
    dphipfcL1Fj1_b2b_  = -999.;
    npfcL1Fj1_         = 0;
    btagpfcL1F_        = false;

    // JPT L2L3 Corrected jets
    ptjptcj1_        = 0.;
    ptjptcj1_b2b_    = -999.;
    dphijptcj1_b2b_  = -999.;
    njptcj1_         = 0;
    btagjptc_        = false;

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

    babyTree_->Branch("run"    , &run_    );
    babyTree_->Branch("ls"     , &ls_     );
    babyTree_->Branch("evt"    , &evt_    );
    babyTree_->Branch("weight" , &weight_ );
  
    // Pileup
    babyTree_->Branch("pu_nPUvertices" , &pu_nPUvertices_ );
    babyTree_->Branch("evt_nvtxs"      , &evt_nvtxs_      );
    babyTree_->Branch("evt_ndavtxs"    , &evt_ndavtxs_    );

    /////////////////////////// 
    // End Event Information //
    ///////////////////////////
        
        
        
    //////////////////////////// 
    // Lepton Information     //
    ////////////////////////////

    babyTree_->Branch("pt"                  , &pt_                  );
    babyTree_->Branch("eta"                 , &eta_                 );
    babyTree_->Branch("phi"                 , &phi_                 );
    babyTree_->Branch("scet"                , &scet_                );
    babyTree_->Branch("tcmet"               , &tcmet_               );
    babyTree_->Branch("tcmetphi"            , &tcmetphi_            );
    babyTree_->Branch("pfmet"               , &pfmet_               );
    babyTree_->Branch("pfmetphi"            , &pfmetphi_            );
    babyTree_->Branch("iso"                 , &iso_                 );
    babyTree_->Branch("iso_nps"             , &iso_nps_             );
    babyTree_->Branch("nt_iso"              , &nt_iso_              );
    babyTree_->Branch("nt_iso_nps"          , &nt_iso_nps_          );
    babyTree_->Branch("trck_iso"            , &trck_iso_            );
    babyTree_->Branch("trck_nt_iso"         , &trck_nt_iso_         );
    babyTree_->Branch("ecal_iso"            , &ecal_iso_            );
    babyTree_->Branch("ecal_iso_nps"        , &ecal_iso_nps_        );
    babyTree_->Branch("ecal_nt_iso"         , &ecal_nt_iso_         );
    babyTree_->Branch("ecal_nt_iso_nps"     , &ecal_nt_iso_nps_     );
    babyTree_->Branch("hcal_iso"            , &hcal_iso_            );
    babyTree_->Branch("hcal_nt_iso"         , &hcal_nt_iso_         );
    babyTree_->Branch("id"                  , &id_                  );
    babyTree_->Branch("closestMuon"         , &closestMuon_         );
    babyTree_->Branch("el_id_smurfV3"       , &el_id_smurfV3_       );
    babyTree_->Branch("el_id_vbtf80"        , &el_id_vbtf80_        );
    babyTree_->Branch("el_id_vbtf90"        , &el_id_vbtf90_        );
    babyTree_->Branch("conv0MissHits"       , &conv0MissHits_       );
    babyTree_->Branch("convHitPattern"      , &convHitPattern_      );
    babyTree_->Branch("convPartnerTrack"    , &convPartnerTrack_    );
    babyTree_->Branch("convMIT"             , &convMIT_             );
    babyTree_->Branch("mt"                  , &mt_                  );
    babyTree_->Branch("pfmt"                , &pfmt_                );
    babyTree_->Branch("q3"                  , &q3_                  );
    babyTree_->Branch("els_exp_innerlayers" , &els_exp_innerlayers_ );
    babyTree_->Branch("mcid"                , &mcid_                );
    babyTree_->Branch("mcmotherid"          , &mcmotherid_          );
    babyTree_->Branch("d0PV_wwV1"           , &d0PV_wwV1_           );
    babyTree_->Branch("dzPV_wwV1"           , &dzPV_wwV1_           );
    babyTree_->Branch("ht_calo"             , &ht_calo_             );
    babyTree_->Branch("ht_calo_L2L3"        , &ht_calo_L2L3_        );
    babyTree_->Branch("ht_jpt_L2L3"         , &ht_jpt_L2L3_         );
    babyTree_->Branch("ht_pf"               , &ht_pf_               );
    babyTree_->Branch("ht_pf_L2L3"          , &ht_pf_L2L3_          );
    babyTree_->Branch("ht_pf_L1FastL2L3"    , &ht_pf_L1FastL2L3_    );

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

    // Electrons
    babyTree_->Branch("num_el_ssV3"   , &num_el_ssV3_    );
    babyTree_->Branch("v1_el_ssV3"    , &v1_el_ssV3_     );
    babyTree_->Branch("v2_el_ssV3"    , &v2_el_ssV3_     );
    babyTree_->Branch("v3_el_ssV3"    , &v3_el_ssV3_     );

    // Muons
    babyTree_->Branch("numNomSSv3",   &numNomSSv3_       );
    babyTree_->Branch("fo_mussV3_04", &fo_mussV3_04_     );

    // OS
 
    // WW, HWW

    // Electrons
    babyTree_->Branch("num_el_smurfV3", &num_el_smurfV3_ );
    babyTree_->Branch("v1_el_smurfV1" , &v1_el_smurfV1_  );
    babyTree_->Branch("v2_el_smurfV1" , &v2_el_smurfV1_  );
    babyTree_->Branch("v3_el_smurfV1" , &v3_el_smurfV1_  );
    babyTree_->Branch("v4_el_smurfV1" , &v4_el_smurfV1_  );

    // Muons
    babyTree_->Branch("num_mu_smurfV3",  &num_mu_smurfV3_ );
    babyTree_->Branch("fo_mu_smurf_04",  &fo_mu_smurf_04_ );
    babyTree_->Branch("fo_mu_smurf_10",  &fo_mu_smurf_10_ );

    //////////
    // 2010 //
    //////////

    // ttbar
    babyTree_->Branch("numOct6"       , &numOct6_       );
    babyTree_->Branch("v1Oct6"        , &v1Oct6_        );
    babyTree_->Branch("v2Oct6"        , &v2Oct6_        );
    babyTree_->Branch("v3Oct6"        , &v3Oct6_        );
    babyTree_->Branch("num"           , &num_           );
    babyTree_->Branch("fo_04"         , &fo_04_         );
    babyTree_->Branch("fo_10"         , &fo_10_         );

    // SS
    babyTree_->Branch("v1SSV2"        , &v1SSV2_        );
    babyTree_->Branch("v2SSV2"        , &v2SSV2_        );
    babyTree_->Branch("v3SSV2"        , &v3SSV2_        );
    babyTree_->Branch("numSSV2"       , &numSSV2_       );
    babyTree_->Branch("numNomSSv2"    , &numNomSSv2_    );
    babyTree_->Branch("fo_mussV2_04"  , &fo_mussV2_04_  );
    babyTree_->Branch("fo_mussV2_10"  , &fo_mussV2_10_  );

    // OS
    babyTree_->Branch("num_OSGv1"     , &num_OSGv1_     );
    babyTree_->Branch("num_OSZv1"     , &num_OSZv1_     );
    babyTree_->Branch("numOSOct18"    , &numOSOct18_    );
    babyTree_->Branch("v1OSOct18"     , &v1OSOct18_     );
    babyTree_->Branch("v2OSOct18"     , &v2OSOct18_     );
    babyTree_->Branch("v3OSOct18"     , &v3OSOct18_     );

    // WW
    babyTree_->Branch("num_wwV1"      , &num_wwV1_      );
    babyTree_->Branch("v1_wwV1"       , &v1_wwV1_       );
    babyTree_->Branch("v2_wwV1"       , &v2_wwV1_       );
    babyTree_->Branch("v3_wwV1"       , &v3_wwV1_       );
    babyTree_->Branch("v4_wwV1"       , &v4_wwV1_       );
    babyTree_->Branch("fo_wwV1_04"    , &fo_wwV1_04_    );
    babyTree_->Branch("fo_wwV1_10"    , &fo_wwV1_10_    );
    babyTree_->Branch("fo_wwV1_10_d0" , &fo_wwV1_10_d0_ );

    //////////////////////////////////////////////////////
    // End Fake Rate Numerator & Denominator Selections //
    //////////////////////////////////////////////////////

    // Electrons
    babyTree_->Branch("ele8_vstar"                                            , &ele8_vstar_                                             );
    babyTree_->Branch("ele8_CaloIdL_TrkIdVL_vstar"                            , &ele8_CaloIdL_TrkIdVL_vstar_                             );
    babyTree_->Branch("ele8_CaloIdL_CaloIsoVL_Jet40_vstar"                    , &ele8_CaloIdL_CaloIsoVL_Jet40_vstar_                     );
    babyTree_->Branch("ele8_CaloIdL_CaloIsoVL_vstar"                          , &ele8_CaloIdL_CaloIsoVL_vstar_                           );
    babyTree_->Branch("ele17_CaloIdL_CaloIsoVL_vstar"                         , &ele17_CaloIdL_CaloIsoVL_vstar_                          );
    babyTree_->Branch("photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_vstar"   , &photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_vstar_    );

    babyTree_->Branch("ele8_version"                                            , &ele8_version_                                             );
    babyTree_->Branch("ele8_CaloIdL_TrkIdVL_version"                            , &ele8_CaloIdL_TrkIdVL_version_                             );
    babyTree_->Branch("ele8_CaloIdL_CaloIsoVL_Jet40_version"                    , &ele8_CaloIdL_CaloIsoVL_Jet40_version_                     );
    babyTree_->Branch("ele8_CaloIdL_CaloIsoVL_version"                          , &ele8_CaloIdL_CaloIsoVL_version_                           );
    babyTree_->Branch("ele17_CaloIdL_CaloIsoVL_version"                         , &ele17_CaloIdL_CaloIsoVL_version_                          );
    babyTree_->Branch("photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_version"   , &photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_version_    );

    babyTree_->Branch("dr_ele8_vstar"                                         , &dr_ele8_vstar_                                          );
    babyTree_->Branch("dr_ele8_CaloIdL_TrkIdVL_vstar"                         , &dr_ele8_CaloIdL_TrkIdVL_vstar_                          );
    babyTree_->Branch("dr_ele8_CaloIdL_CaloIsoVL_Jet40_vstar"                 , &dr_ele8_CaloIdL_CaloIsoVL_Jet40_vstar_                  );
    babyTree_->Branch("dr_ele8_CaloIdL_CaloIsoVL_vstar"                       , &dr_ele8_CaloIdL_CaloIsoVL_vstar_                        );
    babyTree_->Branch("dr_ele17_CaloIdL_CaloIsoVL_vstar"                      , &dr_ele17_CaloIdL_CaloIsoVL_vstar_                       );
    babyTree_->Branch("dr_photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_vstar", &dr_photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_vstar_ );

    babyTree_->Branch("mu3_vstar"          , &mu3_vstar_          );
    babyTree_->Branch("mu5_vstar"          , &mu5_vstar_          );
    babyTree_->Branch("mu8_vstar"          , &mu8_vstar_          );
    babyTree_->Branch("mu12_vstar"         , &mu12_vstar_         );
    babyTree_->Branch("mu15_vstar"         , &mu15_vstar_         );
    babyTree_->Branch("mu20_vstar"         , &mu20_vstar_         );
    babyTree_->Branch("mu24_vstar"         , &mu24_vstar_         );
    babyTree_->Branch("mu30_vstar"         , &mu30_vstar_         );
    babyTree_->Branch("mu8_Jet40_vstar"    , &mu8_Jet40_vstar_    );

    babyTree_->Branch("mu3_version"          , &mu3_version_          );
    babyTree_->Branch("mu5_version"          , &mu5_version_          );
    babyTree_->Branch("mu8_version"          , &mu8_version_          );
    babyTree_->Branch("mu12_version"         , &mu12_version_         );
    babyTree_->Branch("mu15_version"         , &mu15_version_         );
    babyTree_->Branch("mu20_version"         , &mu20_version_         );
    babyTree_->Branch("mu24_version"         , &mu24_version_         );
    babyTree_->Branch("mu30_version"         , &mu30_version_         );
    babyTree_->Branch("mu8_Jet40_version"    , &mu8_Jet40_version_    );

    babyTree_->Branch("dr_mu3_vstar"          , &dr_mu3_vstar_          );
    babyTree_->Branch("dr_mu5_vstar"          , &dr_mu5_vstar_          );
    babyTree_->Branch("dr_mu8_vstar"          , &dr_mu8_vstar_          );
    babyTree_->Branch("dr_mu12_vstar"         , &dr_mu12_vstar_         );
    babyTree_->Branch("dr_mu15_vstar"         , &dr_mu15_vstar_         );
    babyTree_->Branch("dr_mu20_vstar"         , &dr_mu20_vstar_         );
    babyTree_->Branch("dr_mu24_vstar"         , &dr_mu24_vstar_         );
    babyTree_->Branch("dr_mu30_vstar"         , &dr_mu30_vstar_         );
    babyTree_->Branch("dr_mu8_Jet40_vstar"    , &dr_mu8_Jet40_vstar_    );
  
    ///////////////////////  
    // End 2010 Triggers //
    ///////////////////////

    // Jets
    babyTree_->Branch("hlt15u"          , &hlt15u_         );
    babyTree_->Branch("hlt30u"          , &hlt30u_         );
    babyTree_->Branch("hlt50u"          , &hlt50u_         );
    babyTree_->Branch("l16u"            , &l16u_           );
    babyTree_->Branch("l110"            , &l110u_          );
  
    // Electrons
    babyTree_->Branch("ph10"            , &ph10_           );
    babyTree_->Branch("ph15"            , &ph15_           );
    babyTree_->Branch("ph20"            , &ph20_           );
    babyTree_->Branch("el10_lw"         , &el10_lw_        );
    babyTree_->Branch("el10_sw"         , &el10_sw_        );
    babyTree_->Branch("el10_sw_v2"      , &el10_sw_v2_     );
    babyTree_->Branch("el10_lw_id"      , &el10_lw_id_     );
    babyTree_->Branch("el10_sw_id"      , &el10_sw_id_     );
    babyTree_->Branch("el15_lw"         , &el15_lw_        );
    babyTree_->Branch("el15_sw"         , &el15_sw_        );
    babyTree_->Branch("el15_lw_id"      , &el15_lw_id_     );
    babyTree_->Branch("el15_sw_id"      , &el15_sw_id_     );
    babyTree_->Branch("el15_sw_cid"     , &el15_sw_cid_    );
    babyTree_->Branch("el20_sw"         , &el20_sw_        );
    babyTree_->Branch("el25_sw"         , &el25_sw_        );
    babyTree_->Branch("Del10_sw"        , &Del10_sw_       );
    babyTree_->Branch("el17_sw"         , &el17_sw_        );
    babyTree_->Branch("el17_sw_v2"      , &el17_sw_v2_     );
    babyTree_->Branch("el17_iso"        , &el17_iso_       );
    babyTree_->Branch("el17_loose"      , &el17_loose_     );
    babyTree_->Branch("el17_sw_cid"     , &el17_sw_cid_    );
    babyTree_->Branch("el17_sw_id"      , &el17_sw_id_     );
    babyTree_->Branch("el17_tiso"       , &el17_tiso_      );
    babyTree_->Branch("el17_tiso_v1"    , &el17_tiso_v1_   );
  
    babyTree_->Branch("drph10"          , &drph10_         );
    babyTree_->Branch("drph15"          , &drph15_         );
    babyTree_->Branch("drph20"          , drph20_          );
    babyTree_->Branch("drel10_lw"       , &drel10_lw_      );
    babyTree_->Branch("drel10_sw"       , &drel10_sw_      );
    babyTree_->Branch("drel10_sw_v2"    , &drel10_sw_v2_   );
    babyTree_->Branch("drel10_lw_id"    , &drel10_lw_id_   );
    babyTree_->Branch("drel10_sw_id"    , &drel10_sw_id_   );
    babyTree_->Branch("drel15_lw"       , &drel15_lw_      );
    babyTree_->Branch("drel15_sw"       , &drel15_sw_      );
    babyTree_->Branch("drel15_lw_id"    , &drel15_lw_id_   );
    babyTree_->Branch("drel15_sw_id"    , &drel15_sw_id_   );
    babyTree_->Branch("drel15_sw_cid"   , &drel15_sw_cid_  );
    babyTree_->Branch("drel20_sw"       , &drel20_sw_      );
    babyTree_->Branch("drel25_sw"       , &drel25_sw_      );
    babyTree_->Branch("drDel10_sw"      , &drDel10_sw_     );
    babyTree_->Branch("drel17_sw"       , &drel17_sw_      );
    babyTree_->Branch("drel17_sw_v2"    , &drel17_sw_v2_   );
    babyTree_->Branch("drel17_iso"      , &drel17_iso_     );
    babyTree_->Branch("drel17_loose"    , &drel17_loose_   );
    babyTree_->Branch("drel17_sw_cid"   , &drel17_sw_cid_  );
    babyTree_->Branch("drel17_sw_id"    , &drel17_sw_id_   );
    babyTree_->Branch("drel17_tiso"     , &drel17_tiso_    );
    babyTree_->Branch("drel17_tiso_v1"  , &drel17_tiso_v1_ );
  
    // Muons
    babyTree_->Branch("mu17"    , &mu17_ );
    babyTree_->Branch("mu15"    , &mu15_ );
    babyTree_->Branch("mu13"    , &mu13_ );
    babyTree_->Branch("mu11"    , &mu11_ );
    babyTree_->Branch("mu9"     , &mu9_  );
    babyTree_->Branch("mu7"     , &mu7_  );
    babyTree_->Branch("mu5"     , &mu5_  );
  
    babyTree_->Branch("drmu17"  , &drmu17_ );
    babyTree_->Branch("drmu15"  , &drmu15_ );
    babyTree_->Branch("drmu13"  , &drmu13_ );
    babyTree_->Branch("drmu11"  , &drmu11_ );
    babyTree_->Branch("drmu9"   , &drmu9_  );
    babyTree_->Branch("drmu7"   , &drmu7_  );
    babyTree_->Branch("drmu5"   , &drmu5_  );

    ///////////////////////  
    // End 2010 Triggers //
    ///////////////////////


        
    //////////////
    // Jets     //
    //////////////

    //
    babyTree_->Branch("ptj1"         , &ptj1_         );
    babyTree_->Branch("nj1"          , &nj1_          );
    babyTree_->Branch("ptj1_b2b"     , &ptj1_b2b_     );
    babyTree_->Branch("dphij1_b2b"   , &dphij1_b2b_   );
    babyTree_->Branch("ptpfj1"       , &ptpfj1_       );
    babyTree_->Branch("npfj1"        , &npfj1_        );
    babyTree_->Branch("ptpfj1_b2b"   , &ptpfj1_b2b_   );
    babyTree_->Branch("dphipfj1_b2b" , &dphipfj1_b2b_ );
      
    // PF L2L3 Corrected jets
    babyTree_->Branch("ptpfcj1"      , &ptpfcj1_      );
    babyTree_->Branch("npfcj1"       , &npfcj1_       );
    babyTree_->Branch("ptpfcj1_b2b"  , &ptpfcj1_b2b_  );
    babyTree_->Branch("dphipfcj1_b2b", &dphipfcj1_b2b_);
    babyTree_->Branch("btagpfc"      , &btagpfc_      );
      
    // PF L1FastL2L3 Corrected jets         
    babyTree_->Branch("ptpfcL1Fj1"      , &ptpfcL1Fj1_       );       
    babyTree_->Branch("npfcL1Fj1"       , &npfcL1Fj1_        );         
    babyTree_->Branch("ptpfcL1Fj1_b2b"  , &ptpfcL1Fj1_b2b_   );       
    babyTree_->Branch("dphipfcL1Fj1_b2b", &dphipfcL1Fj1_b2b_ );     
    babyTree_->Branch("btagpfcL1F"      , &btagpfcL1F_       );

    // JPT L2L3 Corrected jets
    babyTree_->Branch("ptjptcj1"        , &ptjptcj1_         ); 
    babyTree_->Branch("njptcj1"         , &njptcj1_          ); 
    babyTree_->Branch("ptjptcj1_b2b"    , &ptjptcj1_b2b_     ); 
    babyTree_->Branch("dphijptcj1_b2b"  , &dphijptcj1_b2b_   ); 
    babyTree_->Branch("btagjptc"        , &btagjptc_         );
      
      

    // B Tagging
    babyTree_->Branch("nbjet"        , &nbjet_        );
    babyTree_->Branch("dRNear"       , &dRbNear_      );
    babyTree_->Branch("dRFar"        , &dRbFar_       );
      
    babyTree_->Branch("nbpfcjet"     , &nbpfcjet_     );
    babyTree_->Branch("dRpfcNear"    , &dRbpfcNear_   );
    babyTree_->Branch("dRpfcFar"     , &dRbpfcFar_    );
    
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
