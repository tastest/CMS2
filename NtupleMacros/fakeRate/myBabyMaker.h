#ifndef myBabyMaker_h
#define myBabyMaker_h

// C++ Includes

// ROOT Includes
#include "TFile.h"
#include "TTree.h"
#include "TPRegexp.h"
#include "Math/LorentzVector.h"

// TAS includes
//#include "CMS2.cc"

// lorentz vector of floats 
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

class myBabyMaker {

public:

    myBabyMaker();
    ~myBabyMaker() {
        delete babyFile_;
        delete babyTree_;
    };
    void MakeBabyNtuple (const char *);
    void InitBabyNtuple ();
    void FillBabyNtuple ();
    void CloseBabyNtuple ();
    void ScanChain (TChain *, const char *, bool, int, std::string jetcorrPath="../CORE/jetcorr/data/");
    void SetGoodRunList(const char* fileName, bool goodRunIsJson=false);

private:
      
    // BABY NTUPLE VARIABLES
    TFile *babyFile_;
    TTree *babyTree_;

    // good run list
    Bool_t goodrun_is_json;
     
    /////////////////////////// 
    // Event Information     //
    ///////////////////////////

    // Basic Event Information
    Int_t   run_;
    Int_t   ls_;
    UInt_t   evt_;
    Float_t weight_;

    // Pileup - PUSummaryInfoMaker
    Int_t pu_nPUvertices_;
  
    // Pileup - VertexMaker
    Int_t evt_nvtxs_;
  
    // Pileup - VertexMaker
    Int_t evt_ndavtxs_;

    // event level variables (number of additional objects besides the FO under consideration)
    Int_t nFOels_;
    Int_t nFOmus_;
    Int_t ngsfs_;
    Int_t nmus_;

    /////////////////////////// 
    // End Event Information //
    ///////////////////////////



    //////////////////////////
    // Lepton Variables     //
    //////////////////////////

    // Lepton pt and eta and phi
    Float_t pt_;
    Float_t eta_;
    Float_t sceta_;
    Float_t phi_;
    Float_t scet_;
    Int_t   id_;  // \pm 11 or \pm 13
    Float_t hoe_;
  
    // some useful lepton 4 vectors
    LorentzVector lp4_; // 4-vector of the lepton
    LorentzVector foel_p4_; // 4-vector of the highest additional electron FO in the event
    Int_t foel_id_;
    LorentzVector fomu_p4_; // 4-vector of the highest additional muon FO in the event
    Int_t fomu_id_;
    Float_t foel_mass_;
    Float_t fomu_mass_;

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
    Float_t nt_pfiso03_;      // PF Isolation (not truncated) with a cone size of 0.3
    Float_t nt_pfiso04_;      // PF Isolation (not truncated) with a cone size of 0.4

    // PV
    Float_t d0PV_wwV1_;       // electron_d0PV_wwV1(iEl)
    Float_t dzPV_wwV1_;       // electron_dzPV_wwV1(iEl)


    // Id
    Bool_t closestMuon_;  // true if els_closestMuon().at(index) == -1
    Bool_t el_id_smurfV5_;
    Bool_t el_id_vbtf80_;
    Bool_t el_id_vbtf90_;
    Float_t el_lh_;
    Float_t el_mva_;

    // Z mass variables
    Float_t mz_fo_gsf_;
    Float_t mz_gsf_iso_;
    Float_t mz_fo_ctf_;
    Float_t mz_ctf_iso_;
    Float_t mupsilon_fo_mu_;
    Float_t mupsilon_mu_iso_;

    Bool_t mu_isCosmic_;

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

    // MC truth information
    // Int_t mc1id_;
    // Float_t mc1pt_;
    // Float_t mc1dr_;
    Int_t mc3id_;
    Float_t mc3pt_;
    Float_t mc3dr_;
    Int_t leptonIsFromW_;

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

    // Electrons
    Bool_t num_el_ssV4_;
    Bool_t v1_el_ssV4_;
    Bool_t v2_el_ssV4_;
    Bool_t v3_el_ssV4_;

    // Electrons
    Bool_t num_el_ssV5_;
    Bool_t v1_el_ssV5_;
    Bool_t v2_el_ssV5_;
    Bool_t v3_el_ssV5_;
    Bool_t num_el_ssV5_noIso_;

    // Electrons
    Bool_t num_el_ssV6_;
    Bool_t v1_el_ssV6_;
    Bool_t v2_el_ssV6_;
    Bool_t v3_el_ssV6_;
    Bool_t num_el_ssV6_noIso_;
    
    // Muons
    Bool_t numNomSSv3_;   // NominalSSv3
    Bool_t fo_mussV3_04_; // muonSelectionFO_mu_ssV3

    // Muons
    Bool_t numNomSSv4_;   // NominalSSv4
    Bool_t fo_mussV4_04_; // muonSelectionFO_mu_ssV4
    Bool_t numNomSSv4noIso_; // NominalSSv4 with no isolation applied
    Bool_t fo_mussV4_noIso_; // muonSleectionFO_mu_ssV4 with no isolation applied

    // WW, HWW

    // Electrons
    Bool_t num_el_smurfV6_;
    Bool_t num_el_smurfV6lh_;
    Bool_t v1_el_smurfV1_;
    Bool_t v2_el_smurfV1_;
    Bool_t v3_el_smurfV1_;
    Bool_t v4_el_smurfV1_;

    // Muons
    Bool_t num_mu_smurfV6_;
    Bool_t fo_mu_smurf_04_;
    Bool_t fo_mu_smurf_10_;


    // OS
    Bool_t num_el_OSV2_;   // electronSelection_el_OSV2
    Bool_t num_mu_OSGV2_;  // OSGeneric_v2
    Bool_t num_mu_OSZV2_;  // OSZ_v2
    Bool_t fo_el_OSV2_;    // electronSelection_el_OSV2_FO
    Bool_t fo_mu_OSGV2_;   // OSGeneric_v2_FO

    Bool_t num_el_OSV3_;   // electronSelection_el_OSV3
    Bool_t num_mu_OSGV3_;  // OSGeneric_v3
    Bool_t fo_el_OSV3_;    // electronSelection_el_OSV3_FO
    Bool_t fo_mu_OSGV3_;   // OSGeneric_v3_FO


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
    Int_t ele8_vstar_;                                                          // HLT_Ele8_v*
    Int_t ele8_CaloIdL_TrkIdVL_vstar_;                                          // HLT_Ele8_CaloIdL_TrkIdVL_v*
    Int_t ele8_CaloIdL_CaloIsoVL_Jet40_vstar_;                                  // HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v*
    // moved to 2012 Int_t ele8_CaloIdL_CaloIsoVL_vstar_;                       // HLT_Ele8_CaloIdL_CaloIsoVL_v*
    // moved to 2012 Int_t ele17_CaloIdL_CaloIsoVL_vstar_;                      // HLT_Ele17_CaloIdL_CaloIsoVL_v*
    Int_t ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_vstar_;                       // HLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v*
    // moved to 2012 Int_t ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar_;      // HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*
    Int_t photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_vstar_;                 // HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v*

    Int_t ele8_version_;                                                        // HLT_Ele8_v*
    Int_t ele8_CaloIdL_TrkIdVL_version_;                                        // HLT_Ele8_CaloIdL_TrkIdVL_v*
    Int_t ele8_CaloIdL_CaloIsoVL_Jet40_version_;                                // HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v*
    // moved to 2012 Int_t ele8_CaloIdL_CaloIsoVL_version_;                     // HLT_Ele8_CaloIdL_CaloIsoVL_v*
    // moved to 2012 Int_t ele17_CaloIdL_CaloIsoVL_version_;                    // HLT_Ele17_CaloIdL_CaloIsoVL_v*  
    Int_t ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_version_;                     // HLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v*
    // moved to 2012 Int_t ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version_;    // HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*
    Int_t photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_version_;               // HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v*

    Float_t dr_ele8_vstar_;                                                     // HLT_Ele8_v*
    Float_t dr_ele8_CaloIdL_TrkIdVL_vstar_;                                     // HLT_Ele8_CaloIdL_TrkIdVL_v*
    Float_t dr_ele8_CaloIdL_CaloIsoVL_Jet40_vstar_;                             // HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v*
    // moved to 2012 Float_t dr_ele8_CaloIdL_CaloIsoVL_vstar_;                  // HLT_Ele8_CaloIdL_CaloIsoVL_v*
    // moved to 2012 Float_t dr_ele17_CaloIdL_CaloIsoVL_vstar_;                 // HLT_Ele17_CaloIdL_CaloIsoVL_v*  
    Float_t dr_ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_vstar_;                  // HLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v*
    // moved to 2012 Float_t dr_ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar_; // HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*
    Float_t dr_photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_vstar_;            // HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v*

    // Muons
    Int_t mu3_vstar_;            // HLT_Mu3_v*
    // moved to 2012 Int_t mu5_vstar_;            // HLT_Mu5_v*
    // moved to 2012 Int_t mu8_vstar_;            // HLT_Mu8_v*
    // moved to 2012 Int_t mu12_vstar_;           // HLT_Mu12_v*
    Int_t mu15_vstar_;           // HLT_Mu15_v*
    Int_t mu20_vstar_;           // HLT_Mu20_v*
    Int_t mu24_vstar_;           // HLT_Mu24_v*
    Int_t mu30_vstar_;           // HLT_Mu30_v*
    Int_t mu8_Jet40_vstar_;      // HLT_Mu8_Jet40_v*

    Int_t mu3_version_;          // HLT_Mu3_v*
    // moved to 2012 Int_t mu5_version_;          // HLT_Mu5_v*
    // moved to 2012 Int_t mu8_version_;          // HLT_Mu8_v*
    // moved to 2012 Int_t mu12_version_;         // HLT_Mu12_v*
    Int_t mu15_version_;         // HLT_Mu15_v*
    Int_t mu20_version_;         // HLT_Mu20_v*
    Int_t mu24_version_;         // HLT_Mu24_v*
    Int_t mu30_version_;         // HLT_Mu30_v*
    Int_t mu8_Jet40_version_;    // HLT_Mu8_Jet40_v*

    Float_t dr_mu3_vstar_;       // HLT_Mu5_v*
    // moved to 2012 Float_t dr_mu5_vstar_;       // HLT_Mu5_v*
    // moved to 2012 Float_t dr_mu8_vstar_;       // HLT_Mu8_v*
    // moved to 2012 Float_t dr_mu12_vstar_;      // HLT_Mu12_v*
    Float_t dr_mu15_vstar_;      // HLT_Mu15_v*
    Float_t dr_mu20_vstar_;      // HLT_Mu20_v*
    Float_t dr_mu24_vstar_;      // HLT_Mu24_v*
    Float_t dr_mu30_vstar_;      // HLT_Mu30_v*
    Float_t dr_mu8_Jet40_vstar_; // HLT_Mu8_Jet40_v*

    ///////////////////////  
    // End 2011 Triggers //
    ///////////////////////


    ///////////////////////  
    // 2012 Triggers     //
    ///////////////////////

/*
HLT_Ele17_CaloIdL_CaloIsoVL_v15
HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v4 
HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4
HLT_Ele8_CaloIdL_CaloIsoVL_v15
HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v4
HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v13
HLT_Ele8_CaloIdT_TrkIdVL_v3
HLT_Mu5_v18    
HLT_Mu8_v16
HLT_Mu12_v16
HLT_Mu17_v3
HLT_Mu15_eta2p1_v3
HLT_Mu24_eta2p1_v3
HLT_Mu30_eta2p1_v3
*/

    // Electrons
    Int_t ele17_CaloIdL_CaloIsoVL_vstar_;                              // HLT_Ele17_CaloIdL_CaloIsoVL_v*
    Int_t ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar_;             // HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*
    Int_t ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar_;       // HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v*
    Int_t ele8_CaloIdL_CaloIsoVL_vstar_;                               // HLT_Ele8_CaloIdL_CaloIsoVL_v*
    Int_t ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar_;              // HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*
    Int_t ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar_;        // HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v*
    Int_t ele8_CaloIdT_TrkIdVL_vstar_;                                 // HLT_Ele8_CaloIdT_TrkIdVL_v*

    Int_t ele17_CaloIdL_CaloIsoVL_version_;                            // HLT_Ele17_CaloIdL_CaloIsoVL_v*
    Int_t ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version_;           // HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*
    Int_t ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version_;     // HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v*
    Int_t ele8_CaloIdL_CaloIsoVL_version_;                             // HLT_Ele8_CaloIdL_CaloIsoVL_v*
    Int_t ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version_;            // HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*
    Int_t ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version_;      // HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v*
    Int_t ele8_CaloIdT_TrkIdVL_version_;                               // HLT_Ele8_CaloIdT_TrkIdVL_v*

    Float_t dr_ele17_CaloIdL_CaloIsoVL_vstar_;                         // HLT_Ele17_CaloIdL_CaloIsoVL_v*
    Float_t dr_ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar_;        // HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*
    Float_t dr_ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar_;  // HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v*
    Float_t dr_ele8_CaloIdL_CaloIsoVL_vstar_;                          // HLT_Ele8_CaloIdL_CaloIsoVL_v*
    Float_t dr_ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar_;         // HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*
    Float_t dr_ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar_;   // HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v*
    Float_t dr_ele8_CaloIdT_TrkIdVL_vstar_;                            // HLT_Ele8_CaloIdT_TrkIdVL_v*

    // Muons
    Int_t mu5_vstar_;               // HLT_Mu5_v*  // also in 2011
    Int_t mu8_vstar_;               // HLT_Mu8_v*  // also in 2011
    Int_t mu12_vstar_;              // HLT_Mu12_v* // also in 2011
    Int_t mu17_vstar_;              // HLT_Mu17_v* 
    Int_t mu15_eta2p1_vstar_;       // HLT_Mu15_eta2p1_v* 
    Int_t mu24_eta2p1_vstar_;       // HLT_Mu24_eta2p1_v* 
    Int_t mu30_eta2p1_vstar_;       // HLT_Mu30_eta2p1_v* 

    Int_t mu5_version_;             // HLT_Mu5_v*  // also in 2011
    Int_t mu8_version_;             // HLT_Mu8_v*  // also in 2011
    Int_t mu12_version_;            // HLT_Mu12_v* // also in 2011
    Int_t mu17_version_;            // HLT_Mu17_v* 
    Int_t mu15_eta2p1_version_;     // HLT_Mu15_eta2p1_v* 
    Int_t mu24_eta2p1_version_;     // HLT_Mu24_eta2p1_v* 
    Int_t mu30_eta2p1_version_;     // HLT_Mu30_eta2p1_v* 

    Float_t dr_mu5_vstar_;          // HLT_Mu5_v*  // also in 2011
    Float_t dr_mu8_vstar_;          // HLT_Mu8_v*  // also in 2011
    Float_t dr_mu12_vstar_;         // HLT_Mu12_v* // also in 2011
    Float_t dr_mu17_vstar_;         // HLT_Mu17_v* 
    Float_t dr_mu15_eta2p1_vstar_;  // HLT_Mu15_eta2p1_v* 
    Float_t dr_mu24_eta2p1_vstar_;  // HLT_Mu24_eta2p1_v* 
    Float_t dr_mu30_eta2p1_vstar_;  // HLT_Mu30_eta2p1_v* 

    ///////////////////////  
    // End 2012 Triggers //
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
  
    // Same for PF Corrected jets L2L3
    Float_t ptpfcj1_;       // highest pt jet well separated from the lepton
    Float_t ptpfcj1_b2b_;   // highest pt jet away frmo lepton by dR >= 1.0 and dPhi > 2.5
    Float_t dphipfcj1_b2b_; // dphi between lepton and jet for jets away from lepton by dR >= 1.0
    Int_t   npfcj1_;        // number of jets above 10 GeV and away from lepton by dR >= 1.0
    Bool_t  btagpfc_; 

    // Same for PF Corrected jets L1FastL2L3
    Float_t emfpfcL1Fj1_;      // EMF of hight pt PF jet well separated from lepton
    Float_t ptpfcL1Fj1_;       // highest pt jet well separated from the lepton
    Float_t dphipfcL1Fj1_;     // dphi between highest pt jet well separated from the lepton and lepton
    Float_t ptpfcL1Fj1_b2b_;   // highest pt jet away frmo lepton by dR >= 1.0 and dPhi > 2.5
    Float_t dphipfcL1Fj1_b2b_; // dphi between lepton and jet for jets away from lepton by dR >= 1.0
    Int_t   npfcL1Fj1_;        // number of jets above 10 GeV and away from lepton by dR >= 1.0
    Bool_t  btagpfcL1F_;
    Int_t   npfc30L1Fj1_;      // number of jets above 30 GeV and away from lepton by dR >= 1.0
    Int_t   npfc40L1Fj1_;      // number of jets above 40 GeV and away from lepton by dR >= 1.0

    // Same for btagged PF Corrected jets L1FastL2L3
    Float_t ptbtagpfcL1Fj1_;       // highest pt btagged jet well separated from the lepton
    Float_t dphibtagpfcL1Fj1_;     // dphi between highest pt btagged jet well separated from the lepton and lepton

    // Same for corrected JPT jets
    Float_t ptjptcj1_;       // highest pt jet well separated from the lepton
    Float_t ptjptcj1_b2b_;   // highest pt jet away frmo lepton by dR >= 1.0 and dPhi > 2.5
    Float_t dphijptcj1_b2b_; // dphi between lepton and jet for jets away from lepton by dR >= 1.0
    Int_t   njptcj1_;        // number of jets above 10 GeV and away from lepton by dR >= 1.0
    Bool_t  btagjptc_; 
    
    Float_t rho_;

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

    // regular expressions for triggers
    TPMERegexp ele8_regexp;
    TPMERegexp ele8_CaloIdL_TrkIdVL_regexp;
    TPMERegexp ele8_CaloIdL_CaloIsoVL_regexp;
    TPMERegexp ele8_CaloIdL_CaloIsoVL_Jet40_regexp;
    TPMERegexp ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_regexp;
    TPMERegexp ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_regexp;
    TPMERegexp ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_regexp;
    TPMERegexp ele8_CaloIdT_TrkIdVL_regexp;
    TPMERegexp ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_vstar_regexp;
    TPMERegexp ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar_regexp;
    TPMERegexp ele17_CaloIdL_CaloIsoVL_regexp;
    TPMERegexp ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_regexp;
    TPMERegexp ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_regexp;
    TPMERegexp photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_regexp;
    TPMERegexp mu3_regexp;      
    TPMERegexp mu5_regexp;          
    TPMERegexp mu8_regexp;      
    TPMERegexp mu12_regexp;     
    TPMERegexp mu15_regexp;     
    TPMERegexp mu17_regexp;     
    TPMERegexp mu20_regexp;     
    TPMERegexp mu24_regexp;     
    TPMERegexp mu30_regexp;     
    TPMERegexp mu15_eta2p1_regexp;     
    TPMERegexp mu24_eta2p1_regexp;     
    TPMERegexp mu30_eta2p1_regexp;     
    TPMERegexp mu8_Jet40_regexp;


    // electron ID MVA
    class ElectronIDMVA* electronIdMVA;
};

#endif // myBabyMaker_h
