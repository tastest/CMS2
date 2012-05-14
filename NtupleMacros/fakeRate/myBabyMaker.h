#ifndef myBabyMaker_h
#define myBabyMaker_h

// C++ Includes

// ROOT Includes
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TPRegexp.h"
#include "Math/LorentzVector.h"

// lorentz vector of floats 
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

class myBabyMaker {

public:

    myBabyMaker();
    ~myBabyMaker()
    {
        if(babyFile_) delete babyFile_;
    };
    void MakeBabyNtuple (const char *);
    void InitBabyNtuple ();
    void FillBabyNtuple ();
    void CloseBabyNtuple ();
    void SetNumEvents(int nevt) {nEvents_ = nevt;}
    void SetVerbose(bool verbose) {verbose_ = verbose;}
    void ScanChain (TChain *chain, const char *babyFileName, bool isData, int eormu, bool applyFOfilter = true, const std::string& jetcorrPath="../CORE/jetcorr/data/");
    void SetGoodRunList(const char* fileName, bool goodRunIsJson=false);

private:
      
    // BABY NTUPLE VARIABLES
    TFile *babyFile_;
    TTree *babyTree_;
    int nEvents_;
    bool verbose_;

    // good run list
    Bool_t goodrun_is_json;
     
    /////////////////////////// 
    // Event Information     //
    ///////////////////////////

    // Basic Event Information
    Int_t   run_;
    Int_t   ls_;
    UInt_t  evt_;
    Float_t weight_;

    // Pileup - PUSummaryInfoMaker
    Int_t pu_nPUvertices_;
  
    // Pileup - VertexMaker
    Int_t evt_nvtxs_;

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
    LorentzVector lp4_;     // 4-vector of the lepton
    LorentzVector foel_p4_; // 4-vector of the highest additional electron FO in the event
    Int_t foel_id_;
    LorentzVector fomu_p4_; // 4-vector of the highest additional muon FO in the event
    Int_t fomu_id_;
    Float_t foel_mass_;
    Float_t fomu_mass_;

    // pfmet
    Float_t pfmet_;
    Float_t pfmetphi_;
  
    // isolation
    Float_t iso_;                // Relative Isolation
    Float_t iso_nps_;            // Relative Isolation ( 1 GeV pedestal subtraction in ecal barrel )
    Float_t trck_iso_;           // TRK Isolation (not relative)
    Float_t ecal_iso_;           // ECAL Isolation ( not relative)
    Float_t ecal_iso_nps_;       // ECAL Isolation ( not relaive, 1 GeV pedestal subtraction in ecal barrel )
    Float_t hcal_iso_;           // HCAL Isolation ( not relative )
    Float_t pfiso03_;            // PF Isolation with a cone size of 0.3
    Float_t ch_pfiso03_;         // Charged Hadron PF Isolation with a cone size of 0.3
    Float_t nh_pfiso03_;         // Neutral Hadron PF Isolation with a cone size of 0.3
    Float_t em_pfiso03_;         // E&M PF Isolation with a cone size of 0.3
    Float_t pfiso03_bv_;         // PF Isolation with a cone size of 0.3 (barrel veto)
    Float_t ch_pfiso03_bv_;      // Charged Hadron PF Isolation with a cone size of 0.3 (barrel veto)
    Float_t nh_pfiso03_bv_;      // Neutral Hadron PF Isolation with a cone size of 0.3 (barrel veto)
    Float_t em_pfiso03_bv_;      // E&M PF Isolation with a cone size of 0.3 (barrel veto)
    Float_t pfiso04_;            // PF Isolation with a cone size of 0.4
    Float_t ch_pfiso04_;         // Charged Hadron PF Isolation with a cone size of 0.4
    Float_t nh_pfiso04_;         // Neutral Hadron  PF Isolation with a cone size of 0.4
    Float_t em_pfiso04_;         // E&M PF Isolation with a cone size of 0.4
    Float_t pfiso04_bv_;         // PF Isolation with a cone size of 0.4 (barrel veto)
    Float_t ch_pfiso04_bv_;      // Charged Hadron PF Isolation with a cone size of 0.4 (barrel veto)
    Float_t nh_pfiso04_bv_;      // Neutral Hadron  PF Isolation with a cone size of 0.4 (barrel veto)
    Float_t em_pfiso04_bv_;      // E&M PF Isolation with a cone size of 0.4 (barrel veto)
    Float_t radiso_et1p0_;       // Radial Isolation with a cone size of 0.3, neutral ET threashold of 1.0
    Float_t ch_radiso_et1p0_;    // Charged Hadron Radial Isolation with a cone size of 0.3, neutral ET threashold of 1.0
    Float_t nh_radiso_et1p0_;    // Neutral Hadron  Radial Isolation with a cone size of 0.3, neutral ET threashold of 1.0
    Float_t em_radiso_et1p0_;    // E&M Radial Isolation with a cone size of 0.3, neutral ET threashold of 1.0
    Float_t radiso_et1p0_bv_;    // Radial Isolation with a cone size of 0.3, neutral ET threashold of 1.0 (barrel veto)
    Float_t ch_radiso_et1p0_bv_; // Charged Hadron Radial Isolation with a cone size of 0.3, neutral ET threashold of 1.0 (barrel veto)
    Float_t nh_radiso_et1p0_bv_; // Neutral Hadron  Radial Isolation with a cone size of 0.3, neutral ET threashold of 1.0 (barrel veto)
    Float_t em_radiso_et1p0_bv_; // E&M Radial Isolation with a cone size of 0.3, neutral ET threashold of 1.0 (barrel veto)
    Float_t radiso_et0p5_;       // Radial Isolation with a cone size of 0.3, neutral ET threashold of 0.5
    Float_t ch_radiso_et0p5_;    // Charged Hadron Radial Isolation with a cone size of 0.3, neutral ET threashold of 0.5
    Float_t nh_radiso_et0p5_;    // Neutral Hadron  Radial Isolation with a cone size of 0.3, neutral ET threashold of 0.5
    Float_t em_radiso_et0p5_;    // E&M Radial Isolation with a cone size of 0.3, neutral ET threashold of 0.5
    Float_t radiso_et0p5_bv_;    // Radial Isolation with a cone size of 0.3, neutral ET threashold of 0.5 (barrel veto)
    Float_t ch_radiso_et0p5_bv_; // Charged Hadron Radial Isolation with a cone size of 0.3, neutral ET threashold of 0.5 (barrel veto)
    Float_t nh_radiso_et0p5_bv_; // Neutral Hadron  Radial Isolation with a cone size of 0.3, neutral ET threashold of 0.5 (barrel veto)
    Float_t em_radiso_et0p5_bv_; // E&M Radial Isolation with a cone size of 0.3, neutral ET threashold of 0.5 (barrel veto)
    Float_t pfpupt03_;           // PF Pile Up sum pT (cone 0.3) (a.k.a. DeltaBeta)
    Float_t pfpupt04_;           // PF Pile Up sum pT (cone 0.4) (a.k.a. DeltaBeta)

    // Corrected Isolation
    Float_t cpfiso03_rho_;       // PF Isolation with a cone size of 0.3 (corrected using rho*area_eff -- only filled for electrons)
    Float_t cpfiso03_db_;        // PF Isolation with a cone size of 0.3 (corrected using #DeltaBeta -- only filled for muons)

    // PV
    Float_t d0PV_wwV1_;       // electron_d0PV_wwV1(iEl)
    Float_t dzPV_wwV1_;       // electron_dzPV_wwV1(iEl)

    // Id
    Bool_t  closestMuon_;      // true if els_closestMuon().at(index) == -1
	Float_t el_id_sieie_;
	Float_t el_id_detain_;
	Float_t el_id_dphiin_;
    Bool_t  el_id_smurfV5_;
    Bool_t  el_id_vbtf80_;
    Bool_t  el_id_vbtf90_;

    // effective area
    Float_t el_effarea_;            // 2012 working point effective area (cone = 0.3)
    Float_t mu_effarea03_;          // 2012 working point effective area (combined       , cone = 0.3 , loose)
    Float_t mu_nh_effarea03_;       // 2012 working point effective area (neutral hadron , cone = 0.3 , loose)
    Float_t mu_em_effarea03_;       // 2012 working point effective area (E&M            , cone = 0.3 , loose)
    Float_t mu_effarea03_tight_;    // 2012 working point effective area (combined       , cone = 0.3 , tight)
    Float_t mu_nh_effarea03_tight_; // 2012 working point effective area (neutral hadron , cone = 0.3 , tight)
    Float_t mu_em_effarea03_tight_; // 2012 working point effective area (E&M            , cone = 0.3 , tight)
    Float_t mu_effarea04_;          // 2012 working point effective area (combined       , cone = 0.4 , loose)
    Float_t mu_nh_effarea04_;       // 2012 working point effective area (neutral hadron , cone = 0.4 , loose)
    Float_t mu_em_effarea04_;       // 2012 working point effective area (E&M            , cone = 0.4 , loose)
    Float_t mu_effarea04_tight_;    // 2012 working point effective area (combined       , cone = 0.4 , tight)
    Float_t mu_nh_effarea04_tight_; // 2012 working point effective area (neutral hadron , cone = 0.4 , tight)
    Float_t mu_em_effarea04_tight_; // 2012 working point effective area (E&M            , cone = 0.4 , tight)

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
    float ht_pf_;            
    float ht_pf_L2L3_;       
    float ht_pf_L1FastL2L3_;  

    // MC truth information
    Int_t mcid_;                // els_mc_id or mus_mc_id
    Int_t mcmotherid_;          // els_mc_motherid or mus_mc_motherid
    Int_t mc3id_;
    Float_t mc3pt_;
    Float_t mc3dr_;
	LorentzVector mc3p4_;
    Int_t leptonIsFromW_;

    //////////////////////////
    // End Lepton Variables //
    //////////////////////////

    //////////////////////////////////////////////////////
    // Fake Rate Numerator & Denominator Selections     //
    //////////////////////////////////////////////////////

    //////////
    // 2012 //
    //////////
  
    // SS
      
    // Electrons
    Bool_t num_el_ssV7_;
    Bool_t num_el_ssV7_noIso_;
    Bool_t v1_el_ssV7_;
    Bool_t v2_el_ssV7_;
    Bool_t v3_el_ssV7_;
    
    // Muons
    Bool_t num_mu_ssV5_;        // NominalSSv5
    Bool_t num_mu_ssV5_noIso_;  // NominalSSv5 with no isolation applied
    Bool_t fo_mu_ssV5_;         // muonSelectionFO_mu_ssV5
    Bool_t fo_mu_ssV5_noIso_;   // muonSleectionFO_mu_ssV5 with no isolation applied

    //////////
    // 2011 //
    //////////
  
    // SS
      
    // Electrons
    Bool_t num_el_ssV6_;
    Bool_t v1_el_ssV6_;
    Bool_t v2_el_ssV6_;
    Bool_t v3_el_ssV6_;
    Bool_t num_el_ssV6_noIso_;
    
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

    //////////////////////////////////////////////////////
    // End Fake Rate Numerator & Denominator Selections //
    //////////////////////////////////////////////////////

    ///////////////////////  
    // 2012 Triggers     //
    ///////////////////////

    // Electrons
    Int_t ele8_CaloIdL_CaloIsoVL_vstar_;                                      // HLT_Ele8_CaloIdL_CaloIsoVL_v*
    Int_t ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar_;                     // HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*
    Int_t ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar_;               // HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v*
    Int_t ele8_CaloIdT_TrkIdVL_vstar_;                                        // HLT_Ele8_CaloIdT_TrkIdVL_v*
    Int_t ele8_CaloIdT_TrkIdVL_Jet30_vstar_;                                  // HLT_Ele8_CaloIdT_TrkIdVL_Jet30_v*
    Int_t ele17_CaloIdL_CaloIsoVL_vstar_;                                     // HLT_Ele17_CaloIdL_CaloIsoVL_v*
    Int_t ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar_;                    // HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*
    Int_t ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar_;              // HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v*
    Int_t ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFJet30_vstar_;       // HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFJet30_v*
    Int_t ele27_WP80_vstar_;                                                  // HLT_Ele27_WP80_v*

    Int_t ele8_CaloIdL_CaloIsoVL_version_;                                    // HLT_Ele8_CaloIdL_CaloIsoVL_v*
    Int_t ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version_;                   // HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*
    Int_t ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version_;             // HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v*
    Int_t ele8_CaloIdT_TrkIdVL_version_;                                      // HLT_Ele8_CaloIdT_TrkIdVL_v*
    Int_t ele8_CaloIdT_TrkIdVL_Jet30_version_;                                // HLT_Ele8_CaloIdT_TrkIdVL_Jet30_v*
    Int_t ele17_CaloIdL_CaloIsoVL_version_;                                   // HLT_Ele17_CaloIdL_CaloIsoVL_v*
    Int_t ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version_;                  // HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*
    Int_t ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version_;            // HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v*
    Int_t ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFJet30_version_;     // HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFJet30_v*
    Int_t ele27_WP80_version_;                                                // HLT_Ele27_WP80_v*

    Float_t dr_ele8_CaloIdL_CaloIsoVL_vstar_;                                 // HLT_Ele8_CaloIdL_CaloIsoVL_v*
    Float_t dr_ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar_;                // HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*
    Float_t dr_ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar_;          // HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v*
    Float_t dr_ele8_CaloIdT_TrkIdVL_vstar_;                                   // HLT_Ele8_CaloIdT_TrkIdVL_v*
    Float_t dr_ele8_CaloIdT_TrkIdVL_Jet30_vstar_;                             // HLT_Ele8_CaloIdT_TrkIdVL_Jet30_v*
    Float_t dr_ele17_CaloIdL_CaloIsoVL_vstar_;                                // HLT_Ele17_CaloIdL_CaloIsoVL_v*
    Float_t dr_ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar_;               // HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*
    Float_t dr_ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar_;         // HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v*
    Float_t dr_ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFJet30_vstar_;  // HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFJet30_v*
    Float_t dr_ele27_WP80_vstar_;                                             // HLT_Ele27_WP80_v*

    Int_t hltps_ele8_CaloIdL_CaloIsoVL_vstar_;                                // HLT_Ele8_CaloIdL_CaloIsoVL_v*
    Int_t hltps_ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar_;               // HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*
    Int_t hltps_ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar_;         // HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v*
    Int_t hltps_ele8_CaloIdT_TrkIdVL_vstar_;                                  // HLT_Ele8_CaloIdT_TrkIdVL_v*
    Int_t hltps_ele8_CaloIdT_TrkIdVL_Jet30_vstar_;                             // HLT_Ele8_CaloIdT_TrkIdVL_Jet30_v*
    Int_t hltps_ele17_CaloIdL_CaloIsoVL_vstar_;                               // HLT_Ele17_CaloIdL_CaloIsoVL_v*
    Int_t hltps_ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar_;              // HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*
    Int_t hltps_ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar_;        // HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v*
    Int_t hltps_ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFJet30_vstar_; // HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFJet30_v*
    Int_t hltps_ele27_WP80_vstar_;                                            // HLT_Ele27_WP80_v*

    Int_t l1ps_ele8_CaloIdL_CaloIsoVL_vstar_;                                 // HLT_Ele8_CaloIdL_CaloIsoVL_v*
    Int_t l1ps_ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar_;                // HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*
    Int_t l1ps_ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar_;          // HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v*
    Int_t l1ps_ele8_CaloIdT_TrkIdVL_vstar_;                                   // HLT_Ele8_CaloIdT_TrkIdVL_v*
    Int_t l1ps_ele8_CaloIdT_TrkIdVL_Jet30_vstar_;                             // HLT_Ele8_CaloIdT_TrkIdVL_Jet30_v*
    Int_t l1ps_ele17_CaloIdL_CaloIsoVL_vstar_;                                // HLT_Ele17_CaloIdL_CaloIsoVL_v*
    Int_t l1ps_ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar_;               // HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*
    Int_t l1ps_ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar_;         // HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v*
    Int_t l1ps_ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFJet30_vstar_;  // HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFJet30_v*
    Int_t l1ps_ele27_WP80_vstar_;                                             // HLT_Ele27_WP80_v*

    // Muons
    Int_t mu5_vstar_;                  // HLT_Mu5_v*            // also in 2011
    Int_t mu8_vstar_;                  // HLT_Mu8_v*            // also in 2011
    Int_t mu12_vstar_;                 // HLT_Mu12_v*           // also in 2011
    Int_t mu17_vstar_;                 // HLT_Mu17_v*
    Int_t mu15_eta2p1_vstar_;          // HLT_Mu15_eta2p1_v*
    Int_t mu24_eta2p1_vstar_;          // HLT_Mu24_eta2p1_v*
    Int_t mu30_eta2p1_vstar_;          // HLT_Mu30_eta2p1_v*
    Int_t isoMu20_eta2p1_vstar_;       // HLT_IsoMu20_eta2p1_v*
    Int_t isoMu24_eta2p1_vstar_;       // HLT_IsoMu24_eta2p1_v*
    Int_t isoMu30_eta2p1_vstar_;       // HLT_IsoMu30_eta2p1_v*
    Int_t relIso1p0Mu17_vstar_;        // HLT_RelIso1p0Mu17_v*
    Int_t relIso1p0Mu5_vstar_;         // HLT_RelIso1p0Mu5_v*

    Int_t mu5_version_;                // HLT_Mu5_v*            // also in 2011
    Int_t mu8_version_;                // HLT_Mu8_v*            // also in 2011
    Int_t mu12_version_;               // HLT_Mu12_v*           // also in 2011
    Int_t mu17_version_;               // HLT_Mu17_v*
    Int_t mu15_eta2p1_version_;        // HLT_Mu15_eta2p1_v*
    Int_t mu24_eta2p1_version_;        // HLT_Mu24_eta2p1_v*
    Int_t mu30_eta2p1_version_;        // HLT_Mu30_eta2p1_v*
    Int_t isoMu20_eta2p1_version_;     // HLT_IsoMu20_eta2p1_v*
    Int_t isoMu24_eta2p1_version_;     // HLT_IsoMu24_eta2p1_v*
    Int_t isoMu30_eta2p1_version_;     // HLT_IsoMu30_eta2p1_v*
    Int_t relIso1p0Mu17_version_;      // HLT_RelIso1p0Mu17_v*
    Int_t relIso1p0Mu5_version_;       // HLT_RelIso1p0Mu5_v*

    Float_t dr_mu5_vstar_;             // HLT_Mu5_v*            // also in 2011
    Float_t dr_mu8_vstar_;             // HLT_Mu8_v*            // also in 2011
    Float_t dr_mu12_vstar_;            // HLT_Mu12_v*           // also in 2011
    Float_t dr_mu17_vstar_;            // HLT_Mu17_v*
    Float_t dr_mu15_eta2p1_vstar_;     // HLT_Mu15_eta2p1_v*
    Float_t dr_mu24_eta2p1_vstar_;     // HLT_Mu24_eta2p1_v*
    Float_t dr_mu30_eta2p1_vstar_;     // HLT_Mu30_eta2p1_v*
    Float_t dr_isoMu20_eta2p1_vstar_;  // HLT_IsoMu20_eta2p1_v*
    Float_t dr_isoMu24_eta2p1_vstar_;  // HLT_IsoMu24_eta2p1_v*
    Float_t dr_isoMu30_eta2p1_vstar_;  // HLT_IsoMu30_eta2p1_v*
    Float_t dr_relIso1p0Mu17_vstar_;   // HLT_RelIso1p0Mu17_v*
    Float_t dr_relIso1p0Mu5_vstar_;    // HLT_RelIso1p0Mu5_v*

    Int_t hltps_mu5_vstar_;            // HLT_Mu5_v*            // also in 2011
    Int_t hltps_mu8_vstar_;            // HLT_Mu8_v*            // also in 2011
    Int_t hltps_mu12_vstar_;           // HLT_Mu12_v*           // also in 2011
    Int_t hltps_mu17_vstar_;           // HLT_Mu17_v*
    Int_t hltps_mu15_eta2p1_vstar_;    // HLT_Mu15_eta2p1_v*
    Int_t hltps_mu24_eta2p1_vstar_;    // HLT_Mu24_eta2p1_v*
    Int_t hltps_mu30_eta2p1_vstar_;    // HLT_Mu30_eta2p1_v*
    Int_t hltps_isoMu20_eta2p1_vstar_; // HLT_IsoMu20_eta2p1_v*
    Int_t hltps_isoMu24_eta2p1_vstar_; // HLT_IsoMu24_eta2p1_v*
    Int_t hltps_isoMu30_eta2p1_vstar_; // HLT_IsoMu30_eta2p1_v*
    Int_t hltps_relIso1p0Mu17_vstar_;  // HLT_RelIso1p0Mu17_v*
    Int_t hltps_relIso1p0Mu5_vstar_;   // HLT_RelIso1p0Mu5_v*

    Int_t l1ps_mu5_vstar_;             // HLT_Mu5_v*            // also in 2011
    Int_t l1ps_mu8_vstar_;             // HLT_Mu8_v*            // also in 2011
    Int_t l1ps_mu12_vstar_;            // HLT_Mu12_v*           // also in 2011
    Int_t l1ps_mu17_vstar_;            // HLT_Mu17_v*
    Int_t l1ps_mu15_eta2p1_vstar_;     // HLT_Mu15_eta2p1_v*
    Int_t l1ps_mu24_eta2p1_vstar_;     // HLT_Mu24_eta2p1_v*
    Int_t l1ps_mu30_eta2p1_vstar_;     // HLT_Mu30_eta2p1_v*
    Int_t l1ps_isoMu20_eta2p1_vstar_;  // HLT_IsoMu20_eta2p1_v*
    Int_t l1ps_isoMu24_eta2p1_vstar_;  // HLT_IsoMu24_eta2p1_v*
    Int_t l1ps_isoMu30_eta2p1_vstar_;  // HLT_IsoMu30_eta2p1_v*
    Int_t l1ps_relIso1p0Mu17_vstar_;   // HLT_RelIso1p0Mu17_v*
    Int_t l1ps_relIso1p0Mu5_vstar_;    // HLT_RelIso1p0Mu5_v*

    ///////////////////////  
    // End 2012 Triggers //
    ///////////////////////


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
    //  4. HLT prescale of the trigger object whould be stored for each trigger by prepending "hltps_" to the trigger variable name
    //  4. L1 prescale of the trigger object whould be stored for each trigger by prepending "l1ps_" to the trigger variable name

    // Electrons
    Int_t ele8_vstar_;                                               // HLT_Ele8_v*
    Int_t ele8_CaloIdL_TrkIdVL_vstar_;                               // HLT_Ele8_CaloIdL_TrkIdVL_v*
    Int_t ele8_CaloIdL_CaloIsoVL_Jet40_vstar_;                       // HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v*
    Int_t ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_vstar_;            // HLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v*
    Int_t photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_vstar_;      // HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v*

    Int_t ele8_version_;                                             // HLT_Ele8_v*
    Int_t ele8_CaloIdL_TrkIdVL_version_;                             // HLT_Ele8_CaloIdL_TrkIdVL_v*
    Int_t ele8_CaloIdL_CaloIsoVL_Jet40_version_;                     // HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v*
    Int_t ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_version_;          // HLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v*
    Int_t photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_version_;    // HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v*

    Float_t dr_ele8_vstar_;                                          // HLT_Ele8_v*
    Float_t dr_ele8_CaloIdL_TrkIdVL_vstar_;                          // HLT_Ele8_CaloIdL_TrkIdVL_v*
    Float_t dr_ele8_CaloIdL_CaloIsoVL_Jet40_vstar_;                  // HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v*
    Float_t dr_ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_vstar_;       // HLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v*
    Float_t dr_photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_vstar_; // HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v*

    Int_t hltps_ele8_vstar_;                                            // HLT_Ele8_v*
    Int_t hltps_ele8_CaloIdL_TrkIdVL_vstar_;                            // HLT_Ele8_CaloIdL_TrkIdVL_v*
    Int_t hltps_ele8_CaloIdL_CaloIsoVL_Jet40_vstar_;                    // HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v*
    Int_t hltps_ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_vstar_;         // HLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v*
    Int_t hltps_photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_vstar_;   // HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v*

    // Muons
    Int_t mu3_vstar_;            // HLT_Mu3_v*
    Int_t mu15_vstar_;           // HLT_Mu15_v*
    Int_t mu20_vstar_;           // HLT_Mu20_v*
    Int_t mu24_vstar_;           // HLT_Mu24_v*
    Int_t mu30_vstar_;           // HLT_Mu30_v*
    Int_t mu8_Jet40_vstar_;      // HLT_Mu8_Jet40_v*

    Int_t mu3_version_;          // HLT_Mu3_v*
    Int_t mu15_version_;         // HLT_Mu15_v*
    Int_t mu20_version_;         // HLT_Mu20_v*
    Int_t mu24_version_;         // HLT_Mu24_v*
    Int_t mu30_version_;         // HLT_Mu30_v*
    Int_t mu8_Jet40_version_;    // HLT_Mu8_Jet40_v*

    Float_t dr_mu3_vstar_;       // HLT_Mu5_v*
    Float_t dr_mu15_vstar_;      // HLT_Mu15_v*
    Float_t dr_mu20_vstar_;      // HLT_Mu20_v*
    Float_t dr_mu24_vstar_;      // HLT_Mu24_v*
    Float_t dr_mu30_vstar_;      // HLT_Mu30_v*
    Float_t dr_mu8_Jet40_vstar_; // HLT_Mu8_Jet40_v*

    Int_t hltps_mu3_vstar_;         // HLT_Mu5_v*
    Int_t hltps_mu15_vstar_;        // HLT_Mu15_v*
    Int_t hltps_mu20_vstar_;        // HLT_Mu20_v*
    Int_t hltps_mu24_vstar_;        // HLT_Mu24_v*
    Int_t hltps_mu30_vstar_;        // HLT_Mu30_v*
    Int_t hltps_mu8_Jet40_vstar_;   // HLT_Mu8_Jet40_v*

    ///////////////////////  
    // End 2011 Triggers //
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

    // regular expressions for triggers
    TPMERegexp ele8_regexp;
    TPMERegexp ele8_CaloIdL_TrkIdVL_regexp;
    TPMERegexp ele8_CaloIdL_CaloIsoVL_regexp;
    TPMERegexp ele8_CaloIdL_CaloIsoVL_Jet40_regexp;
    TPMERegexp ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_regexp;
    TPMERegexp ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_regexp;
    TPMERegexp ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_regexp;
    TPMERegexp ele8_CaloIdT_TrkIdVL_regexp;
    TPMERegexp ele8_CaloIdT_TrkIdVL_Jet30_regexp;
    TPMERegexp ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_vstar_regexp;
    TPMERegexp ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar_regexp;
    TPMERegexp ele17_CaloIdL_CaloIsoVL_regexp;
    TPMERegexp ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_regexp;
    TPMERegexp ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_regexp;
    TPMERegexp ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFJet30_rexexp;
    TPMERegexp ele27_WP80_rexexp;
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
    TPMERegexp isoMu20_eta2p1_regexp;
    TPMERegexp isoMu24_eta2p1_regexp;
    TPMERegexp isoMu30_eta2p1_regexp;
    TPMERegexp relIso1p0Mu17_regexp;
    TPMERegexp relIso1p0Mu5_regexp;


    // electron ID MVA
    class ElectronIDMVA* electronIdMVA;
};

#endif // myBabyMaker_h
