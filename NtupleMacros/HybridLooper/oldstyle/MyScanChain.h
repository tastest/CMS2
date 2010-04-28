
#ifndef MYSCANCHAIN_H
#define MYSCANCHAIN_H

// C++ includes
#include <iostream>
#include <vector>
#include<fstream>

// looper includes
#include "EffMulti.h"
#include "Math/VectorUtil.h"


// typedefs
typedef UInt_t      uint32;
typedef ULong64_t   uint64;
typedef uint64  	cuts_t;

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

//
// for det types
//
enum DetType {
    DET_EE = 0,
    DET_EB = 1,
    DET_ALL = 2,
};
static const char det_names[][128] = { "ee", "eb", "all" };

// forward definitions
class TH1F;
class TH2F;
class TChain;
class TDirectory;

class MyScanChain {

    public:

        MyScanChain() {};
        ~MyScanChain() {};

        int ScanChain(bool isData, std::string sampleName, TChain *chain, int nEvents = -1, std::string skimFilePrefix="");

    private:

        // Electrons
        void AnalyseElectrons(const float &weight);

        // Muons
        void AnalyseMuons(const float &weight);

        // get subdetector for histogram filling
        enum DetType getSubdet(int eleIndex);

        // N-1
        // return true if the cuts to apply - the cuts to remove were passed
        // in the cuts that "passed"
        bool CheckCutsNM1(cuts_t apply, cuts_t remove, cuts_t passed);

        // Simple check if the desired cuts to apply are set in the 
        // cuts that "passed"
        bool CheckCuts(cuts_t apply, cuts_t passed);

        // do stuff with histogram
        void FillHist(TH1F** hist, DetType det, const float value, const float weight);
        void Format2DHist(TH2F** hist, std::string name, Int_t nx, Float_t minx, Float_t maxx, Int_t ny, Float_t miny, Float_t maxy);
        void FormatHist(TH1F** hist, std::string name, Int_t n, Float_t min, Float_t max);
        void FormatEffHist(EffMulti** hist, bool lessThan, float thresholdEB, float ThresholdEE, std::string name);

        // sample name
        std::string sampleName_;
        bool isData_;
        // ascifile
        std::ofstream _asciifile; 
        //
        // electrons
        //
        // general  
        TH1F    *h1_ele_pt_[3];
        TH1F    *h1_ele_eta_[3];
        TH1F    *h1_ele_phi_[3];
        TH1F    *h1_ele_tcmet_[3];
        TH1F    *h1_ele_pfmet_[3];
        TH1F    *h1_ele_tcmetdphi_[3];
        TH1F    *h1_ele_pfmetdphi_[3];
        TH1F    *h1_ele_tcmetratio_[3];
        TH1F    *h1_ele_pfmetratio_[3];

	// data-mc comparison for FO
        TH1F    *h1_ele_FO_pt_[3];
        TH1F    *h1_ele_FO_eta_[3];
        TH1F    *h1_ele_FO_iso_[3];

        // inclusive comparisons of bg dists
        TH1F    *h1_ele_incl_r19_[3];
        TH1F    *h1_ele_incl_iso_[3];
        TH1F    *h1_ele_incl_tkIso_[3];
        TH1F    *h1_ele_incl_ecalIso_[3];
        TH1F    *h1_ele_incl_hcalIso_[3];
        TH1F    *h1_ele_incl_pt_[3];
        TH1F    *h1_ele_incl_eta_[3];
        TH1F    *h1_ele_incl_tcmet_[3];
        TH1F    *h1_ele_incl_pfmet_[3];
        TH1F    *h1_ele_incl_pthat_[3];
        TH1F    *h1_ele_incliso_pt_[3];
        TH1F    *h1_ele_incliso_eta_[3];
        TH1F    *h1_ele_incliso_tcmet_[3];
        TH1F    *h1_ele_incliso_pfmet_[3];
        TH1F    *h1_ele_inclnoniso_pt_[3];
        TH1F    *h1_ele_inclnoniso_eta_[3];
        TH1F    *h1_ele_inclnoniso_tcmet_[3];
        TH1F    *h1_ele_inclnoniso_pfmet_[3];

        // nm1
        TH1F    *h1_ele_nm1_pt_[3];
        TH1F 	*h1_ele_nm1_tcmet_pthat_[3];
        TH1F    *h1_ele_nm1_tcmet_pdgidCatagory_[3];
        TH1F 	*h1_ele_nm1_tcmet_[3];
        TH1F    *h1_ele_nm1_pfmet_[3];
        TH1F    *h1_ele_nm1_tcmetdphi_[3];
        TH1F    *h1_ele_nm1_pfmetdphi_[3];
        TH1F    *h1_ele_nm1_tcmetratio_[3];
        TH1F    *h1_ele_nm1_pfmetratio_[3];
        TH1F    *h1_ele_nm1_jetveto_[3];
        TH1F    *h1_ele_nm1_iso_[3];
        TH1F    *h1_ele_nm1_secondpt_[3];
        TH1F    *h1_ele_nm1_r19_[3];

        TH1F    *h1_ele_nm1nor19_tcmet_[3];
        TH1F    *h1_ele_nm1nor19_pfmet_[3];
        TH1F    *h1_ele_nm1nor19_tcmetratio_[3];
        TH1F    *h1_ele_nm1nor19_pfmetratio_[3];

        // after all selections
        TH1F    *h1_ele_selected_pt_[3];
        TH1F    *h1_ele_selected_eta_[3];
        TH1F    *h1_ele_selected_phi_[3];
        TH1F    *h1_ele_selected_tcmet_[3];
        TH1F    *h1_ele_selected_pfmet_[3];
        TH1F    *h1_ele_selected_tcmetdphi_[3];
        TH1F    *h1_ele_selected_pfmetdphi_[3];
        TH1F    *h1_ele_selected_tcmetratio_[3];
        TH1F    *h1_ele_selected_pfmetratio_[3];
        TH1F    *h1_ele_selected_tcmetsignificance_[3];
        TH1F    *h1_ele_selected_pfmetsignificance_[3];
        TH1F    *h1_ele_selected_ptcmet_[3];
        TH1F    *h1_ele_selected_ppfmet_[3];
        TH1F    *h1_ele_selected_pftransmass_[3];
        TH1F    *h1_ele_selected_tctransmass_[3];
        TH1F    *h1_ele_selected_d0corr_[3];
        TH1F    *h1_ele_selected_nmhits_[3];

        TH1F    *h1_ele_selected_sigmaIEtaIEta_[3];
        TH1F    *h1_ele_selected_e2x5MaxOver5x5_[3];
        TH1F    *h1_ele_selected_dEtaIn_[3];
        TH1F    *h1_ele_selected_dPhiIn_[3];
        TH1F    *h1_ele_selected_hOverE_[3];
        TH1F    *h1_ele_selected_fbrem_[3];
        TH1F    *h1_ele_selected_eOverPIn_[3];

        // antiselection on met
        TH1F *h1_ele_antiselected_sigmaIEtaIEta_[3];
        TH1F *h1_ele_antiselected_e2x5MaxOver5x5_[3];
        TH1F *h1_ele_antiselected_dEtaIn_[3];
        TH1F *h1_ele_antiselected_dPhiIn_[3];
        TH1F *h1_ele_antiselected_hOverE_[3];
        TH1F *h1_ele_antiselected_fbrem_[3];
        TH1F *h1_ele_antiselected_eOverPIn_[3];

        // after all selections + electron ID
        TH1F    *h1_ele_selected_cand01_pt_[3];
        TH1F    *h1_ele_selected_cand01_tcmet_[3];
        TH1F    *h1_ele_selected_cand01_pfmet_[3];
        // N-1 plots for selections with electron ID
        TH1F    *h1_ele_nm1_cand01_pt_[3];
        TH1F    *h1_ele_nm1_cand01_tcmet_[3];
        TH1F    *h1_ele_nm1_cand01_pfmet_[3];



        //
        // muons
        //
        // general
        TH1F    *h1_mu_pt_[3];
        TH1F    *h1_mu_eta_[3];
        TH1F    *h1_mu_phi_[3];
        TH1F    *h1_mu_tcmet_[3];
        TH1F    *h1_mu_pfmet_[3];
        TH1F    *h1_mu_tcmetdphi_[3];
        TH1F    *h1_mu_pfmetdphi_[3];
        TH1F    *h1_mu_tcmetratio_[3];
        TH1F    *h1_mu_pfmetratio_[3];

	// data-mc comparison for FO
        TH1F    *h1_mu_FO_pt_[3];
        TH1F    *h1_mu_FO_eta_[3];
        TH1F    *h1_mu_FO_iso_[3];

        // inclusive comparisons of bg dists
        TH1F    *h1_mu_incl_iso_[3];
        TH1F    *h1_mu_incl_pt_[3];
        TH1F    *h1_mu_incl_eta_[3];
        TH1F    *h1_mu_incl_tcmet_[3];
        TH1F    *h1_mu_incl_pfmet_[3];
        TH1F    *h1_mu_incliso_pt_[3];
        TH1F    *h1_mu_incliso_eta_[3];
        TH1F    *h1_mu_incliso_tcmet_[3];
        TH1F    *h1_mu_incliso_pfmet_[3];
        TH1F    *h1_mu_inclnoniso_pt_[3];
        TH1F    *h1_mu_inclnoniso_eta_[3];
        TH1F    *h1_mu_inclnoniso_tcmet_[3];
        TH1F    *h1_mu_inclnoniso_pfmet_[3];

        // nm1
        TH1F    *h1_mu_nm1_pt_[3];
        TH1F    *h1_mu_nm1_secondpt_[3];
        TH1F    *h1_mu_nm1_tcmet_[3];
        TH1F    *h1_mu_nm1_pfmet_[3];
        TH1F    *h1_mu_nm1_iso_[3];

        // after all selections
        TH1F    *h1_mu_selected_pt_[3];
        TH1F    *h1_mu_selected_eta_[3];
        TH1F    *h1_mu_selected_phi_[3];
        TH1F    *h1_mu_selected_tcmet_[3];
        TH1F    *h1_mu_selected_pfmet_[3];
        TH1F    *h1_mu_selected_tcmetdphi_[3];
        TH1F    *h1_mu_selected_pfmetdphi_[3];
        TH1F    *h1_mu_selected_tcmetratio_[3];
        TH1F    *h1_mu_selected_pfmetratio_[3];
        TH1F    *h1_mu_selected_tcmetsignificance_[3];
        TH1F    *h1_mu_selected_pfmetsignificance_[3];
        TH1F    *h1_mu_selected_ptcmet_[3];
        TH1F    *h1_mu_selected_ppfmet_[3];
        TH1F    *h1_mu_selected_pftransmass_[3];
        TH1F    *h1_mu_selected_tctransmass_[3];
        TH1F    *h1_mu_selected_d0corr_[3];

        TH1F    * h1_mu_d0corr_[3];
        TH1F    * h1_mu_nChi2_[3];
        TH1F    * h1_mu_type_[3];
        TH1F    * h1_mu_validHits_[3];
        TH1F    * h1_mu_ecalvetoDep_[3];
        TH1F    * h1_mu_hcalvetoDep_[3];
        TH1F    * h1_mu_validSTAHits_[3];
        TH1F    * h1_mu_muonIsoValue_[3];
        TH1F    * h1_mu_isCosmics_[3];

        TH1F    * h1_mu_selected_nChi2_[3];
        TH1F    * h1_mu_selected_type_[3];
        TH1F    * h1_mu_selected_validHits_[3];
        TH1F    * h1_mu_selected_ecalvetoDep_[3];
        TH1F    * h1_mu_selected_hcalvetoDep_[3];
        TH1F    * h1_mu_selected_validSTAHits_[3];
        TH1F    * h1_mu_selected_muonIsoValue_[3];
        TH1F    * h1_mu_selected_isCosmics_[3];
        TH1F    * h1_mu_selected_caloCompatibility_[3];
        TH1F    * h1_mu_selected_pid_[3];


        TH1F    *h1_mu_antiselected_pt_[3];
        TH1F    *h1_mu_antiselected_eta_[3];
        TH1F    *h1_mu_antiselected_phi_[3];
        TH1F    *h1_mu_antiselected_tcmet_[3];
        TH1F    *h1_mu_antiselected_pfmet_[3];
        TH1F    *h1_mu_antiselected_tcmetdphi_[3];
        TH1F    *h1_mu_antiselected_pfmetdphi_[3];
        TH1F    *h1_mu_antiselected_tcmetratio_[3];
        TH1F    *h1_mu_antiselected_pfmetratio_[3];
        TH1F    *h1_mu_antiselected_tcmetsignificance_[3];
        TH1F    *h1_mu_antiselected_pfmetsignificance_[3];
        TH1F    *h1_mu_antiselected_ptcmet_[3];
        TH1F    *h1_mu_antiselected_ppfmet_[3];
        TH1F    *h1_mu_antiselected_pftransmass_[3];
        TH1F    *h1_mu_antiselected_tctransmass_[3];

        TH1F    * h1_mu_antiselected_d0corr_[3];
        TH1F    * h1_mu_antiselected_nChi2_[3];
        TH1F    * h1_mu_antiselected_type_[3];
        TH1F    * h1_mu_antiselected_validHits_[3];
        TH1F    * h1_mu_antiselected_ecalvetoDep_[3];
        TH1F    * h1_mu_antiselected_hcalvetoDep_[3];
        TH1F    * h1_mu_antiselected_validSTAHits_[3];
        TH1F    * h1_mu_antiselected_muonIsoValue_[3];
        TH1F    * h1_mu_antiselected_isCosmics_[3];
        TH1F    * h1_mu_antiselected_caloCompatibility_[3];
        TH1F    * h1_mu_antiselected_pid_[3];

};

#endif
