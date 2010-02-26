
#ifndef MYSCANCHAIN_H
#define MYSCANCHAIN_H

// C++ includes
#include <iostream>
#include <vector>

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


		void Fill(TH1F** hist, const unsigned int hyp, const float &val, const float &weight);

		void FillAllEleIdHistogramsNoHyp(const float &weight, const TString &sampleName);
		void FillAllEleIdHistogramsHyp(const unsigned int h, const float &weight, const TString &sampleName);
		void FillAllEleIdHistograms(const unsigned int index, const float &weight, const TString &sampleName);

		void FormatHist(TH1F** hist, std::string sampleName, std::string name, int n, float min, float max);
		void FormatAllEleIdHistograms(std::string sampleName);

		// for 2D
		void Fill2D(TH2F** hist, const unsigned int hyp, const float &valx, const float &valy, const float &weight);
		void FormatHist2D(TH2F** hist, std::string sampleName, std::string name, int nx, float minx, float maxx, int ny, float miny, float maxy);

		//
		// ele ID plots


        // NEW
        TH1F *h1_hyp_ltid_sigmaIEtaIEta_[2][4];



		//
		// basic selection quantities
		//

        TH1F *h1_hyp_lt_eb_pt_[4];
        TH1F *h1_hyp_lt_ee_pt_[4];

		TH1F *h1_hyp_lt_eb_hoe_[4];
		TH1F *h1_hyp_lt_eb_sigmaIEtaIEta_[4];
		TH1F *h1_hyp_lt_eb_dEtaIn_[4];
		TH1F *h1_hyp_lt_eb_dPhiIn_[4];
		TH1F *h1_hyp_lt_eb_d0_[4];
		TH1F *h1_hyp_lt_eb_E2x5MaxOver5x5_[4];
		TH1F *h1_hyp_lt_eb_ecalIso_[4];
		TH1F *h1_hyp_lt_eb_hcalIso_[4];
		TH1F *h1_hyp_lt_eb_tkIso_[4];
        TH1F *h1_hyp_lt_eb_relsusy_[4];

		TH1F *h1_hyp_lt_ee_hoe_[4];
		TH1F *h1_hyp_lt_ee_sigmaIEtaIEta_[4];
		TH1F *h1_hyp_lt_ee_dEtaIn_[4];
		TH1F *h1_hyp_lt_ee_dPhiIn_[4];
		TH1F *h1_hyp_lt_ee_d0_[4];
		TH1F *h1_hyp_lt_ee_E2x5MaxOver5x5_[4];
		TH1F *h1_hyp_lt_ee_ecalIso_[4];
		TH1F *h1_hyp_lt_ee_hcalIso_[4];
		TH1F *h1_hyp_lt_ee_tkIso_[4];
		TH1F *h1_hyp_lt_ee_relsusy_[4];

		// N-1
        TH1F *h1_hyp_lt_eb_nm1_hoe_[4];
        TH1F *h1_hyp_lt_eb_nm1_sigmaIEtaIEta_[4];
        TH1F *h1_hyp_lt_eb_nm1_dEtaIn_[4];
        TH1F *h1_hyp_lt_eb_nm1_dPhiIn_[4];
        TH1F *h1_hyp_lt_eb_nm1_d0_[4];
        TH1F *h1_hyp_lt_eb_nm1_E2x5MaxOver5x5_[4];
		TH2F *h1_hyp_lt_eb_nm1_lateral_[4];

        TH1F *h1_hyp_lt_ee_nm1_hoe_[4];
        TH1F *h1_hyp_lt_ee_nm1_sigmaIEtaIEta_[4];
        TH1F *h1_hyp_lt_ee_nm1_dEtaIn_[4];
        TH1F *h1_hyp_lt_ee_nm1_dPhiIn_[4];
        TH1F *h1_hyp_lt_ee_nm1_d0_[4];
        TH1F *h1_hyp_lt_ee_nm1_E2x5MaxOver5x5_[4];
        TH2F *h1_hyp_lt_ee_nm1_lateral_[4];

        TH1F *h1_hyp_lt_eb_afterid_relsusy_[4];
        TH1F *h1_hyp_lt_ee_afterid_relsusy_[4];

        TH1F *h1_hyp_lt_eb_afterid_fbrem_[4];
        TH1F *h1_hyp_lt_ee_afterid_fbrem_[4];

        TH1F *h1_hyp_lt_eb_afterid_eopin_[4];
        TH1F *h1_hyp_lt_ee_afterid_eopin_[4];

        TH1F *h1_hyp_lt_eb_afterid_relsusy_lowfbrem_[4];
        TH1F *h1_hyp_lt_ee_afterid_relsusy_lowfbrem_[4];
		
        TH1F *h1_hyp_lt_eb_afterid_relsusy_highfbrem_[4];
        TH1F *h1_hyp_lt_ee_afterid_relsusy_highfbrem_[4];

        TH1F *h1_hyp_lt_eb_afterid_eopin_lowfbrem_[4];
        TH1F *h1_hyp_lt_ee_afterid_eopin_lowfbrem_[4];
        TH1F *h1_hyp_lt_eb_afterid_eopin_highfbrem_[4];
        TH1F *h1_hyp_lt_ee_afterid_eopin_highfbrem_[4];

		// dPhiIn after id 
        TH1F *h1_hyp_lt_eb_afterid_dPhiIn_lowfbrem_[4];
        TH1F *h1_hyp_lt_ee_afterid_dPhiIn_lowfbrem_[4];
        TH1F *h1_hyp_lt_eb_afterid_dPhiIn_highfbrem_[4];
        TH1F *h1_hyp_lt_ee_afterid_dPhiIn_highfbrem_[4];

        // dEtaIn after id 
        TH1F *h1_hyp_lt_eb_afterid_dEtaIn_lowfbrem_[4];
        TH1F *h1_hyp_lt_ee_afterid_dEtaIn_lowfbrem_[4];
        TH1F *h1_hyp_lt_eb_afterid_dEtaIn_highfbrem_[4];
        TH1F *h1_hyp_lt_ee_afterid_dEtaIn_highfbrem_[4];

		// id
		//
		//
		TH1F *h1_hyp_lt_eb_pt_cand01_[4];
		TH1F *h1_hyp_lt_eb_pt_cand02_[4];
        TH1F *h1_hyp_lt_eb_pt_cand02_extra_[4];
        TH1F *h1_hyp_lt_eb_pt_cand02_extra_v2_[4];
		TH1F *h1_hyp_lt_eb_eta_cand02_extra_v2_[4];

        TH1F *h1_hyp_lt_ee_pt_cand01_[4];
        TH1F *h1_hyp_lt_ee_pt_cand02_[4];
        TH1F *h1_hyp_lt_ee_pt_cand02_extra_[4];
        TH1F *h1_hyp_lt_ee_pt_cand02_extra_v2_[4];
		TH1F *h1_hyp_lt_ee_eta_cand02_extra_v2_[4];
		//
		//

		TH1F *h1_hyp_lt_eb_pt_idnew_[4];
		TH1F *h1_hyp_lt_ee_pt_idnew_[4];
		TH1F *h1_hyp_lt_eb_pt_idold_[4];
		TH1F *h1_hyp_lt_ee_pt_idold_[4];

		// conv
		//
		TH1F *h1_hyp_lt_eb_pt_conv_[4];
		TH1F *h1_hyp_lt_ee_pt_conv_[4];


		// iso
		//
		TH1F *h1_hyp_lt_eb_pt_isonew_cand1_[4];
		TH1F *h1_hyp_lt_ee_pt_isonew_cand1_[4];
		TH1F *h1_hyp_lt_eb_pt_isonew_[4];
		TH1F *h1_hyp_lt_ee_pt_isonew_[4];
		TH1F *h1_hyp_lt_eb_pt_isoold_[4];
		TH1F *h1_hyp_lt_ee_pt_isoold_[4];

		// iso and conv
		//
		TH1F *h1_hyp_lt_eb_pt_id1_iso1_conv_[4];
		TH1F *h1_hyp_lt_ee_pt_id1_iso1_conv_[4];

		//
		// validation plots
		//


};

#endif

