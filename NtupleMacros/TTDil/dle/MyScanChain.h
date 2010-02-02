
#ifndef MYSCANCHAIN_H
#define MYSCANCHAIN_H

// C++ includes
#include <iostream>
#include <vector>



class TH1F;
class TChain;
class TDirectory;

class MyScanChain {

	public:

		MyScanChain() {};
		~MyScanChain() {};

		int ScanChain(bool isData, std::string sampleName, TChain *chain, int nEvents = -1, std::string skimFilePrefix="");

	private:


		void Fill(TH1F** hist, unsigned int hyp, float val, float weight);
		void FillAllEleIdHistograms(const unsigned int h, const float &weight);

		void FormatHist(TH1F** hist, std::string sampleName, std::string name, int n, float min, float max);
		void FormatAllEleIdHistograms(std::string sampleName);
		void FormatAllAnaHistograms(std::string sampleName);

		//
		// ttbar analysis plots
		//

		TH1F *h1_hyp_njets_[4];

		//
		// ele ID plots
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

		TH1F *h1_hyp_lt_ee_hoe_[4];
		TH1F *h1_hyp_lt_ee_sigmaIEtaIEta_[4];
		TH1F *h1_hyp_lt_ee_dEtaIn_[4];
		TH1F *h1_hyp_lt_ee_dPhiIn_[4];
		TH1F *h1_hyp_lt_ee_d0_[4];
		TH1F *h1_hyp_lt_ee_E2x5MaxOver5x5_[4];
		TH1F *h1_hyp_lt_ee_ecalIso_[4];
		TH1F *h1_hyp_lt_ee_hcalIso_[4];
		TH1F *h1_hyp_lt_ee_tkIso_[4];

		// id
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


};

#endif

