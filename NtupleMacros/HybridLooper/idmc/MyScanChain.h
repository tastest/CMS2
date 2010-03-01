
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

        // id
        TH1F *h1_hyp_debug_after_idcand01_pt_[2][4];
        TH1F *h1_hyp_debug_after_idcand02_pt_[2][4];

        // iso
        TH1F *h1_hyp_debug_after_isocand01_pt_[2][4];
        TH1F *h1_hyp_debug_after_isocand02_pt_[2][4];

        // conv
        TH1F *h1_hyp_debug_after_convcand01_pt_[2][4];
        TH1F *h1_hyp_debug_after_convcand02_pt_[2][4];

        // iso and id
        TH1F *h1_hyp_debug_after_idisocand01_pt_[2][4];
        TH1F *h1_hyp_debug_after_idisocand02_pt_[2][4];

        // iso and id and conv
        TH1F *h1_hyp_debug_after_idisoconvcand01_pt_[2][4];
        TH1F *h1_hyp_debug_after_idisoconvcand02_pt_[2][4];


        // pdg id
        TH1F *h1_hyp_debug_pdgid_[2][4];
        TH1F *h1_hyp_debug_after_cand01_pdgid_[2][4];
        TH1F *h1_hyp_debug_after_cand02_pdgid_[2][4];

        // pt
        TH1F *h1_hyp_debug_pt_[2][4];
        TH1F *h1_hyp_debug_after_cand01_pt_[2][4];
        TH1F *h1_hyp_debug_after_cand02_pt_[2][4];

        // eta
        TH1F *h1_hyp_debug_eta_[2][4];
        TH1F *h1_hyp_debug_after_cand01_eta_[2][4];
        TH1F *h1_hyp_debug_after_cand02_eta_[2][4];

        // sigmaIEtaIEta
        TH1F *h1_hyp_debug_sigmaIEtaIEta_[2][4];
        TH1F *h1_hyp_debug_after_cand01_sigmaIEtaIEta_[2][4];
        TH1F *h1_hyp_debug_after_cand02_sigmaIEtaIEta_[2][4];

        // hoe
        TH1F *h1_hyp_debug_hoe_[2][4];
        TH1F *h1_hyp_debug_after_cand01_hoe_[2][4];
        TH1F *h1_hyp_debug_after_cand02_hoe_[2][4];

        // dPhiIn
        TH1F *h1_hyp_debug_dPhiIn_[2][4];
        TH1F *h1_hyp_debug_after_cand01_dPhiIn_[2][4];
        TH1F *h1_hyp_debug_after_cand02_dPhiIn_[2][4];

        // dEtaIn
        TH1F *h1_hyp_debug_dEtaIn_[2][4];
        TH1F *h1_hyp_debug_after_cand01_dEtaIn_[2][4];
        TH1F *h1_hyp_debug_after_cand02_dEtaIn_[2][4];

        // d0
        TH1F *h1_hyp_debug_d0_[2][4];
        TH1F *h1_hyp_debug_after_cand01_d0_[2][4];
        TH1F *h1_hyp_debug_after_cand02_d0_[2][4];

        // E2x5MaxOver5x5
        TH1F *h1_hyp_debug_E2x5MaxOver5x5_[2][4];
        TH1F *h1_hyp_debug_after_cand01_E2x5MaxOver5x5_[2][4];
        TH1F *h1_hyp_debug_after_cand02_E2x5MaxOver5x5_[2][4];

        // reliso
        TH1F *h1_hyp_debug_reliso_[2][4];
        TH1F *h1_hyp_debug_after_cand01_reliso_[2][4];
        TH1F *h1_hyp_debug_after_cand02_reliso_[2][4];


};

#endif

