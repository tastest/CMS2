
#ifndef MYSCANCHAIN_H
#define MYSCANCHAIN_H

// C++ includes
#include <iostream>
#include <vector>

class TH1F;
class TH2F;
class TChain;
class TDirectory;
class TProfile;

class MyScanChain {

	public:

		MyScanChain(unsigned int electronId) : electronId_(electronId) {};
		~MyScanChain() {};

		int ScanChain(bool isData, std::string sampleName, TChain *chain, int nEvents = -1, std::string skimFilePrefix="");

	private:


		void Fill(TH1F** hist, const unsigned int hyp, const float &val, const float &weight);
		void FillAllDYEstHistograms(const unsigned int h, const float &weight, const unsigned int njet);

		void FormatHist(TH1F** hist, std::string sampleName, std::string name, int n, float min, float max);
		void FormatAllDYEstHistograms(std::string sampleName);
        void FormatAllAnaHistograms(std::string sampleName);

		// for 2D
		void Fill2D(TH2F** hist, const unsigned int hyp, const float &valx, const float &valy, const float &weight);
		void FormatHist2D(TH2F** hist, std::string sampleName, std::string name, int nx, float minx, float maxx, int ny, float miny, float maxy);

        //
        // dealing with selections
        //
        unsigned int leptonSelect(const int id, const unsigned int lepIdx);
        unsigned int electronId_;
        float leptonIsolation(const unsigned int index, const int id);


        // early data studies
        void TopNotTopDown(const float &weight);
        unsigned int leptonSelectTopNotTopDown(const int id, const unsigned int lepIdx);

        // debugging selection functions
        void Debug(const float &weight);

        // spike check
        void SpikeCheck(const float &weight);
        TH2F *h2_spike_scattermet_;
        TH2F *h2_spike_scattermet_fixed_;
        TH2F *h2_spike_scatteret_;
        TH2F *h2_spike_scatteret_eid_;

        // checking the deta
        void SuperClusterCheck(const float &weight); 
        TProfile *p1_detain_eeplus_;
        TProfile *p1_detain_corrected_eeplus_;
        TProfile *p1_detain_correctedfit_eeplus_;

        TProfile *p1_detain_eeminus_;
        TProfile *p1_detain_corrected_eeminus_;
        TProfile *p1_detain_correctedfit_eeminus_;

        TH1F *h1_dx_eeminus_;
        TH1F *h1_dy_eeminus_;
        TH1F *h1_dx_eeplus_;
        TH1F *h1_dy_eeplus_;




		//
		// ttbar analysis plots
		//

		TH1F *h1_hyp_njets_[4];

        // histograms after basic selection
        TH1F *h1_hyp_tcmet_[4][4];
        TH1F *h1_hyp_lt_pt_[4][4];
        TH1F *h1_hyp_ll_pt_[4][4];
        TH1F *h1_hyp_lt_d0corr_[4][4];
        TH1F *h1_hyp_ll_d0corr_[4][4];

        TH1F *h1_hyp_allmuon_d0corr_[4][4];
        TH1F *h1_hyp_glbtrkmuon_d0corr_[4][4];

        // histograms after basic selection + preselection
        // e.g. d0 < 0.04 (electrons and muons)
        // |eta| < 2.5 (electrons and muons)
        // ecal driven (electrons)
        // no close muon (electrons)
        TH1F *h1_hyp_presel_njets_[4];
        TH1F *h1_hyp_presel_mll_[4][4];
        TH1F *h1_hyp_presel_pt_[4][4];
        TH1F *h1_hyp_presel_tcmet_[4][4];
        TH1F *h1_hyp_presel_lt_iso_[4][4];
        TH1F *h1_hyp_presel_ll_iso_[4][4];
        TH1F *h1_hyp_presel_lt_pt_[4][4];
        TH1F *h1_hyp_presel_ll_pt_[4][4];
        TH1F *h1_hyp_presel_lt_d0corr_[4][4];
        TH1F *h1_hyp_presel_ll_d0corr_[4][4];

        // histograms after basic selection + preselection
        // + loose iso
        // e.g. d0 < 0.04 (electrons and muons)
        // |eta| < 2.5 (electrons and muons)
        // ecal driven (electrons)
        // no close muon (electrons)
        // highest pt leg has iso < 0.4
        TH1F *h1_hyp_presel_looseiso_njets_[4];
        TH1F *h1_hyp_presel_looseiso_mll_[4][4];
        TH1F *h1_hyp_presel_looseiso_pt_[4][4];
        TH1F *h1_hyp_presel_looseiso_tcmet_[4][4];

        TH1F *h1_hyp_presel_looseiso_lt_passEle10_[4];
        TH1F *h1_hyp_presel_looseiso_ll_pt_reducedBias_[4];
        TH1F *h1_hyp_presel_looseiso_ll_pt_reducedBias_passEle10_[4];
        TH1F *h1_hyp_presel_looseiso_ll_eta_reducedBias_[4];
        TH1F *h1_hyp_presel_looseiso_ll_eta_reducedBias_passEle10_[4];


		//
		// DY estimate plots
		//

		TH1F *h1_dyest_mll_met_[4][4];
        TH1F *h1_dyest_mll_nomet_[4][4];
		TH1F *h1_dyest_met_in_[4][4];
		TH1F *h1_dyest_met_out_[4][4];

};

#endif

