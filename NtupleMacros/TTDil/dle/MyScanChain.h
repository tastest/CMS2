
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

		//
		// ttbar analysis plots
		//

		TH1F *h1_hyp_njets_[4];

        TH1F *h1_hyp_mll_[4][4];
        TH1F *h1_hyp_pt_[4][4];
        TH1F *h1_hyp_pfmet_[4][4];
        TH1F *h1_hyp_tcmet_[4][4];


		//
		// DY estimate plots
		//

		TH1F *h1_dyest_mll_met_[4][4];
        TH1F *h1_dyest_mll_nomet_[4][4];
		TH1F *h1_dyest_met_in_[4][4];
		TH1F *h1_dyest_met_out_[4][4];

};

#endif

