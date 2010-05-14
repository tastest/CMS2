
#ifndef MYSCANCHAIN_H
#define MYSCANCHAIN_H

// C++ includes
#include <iostream>
#include <vector>

class TH1F;
class TH2F;
class TChain;
class TDirectory;
class TString;

#include "../../CORE/electronSelections.h"
#include "../../Tools/DileptonHypType.h"

class MyScanChain {

	public:

		MyScanChain() {};
		MyScanChain(cuts_t configured_cuts) : configured_cuts_(configured_cuts) { };
		~MyScanChain() {};

		int ScanChain(bool isData, std::string sampleName, TChain *chain, int nEvents = -1, std::string skimFilePrefix="");

	private:


		void Fill(TH1F** hist, const unsigned int hyp, const float &val, const float &weight);

		void FillAllEleIdHistograms(const unsigned int index, const float &weight, const TString &sampleName, const unsigned int hyp);

		void FormatHist(TH1F** hist, std::string sampleName, std::string name, int n, float min, float max);
		void FormatAllEleIdHistograms(std::string sampleName);

		// for 2D
		void Fill2D(TH2F** hist, const unsigned int hyp, const float &valx, const float &valy, const float &weight);
		void FormatHist2D(TH2F** hist, std::string sampleName, std::string name, int nx, float minx, float maxx, int ny, float miny, float maxy);

		// dealing with cuts
		bool CheckCutsNM1(cuts_t apply, cuts_t remove, cuts_t passed);
		bool CheckCuts(cuts_t apply, cuts_t passed);

		// misc
		enum DileptonHypType hyp_typeToHypType (int hyp_type);

		// configured cuts
		//
		cuts_t configured_cuts_;

		//
		// plots
		//
        TH1F *h1_hyp_pt_[2][4];


};

#endif

