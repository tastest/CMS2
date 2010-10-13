
#ifndef HISTOGRAMUTILITIES_H
#define HISTOGRAMUTILITIES_H

#include "TFile.h"

// system includes
#include <iostream>
#include <vector>

// root includes
#include "TH1F.h"
#include "THStack.h"
#include "TLegend.h"
#include "TString.h"

// local utilities
#include "Utilities.h"
#include "DataSource.h"

class HistogramUtilities {

	public:
		HistogramUtilities(TString fileName, Double_t metCut, bool doWZ, Double_t lumiNorm);
		~HistogramUtilities() { delete file_; }

		Double_t 	getLumiNorm() 	{ return lumiNorm_; }
		Int_t 		getZLow() 	{ return zLow_; }
		Int_t 		getZHigh() 	{ return zHigh_; }
		Int_t 		getMetCutBin()	{ return metCutBin_; }
		bool 		getDoWZ() 	{ return doWZ_; }

		void printInOutTruth(TString nJets, TString hyp_type); 
		Float_t getN(bool in, const sources_t &theSources, 
			TString nJets, TString hyp_type, Float_t &err2);

		TH1F* getHistogram(sources_t theSources, TString var, TString nJets, TString hyp_type);
		THStack* getStack(sources_t theSources, TString var, TString nJets, TString hyp_type);
		TLegend* getLegend(sources_t theSources);

		TH1F* getRHist(DataSource theSource, TString nJets, TString hyp_type);


	private:

		sources_t makeBit(sources_t source) {
			return 1ll << source;
		}

		TFile *file_;
		std::vector<DataSource> sources_;

		bool doWZ_;
                Double_t lumiNorm_;
		Int_t zLow_;
		Int_t zHigh_;
		Int_t metCutBin_;
};

#endif

