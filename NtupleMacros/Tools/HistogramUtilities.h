
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
#include "Tools/Utilities.h"
#include "Tools/DataSource.h"

class HistogramUtilities {

	public:
		HistogramUtilities(TString fileName, Double_t lumiNorm);
		~HistogramUtilities() { delete file_; }

		TH1F* getHistogram(sources_t theSources, TString var, TString nJets, TString hyp_type, Int_t rebin = 1);
		THStack* getStack(sources_t theSources, TString var, TString nJets, TString hyp_type, Int_t rebin = 1);
		TLegend* getLegend(sources_t theSources, TString var, TString nJets, TString hyp_type);

		void setOrder(std::vector<DataSource> potentialSources);


	private:

		sources_t makeBit(sources_t source) {
			return 1ll << source;
		}

		TFile *file_;
		std::vector<DataSource> sources_;
		Double_t lumiNorm_;
};

#endif

