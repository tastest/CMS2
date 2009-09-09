
#ifndef HISTOGRAMUTILITIES_H
#define HISTOGRAMUTILITIES_H

#include "TFile.h"

// system includes
#include <iostream>
#include <vector>

// root includes
#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "TLegend.h"
#include "TString.h"

// local utilities
#include "Tools/Utilities.h"
#include "Tools/DataSource.h"

class HistogramUtilities {

public:
  HistogramUtilities(TString fileName, Double_t lumiNorm = 1.0);
  HistogramUtilities(TString fileName, TString fileName2, Double_t lumiNorm = 1.0);
  HistogramUtilities(TString fileName, TString fileName2, TString fileName3, Double_t lumiNorm = 1.0);
  ~HistogramUtilities() { delete file_; }

  TH1F* getHistogram(sources_t theSources, TString var, TString nJets, TString hyp_type, Int_t rebin = 1, TString nameprefix="");
  //TH1F* getHistogramSum(sources_t theSources, TString var, TString nJets, TString hyp_type, Int_t rebin = 1);
  TH2F* get2dHistogram(sources_t theSources, TString var, TString nJets, TString hyp_type, Int_t rebin = 1);

  THStack* getStack(sources_t theSources, TString var, TString nJets, TString hyp_type, Int_t rebin = 1);
  THStack* getSumStack(sources_t theSources, TString var, TString nJets, TString hyp1, TString hyp2, Int_t rebin = 1, TString var2="", double scale=1);
  THStack* get2fileStack(sources_t theSources, TString var, TString nJets, TString hyp1, Int_t rebin = 1);
  THStack* getSumDifStack(sources_t theSources, TString var, TString nJets, TString hyp1, TString hyp2, TString hyp3, Int_t rebin = 1);

  TLegend* getLegend(sources_t theSources, TString var, TString nJets, TString hyp_type);

  void setOrder(std::vector<DataSource> potentialSources);

  void setVerbose( bool set = false ) { verbose_ = set; }

private:

  void setSources();

  sources_t makeBit(sources_t source) {
	return 1ll << source;
  }

  TFile *file_;
  TFile *file2_;
  TFile *file3_;
  
  std::vector<DataSource> sources_;
  Double_t lumiNorm_;
  bool verbose_;
};

#endif

