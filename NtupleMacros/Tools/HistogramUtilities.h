
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
#include "Utilities.h"
#include "DataSource.h"

class HistogramUtilities {

public:
  HistogramUtilities(TString fileName, std::vector<DataSource> potentialSources, Double_t lumiNorm = 1.0);
  HistogramUtilities(TString fileName, TString fileName2, std::vector<DataSource> potentialSources, Double_t lumiNorm = 1.0);
  ~HistogramUtilities() { delete file_; }

  TH1F* getHistogram(sources_t theSources, TString var, TString nJets, TString hyp_type, Int_t rebin = 1, TString nameprefix="");
  TH1F* getHistogramSum(sources_t theSources, TString var, TString nJets, TString hyp1, TString hyp2, Int_t rebin = 1);
  TH2F* get2dHistogram(sources_t theSources, TString var, TString nJets, TString hyp_type, Int_t rebin = 1);

  THStack* getStack(sources_t theSources, TString var, TString nJets, TString hyp_type, Int_t rebin = 1, bool plotCumulated = false, bool cumulateAscending = false);
  THStack* getSumStack(sources_t theSources, TString var, TString nJets, TString hyp1, TString hyp2, Int_t rebin = 1, TString var2="", double scale=1);
  THStack* get2fileStack(sources_t theSources, TString var, TString nJets, TString hyp1, Int_t rebin = 1);
  THStack* getSumDifStack(sources_t theSources, TString var, TString nJets, TString hyp1, TString hyp2, TString hyp3, Int_t rebin = 1, Int_t posneg=0);

  TLegend* getLegend(sources_t theSources, TString var, TString nJets, TString hyp_type);

  void setOrder(std::vector<DataSource> potentialSources);

  void setVerbose( bool set = false ) { verbose_ = set; }

private:

  sources_t makeBit(sources_t source) {
	return (1ll << source);
  }

  TFile *file_;
  TFile *file2_;
  //TFile *file3_;
  
  std::vector<DataSource> sources_;
  Double_t lumiNorm_;
  bool verbose_;
};


//This standalone function is utility for separating positive or negative bin contents
TH1F* GetPosNeg( TH1F* h1_temp, Int_t posneg );


#endif

