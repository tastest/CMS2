
#ifndef RESULTSHISTOGRAMS_H
#define RESULTSHISTOGRAMS_H

#include "TCanvas.h"
#include "TLegend.h"
#include "TH1F.h"

enum JetBins_t {
        J0,
        J1,
        J2
};

class ResultsHistograms {
        public:
                ResultsHistograms(Color_t bandColor, TString title, TString titleX, TString titleY);
                ~ResultsHistograms() {}

                void add(JetBins_t bin, Float_t truth, Float_t errTruth,
                                Float_t est, Float_t errEst);
                TCanvas* results(Float_t min, Float_t max);

        private:
                TH1F* h1_true_;
                TH1F* h1_true_errLow_;
                TH1F* h1_true_errHigh_;
                TH1F* h1_est_;
        	TLegend *lg_estimate_;
};

#endif

