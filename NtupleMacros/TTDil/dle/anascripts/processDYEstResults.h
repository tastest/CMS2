#ifndef PROCESSDYESTRESULTS_H
#define PROCESSDYESTRESULTS_H

#include "../../../Tools/DataSource.h"
#include "../../../Tools/HistogramUtilities.h"
#include "../../../Tools/Utilities.h"

#include "TH1F.h"
#include "TCanvas.h"

#include "ResultsHistograms.h"


        // ****** R_{out/in} for this jet bin and met cuts
void get_R(HistogramUtilities *hUtil, TString hyp_type, TString nJets, Float_t &R, Float_t &err2);

        // ******* k including correction for e/mu eff
//void get_k(TString nJets, TString hyp_type, Float_t &k, Float_t &err2);

void estimate(HistogramUtilities *hUtil, JetBins_t nJets, TString hyp_type, 
		ResultsHistograms &estResults, ResultsHistograms &nonpeakResults);

void processDYEstResults(TString fileName);

#endif

