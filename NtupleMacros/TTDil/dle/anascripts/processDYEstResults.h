#ifndef PROCESSDYESTRESULTS_H
#define PROCESSDYESTRESULTS_H

#include "../../../Tools/DataSource.h"
#include "../../../Tools/HistogramUtilities.h"
#include "../../../Tools/Utilities.h"

#include "TH1F.h"
#include "TCanvas.h"

#include "ResultsHistograms.h"

#include "../../../Tools/AllDataSources.h"

#include "TROOT.h"
#include "TStyle.h"

const static sources_t sources_dy =
(1ll << H_DYMM)        |
(1ll << H_DYEE);

const static sources_t sources_peaking =
(1ll << H_DYMM)        |
(1ll << H_DYEE)        |
(1ll << H_ZZ)          |
(1ll << H_WZ);

const static sources_t sources_nonpeak =
(1ll << H_WW)          |
(1ll << H_TTBAR)       |
(1ll << H_DYTT)	|
(1ll << H_WJETS);

const static sources_t sources_all =
(1ll << H_TTBAR) |
(1ll << H_WW) |
(1ll << H_WZ) |
(1ll << H_ZZ) |
(1ll << H_DYMM) |
(1ll << H_DYEE)	|
(1ll << H_DYTT) |
(1ll << H_WJETS);

const static int metCutBin = 7;

		// ****** Get the histogram of R as a function of MET for the
		// ****** desired datasource
TH1F* getRHist(HistogramUtilities *hUtil, sources_t theSource, TString nJets, TString hyp_type);

        // ****** R_{out/in} for this jet bin and met cuts
void get_R(HistogramUtilities *hUtil, TString hyp_type, TString nJets, Float_t &R, Float_t &err2);

        // ******* k including correction for e/mu eff
void get_k(HistogramUtilities *hUtil, TString nJets, TString hyp_type, Float_t &k, Float_t &err2);

void estimate(HistogramUtilities *hUtil, JetBins_t nJets, TString hyp_type, 
		ResultsHistograms &estResults, ResultsHistograms &nonpeakResults);

//void processDYEstResults(TString fileName);
void processDYEstResults(TString det, TString fileStamp, TString refFileStamp, TString norm, const float &luminorm);

#endif

