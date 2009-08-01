
#include "plotResults.h"
#include "Tools/HistogramUtilities.h"
#include "Tools/Utilities.h"
#include "Tools/histtools.cc"

#include "TROOT.h"
#include "TGraphAsymmErrors.h"

        const static sources_t theSources_31X =
                (1ll << H_WENU)          |
                (1ll << H_EM30_80)       |
                (1ll << H_BC30_80);

        const static sources_t theSignal_31X =
                (1ll << H_WENU);

        const static sources_t theBackground_31X =
                (1ll << H_EM30_80)       |
                (1ll << H_BC30_80);

        const static sources_t theSources_22X =
                (1ll << H_QCD30) 	|
                (1ll << H_QCD80)        |
		(1ll << H_WJET_ALP);

	const static sources_t theSignal_22X = 
                (1ll << H_WJET_ALP);

	const static sources_t theBackground_22X = 
                (1ll << H_QCD30);// 	|
//                (1ll << H_QCD80);

	// for 2_2_X
	const static sources_t &theSignal = theSignal_22X;
	const static sources_t &theBackground = theBackground_22X;
	const static sources_t &theSources = theSources_22X;


void plotResults(TString det)
{

	gROOT->ProcessLine(".L ~/tdrStyle.C");
	gROOT->ProcessLine("setTDRStyle()");

	// luminorm for 1pb-1
	HistogramUtilities h1("Results.root", 0.001);
	THStack *st_pt = h1.getStack(theSources, "h1_pt", "", det, 2);
	TLegend *lg_all = h1.getLegend(theSources, "h1_pt", "", det);

	TCanvas *c1 = new TCanvas();
	c1->cd();
	st_pt->Draw();
	lg_all->Draw();
	Utilities::saveCanvas(c1, "results/pt_" + det);	
}

