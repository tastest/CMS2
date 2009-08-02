
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
                //(1ll << H_QCD80)        |
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


void plotResults(TString hyp)
{

	gROOT->ProcessLine(".L ~/tdrStyle.C");
	gROOT->ProcessLine("setTDRStyle()");

	// luminorm for 1pb-1
	HistogramUtilities h1("Results.root", 0.001);
	TLegend *lg_all = h1.getLegend(theSources, "lep_pt", "", hyp);

        THStack *st_pt = h1.getStack(theSources, "lep_pt", "", hyp);
        THStack *st_met = h1.getStack(theSources, "lep_met", "", hyp);
        THStack *st_tkIso = h1.getStack(theSources, "lep_tkIso", "", hyp);


	TCanvas *c1 = new TCanvas();
	c1->cd();
	st_pt->Draw();
	lg_all->Draw();
	Utilities::saveCanvas(c1, "results/lep_pt_" + hyp);	

        TCanvas *c2 = new TCanvas();
        c2->cd();
        st_met->Draw();
        lg_all->Draw();
        Utilities::saveCanvas(c2, "results/lep_met_" + hyp);

        TCanvas *c3 = new TCanvas();
        c3->cd();
        st_tkIso->Draw();
        lg_all->Draw();
        Utilities::saveCanvas(c3, "results/lep_tkIso_" + hyp);


}

