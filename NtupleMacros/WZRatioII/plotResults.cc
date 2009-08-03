
#include "plotResults.h"
#include "Tools/HistogramUtilities.h"
#include "Tools/Utilities.h"
#include "Tools/histtools.cc"

#include "TROOT.h"
#include "TGraphAsymmErrors.h"

        const static sources_t theSources_22X =
                (1ll << H_ZEEJET_ALP) 	|
             	(1ll << H_ZMMJET_ALP)   |
                (1ll << H_ZTTJET_ALP)   |
		(1ll << H_WJET_ALP);

	// for 2_2_X
	const static sources_t &theSources = theSources_22X;


void plotResultsLep(TString hyp)
{

        gROOT->ProcessLine(".L ~/tdrStyle.C");
        gROOT->ProcessLine("setTDRStyle()");

        // sources/ordering for Z stack plots
        std::vector<DataSource> zSources;

        if (hyp == "e") {
                zSources.push_back(     fH_WJET_ALP()   );
                zSources.push_back(     fH_ZEEJET_ALP()   );
                zSources.push_back(     fH_ZTTJET_ALP()   );
                zSources.push_back(     fH_ZMMJET_ALP()   );

        } else if (hyp == "m") {
                zSources.push_back(     fH_WJET_ALP()   );
                zSources.push_back(     fH_ZMMJET_ALP()   );
                zSources.push_back(     fH_ZTTJET_ALP()   );
                zSources.push_back(     fH_ZEEJET_ALP()   );
        }

        // luminorm for 1pb-1
        HistogramUtilities h1("Results.root", 0.001);
        h1.setOrder(zSources);

        TLegend *lg_all = h1.getLegend(theSources, "lep_met", "", hyp);
        THStack *st_lep_met = h1.getStack(theSources, "lep_met", "", hyp);

        TCanvas *c1 = new TCanvas();
        c1->cd();
        st_lep_met->Draw();
        lg_all->Draw();
        Utilities::saveCanvas(c1, "results/lep_met_" + hyp);

}


void plotResultsDilep(TString hyp)
{

	gROOT->ProcessLine(".L ~/tdrStyle.C");
	gROOT->ProcessLine("setTDRStyle()");

	// sources/ordering for Z stack plots
	std::vector<DataSource> zSources;

	if (hyp == "ee") {
        	zSources.push_back(     fH_ZEEJET_ALP()   );
	        zSources.push_back(     fH_WJET_ALP()   );
                zSources.push_back(     fH_ZTTJET_ALP()   );
                zSources.push_back(     fH_ZMMJET_ALP()   );

	} else if (hyp == "mm") {
                zSources.push_back(     fH_ZMMJET_ALP()   );
                zSources.push_back(     fH_WJET_ALP()   );
                zSources.push_back(     fH_ZTTJET_ALP()   );
                zSources.push_back(     fH_ZEEJET_ALP()   );
	}

	// luminorm for 1pb-1
	HistogramUtilities h1("Results.root", 0.001);
	h1.setOrder(zSources);

	TLegend *lg_all = h1.getLegend(theSources, "dilep_mass", "", hyp);
        THStack *st_dilep_mass = h1.getStack(theSources, "dilep_mass", "", hyp);
        THStack *st_dilep_met = h1.getStack(theSources, "dilep_met", "", hyp);

	TCanvas *c1 = new TCanvas();
	c1->cd();
	st_dilep_mass->Draw();
	lg_all->Draw();
	Utilities::saveCanvas(c1, "results/dilep_mass_" + hyp);	

        TCanvas *c2 = new TCanvas();
        c2->cd();
        st_dilep_met->Draw();
        lg_all->Draw();
        Utilities::saveCanvas(c2, "results/dilep_met_" + hyp);


}

