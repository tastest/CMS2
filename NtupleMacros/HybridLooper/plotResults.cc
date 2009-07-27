
#include "plotResults.h"
#include "Tools/HistogramUtilities.h"
#include "Tools/Utilities.h"
#include "Tools/histtools.cc"

        const static sources_t theSources =
                (1ll << H_WENU)          |
                (1ll << H_EM30_80)       |
                (1ll << H_BC30_80);

        const static sources_t theSignal =
                (1ll << H_WENU);

        const static sources_t theBackground =
                (1ll << H_EM30_80)       |
                (1ll << H_BC30_80);


void plotEff(HistogramUtilities &h1, TString name, TString det)
{
	TH1F *h1_signal = h1.getHistogram(theSignal, name, "", det);
        TH1F *h1_background = h1.getHistogram(theBackground, name, "", det);

	// latter bool is "ascending"
        TGraph *gr = (TGraph*)(eff_rej(*h1_signal, *h1_background, true, true).Clone());
        gr->SetTitle(name + ";Signal;Background");
	gr->SetMarkerStyle(23);

        TCanvas *c = new TCanvas();
        c->cd();
        gr->Draw("AP");
	Utilities::saveCanvas(c, "eff_" + name + "_" + det);

	delete gr;
	delete h1_signal;
	delete h1_background;
	delete c;

}

void plotResults(TString det)
{

	// luminorm for 1pb-1
	HistogramUtilities h1("Results.root", 0.001);
	
	THStack *st_pt = h1.getStack(theSources, "h1_pt", "", det);
        THStack *st_eta = h1.getStack(theSources, "h1_eta", "", det);
	TLegend *lg_all = h1.getLegend(theSources, "h1_pt", "", det);

        THStack *st_ecalIso03 = h1.getStack(theSources, "h1_ecalIso03", "", det, 2);
        THStack *st_hcalIso03 = h1.getStack(theSources, "h1_hcalIso03", "", det, 2);
        THStack *st_tkIso03 = h1.getStack(theSources, "h1_tkIso03", "", det, 2);
        THStack *st_esJuraIso03 = h1.getStack(theSources, "h1_esJuraIso03", "", det);
        THStack *st_wwIso = h1.getStack(theSources, "h1_wwIso", "", det, 2);
	
	plotEff(h1, "h1_esJuraIso03", det);
        plotEff(h1, "h1_ecalIso03", det);
        plotEff(h1, "h1_hcalIso03", det);
        plotEff(h1, "h1_tkIso03", det);

	TCanvas *c1 = new TCanvas();
	c1->cd();
	st_pt->Draw();
	lg_all->Draw();
	Utilities::saveCanvas(c1, "pt_" + det);	

        TCanvas *c2 = new TCanvas();
        c2->cd();
        st_eta->Draw();
        lg_all->Draw();
        Utilities::saveCanvas(c2, "eta_" + det);

        TCanvas *c3 = new TCanvas();
        c3->cd();
        st_ecalIso03->Draw();
        lg_all->Draw();
        Utilities::saveCanvas(c3, "ecalIso03_" + det);

        TCanvas *c4 = new TCanvas();
        c4->cd();
        st_hcalIso03->Draw();
        lg_all->Draw();
        Utilities::saveCanvas(c4, "hcalIso03_" + det);

        TCanvas *c5 = new TCanvas();
        c5->cd();
        st_tkIso03->Draw();
        lg_all->Draw();
        Utilities::saveCanvas(c5, "tkIso03_" + det);

        TCanvas *c6 = new TCanvas();
        c6->cd();
        st_wwIso->Draw();
        lg_all->Draw();
        Utilities::saveCanvas(c6, "wwIso_" + det);

        TCanvas *c7 = new TCanvas();
        c7->cd();
        st_esJuraIso03->Draw();
        lg_all->Draw();
        Utilities::saveCanvas(c7, "esJuraIso_" + det);

}


