
#include "plotResults.h"
#include "Tools/HistogramUtilities.h"
#include "Tools/Utilities.h"
#include "Tools/histtools.cc"

#include "TROOT.h"

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
                (1ll << H_QCD30);

	const static sources_t theSignal_22X = 
                (1ll << H_QCD30);

	const static sources_t theBackground_22X = 
                (1ll << H_QCD30);

	// for 2_2_X
	const static sources_t &theSignal = theSignal_22X;
	const static sources_t &theBackground = theBackground_22X;
	const static sources_t &theSources = theSources_22X;

void plotEff(HistogramUtilities &h1, TString name, TString det, bool ascending)
{

	TH1F *h1_signal = h1.getHistogram(theSignal, name, "", det);
        TH1F *h1_background = h1.getHistogram(theBackground, name, "", det);

	// latter bool is "ascending"
        TGraph *gr = (TGraph*)(eff_rej(*h1_signal, *h1_background, true, ascending).Clone());
        gr->SetTitle(name + ";Signal;Background");
	gr->SetMarkerStyle(23);

        TCanvas *c = new TCanvas();
        c->cd();
        gr->Draw("AP");
	Utilities::saveCanvas(c, "results/eff_" + name + "_" + det);

	delete gr;
	delete h1_signal;
	delete h1_background;
	delete c;

}

void plotResults(TString det)
{

	gROOT->ProcessLine(".L ~/tdrStyle.C");
	gROOT->ProcessLine("setTDRStyle()");

	// luminorm for 1pb-1
	HistogramUtilities h1("Results.root", 0.001);
	THStack *st_pt = h1.getStack(theSources, "h1_pt", "", det, 2);
        THStack *st_eta = h1.getStack(theSources, "h1_eta", "", det);
	TLegend *lg_all = h1.getLegend(theSources, "h1_pt", "", det);

        THStack *st_ecalIso03 = h1.getStack(theSources, "h1_ecalIso03", "", det);
        THStack *st_hcalIso03 = h1.getStack(theSources, "h1_hcalIso03", "", det);
        THStack *st_tkIso03 = h1.getStack(theSources, "h1_tkIso03", "", det);
        //THStack *st_esJuraIso03 = h1.getStack(theSources, "h1_esJuraIso03", "", det);
        THStack *st_wwIso = h1.getStack(theSources, "h1_wwIso", "", det, 4);

	//plotEff(h1, "h1_esJuraIso03", det, true);
        plotEff(h1, "h1_ecalIso03", det, true);
        plotEff(h1, "h1_hcalIso03", det, true);
        plotEff(h1, "h1_tkIso03", det, true);
	plotEff(h1, "h1_wwIso", det, false);

	TCanvas *c1 = new TCanvas();
	c1->cd();
	st_pt->Draw();
	lg_all->Draw();
	Utilities::saveCanvas(c1, "results/pt_" + det);	

        TCanvas *c2 = new TCanvas();
        c2->cd();
        st_eta->Draw();
        lg_all->Draw();
        Utilities::saveCanvas(c2, "results/eta_" + det);

        TCanvas *c3 = new TCanvas();
        c3->cd();
        st_ecalIso03->Draw();
        lg_all->Draw();
        Utilities::saveCanvas(c3, "results/ecalIso03_" + det);

        TCanvas *c4 = new TCanvas();
        c4->cd();
        st_hcalIso03->Draw();
        lg_all->Draw();
        Utilities::saveCanvas(c4, "results/hcalIso03_" + det);

        TCanvas *c5 = new TCanvas();
        c5->cd();
        st_tkIso03->Draw();
        lg_all->Draw();
        Utilities::saveCanvas(c5, "results/tkIso03_" + det);

        TCanvas *c6 = new TCanvas();
        c6->cd();
        st_wwIso->Draw();
        lg_all->Draw();
        Utilities::saveCanvas(c6, "results/wwIso_" + det);

        TCanvas *c7 = new TCanvas();
        c7->cd();
        st_esJuraIso03->Draw();
        lg_all->Draw();
        Utilities::saveCanvas(c7, "results/esJuraIso_" + det);

	// electron id related
	//
        THStack *st_dEtaIn = h1.getStack(theSources, "h1_dEtaIn", "", det, 2);
        THStack *st_dPhiIn = h1.getStack(theSources, "h1_dPhiIn", "", det, 2);
        THStack *st_hoe = h1.getStack(theSources, "h1_hoe", "", det, 4);
        THStack *st_sigmaIEtaIEta = h1.getStack(theSources, "h1_sigmaIEtaIEta", "", det);

        plotEff(h1, "h1_hoe", det, true);
	
        TCanvas *c_id1 = new TCanvas();
        c_id1->cd();
        st_dEtaIn->Draw();
        lg_all->Draw();
        Utilities::saveCanvas(c_id1, "results/dEtaIn_" + det);

        TCanvas *c_id2 = new TCanvas();
        c_id2->cd();
        st_dPhiIn->Draw();
        lg_all->Draw();
        Utilities::saveCanvas(c_id2, "results/dPhiIn_" + det);

        TCanvas *c_id3 = new TCanvas();
        c_id3->cd();
        st_hoe->Draw();
        lg_all->Draw();
        Utilities::saveCanvas(c_id3, "results/hoe_" + det);

        TCanvas *c_id4 = new TCanvas();
        c_id4->cd();
        st_sigmaIEtaIEta->Draw();
        lg_all->Draw();
        Utilities::saveCanvas(c_id4, "results/sigmaIEtaIEta_" + det);



}


