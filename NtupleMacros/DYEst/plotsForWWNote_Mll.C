
#include "DataSource.h"

void setLabels(THStack *st, TString titleX, TString titleY, Float_t max)
{
	st->GetXaxis()->SetTitle(titleX);
	st->GetYaxis()->SetTitle(titleY);
	st->SetMaximum(max);
}
	
void plotsForWWNote_Mll(TString hyp_type)
{

        gROOT->ProcessLine(".L tdrStyle.C"); 
        gROOT->ProcessLine("setTDRStyle()");
        gStyle->SetOptTitle(0);

	gROOT->ProcessLine(".L DataSource.cc+");
        gROOT->ProcessLine(".L Utilities.cc+");
	gROOT->ProcessLine(".L HistogramUtilities.cc+");

	// 1fb-1
	HistogramUtilities *h1 = new HistogramUtilities("DYEstResults_ForWW_MET45_INCL.root", 20.0, true, 0.1);
        gROOT->cd();

	// test get a stack and legend
	THStack *st_peaking = h1->getStack(sources_peaking, "mll", "0j", hyp_type);
        THStack *st_nonpeaking = h1->getStack(sources_nonpeaking, "mll", "0j", hyp_type);
        THStack *st_all = h1->getStack(sources_all, "mll", "0j", hyp_type);

	TLegend *lg_peaking = h1->getLegend(sources_peaking);
        TLegend *lg_nonpeaking = h1->getLegend(sources_nonpeaking);
        TLegend *lg_all = h1->getLegend(sources_all);

	TString hyp_type_tex = "#mu#mu";
	if (hyp_type == "ee") hyp_type_tex = hyp_type;
	if (hyp_type == "em") hyp_type_tex = "e#mu";

	TLatex *l1_lumi = new TLatex(0.2, 0.8, "#int Ldt = 100 pb^{-1}");
	l1_lumi->SetNDC(true);
	l1_lumi->SetTextSize(0.04);

	// draw the results
	TCanvas *c1 = new TCanvas();
	c1->cd();
	st_peaking->Draw("HIST");
	setLabels(st_peaking, "M_{" + hyp_type_tex + "} GeV", "#Events/5 GeV", 3.0);
	lg_peaking->Draw("SAME");
        l1_lumi->Draw();

        TCanvas *c2 = new TCanvas();
        c2->cd();
        st_nonpeaking->Draw("HIST");
        setLabels(st_nonpeaking, "M_{" + hyp_type_tex + "} GeV", "#Events/5 GeV", 3.0);
        lg_nonpeaking->Draw("SAME");
        l1_lumi->Draw();

        TCanvas *c3 = new TCanvas();
        c3->cd();
        st_all->Draw("HIST");
        setLabels(st_all, "M_{" + hyp_type_tex + "} GeV", "#Events/5 GeV", 3.0);
        lg_all->Draw("SAME");
        l1_lumi->Draw();

	// tidy up
	delete h1;

}

