
#include "DataSource.h"
	
void plotsForWWNote_R(TString source)
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

	// R hists
	DataSource *s_mm;
	DataSource *s_ee;
	if (source == "DY") {
		s_mm = new DataSource("dymm", H_DYMM);
        	s_ee = new DataSource("dyee", H_DYEE);
	}
	else if (source == "ZZ") {
		s_mm = new DataSource("zz", H_ZZ);
                s_ee = new DataSource("zz", H_ZZ);
 	}
	else if (source == "WZ") {
                s_mm = new DataSource("wz", H_ZZ);
                s_ee = new DataSource("wz", H_ZZ);
	}
	else break;

        TH1F *h1_Rdymm   = h1->getRHist(*s_mm, "0j", "mm");
		h1_Rdymm->SetName("h1_Rdymm");
		h1_Rdymm->SetLineColor(kBlue);
		h1_Rdymm->SetFillColor(kWhite);
		h1_Rdymm->SetMarkerStyle(20);
		h1_Rdymm->SetMarkerColor(kBlue);
		h1_Rdymm->GetXaxis()->SetTitle("MET (GeV)");
		h1_Rdymm->GetYaxis()->SetTitle("R_{out/in}");
		h1_Rdymm->GetYaxis()->SetTitleOffset(1.7);
        TH1F *h1_Rdyee   = h1->getRHist(*s_ee, "0j", "ee");
                h1_Rdyee->SetName("h1_Rdyee");
		h1_Rdyee->SetLineColor(kRed);
		h1_Rdyee->SetFillColor(kWhite);
                h1_Rdyee->SetMarkerStyle(26);
                h1_Rdyee->SetMarkerColor(kRed);

	TLegend *lg_Rdy = new TLegend(0.55, 0.75, 0.90, 0.90);
        lg_Rdy->SetFillColor(kWhite);
        lg_Rdy->SetLineColor(kWhite);
	lg_Rdy->SetShadowColor(kWhite);
	lg_Rdy->AddEntry(h1_Rdymm, source + "#rightarrow #mu#mu", "lp");
        lg_Rdy->AddEntry(h1_Rdyee, source + "#rightarrow ee", "lp");	

	TCanvas *c4 = new TCanvas();
	c4->cd();
	h1_Rdymm->Draw("HIST E1");
	h1_Rdyee->Draw("HIST E1 SAME");
	lg_Rdy->Draw();
	h1_Rdymm->GetYaxis()->SetRangeUser(0, 0.5);
	h1_Rdymm->GetXaxis()->SetRangeUser(20, 60);
	c1->RedrawAxis();

	// tidy up
	delete h1;

}

