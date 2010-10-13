
void plot(TString fileName, TString det)
{

	gROOT->ProcessLine(".L tdrStyle.C");
	gROOT->ProcessLine("setTDRStyle()");

	TFile f(fileName, "READ");
	TH1F *h1_numerator = (TH1F*)f.Get("ttbar_hyp_lt_" + det + "_pt_ee");
	TH1F *h1_denom_old = (TH1F*)f.Get("ttbar_hyp_lt_" + det + "_pt_idold_ee");
        TH1F *h1_denom_new = (TH1F*)f.Get("ttbar_hyp_lt_" + det + "_pt_idnew_ee");

	TGraphAsymmErrors *gr_eff_old = new TGraphAsymmErrors();
	gr_eff_old->BayesDivide(h1_denom_old, h1_numerator);
	gr_eff_old->GetYaxis()->SetRangeUser(0, 1.1);
	gr_eff_old->GetYaxis()->SetTitle("Efficiency");
	gr_eff_old->GetXaxis()->SetTitle("LT p_{T} (GeV)");
	gr_eff_old->SetMarkerColor(kRed);
	gr_eff_old->SetLineColor(kRed);

        TGraphAsymmErrors *gr_eff_new = new TGraphAsymmErrors();
        gr_eff_new->BayesDivide(h1_denom_new, h1_numerator);
        gr_eff_new->GetYaxis()->SetRangeUser(0, 1.1);
        gr_eff_new->GetYaxis()->SetTitle("Efficiency");
        gr_eff_new->GetXaxis()->SetTitle("LT p_{T} (GeV)");
	gr_eff_new->SetMarkerColor(kBlue);
	gr_eff_new->SetLineColor(kBlue);

	TLegend *lg = new TLegend(0.5, 0.2, 0.9, 0.4);
	lg->SetFillColor(kWhite);
	lg->SetLineColor(kWhite);
	lg->SetShadowColor(kWhite);
	lg->AddEntry(gr_eff_old, "egamma_looseId", "lp");
	lg->AddEntry(gr_eff_new, "cand01", "lp");

	TCanvas *c1 = new TCanvas();
	c1->cd();
	gr_eff_old->Draw("AP");
	gr_eff_new->Draw("P");
	lg->Draw();
	c1->SaveAs("eff_mc_" + det + ".png");

}

