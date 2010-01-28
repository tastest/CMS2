
{

	TFile f("histos_mc_2x.root", "READ");
	TH1F *h1_numerator = (TH1F*)f.Get("ttbar_hyp_lt_pt_ee");
	TH1F *h1_denom_old = (TH1F*)f.Get("ttbar_hyp_lt_pt_idold_ee");

	TGraphAsymmErrors *gr_eff_old = new TGraphAsymmErrors();
	gr_eff_old->BayesDivide(h1_denom_old, h1_numerator);
	gr_eff_old->GetYaxis()->SetRangeUser(0, 1.1);
	gr_eff_old->GetYaxis()->SetTitle("Efficiency");
	gr_eff_old->GetXaxis()->SetTitle("LT p_{T} (GeV)");

	TCanvas *c1 = new TCanvas();
	c1->cd();
	gr_eff_old->Draw("AP");

}

