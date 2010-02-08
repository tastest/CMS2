
void plot(TString fileName, TString det, TString version)
{

	gROOT->ProcessLine(".L tdrStyle.C");
	gROOT->ProcessLine("setTDRStyle()");

	TFile f(fileName, "READ");
	TH1F *h1_numerator = (TH1F*)f.Get("ttbar_hyp_lt_" + det + "_pt_ee");
	TH1F *h1_denom_old = (TH1F*)f.Get("ttbar_hyp_lt_" + det + "_pt_idold_ee");
	TH1F *h1_denom_new = (TH1F*)f.Get("ttbar_hyp_lt_" + det + "_pt_idnew_ee");

	TH1F *h1_denom_iso_old = (TH1F*)f.Get("ttbar_hyp_lt_" + det + "_pt_isoold_ee");
	TH1F *h1_denom_iso_new = (TH1F*)f.Get("ttbar_hyp_lt_" + det + "_pt_isonew_ee");

	TH1F *h1_denom_iso_new_cand1 = (TH1F*)f.Get("ttbar_hyp_lt_" + det + "_pt_isonew_cand1_ee");

	TH1F *h1_denom_id1_iso1_conv = (TH1F*)f.Get("ttbar_hyp_lt_" + det + "_pt_id1_iso1_conv_ee");

	TH1F *h1_denom_conv = (TH1F*)f.Get("ttbar_hyp_lt_" + det + "_pt_conv_ee");

	//
	// id
	//

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

	//
	// iso
	//

	TGraphAsymmErrors *gr_eff_iso_old = new TGraphAsymmErrors();
	gr_eff_iso_old->BayesDivide(h1_denom_iso_old, h1_numerator);
	gr_eff_iso_old->GetYaxis()->SetRangeUser(0, 1.1);
	gr_eff_iso_old->GetYaxis()->SetTitle("Efficiency");
	gr_eff_iso_old->GetXaxis()->SetTitle("LT p_{T} (GeV)");
	gr_eff_iso_old->SetMarkerColor(kRed);
	gr_eff_iso_old->SetLineColor(kRed);

	TGraphAsymmErrors *gr_eff_iso_new = new TGraphAsymmErrors();
	gr_eff_iso_new->BayesDivide(h1_denom_iso_new, h1_numerator);
	gr_eff_iso_new->GetYaxis()->SetRangeUser(0, 1.1);
	gr_eff_iso_new->GetYaxis()->SetTitle("Efficiency");
	gr_eff_iso_new->GetXaxis()->SetTitle("LT p_{T} (GeV)");
	gr_eff_iso_new->SetMarkerColor(kBlue);
	gr_eff_iso_new->SetLineColor(kBlue);

	TGraphAsymmErrors *gr_eff_iso_new_cand1 = new TGraphAsymmErrors();
	gr_eff_iso_new_cand1->BayesDivide(h1_denom_iso_new_cand1, h1_numerator);
	gr_eff_iso_new_cand1->GetYaxis()->SetRangeUser(0, 1.1);
	gr_eff_iso_new_cand1->GetYaxis()->SetTitle("Efficiency");
	gr_eff_iso_new_cand1->GetXaxis()->SetTitle("LT p_{T} (GeV)");
	gr_eff_iso_new_cand1->SetMarkerColor(kBlack);
	gr_eff_iso_new_cand1->SetLineColor(kBlack);

	TGraphAsymmErrors *gr_eff_id1_iso1_conv = new TGraphAsymmErrors();
	gr_eff_id1_iso1_conv->BayesDivide(h1_denom_id1_iso1_conv, h1_numerator);
	gr_eff_id1_iso1_conv->GetYaxis()->SetRangeUser(0, 1.1);
	gr_eff_id1_iso1_conv->GetYaxis()->SetTitle("Efficiency");
	gr_eff_id1_iso1_conv->GetXaxis()->SetTitle("LT p_{T} (GeV)");
	gr_eff_id1_iso1_conv->SetMarkerColor(kMagenta);
	gr_eff_id1_iso1_conv->SetLineColor(kMagenta);
	gr_eff_id1_iso1_conv->SetMarkerStyle(25);

	TGraphAsymmErrors *gr_eff_conv = new TGraphAsymmErrors();
	gr_eff_conv->BayesDivide(h1_denom_conv, h1_numerator);
	gr_eff_conv->GetYaxis()->SetRangeUser(0, 1.1);
	gr_eff_conv->GetYaxis()->SetTitle("Efficiency");
	gr_eff_conv->GetXaxis()->SetTitle("LT p_{T} (GeV)");
	gr_eff_conv->SetMarkerColor(kMagenta);
	gr_eff_conv->SetLineColor(kMagenta);
	gr_eff_conv->SetMarkerStyle(22);


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
	c1->SaveAs("results/eff_mc" + version + "_" + det + ".png");


	lg->Clear();
	lg->AddEntry(gr_eff_iso_old, "relsusy_iso", "lp");
	lg->AddEntry(gr_eff_iso_new, "relsusy_iso_cand0", "lp");
	lg->AddEntry(gr_eff_iso_new_cand1, "relsusy_iso_cand1", "lp");

	c1->cd();
	gr_eff_iso_old->Draw("AP");
	gr_eff_iso_new->Draw("P");
	gr_eff_iso_new_cand1->Draw("P");
	lg->Draw();
	c1->SaveAs("results/effiso_mc" + version + "_" + det + ".png");

	lg->Clear();
	lg->AddEntry(gr_eff_id1_iso1_conv, "id1_iso1_conv", "lp");

	c1->cd();
	gr_eff_id1_iso1_conv->Draw("AP");
	lg->Draw();
	c1->SaveAs("results/eff_id1_iso1_conv_mc" + version + "_" + det + ".png");

	lg->Clear();
	lg->AddEntry(gr_eff_conv, "conv", "lp");

	c1->cd();
	gr_eff_conv->Draw("AP");
	lg->Draw();
	c1->SaveAs("results/eff_conv_mc" + version + "_" + det + ".png");


}

