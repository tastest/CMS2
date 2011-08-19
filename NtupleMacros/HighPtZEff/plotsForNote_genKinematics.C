
{

	gROOT->ProcessLine(".L ~/tdrStyle.C");
	setTDRStyle();
	gROOT->ForceStyle();
	gStyle->SetOptTitle(1);
	
	TFile file("Results_v3.root");

	// Gen Z mass
	//
	/*
	TH1F *h1_mll = (TH1F*)file.Get("dyee_nofilter_Gen_Z_Mass");
	h1_mll->GetXaxis()->SetRangeUser(20, 200);
	h1_mll->GetXaxis()->SetTitle("M_{ll} (GeV/c^{2})");
	h1_mll->GetYaxis()->SetTitle("#Events / 2 GeV/c^{2}");
	
	TCanvas *c1 = new TCanvas();
	c1->cd();
	c1->SetLogy();
	h1_mll->Draw("HIST E1");
	*/
	
	// Gen Z Pt
	
	TH1F *h1_zpt = (TH1F*)file.Get("dyee_nofilter_Gen_Z_Pt");
	h1_zpt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	h1_zpt->GetYaxis()->SetTitle("#Events / 10 GeV/c");
	
	TCanvas *c2 = new TCanvas();
	c2->cd();
	c2->SetLogy();
	h1_zpt->Draw("HIST E1");	
	
	// Gen Z P
	
	TH1F *h1_zp = (TH1F*)file.Get("dyee_nofilter_Gen_Z_Momentum");
	h1_zp->GetXaxis()->SetTitle("p (GeV/c)");
	h1_zp->GetYaxis()->SetTitle("#Events / 10 GeV/c");
	
	TCanvas *c20 = new TCanvas();
	c20->cd();
	c20->SetLogy();
	h1_zp->Draw("HIST E1");	
	
	// Gen Z Eta
	
	TH1F *h1_zetas3 = (TH1F*)file.Get("dyee_nofilter_Gen_Z_Eta");
	h1_zetas3->GetXaxis()->SetTitle("#eta");
	h1_zetas3->GetYaxis()->SetTitle("#Events / 0.2");
	h1_zetas3->SetMarkerColor(kBlack);
	h1_zetas3->SetLineColor(kBlack);
	
	TH1F *h1_zetas1 = (TH1F*)file.Get("dyee_nofilter_Gen1_Z_Eta");
	h1_zetas1->GetXaxis()->SetTitle("#eta");
	h1_zetas1->GetYaxis()->SetTitle("#Events / 0.2");
	h1_zetas3->SetMarkerColor(kRed);
	h1_zetas3->SetLineColor(kRed);
	
	TLegend *lg_genz = new TLegend(0.6, 0.75, 0.9, 0.9);	
	lg_genz->SetFillColor(kWhite);
	lg_genz->SetLineColor(kWhite);
	lg_genz->SetBorderSize(0);
	lg_genz->AddEntry(h1_zetas3, "Status = 3", "lp");
	lg_genz->AddEntry(h1_zetas1, "Status = 1 Daughters", "lp");		
	
	TCanvas *c3 = new TCanvas();
	c3->cd();
	c3->SetLogy();
	h1_zetas3->Draw("HIST E1");	
	h1_zetas1->Draw("SAME HIST E1");	
	lg_genz->Draw();
	c3->SaveAs("plot_genZEta.eps");
	
	// Gen Z mass (from Status = 1 leptons)
	TH1F *h1_mll = (TH1F*)file.Get("dyee_nofilter_Gen_Z_Mass");
	h1_mll->GetXaxis()->SetRangeUser(20, 200);
	h1_mll->GetXaxis()->SetTitle("M_{ll} (GeV/c^{2})");
	h1_mll->GetYaxis()->SetTitle("#Events / 2 GeV/c^{2}");
	
	TH1F *h1_mlls1 = (TH1F*)file.Get("dyee_nofilter_Gen1_Lep_Z_Mass");
	h1_mlls1->SetLineColor(kRed);
	h1_mlls1->SetMarkerColor(kRed);
	h1_mlls1->GetXaxis()->SetTitle("M_{ll} (GeV/c^{2})");
	h1_mlls1->GetYaxis()->SetTitle("#Events / 2 GeV/c^{2}");
	
	TCanvas *c4 = new TCanvas();
	c4->cd();
	c4->SetLogy();
	h1_mll->Draw("HIST E1");
	h1_mlls1->Draw("SAME HIST E1");	
	lg_genz->Draw();
	c4->SaveAs("plot_genZMass.eps");
	
	// Status = 1 Leptons	
	// Pt
	TH1F *h1_l0_pt = (TH1F*)file.Get("dyee_nofilter_Gen1_Lep_Z_Pt0");
	h1_l0_pt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	h1_l0_pt->GetYaxis()->SetTitle("#Events / 2 GeV/c^{2}");
	h1_l0_pt->SetLineColor(kRed);
	h1_l0_pt->SetMarkerColor(kRed);
	TH1F *h1_l1_pt = (TH1F*)file.Get("dyee_nofilter_Gen1_Lep_Z_Pt1");
	h1_l1_pt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	h1_l1_pt->GetYaxis()->SetTitle("#Events / 2 GeV/c^{2}");
	h1_l1_pt->SetLineColor(kBlue);
	h1_l1_pt->SetMarkerColor(kBlue);
	
	/*
	TH1F *h1_l0s3_pt = (TH1F*)file.Get("dyee_nofilter_Gen3_Lep_Z_Pt0");
	h1_l0s3_pt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	h1_l0s3_pt->GetYaxis()->SetTitle("#Events / 2 GeV/c^{2}");
	h1_l0s3_pt->SetLineColor(kBlack);
	h1_l0s3_pt->SetMarkerColor(kBlack);
	h1_l0s3_pt->SetMarkerStyle(23);
	TH1F *h1_l1s3_pt = (TH1F*)file.Get("dyee_nofilter_Gen3_Lep_Z_Pt1");
	h1_l1s3_pt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	h1_l1s3_pt->GetYaxis()->SetTitle("#Events / 2 GeV/c^{2}");
	h1_l1s3_pt->SetMarkerStyle(23);
	h1_l1s3_pt->SetLineColor(kBlack);
	h1_l1s3_pt->SetMarkerColor(kBlack);
	*/
	  
	TLegend *lg_leptons = new TLegend(0.7, 0.75, 0.9, 0.9);	
	lg_leptons->SetFillColor(kWhite);
	lg_leptons->SetLineColor(kWhite);
	lg_leptons->SetBorderSize(0);
	lg_leptons->AddEntry(h1_l0_pt, "Leading lepton (S1)", "lp");
	lg_leptons->AddEntry(h1_l1_pt, "Second lepton (S1)", "lp");
	//lg_leptons->AddEntry(h1_l0s3_pt, "Leading lepton (S3)", "lp");
	//lg_leptons->AddEntry(h1_l1s3_pt, "Second lepton (S3)", "lp");
	
	TCanvas *c5 = new TCanvas();
	c5->cd();
	c5->SetLogy();
	h1_l0_pt->Draw("HIST E1");		
	h1_l1_pt->Draw("SAME HIST E1");	
	//h1_l0s3_pt->Draw("SAME HIST E1");
	//h1_l1s3_pt->Draw("SAME HIST E1");	
	lg_leptons->Draw();
	
	// Eta
	TH1F *h1_l0_eta = (TH1F*)file.Get("dyee_nofilter_Gen1_Lep_Z_Eta0");
	h1_l0_eta->GetXaxis()->SetTitle("#eta");
	h1_l0_eta->GetYaxis()->SetTitle("#Events / 0.2");
	h1_l0_eta->SetLineColor(kRed);
	h1_l0_eta->SetMarkerColor(kRed);
	TH1F *h1_l1_eta = (TH1F*)file.Get("dyee_nofilter_Gen1_Lep_Z_Eta1");
	h1_l1_eta->GetXaxis()->SetTitle("#eta");
	h1_l1_eta->GetYaxis()->SetTitle("#Events / 0.2");
	h1_l1_eta->SetLineColor(kBlue);
	h1_l1_eta->SetMarkerColor(kBlue);
	
	TCanvas *c6 = new TCanvas();
	c6->cd();
	c5->SetLogy();
	h1_l0_eta->Draw("HIST E1");		
	h1_l1_eta->Draw("SAME HIST E1");			
	lg_leptons->Draw();
	
	
}