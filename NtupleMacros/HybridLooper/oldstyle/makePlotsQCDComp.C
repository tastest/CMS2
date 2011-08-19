{

        TFile f("histos_data.root", "READ");
        TFile f_ref("histos_reference_qcdonly.root");
        
        TCanvas *c1 = new TCanvas();
        c1->cd();

	TH1F *h1_whunt_tcmet = (TH1F*)f.Get("whunt_ele_nm1_tcmet_all");
	h1_whunt_tcmet->SetMarkerStyle(20);
	h1_whunt_tcmet->SetDirectory(gDirectory);
	
	TH1F *h1_qcd15_tcmet = (TH1F*)f_ref.Get("QCDpt15_ele_nm1_tcmet_all");
	h1_qcd15_tcmet->SetLineColor(kRed);
	h1_qcd15_tcmet->Scale(0.201/1e+06);
	h1_qcd15_tcmet->SetDirectory(gDirectory);

	TH1F *h1_minbias_tcmet = (TH1F*)f_ref.Get("minbias_ele_nm1_tcmet_all");
	h1_minbias_tcmet->SetLineColor(kBlue);
	h1_minbias_tcmet->Scale(0.201/1e+03);
	h1_minbias_tcmet->SetDirectory(gDirectory);


	h1_whunt_tcmet->Draw("ep");
	h1_qcd15_tcmet->Draw("SAMES");
	h1_minbias_tcmet->Draw("SAMES");


	TH1F *h1_qcd15_tcmet_pthat = (TH1F*)f_ref.Get("QCDpt15_ele_nm1_tcmet_pthat_all");
	h1_qcd15_tcmet_pthat->SetLineColor(kRed);
	h1_qcd15_tcmet_pthat->Scale(0.201/1e+06);
	h1_qcd15_tcmet_pthat->SetDirectory(gDirectory);

	TH1F *h1_minbias_tcmet_pthat = (TH1F*)f_ref.Get("minbias_ele_nm1_tcmet_pthat_all");
	h1_minbias_tcmet_pthat->SetLineColor(kBlue);
	h1_minbias_tcmet_pthat->Scale(0.201/1e+03);
	h1_minbias_tcmet_pthat->SetDirectory(gDirectory);


	TCanvas *c2 = new TCanvas();
	c2->cd();
	h1_minbias_tcmet_pthat->Draw();
	h1_qcd15_tcmet_pthat->Draw("sames");
	

}

