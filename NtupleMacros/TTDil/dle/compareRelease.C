
void formatHist(TH1F *hist, bool data)
{
        hist->Sumw2();
        if (data) {
                hist->SetMarkerStyle(20);
                hist->SetLineColor(kRed);
                hist->SetMarkerColor(kRed);
        }
        else {
                hist->SetLineColor(kBlack);
                hist->SetFillColor(kYellow);
        }
}

double histMax(TH1F *hist1, TH1F *hist2)
{

        if (hist1->GetMaximum() > hist2->GetMaximum()) return hist1->GetMaximum()*1.1;
        return hist2->GetMaximum()*1.1;
}


double histMaxLimited(TH1F *hist1, TH1F *hist2)
{

        if (hist1->GetMaximum(hist1->GetMaximum()) > hist2->GetMaximum(hist2->GetMaximum())) return hist1->GetMaximum(hist1->GetMaximum())*1.1;
        return hist2->GetMaximum(hist2->GetMaximum())*1.1;
}

void saveCanvas(TCanvas *c1, const TString name, bool log)
{

        TString dir = "results/";
        c1->SetLogy(log);
        c1->RedrawAxis();
        if (log) {
                c1->SaveAs(dir + name + "_log.png");
        }
        else {
                c1->SaveAs(dir + name + "_lin.png");
        }
}

void compareRelease(TString fileName2X, TString fileName3X, TString det)
{

	gROOT->ProcessLine(".L tdrStyle.C");
	gROOT->ProcessLine("setTDRStyle()");

	//
	// get files to compare
	//
	TFile f2X(fileName2X, "READ");
	TFile f3X(fileName3X, "READ");

	//
	// get histograms of each variable
	//
	TH1F *h1_hoe_2x = (TH1F*)f2X.Get("ttbar_hyp_lt_" + det + "_hoe_ee");
        TH1F *h1_hoe_3x = (TH1F*)f3X.Get("ttbar_hyp_lt_" + det + "_hoe_ee");
	formatHist(h1_hoe_2x, false);
        formatHist(h1_hoe_3x, true);

        TH1F *h1_dEtaIn_2x = (TH1F*)f2X.Get("ttbar_hyp_lt_" + det + "_dEtaIn_ee");
        TH1F *h1_dEtaIn_3x = (TH1F*)f3X.Get("ttbar_hyp_lt_" + det + "_dEtaIn_ee");
        formatHist(h1_dEtaIn_2x, false);
        formatHist(h1_dEtaIn_3x, true);

        TH1F *h1_dPhiIn_2x = (TH1F*)f2X.Get("ttbar_hyp_lt_" + det + "_dPhiIn_ee");
        TH1F *h1_dPhiIn_3x = (TH1F*)f3X.Get("ttbar_hyp_lt_" + det + "_dPhiIn_ee");
        formatHist(h1_dPhiIn_2x, false);
        formatHist(h1_dPhiIn_3x, true);

        TH1F *h1_sigmaIEtaIEta_2x = (TH1F*)f2X.Get("ttbar_hyp_lt_" + det + "_sigmaIEtaIEta_ee");
        TH1F *h1_sigmaIEtaIEta_3x = (TH1F*)f3X.Get("ttbar_hyp_lt_" + det + "_sigmaIEtaIEta_ee");
        formatHist(h1_sigmaIEtaIEta_2x, false);
        formatHist(h1_sigmaIEtaIEta_3x, true);

        TH1F *h1_d0_2x = (TH1F*)f2X.Get("ttbar_hyp_lt_" + det + "_d0_ee");
        TH1F *h1_d0_3x = (TH1F*)f3X.Get("ttbar_hyp_lt_" + det + "_d0_ee");
        formatHist(h1_d0_2x, false);
        formatHist(h1_d0_3x, true);

        TH1F *h1_E2x5MaxOver5x5_2x = (TH1F*)f2X.Get("ttbar_hyp_lt_" + det + "_E2x5MaxOver5x5_ee"); 
        TH1F *h1_E2x5MaxOver5x5_3x = (TH1F*)f3X.Get("ttbar_hyp_lt_" + det + "_E2x5MaxOver5x5_ee"); 
        formatHist(h1_E2x5MaxOver5x5_2x, false);
        formatHist(h1_E2x5MaxOver5x5_3x, true);

        TH1F *h1_ecalIso_2x = (TH1F*)f2X.Get("ttbar_hyp_lt_" + det + "_ecalIso_ee");
        TH1F *h1_ecalIso_3x = (TH1F*)f3X.Get("ttbar_hyp_lt_" + det + "_ecalIso_ee");
        formatHist(h1_ecalIso_2x, false);
        formatHist(h1_ecalIso_3x, true);

        TH1F *h1_hcalIso_2x = (TH1F*)f2X.Get("ttbar_hyp_lt_" + det + "_hcalIso_ee");
        TH1F *h1_hcalIso_3x = (TH1F*)f3X.Get("ttbar_hyp_lt_" + det + "_hcalIso_ee");
        formatHist(h1_hcalIso_2x, false);
        formatHist(h1_hcalIso_3x, true);

        TH1F *h1_tkIso_2x = (TH1F*)f2X.Get("ttbar_hyp_lt_" + det + "_tkIso_ee");
        TH1F *h1_tkIso_3x = (TH1F*)f3X.Get("ttbar_hyp_lt_" + det + "_tkIso_ee");
        formatHist(h1_tkIso_2x, false);
        formatHist(h1_tkIso_3x, true);


	//
	// legend
	//
        TLegend *lg = new TLegend(0.6, 0.8, 0.95, 0.94);
        lg->SetFillColor(kWhite);
        lg->SetLineColor(kWhite);
        lg->SetShadowColor(kWhite);
        TString upperDet;
        if (det == "ee") upperDet = "|#eta| > 1.5";
        else upperDet = "|#eta| < 1.5";
        lg->AddEntry(h1_hoe_2x, "2X MC (" + upperDet + ")", "f");
        lg->AddEntry(h1_hoe_3x, "3X MC (" + upperDet + ")", "lp");

	//
	// plot
	//

	TCanvas *c1 = new TCanvas();
	c1->cd();
	h1_hoe_2x->Scale(h1_hoe_3x->Integral()/h1_hoe_2x->Integral());
	h1_hoe_2x->GetYaxis()->SetRangeUser(0.01, histMax(h1_hoe_2x, h1_hoe_3x));
	h1_hoe_2x->Draw("HIST");
	h1_hoe_3x->Draw("SAME E1");
	lg->Draw();
	saveCanvas(c1, "comparison_hoe_" + det, false);
        saveCanvas(c1, "comparison_hoe_" + det, true);

        h1_dEtaIn_2x->Scale(h1_dEtaIn_3x->Integral()/h1_dEtaIn_2x->Integral());
        h1_dEtaIn_2x->GetYaxis()->SetRangeUser(0.01, histMax(h1_dEtaIn_2x, h1_dEtaIn_3x));
        h1_dEtaIn_2x->Draw("HIST");
        h1_dEtaIn_3x->Draw("SAME E1");
        lg->Draw();
        saveCanvas(c1, "comparison_dEtaIn_" + det, false);
        saveCanvas(c1, "comparison_dEtaIn_" + det, true);

        h1_dPhiIn_2x->Scale(h1_dPhiIn_3x->Integral()/h1_dPhiIn_2x->Integral());
        h1_dPhiIn_2x->GetYaxis()->SetRangeUser(0.01, histMax(h1_dPhiIn_2x, h1_dPhiIn_3x));
        h1_dPhiIn_2x->Draw("HIST");
        h1_dPhiIn_3x->Draw("SAME E1");
        lg->Draw();
        saveCanvas(c1, "comparison_dPhiIn_" + det, false);
        saveCanvas(c1, "comparison_dPhiIn_" + det, true);

        h1_sigmaIEtaIEta_2x->Scale(h1_sigmaIEtaIEta_3x->Integral()/h1_sigmaIEtaIEta_2x->Integral());
        h1_sigmaIEtaIEta_2x->GetYaxis()->SetRangeUser(0.01, histMax(h1_sigmaIEtaIEta_2x, h1_sigmaIEtaIEta_3x));
        h1_sigmaIEtaIEta_2x->Draw("HIST");
        h1_sigmaIEtaIEta_3x->Draw("SAME E1");
        lg->Draw();
        saveCanvas(c1, "comparison_sigmaIEtaIEta_" + det, false);
        saveCanvas(c1, "comparison_sigmaIEtaIEta_" + det, true);

        h1_d0_2x->Scale(h1_d0_3x->Integral()/h1_d0_2x->Integral());
        h1_d0_2x->GetYaxis()->SetRangeUser(0.01, histMax(h1_d0_2x, h1_d0_3x));
        h1_d0_2x->Draw("HIST");
        h1_d0_3x->Draw("SAME E1");
        lg->Draw();
        saveCanvas(c1, "comparison_d0_" + det, false);
        saveCanvas(c1, "comparison_d0_" + det, true);

        h1_E2x5MaxOver5x5_2x->Scale(h1_E2x5MaxOver5x5_3x->Integral()/h1_E2x5MaxOver5x5_2x->Integral());
        h1_E2x5MaxOver5x5_2x->GetYaxis()->SetRangeUser(0.01, histMax(h1_E2x5MaxOver5x5_2x, h1_E2x5MaxOver5x5_3x));
        h1_E2x5MaxOver5x5_2x->Draw("HIST");
        h1_E2x5MaxOver5x5_3x->Draw("SAME E1");
        lg->Draw();
        saveCanvas(c1, "comparison_E2x5MaxOver5x5_" + det, false);
        saveCanvas(c1, "comparison_E2x5MaxOver5x5_" + det, true);

        h1_ecalIso_2x->Scale(h1_ecalIso_3x->Integral()/h1_ecalIso_2x->Integral());
        h1_ecalIso_2x->GetYaxis()->SetRangeUser(0.01, histMax(h1_ecalIso_2x, h1_ecalIso_3x));
        h1_ecalIso_2x->Draw("HIST");
        h1_ecalIso_3x->Draw("SAME E1");
        lg->Draw();
        saveCanvas(c1, "comparison_ecalIso_" + det, false);
        saveCanvas(c1, "comparison_ecalIso_" + det, true);

        h1_hcalIso_2x->Scale(h1_hcalIso_3x->Integral()/h1_hcalIso_2x->Integral());
        h1_hcalIso_2x->GetYaxis()->SetRangeUser(0.01, histMax(h1_hcalIso_2x, h1_hcalIso_3x));
        h1_hcalIso_2x->Draw("HIST");
        h1_hcalIso_3x->Draw("SAME E1");
        lg->Draw();
        saveCanvas(c1, "comparison_hcalIso_" + det, false);
        saveCanvas(c1, "comparison_hcalIso_" + det, true); 

        h1_tkIso_2x->Scale(h1_tkIso_3x->Integral()/h1_tkIso_2x->Integral());
        h1_tkIso_2x->GetYaxis()->SetRangeUser(0.01, histMax(h1_tkIso_2x, h1_tkIso_3x));
        h1_tkIso_2x->Draw("HIST");
        h1_tkIso_3x->Draw("SAME E1");
        lg->Draw();
        saveCanvas(c1, "comparison_tkIso_" + det, false);
        saveCanvas(c1, "comparison_tkIso_" + det, true); 


}

