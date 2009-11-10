
{

	gROOT->ProcessLine(".L ~/tdrStyle.C");
	gROOT->ProcessLine("setTDRStyle()");

	Int_t rebin = 2;

	gROOT->ProcessLine(".x results/IDStudy_s_lintkIso03AllNM1_ee.C");
	s_tkIso03AllNM1_ee->SetFillColor(kBlue);
        s_tkIso03AllNM1_ee->SetLineStyle(1);
        s_tkIso03AllNM1_ee->Rebin(rebin);
//        gROOT->ProcessLine(".x results/IDStudy_s_lintkIso03AllRe_ee.C");
//	s_tkIso03AllRe_ee->SetMarkerColor(kBlack);
//	s_tkIso03AllRe_ee->SetMarkerStyle(22);
//	s_tkIso03AllRe_ee->SetLineStyle(1);
 //       s_tkIso03AllRe_ee->Rebin(rebin);

        gROOT->ProcessLine(".x results/IDStudy_s_lintkIso03AllReJura01In015NM1_ee.C");
        s_tkIso03AllReJura01In015NM1_ee->SetMarkerColor(kBlack);
        s_tkIso03AllReJura01In015NM1_ee->SetMarkerStyle(24);
        s_tkIso03AllReJura01In015NM1_ee->SetLineStyle(3);
        s_tkIso03AllReJura01In015NM1_ee->SetLineColor(kRed);
        s_tkIso03AllReJura01In015NM1_ee->Rebin(rebin);

        gROOT->ProcessLine(".x results/IDStudy_s_lintkIso03AllReJura01In000NM1_ee.C");
        s_tkIso03AllReJura01In000NM1_ee->SetMarkerColor(kBlack);
        s_tkIso03AllReJura01In000NM1_ee->SetMarkerStyle(24);
        s_tkIso03AllReJura01In000NM1_ee->SetLineStyle(kDashed);
        s_tkIso03AllReJura01In000NM1_ee->SetLineColor(kBlack);
        s_tkIso03AllReJura01In000NM1_ee->Rebin(rebin);

        gROOT->ProcessLine(".x results/IDStudy_s_lintkIso03AllNM1_eb.C");
        s_tkIso03AllNM1_eb->SetFillColor(kBlue);
        s_tkIso03AllNM1_eb->SetLineStyle(1);
	s_tkIso03AllNM1_eb->Rebin(rebin);
//        gROOT->ProcessLine(".x results/IDStudy_s_lintkIso03AllRe_eb.C");
//        s_tkIso03AllRe_eb->SetMarkerColor(kBlack);
//        s_tkIso03AllRe_eb->SetMarkerStyle(22);
//        s_tkIso03AllRe_eb->SetLineStyle(1);
//        s_tkIso03AllRe_eb->Rebin(rebin);

        gROOT->ProcessLine(".x results/IDStudy_s_lintkIso03AllReJura01In015NM1_eb.C");
        s_tkIso03AllReJura01In015NM1_eb->SetMarkerColor(kBlack);
        s_tkIso03AllReJura01In015NM1_eb->SetMarkerStyle(24);
        s_tkIso03AllReJura01In015NM1_eb->SetLineStyle(3);
        s_tkIso03AllReJura01In015NM1_eb->SetLineColor(kRed);
        s_tkIso03AllReJura01In015NM1_eb->Rebin(rebin);

        gROOT->ProcessLine(".x results/IDStudy_s_lintkIso03AllReJura01In000NM1_eb.C");
        s_tkIso03AllReJura01In000NM1_eb->SetMarkerColor(kBlack);
        s_tkIso03AllReJura01In000NM1_eb->SetMarkerStyle(24);
        s_tkIso03AllReJura01In000NM1_eb->SetLineStyle(kDashed);
        s_tkIso03AllReJura01In000NM1_eb->SetLineColor(kBlack);
        s_tkIso03AllReJura01In000NM1_eb->Rebin(rebin);

	TLegend *lg = new TLegend(0.4, 0.6, 0.9, 0.9);
        lg->SetFillColor(kWhite);
        lg->SetLineColor(kWhite);
        lg->SetShadowColor(kWhite);

	TCanvas *c1 = new TCanvas();
	c1->cd();
	c1->SetLogy();

        lg->AddEntry(s_tkIso03AllNM1_ee, "Standard NM1 (EE)", "f");
//        lg->AddEntry(s_tkIso03AllRe_ee, "Recomputed (EE)", "lp");
        lg->AddEntry(s_tkIso03AllReJura01In015NM1_ee, "#Delta#eta #pm 0.01, #DeltaR > 0.015 NM1 (EE)", "l");
        lg->AddEntry(s_tkIso03AllReJura01In000NM1_ee, "#Delta#eta #pm 0.01, #DeltaR > 0.000 NM1 (EE)", "l");
	s_tkIso03AllNM1_ee->Draw("HIST");
//	s_tkIso03AllRe_ee->Draw("SAME E1");
        s_tkIso03AllReJura01In000NM1_ee->Draw("SAME");
        s_tkIso03AllReJura01In015NM1_ee->Draw("SAME");
	lg->Draw();
	c1->SaveAs("results/tkvalidation000NM1_ee.png");

	lg->Clear();
        lg->AddEntry(s_tkIso03AllNM1_eb, "Standard NM1 (EB)", "f");
//        lg->AddEntry(s_tkIso03AllRe_eb, "Recomputed (EB)", "lp");
        lg->AddEntry(s_tkIso03AllReJura01In015NM1_eb, "#Delta#eta #pm 0.01, #DeltaR > 0.015 NM1 (EB)", "l");
        lg->AddEntry(s_tkIso03AllReJura01In000NM1_ee, "#Delta#eta #pm 0.01, #DeltaR > 0.000 NM1 (EB)", "l");
        s_tkIso03AllNM1_eb->Draw("HIST");
//        s_tkIso03AllRe_eb->Draw("SAME E1");
        s_tkIso03AllReJura01In000NM1_eb->Draw("SAME");
        s_tkIso03AllReJura01In015NM1_eb->Draw("SAME");
        lg->Draw();
        c1->SaveAs("results/tkvalidation000NM1_eb.png");

}

