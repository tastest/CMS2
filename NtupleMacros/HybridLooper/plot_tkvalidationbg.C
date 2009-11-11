
{

	gROOT->ProcessLine(".L ~/tdrStyle.C");
	gROOT->ProcessLine("setTDRStyle()");

	Int_t rebin = 2;

	gROOT->ProcessLine(".x results/IDStudy_sb_lintkIso03AllNM1_ee.C");
	b_tkIso03AllNM1_ee->SetFillColor(kGreen);
        b_tkIso03AllNM1_ee->SetLineStyle(1);
        b_tkIso03AllNM1_ee->Rebin(rebin);
//        gROOT->ProcessLine(".x results/IDStudy_sb_lintkIso03AllRe_ee.C");
//	b_tkIso03AllRe_ee->SetMarkerColor(kBlack);
//	b_tkIso03AllRe_ee->SetMarkerStyle(22);
//	b_tkIso03AllRe_ee->SetLineStyle(1);
 //       b_tkIso03AllRe_ee->Rebin(rebin);

        gROOT->ProcessLine(".x results/IDStudy_sb_lintkIso03AllReJura01In015NM1_ee.C");
        b_tkIso03AllReJura01In015NM1_ee->SetMarkerColor(kBlack);
        b_tkIso03AllReJura01In015NM1_ee->SetMarkerStyle(24);
        b_tkIso03AllReJura01In015NM1_ee->SetLineStyle(3);
        b_tkIso03AllReJura01In015NM1_ee->SetLineColor(kRed);
        b_tkIso03AllReJura01In015NM1_ee->Rebin(rebin);

        gROOT->ProcessLine(".x results/IDStudy_sb_lintkIso03AllReJura01In000NM1_ee.C");
        b_tkIso03AllReJura01In000NM1_ee->SetMarkerColor(kBlack);
        b_tkIso03AllReJura01In000NM1_ee->SetMarkerStyle(24);
        b_tkIso03AllReJura01In000NM1_ee->SetLineStyle(kDashed);
        b_tkIso03AllReJura01In000NM1_ee->SetLineColor(kBlack);
        b_tkIso03AllReJura01In000NM1_ee->Rebin(rebin);

        gROOT->ProcessLine(".x results/IDStudy_sb_lintkIso03AllNM1_eb.C");
        b_tkIso03AllNM1_eb->SetFillColor(kGreen);
        b_tkIso03AllNM1_eb->SetLineStyle(1);
	b_tkIso03AllNM1_eb->Rebin(rebin);
//        gROOT->ProcessLine(".x results/IDStudy_sb_lintkIso03AllRe_eb.C");
//        b_tkIso03AllRe_eb->SetMarkerColor(kBlack);
//        b_tkIso03AllRe_eb->SetMarkerStyle(22);
//        b_tkIso03AllRe_eb->SetLineStyle(1);
//        b_tkIso03AllRe_eb->Rebin(rebin);

        gROOT->ProcessLine(".x results/IDStudy_sb_lintkIso03AllReJura01In015NM1_eb.C");
        b_tkIso03AllReJura01In015NM1_eb->SetMarkerColor(kBlack);
        b_tkIso03AllReJura01In015NM1_eb->SetMarkerStyle(24);
        b_tkIso03AllReJura01In015NM1_eb->SetLineStyle(3);
        b_tkIso03AllReJura01In015NM1_eb->SetLineColor(kRed);
        b_tkIso03AllReJura01In015NM1_eb->Rebin(rebin);

        gROOT->ProcessLine(".x results/IDStudy_sb_lintkIso03AllReJura01In000NM1_eb.C");
        b_tkIso03AllReJura01In000NM1_eb->SetMarkerColor(kBlack);
        b_tkIso03AllReJura01In000NM1_eb->SetMarkerStyle(24);
        b_tkIso03AllReJura01In000NM1_eb->SetLineStyle(kDashed);
        b_tkIso03AllReJura01In000NM1_eb->SetLineColor(kBlack);
        b_tkIso03AllReJura01In000NM1_eb->Rebin(rebin);

	TLegend *lg = new TLegend(0.4, 0.6, 0.9, 0.9);
        lg->SetFillColor(kWhite);
        lg->SetLineColor(kWhite);
        lg->SetShadowColor(kWhite);

	TCanvas *c1 = new TCanvas();
	c1->cd();

        lg->AddEntry(b_tkIso03AllNM1_ee, "Standard NM1 (EE)", "f");
//        lg->AddEntry(b_tkIso03AllRe_ee, "Recomputed (EE)", "lp");
        lg->AddEntry(b_tkIso03AllReJura01In015NM1_ee, "#Delta#eta #pm 0.01, #DeltaR > 0.015 NM1 (EE)", "l");
        lg->AddEntry(b_tkIso03AllReJura01In000NM1_ee, "#Delta#eta #pm 0.01, #DeltaR > 0.000 NM1 (EE)", "l");
	b_tkIso03AllNM1_ee->Draw("HIST");
//	b_tkIso03AllRe_ee->Draw("SAME E1");
        b_tkIso03AllReJura01In000NM1_ee->Draw("SAME");
        b_tkIso03AllReJura01In015NM1_ee->Draw("SAME");
	lg->Draw();
	c1->SaveAs("results/tkvalidationsb000NM1_ee.png");

	lg->Clear();
        lg->AddEntry(b_tkIso03AllNM1_eb, "Standard NM1 (EB)", "f");
//        lg->AddEntry(b_tkIso03AllRe_eb, "Recomputed (EB)", "lp");
        lg->AddEntry(b_tkIso03AllReJura01In015NM1_eb, "#Delta#eta #pm 0.01, #DeltaR > 0.015 NM1 (EB)", "l");
        lg->AddEntry(b_tkIso03AllReJura01In000NM1_ee, "#Delta#eta #pm 0.01, #DeltaR > 0.000 NM1 (EB)", "l");
        b_tkIso03AllNM1_eb->Draw("HIST");
//        b_tkIso03AllRe_eb->Draw("SAME E1");
        b_tkIso03AllReJura01In000NM1_eb->Draw("SAME");
        b_tkIso03AllReJura01In015NM1_eb->Draw("SAME");
        lg->Draw();
        c1->SaveAs("results/tkvalidationsb000NM1_eb.png");

}

