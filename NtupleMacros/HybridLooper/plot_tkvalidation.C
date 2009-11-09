
{

	gROOT->ProcessLine(".L ~/tdrStyle.C");
	gROOT->ProcessLine("setTDRStyle()");

	Int_t rebin = 2;

	gROOT->ProcessLine(".x results/IDStudy_sb_lintkIso03All_ee.C");
	s_tkIso03All_ee->SetFillColor(kBlue);
        s_tkIso03All_ee->SetLineStyle(1);
        s_tkIso03All_ee->Rebin(rebin);
        gROOT->ProcessLine(".x results/IDStudy_sb_lintkIso03AllRe_ee.C");
	s_tkIso03AllRe_ee->SetMarkerColor(kBlack);
	s_tkIso03AllRe_ee->SetMarkerStyle(22);
	s_tkIso03AllRe_ee->SetLineStyle(1);
        s_tkIso03AllRe_ee->Rebin(rebin);
        gROOT->ProcessLine(".x results/IDStudy_sb_lintkIso03AllReJura01In015_ee.C");
        s_tkIso03AllReJura01In015_ee->SetMarkerColor(kBlack);
        s_tkIso03AllReJura01In015_ee->SetMarkerStyle(24);
        s_tkIso03AllReJura01In015_ee->SetLineStyle(1);
        s_tkIso03AllReJura01In015_ee->SetLineColor(kRed);
        s_tkIso03AllReJura01In015_ee->Rebin(rebin);

        gROOT->ProcessLine(".x results/IDStudy_sb_lintkIso03All_eb.C");
        s_tkIso03All_eb->SetFillColor(kBlue);
        s_tkIso03All_eb->SetLineStyle(1);
	s_tkIso03All_eb->Rebin(rebin);
        gROOT->ProcessLine(".x results/IDStudy_sb_lintkIso03AllRe_eb.C");
        s_tkIso03AllRe_eb->SetMarkerColor(kBlack);
        s_tkIso03AllRe_eb->SetMarkerStyle(22);
        s_tkIso03AllRe_eb->SetLineStyle(1);
        s_tkIso03AllRe_eb->Rebin(rebin);
        gROOT->ProcessLine(".x results/IDStudy_sb_lintkIso03AllReJura01In015_eb.C");
        s_tkIso03AllReJura01In015_eb->SetMarkerColor(kBlack);
        s_tkIso03AllReJura01In015_eb->SetMarkerStyle(24);
        s_tkIso03AllReJura01In015_eb->SetLineStyle(1);
        s_tkIso03AllReJura01In015_eb->SetLineColor(kRed);
        s_tkIso03AllReJura01In015_eb->Rebin(rebin);

	TLegend *lg = new TLegend(0.4, 0.6, 0.9, 0.9);
        lg->SetFillColor(kWhite);
        lg->SetLineColor(kWhite);
        lg->SetShadowColor(kWhite);

	TCanvas *c1 = new TCanvas();
	c1->cd();
	c1->SetLogy();

        lg->AddEntry(s_tkIso03All_ee, "Standard (EE)", "f");
        lg->AddEntry(s_tkIso03AllRe_ee, "Recomputed (EE)", "lp");
        lg->AddEntry(s_tkIso03AllReJura01In015_ee, "#Delta#eta #pm 0.1, #DeltaR > 0.015 (EE)", "l");
	s_tkIso03All_ee->Draw("HIST");
	s_tkIso03AllRe_ee->Draw("SAME E1");
        s_tkIso03AllReJura01In015_ee->Draw("SAME");
	lg->Draw();
	c1->SaveAs("results/tkvalidation_ee.png");

	lg->Clear();
        lg->AddEntry(s_tkIso03All_eb, "Standard (EB)", "f");
        lg->AddEntry(s_tkIso03AllRe_eb, "Recomputed (EB)", "lp");
        lg->AddEntry(s_tkIso03AllReJura01In015_ee, "#Delta#eta #pm 0.1, #DeltaR > 0.015 (EB)", "l");
        s_tkIso03All_eb->Draw("HIST");
        s_tkIso03AllRe_eb->Draw("SAME E1");
        s_tkIso03AllReJura01In015_eb->Draw("SAME");
        lg->Draw();
        c1->SaveAs("results/tkvalidation_eb.png");





}


