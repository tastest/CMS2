
plot_tkiso_comparison(TString det)
{

	gROOT->ProcessLine(".L ~/tdrStyle.C");
	gROOT->ProcessLine("setTDRStyle()");
	TString DET = det;
	DET = DET.ToUpper();

	gROOT->ProcessLine(".x results/IDStudy_eff_tkIso03All_" + det + ".C");
	TGraph *gr_tkIso03All = (TGraph*)graph->Clone("gr_tkIso03All");
	gr_tkIso03All->SetMarkerColor(kBlack);
	gr_tkIso03All->SetMarkerStyle(22);

/*
        gROOT->ProcessLine(".x results/IDStudy_eff_tkIso03AllReJura01_" + det + ".C");
        TGraph *gr_tkIso03AllReJura01 = (TGraph*)graph->Clone("gr_tkIso03AllReJura01");
        gr_tkIso03AllReJura01->SetMarkerColor(kRed);
	gr_tkIso03AllReJura01->SetMarkerStyle(23);

        gROOT->ProcessLine(".x results/IDStudy_eff_tkIso03AllReJura02_" + det + ".C");
        TGraph *gr_tkIso03AllReJura02 = (TGraph*)graph->Clone("gr_tkIso03AllReJura02");
        gr_tkIso03AllReJura02->SetMarkerColor(kGreen);
        gr_tkIso03AllReJura02->SetMarkerStyle(23);

        gROOT->ProcessLine(".x results/IDStudy_eff_tkIso03AllReJura03_" + det + ".C");
        TGraph *gr_tkIso03AllReJura03 = (TGraph*)graph->Clone("gr_tkIso03AllReJura03");
        gr_tkIso03AllReJura03->SetMarkerColor(kBlue);
        gr_tkIso03AllReJura03->SetMarkerStyle(23);

        gROOT->ProcessLine(".x results/IDStudy_eff_tkIso03AllReJura01In015_" + det + ".C");
        TGraph *gr_tkIso03AllReJura01In015 = (TGraph*)graph->Clone("gr_tkIso03AllReJura01In015");
        gr_tkIso03AllReJura01In015->SetMarkerColor(kMagenta);
        gr_tkIso03AllReJura01In015->SetMarkerStyle(25);

*/

	//
	//

        gROOT->ProcessLine(".x results/IDStudy_eff_wwIsoAll_" + det + ".C");
        TGraph *gr_wwIsoAll = (TGraph*)graph->Clone("gr_wwIsoAll");
        gr_wwIsoAll->SetMarkerColor(kBlack);
        gr_wwIsoAll->SetMarkerStyle(22);

        gROOT->ProcessLine(".x results/IDStudy_eff_wwIsoV1All_" + det + ".C");
        TGraph *gr_wwIsoV1All = (TGraph*)graph->Clone("gr_wwIsoV1All");
        gr_wwIsoV1All->SetMarkerColor(kMagenta);
        gr_wwIsoV1All->SetMarkerStyle(25);

	//
	//

        gROOT->ProcessLine(".x results/IDStudy_eff_tkIso03AllNM1_" + det + ".C");
        TGraph *gr_tkIso03AllNM1 = (TGraph*)graph->Clone("gr_tkIso03AllNM1");
        gr_tkIso03AllNM1->SetMarkerColor(kBlack);
        gr_tkIso03AllNM1->SetMarkerStyle(22);

        gROOT->ProcessLine(".x results/IDStudy_eff_tkIso03AllConvNM1_" + det + ".C");
        TGraph *gr_tkIso03AllConvNM1 = (TGraph*)graph->Clone("gr_tkIso03AllConvNM1");
        gr_tkIso03AllConvNM1->SetMarkerColor(kGreen);
        gr_tkIso03AllConvNM1->SetMarkerStyle(22);

        gROOT->ProcessLine(".x results/IDStudy_eff_tkIso03AllIDNM1_" + det + ".C");
        TGraph *gr_tkIso03AllIDNM1 = (TGraph*)graph->Clone("gr_tkIso03AllIDNM1");
        gr_tkIso03AllIDNM1->SetMarkerColor(kRed);
        gr_tkIso03AllIDNM1->SetMarkerStyle(22);

        gROOT->ProcessLine(".x results/IDStudy_eff_tkIso03AllConvIDNM1_" + det + ".C");
        TGraph *gr_tkIso03AllConvIDNM1 = (TGraph*)graph->Clone("gr_tkIso03AllConvIDNM1");
        gr_tkIso03AllConvIDNM1->SetMarkerColor(kBlue);
        gr_tkIso03AllConvIDNM1->SetMarkerStyle(22);

        gROOT->ProcessLine(".x results/IDStudy_eff_tkIso03AllReJura01In015NM1_" + det + ".C");
        TGraph *gr_tkIso03AllReJura01In015NM1 = (TGraph*)graph->Clone("gr_tkIso03AllReJura01In015NM1");
        gr_tkIso03AllReJura01In015NM1->SetMarkerColor(kMagenta);
        gr_tkIso03AllReJura01In015NM1->SetMarkerStyle(25);

        gROOT->ProcessLine(".x results/IDStudy_eff_tkIso03AllReJura01In000NM1_" + det + ".C");
        TGraph *gr_tkIso03AllReJura01In000NM1 = (TGraph*)graph->Clone("gr_tkIso03AllReJura01In000NM1");
        gr_tkIso03AllReJura01In000NM1->SetMarkerColor(kBlack);
        gr_tkIso03AllReJura01In000NM1->SetMarkerStyle(26);

        gROOT->ProcessLine(".x results/IDStudy_eff_tkIso03AllReJura01In015IDNM1_" + det + ".C");
        TGraph *gr_tkIso03AllReJura01In015IDNM1 = (TGraph*)graph->Clone("gr_tkIso03AllReJura01In015IDNM1");
        gr_tkIso03AllReJura01In015IDNM1->SetMarkerColor(kRed);
        gr_tkIso03AllReJura01In015IDNM1->SetMarkerStyle(25);

        gROOT->ProcessLine(".x results/IDStudy_eff_tkIso03AllReJura01In015ConvNM1_" + det + ".C");
        TGraph *gr_tkIso03AllReJura01In015ConvNM1 = (TGraph*)graph->Clone("gr_tkIso03AllReJura01In015ConvNM1");
        gr_tkIso03AllReJura01In015ConvNM1->SetMarkerColor(kGreen);
        gr_tkIso03AllReJura01In015ConvNM1->SetMarkerStyle(25);

        gROOT->ProcessLine(".x results/IDStudy_eff_tkIso03AllReJura01In015ConvIDNM1_" + det + ".C");
        TGraph *gr_tkIso03AllReJura01In015ConvIDNM1 = (TGraph*)graph->Clone("gr_tkIso03AllReJura01In015ConvIDNM1");
        gr_tkIso03AllReJura01In015ConvIDNM1->SetMarkerColor(kBlue);
        gr_tkIso03AllReJura01In015ConvIDNM1->SetMarkerStyle(25);

        gROOT->ProcessLine(".x results/IDStudy_eff_tkIso03AllReShCutNM1_" + det + ".C");
        TGraph *gr_tkIso03AllReShCutNM1 = (TGraph*)graph->Clone("gr_tkIso03AllReShCutNM1");
        gr_tkIso03AllReShCutNM1->SetMarkerColor(kCyan);
        gr_tkIso03AllReShCutNM1->SetMarkerStyle(30);

	TCanvas *c1 = new TCanvas();
	c1->cd();

/*
        TLegend *lg = new TLegend(0.2, 0.4, 0.7, 0.9);
        lg->SetFillColor(kWhite);
        lg->SetLineColor(kWhite);
        lg->SetShadowColor(kWhite);
        lg->AddEntry(gr_tkIso03All, "Standard (" + DET + ")", "lp");
        lg->AddEntry(gr_tkIso03AllReJura01, "#Delta#eta #pm 0.01, #DeltaR > 0.04 (" + DET + ")", "lp");
        lg->AddEntry(gr_tkIso03AllReJura02, "#Delta#eta #pm 0.02, #DeltaR > 0.04 (" + DET + ")", "lp");
        lg->AddEntry(gr_tkIso03AllReJura03, "#Delta#eta #pm 0.03, #DeltaR > 0.04 (" + DET + ")", "lp");
        lg->AddEntry(gr_tkIso03AllReJura01In015, "#Delta#eta #pm 0.01, #DeltaR > 0.015 (" + DET + ")", "lp");
	gr_tkIso03All->GetXaxis()->SetRangeUser(0.0, 1.05);
	gr_tkIso03All->GetYaxis()->SetRangeUser(0.0, 1.05);
	gr_tkIso03All->Draw("AP");
	gr_tkIso03AllReJura01->Draw("P");
        gr_tkIso03AllReJura02->Draw("P");
        gr_tkIso03AllReJura03->Draw("P");
        gr_tkIso03AllReJura01In015->Draw("P");
	lg->Draw();

	c1->SaveAs("results/effrej_tkIso03AllReJura_comparison_" + det + ".png");
*/

        TLegend *lg2 = new TLegend(0.2, 0.7, 0.7, 0.9);
        lg2->SetFillColor(kWhite);
        lg2->SetLineColor(kWhite);
        lg2->SetShadowColor(kWhite);
        lg2->AddEntry(gr_tkIso03AllNM1, "Standard N-1 (" + DET + ")", "lp");
        lg2->AddEntry(gr_tkIso03AllReJura01In015NM1, "#Delta#eta #pm 0.01, #DeltaR > 0.015 N-1 (" + DET + ")", "lp");
        gr_tkIso03AllNM1->GetXaxis()->SetRangeUser(0.0, 1.05);
        gr_tkIso03AllNM1->GetYaxis()->SetRangeUser(0.0, 1.05);
        gr_tkIso03AllNM1->Draw("AP");
        gr_tkIso03AllReJura01In015NM1->Draw("P");
        lg2->Draw();

        c1->SaveAs("results/effrej_tkIso03AllNM1_comparison_" + det + ".png");

        TLegend *lg2_2 = new TLegend(0.2, 0.7, 0.7, 0.9);
        lg2_2->SetFillColor(kWhite);
        lg2_2->SetLineColor(kWhite);
        lg2_2->SetShadowColor(kWhite);
        lg2_2->AddEntry(gr_tkIso03AllNM1, "Standard N-1 (" + DET + ")", "lp");
        lg2_2->AddEntry(gr_tkIso03AllReJura01In015NM1, "#Delta#eta #pm 0.01, #DeltaR > 0.015 N-1 (" + DET + ")", "lp");
        lg2_2->AddEntry(gr_tkIso03AllReJura01In000NM1, "#Delta#eta #pm 0.01 N-1 (" + DET + ")", "lp");
        gr_tkIso03AllNM1->GetXaxis()->SetRangeUser(0.0, 1.05);
        gr_tkIso03AllNM1->GetYaxis()->SetRangeUser(0.0, 1.05);
        gr_tkIso03AllNM1->Draw("AP");
        gr_tkIso03AllReJura01In015NM1->Draw("P");
        gr_tkIso03AllReJura01In000NM1->Draw("P");
        lg2_2->Draw();

        c1->SaveAs("results/effrej_tkIso03AllNM1_comparison2_" + det + ".png");
        c1->SaveAs("results/effrej_tkIso03AllNM1_comparison2_" + det + ".root");

	c1->Clear();
	c1->Divide(2, 1);
	c1->cd(1);
        TLegend *lg3 = new TLegend(0.2, 0.2, 0.7, 0.9);
        lg3->SetFillColor(kWhite);
	lg3->SetFillStyle(0);
        lg3->SetLineColor(kWhite);
        lg3->SetShadowColor(kWhite);
        lg3->AddEntry(gr_tkIso03AllNM1, "tk N-1 (" + DET + ")", "lp");
        lg3->AddEntry(gr_tkIso03AllConvNM1, "tk N-1 + !conv (" + DET + ")", "lp");
        lg3->AddEntry(gr_tkIso03AllIDNM1, "tk N-1 + tightID (" + DET + ")", "lp");
        lg3->AddEntry(gr_tkIso03AllConvIDNM1, "tk N-1 + !conv + tightID (" + DET + ")", "lp");
        gr_tkIso03AllNM1->GetXaxis()->SetRangeUser(0.0, 1.05);
        gr_tkIso03AllNM1->GetYaxis()->SetRangeUser(0.0, 1.05);
        gr_tkIso03AllNM1->Draw("AP");
	gr_tkIso03AllConvNM1->Draw("P");
        gr_tkIso03AllIDNM1->Draw("P");
        gr_tkIso03AllConvIDNM1->Draw("P");
        lg3->Draw();

	c1->cd(2);
        TLegend *lg4 = new TLegend(0.2, 0.2, 0.7, 0.9);
        lg4->SetFillColor(kWhite);
        lg4->SetFillStyle(0);
        lg4->SetLineColor(kWhite);
        lg4->SetShadowColor(kWhite);
        lg4->AddEntry(gr_tkIso03AllReJura01In015NM1, "tk015-01 N-1 (" + DET + ")", "lp");
        lg4->AddEntry(gr_tkIso03AllReJura01In015ConvNM1, "tk015-01 N-1 + !conv (" + DET + ")", "lp");
        lg4->AddEntry(gr_tkIso03AllReJura01In015IDNM1, "tk015-01 N-1 + tightID (" + DET + ")", "lp");
        lg4->AddEntry(gr_tkIso03AllReJura01In015ConvIDNM1, "tk015-01 N-1 + tightID && !isconv (" + DET + ")", "lp");
        gr_tkIso03AllReJura01In015NM1->GetXaxis()->SetRangeUser(0.0, 1.05);
        gr_tkIso03AllReJura01In015NM1->GetYaxis()->SetRangeUser(0.0, 1.05);
        gr_tkIso03AllReJura01In015NM1->Draw("AP");
        gr_tkIso03AllReJura01In015ConvNM1->Draw("P");
        gr_tkIso03AllReJura01In015IDNM1->Draw("P");
        gr_tkIso03AllReJura01In015ConvIDNM1->Draw("P");
        lg4->Draw();

        c1->SaveAs("results/effrej_tkIso03AllConvIDTogetherNM1_comparison_" + det + ".png");

	c1->Divide(1,1);
        c1->Clear();
	c1->cd();
        TLegend *lg5 = new TLegend(0.2, 0.7, 0.75, 0.9);
        lg5->SetFillColor(kWhite);
        lg5->SetFillStyle(0);
        lg5->SetLineColor(kWhite);
        lg5->SetShadowColor(kWhite);
        lg5->AddEntry(gr_tkIso03AllConvIDNM1, "tk N-1 + tightID + !conv (" + DET + ")", "lp");
        lg5->AddEntry(gr_tkIso03AllReJura01In015ConvIDNM1, "tk015-01 N-1 + tightID && !isconv (" + DET + ")", "lp");
        gr_tkIso03AllConvIDNM1->GetXaxis()->SetRangeUser(0.0, 1.05);
        gr_tkIso03AllConvIDNM1->GetYaxis()->SetRangeUser(0.0, 1.05);
        gr_tkIso03AllConvIDNM1->Draw("AP");
        gr_tkIso03AllReJura01In015ConvIDNM1->Draw("P");
        lg5->Draw();

        c1->SaveAs("results/effrej_tkIso03AllConvIDNM1_comparison_" + det + ".png");

        c1->Divide(1,1);
        c1->Clear();
        c1->cd();
        TLegend *lg6 = new TLegend(0.2, 0.7, 0.75, 0.9);
        lg6->SetFillColor(kWhite);
        lg6->SetFillStyle(0);
        lg6->SetLineColor(kWhite);
        lg6->SetShadowColor(kWhite);
        lg6->AddEntry(gr_tkIso03AllConvNM1, "tk N-1 + !conv (" + DET + ")", "lp");
        lg6->AddEntry(gr_tkIso03AllReJura01In015ConvNM1, "tk015-01 N-1 +!conv (" + DET + ")", "lp");
        gr_tkIso03AllConvNM1->GetXaxis()->SetRangeUser(0.0, 1.05);
        gr_tkIso03AllConvNM1->GetYaxis()->SetRangeUser(0.0, 1.05);
        gr_tkIso03AllConvNM1->Draw("AP");
        gr_tkIso03AllReJura01In015ConvNM1->Draw("P");
        lg6->Draw();

        c1->SaveAs("results/effrej_tkIso03AllConvNM1_comparison_" + det + ".png");


        c1->Divide(1,1);
        c1->Clear();
        c1->cd();
        TLegend *lg7 = new TLegend(0.2, 0.7, 0.75, 0.9);
        lg7->SetFillColor(kWhite);
        lg7->SetFillStyle(0);
        lg7->SetLineColor(kWhite);
        lg7->SetShadowColor(kWhite);
        lg7->AddEntry(gr_tkIso03AllIDNM1, "tk N-1 + tightID (" + DET + ")", "lp");
        lg7->AddEntry(gr_tkIso03AllReJura01In015IDNM1, "tk015-01 N-1 + tightID (" + DET + ")", "lp");
        gr_tkIso03AllIDNM1->GetXaxis()->SetRangeUser(0.0, 1.05);
        gr_tkIso03AllIDNM1->GetYaxis()->SetRangeUser(0.0, 1.05);
        gr_tkIso03AllIDNM1->Draw("AP");
        gr_tkIso03AllReJura01In015IDNM1->Draw("P");
        lg7->Draw();

        c1->SaveAs("results/effrej_tkIso03AllIDNM1_comparison_" + det + ".png");

        c1->Divide(1,1);
        c1->Clear();
        c1->cd();
        TLegend *lg8 = new TLegend(0.2, 0.7, 0.75, 0.9);
        lg8->SetFillColor(kWhite);
        lg8->SetFillStyle(0);
        lg8->SetLineColor(kWhite);
        lg8->SetShadowColor(kWhite);
        lg8->AddEntry(gr_tkIso03AllReJura01In015NM1, "tk015-01 N-1 (" + DET + ")", "lp");
        lg8->AddEntry(gr_tkIso03AllReShCutNM1, "ShCut45 N-1 (" + DET + ")", "lp");
        gr_tkIso03AllReJura01In015NM1->GetXaxis()->SetRangeUser(0.0, 1.05);
        gr_tkIso03AllReJura01In015NM1->GetYaxis()->SetRangeUser(0.0, 1.05);
        gr_tkIso03AllReJura01In015NM1->Draw("AP");
        gr_tkIso03AllReShCutNM1->Draw("P");
        lg8->Draw();

        c1->SaveAs("results/effrej_tkIso03AllShCut_comparison_" + det + ".png");

	c1->Clear();
	c1->cd();
        TLegend *lg9 = new TLegend(0.2, 0.7, 0.75, 0.9);
        lg9->SetFillColor(kWhite);
        lg9->SetFillStyle(0);
        lg9->SetLineColor(kWhite);
        lg9->SetShadowColor(kWhite);
        lg9->AddEntry(gr_wwIsoAll, "RelIso (" + DET + ")", "lp");
        lg9->AddEntry(gr_wwIsoV1All, "RelIso(Jura 0.015, 0.01) (" + DET + ")", "lp"); 
	gr_wwIsoAll->Draw("AP");
	gr_wwIsoAll->GetYaxis()->SetRangeUser(0, 0.8);
	gr_wwIsoAll->GetXaxis()->SetRangeUser(0.9, 1.0);
	gr_wwIsoV1All->Draw("P");
	lg9->Draw();

        c1->SaveAs("results/effrej_relIsoComparison_" + det + ".png");

}

