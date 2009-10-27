
#include "plotResults.h"
#include "Tools/HistogramUtilities.h"
#include "Tools/Utilities.h"
#include "Tools/histtools.cc"

#include "TROOT.h"

#include "TGraphAsymmErrors.h"
#include "TArrow.h"

        const static sources_t theSources_31X =
//		(1ll << H_ZZ);
                (1ll << H_WENU_7TeV) 	|       
		(1ll << H_QCD30_7TeV) 	|
		(1ll << H_PHOTONJET_7TeV);
//                (1ll << H_EM30_80)       |
//                (1ll << H_BC30_80);

        const static sources_t theSignal_31X =
//		(1ll << H_ZZ);
                (1ll << H_WENU_7TeV);

        const static sources_t theBackground_31X =
		(1ll << H_QCD30_7TeV);
//		(1ll << H_ZZ);
//                (1ll << H_EM30_80)       |
//                (1ll << H_BC30_80);

        const static sources_t theSources_22X =
                (1ll << H_QCD30) 	|
                (1ll << H_QCD80)        |
		(1ll << H_WJET_ALP);

	const static sources_t theSignal_22X = 
                (1ll << H_WJET_ALP);

	const static sources_t theBackground_22X = 
                (1ll << H_QCD30) 	|
                (1ll << H_QCD80);

	// for 2_2_X
	const static sources_t &theSignal = theSignal_31X;
	const static sources_t &theBackground = theBackground_31X;
	const static sources_t &theSources = theSources_31X;

TArrow *getArrow(THStack *st, TString det, float cutValEB, float cutValEE)
{

        float cutVal; 
        if (det == "eb") cutVal = cutValEB;
        else cutVal = cutValEE;
        TArrow *arr_cut = new TArrow(cutVal, st->GetMaximum()/2.0, cutVal, 0, 0.05, "|>");
        arr_cut->SetLineWidth(2.0);
	return arr_cut;
}

void plot2DSB(HistogramUtilities &h1, TString name, TString xTitle, TString yTitle, TString saveName, TString det)
{

        TH2F *h2_signal = h1.get2dHistogram(theSignal, name, "", det, 1);
	h2_signal->SetTitle(";" + xTitle + ";" + yTitle);
        TH2F *h2_background = h1.get2dHistogram(theBackground, name, "", det, 1);
        h2_background->SetTitle(";" + xTitle + ";" + yTitle);
	TCanvas *c1 = new TCanvas();
	c1->SetCanvasSize(600, 300);
	c1->Divide(2, 1);
	c1->cd(1);
	h2_signal->Draw("COLZ");
	c1->cd(2);
	h2_background->Draw("COLZ");
        Utilities::saveCanvas(c1, "results/" + saveName + "_2D_linz_" + name + "_" + det);

        c1->cd(1);
	c1->cd(1)->SetLogz();
        h2_signal->Draw("COLZ");
        c1->cd(2);
	c1->cd(2)->SetLogz();
        h2_background->Draw("COLZ");
        Utilities::saveCanvas(c1, "results/" + saveName + "_2D_logz_" + name + "_" + det);

	delete c1;
	delete h2_signal;
	delete h2_background;

}

void plotEff(HistogramUtilities &h1, TString name, TString saveName, TString det, bool ascending, int rebin, bool legendOnRight, float cutValEB, float cutValEE)
{

	std::cout << "[plotEff] " << name << std::endl;
	TH1F *h1_signal = h1.getHistogram(theSignal, name, "", det, rebin, "s_");
        TH1F *h1_background = h1.getHistogram(theBackground, name, "", det, rebin, "b_");
	TH1F *h1_both = h1.getHistogram(theBackground | theSignal, name, "", det, rebin, "");

	//
	// S/B plots
	//

	TH1F *h1_signal_cumulated = (TH1F*)(cumulate(*h1_signal, ascending).Clone());
        TH1F *h1_background_cumulated = (TH1F*)(cumulate(*h1_background, ascending).Clone());
	TH1F *h1_sob = (TH1F*)h1_signal_cumulated->Clone();
	float sTotal = h1_signal->Integral(0, h1_signal_cumulated->GetNbinsX() + 1);
        float bTotal = h1_background->Integral(0, h1_background_cumulated->GetNbinsX() + 1);

	//int bin_bg50 = 0;
	int bin_eff99 = 0;
	int bin_eff98 = 0;
	int bin_eff95 = 0;
	int bin_eff80 = 0;
	for (Int_t bin = 0; bin < h1_signal_cumulated->GetNbinsX() + 1; ++bin)
	{
		float s = h1_signal_cumulated->GetBinContent(bin);
		float b = h1_background_cumulated->GetBinContent(bin);
		float sob = 0;
		if (b > 0.0) { 
			sob = s/b;
		}
		h1_sob->SetBinContent(bin, sob);
		//std::cout << s/sTotal << std::endl;
		if (ascending) {
			if (s/sTotal >= 0.99 && bin_eff99 == 0) bin_eff99 = bin;
                	if (s/sTotal >= 0.98 && bin_eff98 == 0) bin_eff98 = bin;
                	if (s/sTotal >= 0.95 && bin_eff95 == 0) bin_eff95 = bin;
                	if (s/sTotal >= 0.80 && bin_eff80 == 0) bin_eff80 = bin;
		}
		else {
                        if (s/sTotal <= 0.99 && bin_eff99 == 0) bin_eff99 = bin;
                        if (s/sTotal <= 0.98 && bin_eff98 == 0) bin_eff98 = bin;
                        if (s/sTotal <= 0.95 && bin_eff95 == 0) bin_eff95 = bin;
                        if (s/sTotal <= 0.80 && bin_eff80 == 0) bin_eff80 = bin;
		}

		//if (b/bTotal >= 0.50 && bin_bg50 == 0) bin_bg50 = bin;
	}

        TArrow *arr_eff99 = new TArrow(h1_signal->GetBinCenter(bin_eff99), h1_both->GetMaximum()/2.0, h1_signal->GetBinCenter(bin_eff99), 0, 0.05, "|>");
	arr_eff99->SetLineColor(kGreen);
	arr_eff99->SetFillColor(kGreen);
	arr_eff99->SetLineWidth(2);
        TArrow *arr_99 = new TArrow(0.99, 1.0, 0.99, 0, 0.05, "|>");
        arr_99->SetLineColor(kGreen);
        arr_99->SetFillColor(kGreen);
        arr_99->SetLineWidth(2);

        TArrow *arr_eff98 = new TArrow(h1_signal->GetBinCenter(bin_eff98), h1_both->GetMaximum()/2.0, h1_signal->GetBinCenter(bin_eff98), 0, 0.05, "|>");
        arr_eff98->SetLineColor(kBlue);
        arr_eff98->SetFillColor(kBlue);
        arr_eff98->SetLineWidth(2);
        TArrow *arr_98 = new TArrow(0.98, 1.0, 0.98, 0, 0.05, "|>");
        arr_98->SetLineColor(kBlue);
        arr_98->SetFillColor(kBlue);
        arr_98->SetLineWidth(2);

        TArrow *arr_eff95 = new TArrow(h1_signal->GetBinCenter(bin_eff95), h1_both->GetMaximum()/2.0, h1_signal->GetBinCenter(bin_eff95), 0, 0.05, "|>");
        arr_eff95->SetLineColor(kRed);
        arr_eff95->SetFillColor(kRed);
        arr_eff95->SetLineWidth(2);
        TArrow *arr_95 = new TArrow(0.95, 1.0, 0.95, 0, 0.05, "|>");
        arr_95->SetLineColor(kRed);
        arr_95->SetFillColor(kRed);
        arr_95->SetLineWidth(2);

        TArrow *arr_eff80 = new TArrow(h1_signal->GetBinCenter(bin_eff80), h1_both->GetMaximum()/2.0, h1_signal->GetBinCenter(bin_eff80), 0, 0.05, "|>");
        arr_eff80->SetLineColor(kMagenta);
        arr_eff80->SetFillColor(kMagenta);
        arr_eff80->SetLineWidth(2);
        TArrow *arr_80 = new TArrow(0.80, 1.0, 0.80, 0, 0.05, "|>");
        arr_80->SetLineColor(kMagenta);
        arr_80->SetFillColor(kMagenta);
        arr_80->SetLineWidth(2);

        //TArrow *arr_bg50 = new TArrow(h1_background->GetBinCenter(bin_bg50), h1_both->GetMaximum()/2.0, h1_background->GetBinCenter(bin_bg50), 0, 0.05, "|>");
        //arr_bg50->SetLineColor(kBlack);
        //arr_bg50->SetFillColor(kBlack);
        //arr_bg50->SetLineWidth(2);

        TCanvas *c = new TCanvas();
	c->cd();	
	h1_sob->Draw();
        Utilities::saveCanvas(c, "results/" + saveName + "_sob_" + name + "_" + det);	

        //
        // Acceptance rejection curves
        //

        // latter bool is "ascending"
        bool normalise = true;
        TGraph *gr = (TGraph*)(eff_rej(*h1_signal, *h1_background, normalise, ascending).Clone());
        gr->SetTitle(name + ";Signal;Background");
        gr->SetMarkerStyle(23);

        c->cd();
        gr->Draw("AP");
        gr->GetXaxis()->SetRangeUser(0.80, 1.1);
        gr->GetYaxis()->SetRangeUser(0.00, 1.1); 
        arr_99->Draw();
        arr_98->Draw();
        arr_95->Draw();
        Utilities::saveCanvas(c, "results/" + saveName + "_eff_" + name + "_" + det);

        gr->GetXaxis()->SetRangeUser(0.00, 1.1);
        gr->GetYaxis()->SetRangeUser(0.00, 1.1);
        arr_80->Draw();
        Utilities::saveCanvas(c, "results/" + saveName + "_effzoomout_" + name + "_" + det);


	//
	// S and B overlays
	//

        h1_signal->SetLineWidth(2);
	h1_signal->SetLineColor(kBlue);
	h1_signal->SetLineStyle(kDashed);
        h1_background->SetLineWidth(2);
	h1_background->SetLineColor(kGreen);
	TLegend *lg = 0;
	if (! legendOnRight) lg = new TLegend(0.2, 0.75, 0.5, 0.9);
	else lg = new TLegend(0.55, 0.7, 0.9, 0.9);

	float cutVal;
	if (det == "eb") cutVal = cutValEB;
	else cutVal = cutValEE;
	TArrow *arr_cut = new TArrow(cutVal, h1_both->GetMaximum()/2.0, cutVal, 0, 0.05, "|>");
	arr_cut->SetLineWidth(2.0);

  	lg->SetFillColor(kWhite);
  	lg->SetLineColor(kWhite);
  	lg->SetShadowColor(kWhite);
	TString upperDet = det;
	upperDet.ToUpper();
	lg->AddEntry(h1_signal, "Signal (" + upperDet + ")", "l");
	lg->AddEntry(h1_background, "Background (" + upperDet + ")", "l");
	TCanvas *c_sb = new TCanvas();
	c_sb->cd();
	if (h1_signal->GetMaximum() > h1_background->GetMaximum()) {
		h1_signal->Draw();
		h1_signal->GetYaxis()->SetRangeUser(0, h1_signal->GetMaximum() 
							+ h1_signal->GetMaximum()*0.1);
		h1_background->Draw("SAME");
	}
	else {
        	h1_background->Draw();
                h1_background->GetYaxis()->SetRangeUser(0, h1_background->GetMaximum() 
								+ h1_background->GetMaximum()*0.1);
        	h1_signal->Draw("SAME");
	}
	if (cutVal > 0.0) arr_cut->Draw();
	if (bin_eff99 != 0) arr_eff99->Draw();
        if (bin_eff98 != 0) arr_eff98->Draw();
        if (bin_eff95 != 0) arr_eff95->Draw();
	if (bin_eff80 != 0) arr_eff80->Draw();
	//if (bin_bg50 != 0) arr_bg50->Draw();

	lg->Draw();
        Utilities::saveCanvas(c_sb, "results/" + saveName + "_sb_lin" + name + "_" + det);

	h1_signal->Draw();
        Utilities::saveCanvas(c_sb, "results/" + saveName + "_s_lin" + name + "_" + det);


	TCanvas *c_sb_log = new TCanvas();
	c_sb_log->cd();
	c_sb_log->SetLogy();
        if (h1_signal->GetMaximum() > h1_background->GetMaximum()) {
                h1_signal->Draw();
                h1_signal->GetYaxis()->SetRangeUser(0.1, h1_signal->GetMaximum() 
                                                        + h1_signal->GetMaximum());
                h1_background->Draw("SAME");
        }
        else {
                h1_background->Draw();
                h1_background->GetYaxis()->SetRangeUser(0.1, h1_background->GetMaximum() 
                                                                + h1_background->GetMaximum());
                h1_signal->Draw("SAME");
        }
        if (cutVal > 0.0) arr_cut->Draw();
        if (bin_eff99 != 0) arr_eff99->Draw();
        if (bin_eff98 != 0) arr_eff98->Draw();
        if (bin_eff95 != 0) arr_eff95->Draw();
	if (bin_eff80 != 0) arr_eff80->Draw();
        //if (bin_bg50 != 0) arr_bg50->Draw();

        lg->Draw();
	Utilities::saveCanvas(c_sb_log, "results/" + saveName + "_sb_log" + name + "_" + det);

        h1_signal->Draw();
        Utilities::saveCanvas(c_sb_log, "results/" + saveName + "_s_log" + name + "_" + det);

	delete gr;
	delete h1_signal;
	delete h1_background;
	delete c;
	delete c_sb;
	delete c_sb_log;
	delete lg;
	delete arr_cut;
	delete h1_both;
	delete arr_eff99;
	delete arr_eff98;
	delete arr_eff95;
	//delete arr_bg50;
	std::cout << "[plotEff] Done" << std::endl;
}



void plotEffVar(HistogramUtilities &h1, TString name, TString det, TString saveName, Int_t rebin = 1)
{

        TCanvas *c = new TCanvas();
        TH1F *h1_eff = 0;
        TH1F *h1_total = 0;
        TH1F *h1_pass = 0;
	TGraphAsymmErrors *gr_eff_s = new TGraphAsymmErrors();
        TGraphAsymmErrors *gr_eff_bg = new TGraphAsymmErrors();

        h1_total = h1.getHistogram(theSignal, name, "", det + "_denom", rebin, "");
        h1_pass = h1.getHistogram(theSignal, name, "", det + "_numer", rebin, "");
        //h1_total->Rebin(rebin);
        //h1_pass->Rebin(rebin);
	gr_eff_s->BayesDivide(h1_pass, h1_total);
        gr_eff_s->SetMarkerColor(kBlue);
        gr_eff_s->SetLineColor(kBlue);
        //h1_eff = (TH1F*)h1_pass->Clone();
        //h1_eff->Reset();
        //h1_eff->Divide(h1_pass, h1_total, 1.0, 1.0, "B");
        //h1_eff->SetName(h1_pass->GetName());
        c->SetName(TString("c_s_") + h1_pass->GetName());
        c->cd();
        //h1_eff->SetLineWidth(2);
        //h1_eff->SetMarkerStyle(20);
	//h1_eff->GetYaxis()->SetRangeUser(0, 1.1);
	//h1_eff->GetXaxis()->SetRangeUser(0, 100.0);
        //h1_eff->Draw("E1");
	gr_eff_s->Draw("AP");
	gr_eff_s->GetXaxis()->SetRangeUser(0, 100.0);
        gr_eff_s->GetYaxis()->SetRangeUser(0, 1.1);
	gr_eff_s->GetXaxis()->SetTitle(h1_total->GetXaxis()->GetTitle());
        gr_eff_s->GetYaxis()->SetTitle("Efficiency");
        Utilities::saveCanvas(c, "results/" + saveName + "_effVar_s_" + name + "_" + det);

        h1_total = h1.getHistogram(theBackground, name, "", det + "_denom", rebin, "");
        h1_pass = h1.getHistogram(theBackground, name, "", det + "_numer", rebin, "");
	//h1_total->Rebin(rebin);
	//h1_pass->Rebin(rebin);
        gr_eff_bg->BayesDivide(h1_pass, h1_total);
	gr_eff_bg->SetMarkerColor(kGreen);
	gr_eff_bg->SetLineColor(kGreen);
        //h1_eff->Reset();
        //h1_eff->Divide(h1_pass, h1_total, 1.0, 1.0, "B");
        //h1_eff->SetName(h1_pass->GetName());
        c->SetName(TString("c_bg_") + h1_pass->GetName());
        c->cd();
	//h1_eff->SetLineWidth(2);
	//h1_eff->SetMarkerStyle(20);
        //h1_eff->GetYaxis()->SetRangeUser(0, 1.1);
        //h1_eff->GetXaxis()->SetRangeUser(0, 100.0);
        //h1_eff->Draw("E1");
	gr_eff_bg->Draw("AP");
        gr_eff_bg->GetXaxis()->SetRangeUser(0, 100.0);
        gr_eff_bg->GetYaxis()->SetRangeUser(0, 1.1);
        gr_eff_bg->GetXaxis()->SetTitle(h1_total->GetXaxis()->GetTitle());
	gr_eff_bg->GetYaxis()->SetTitle("Efficiency");
        Utilities::saveCanvas(c, "results/" + saveName + "_effVar_bg_" + name + "_" + det);

        TLegend *lg = new TLegend(0.5, 0.2, 0.9, 0.5);
        lg->SetFillColor(kWhite);
        lg->SetLineColor(kWhite);
        lg->SetShadowColor(kWhite);
        TString upperDet = det;
        upperDet.ToUpper();
        lg->AddEntry(gr_eff_s, "Signal (" + upperDet + ")", "lp");
        lg->AddEntry(gr_eff_bg, "Background (" + upperDet + ")", "lp");

        c->SetName(TString("c_sb_") + h1_pass->GetName());
        c->cd();
	c->Clear();
        gr_eff_bg->Draw("AP");
        gr_eff_s->Draw("P");
	lg->Draw();
        Utilities::saveCanvas(c, "results/" + saveName + "_effVar_sb_" + name + "_" + det);

        delete c;
	delete lg;
        //delete h1_eff;
        delete h1_pass;
        delete h1_total;
	delete gr_eff_s;
	delete gr_eff_bg;
}

void plotStack(HistogramUtilities &h1, TString name, TString titleX, TString saveName, TString det, int rebin, float cutValEB, float cutValEE)
{

        THStack *st = h1.getStack(theSources, name, "", det, rebin);
        TLegend *lg_all = h1.getLegend(theSources, name, "", det);

        TCanvas *c1 = new TCanvas();
        c1->cd();
        st->Draw();
	st->GetXaxis()->SetTitle(titleX); 
        lg_all->Draw();
	if (det == "ee" && cutValEE != -1) getArrow(st, det, cutValEB, cutValEE)->Draw();
        if (det == "eb" && cutValEB != -1) getArrow(st, det, cutValEB, cutValEE)->Draw();
        Utilities::saveCanvas(c1, "results/" + saveName  + "_lin_" + name + "_" + det);

	c1->SetLogy();
	st->SetMinimum(1.0);
        st->Draw();	
	lg_all->Draw();
        if (det == "ee" && cutValEE != -1) getArrow(st, det, cutValEB, cutValEE)->Draw();
        if (det == "eb" && cutValEB != -1) getArrow(st, det, cutValEB, cutValEE)->Draw();
        Utilities::saveCanvas(c1, "results/" + saveName  + "_log_" + name + "_" + det);

	delete c1;
	delete st;
	delete lg_all;

}

void plotAllResultsID()
{
	plotResultsID("ee", "isoV0_studies");
	plotResultsID("eb", "isoV0_studies");
}

void plotAllResultsW()
{
        //plotResultsW("ee", "iso10_jpt25_tcmet30");
        //plotResultsW("eb", "iso10_jpt25_tcmet30");

        //plotResultsW("ee", "iso10_jptphimax110_tcmet30");
        //plotResultsW("eb", "iso10_jptphimax110_tcmet30");


	//
	// jet veto studies
	//
	/*
        plotResultsW("ee", "pt20_isoV1_tcmet30");
        plotResultsW("eb", "pt20_isoV1_tcmet30");

        plotResultsW("ee", "pt30_isoV1_tcmet30");
        plotResultsW("eb", "pt30_isoV1_tcmet30");

        plotResultsW("ee", "pt20_isoV1_phimax130_tcmet30");
        plotResultsW("eb", "pt20_isoV1_phimax130_tcmet30");
	*/

	//
	// try to apply candidate electron id
	// tasElectron_v0
	//
	plotResultsW("ee", "pt20_isoV1_phimax130_tasv1_tcmet30");
        plotResultsW("eb", "pt20_isoV1_phimax130_tasv1_tcmet30");

}

void plotAllResultsAN2009_098()
{
//        plotResultsW("ee", "AN2009_098_studies");
//        plotResultsW("eb", "AN2009_098_studies");

        plotResultsW("ee", "AN2009_098_studies_30v1");
        plotResultsW("eb", "AN2009_098_studies_30v1");

//        plotResultsAN2009_098("ee", "AN2009_098_studies");
//        plotResultsAN2009_098("eb", "AN2009_098_studies");

        plotResultsAN2009_098("ee", "AN2009_098_studies_30v1");
        plotResultsAN2009_098("eb", "AN2009_098_studies_30v1");

}

void plotResultsAN2009_098(TString det, TString fileStamp)
{

        gROOT->ProcessLine(".L ~/tdrStyle.C");
        gROOT->ProcessLine("setTDRStyle()");

        // luminorm for 1pb-1
        // luminosity is already normalised to 1pb-1 in the looper
        HistogramUtilities h1("Results_" + fileStamp + ".root", 1.0);

        plotStack(h1, "AN2009_098_pt2", "Second p_{T} (GeV)", fileStamp, det, 2, 20.0, 20.0);
        plotStack(h1, "AN2009_098_eta1", "Electron #eta", fileStamp, det, 2);
        plotStack(h1, "AN2009_098_ecalIso", "Ecal Iso", fileStamp, det, 2, 4.2, 3.4);
        plotStack(h1, "AN2009_098_hcalIso", "Hcal Iso", fileStamp, det, 2, 2.0, 1.3);
        plotStack(h1, "AN2009_098_tkIso", "Track Iso", fileStamp, det, 2, 2.2, 1.1);
        plotStack(h1, "AN2009_098_tcmet_after_selection", "tcMet (GeV)", fileStamp, det, 2, 30.0, 30.0);

	// clearer versions of plots...
        plotEff(h1, "AN2009_098_pt2", fileStamp, det, true, 2, true, 20.0, 20.0);
        plotEff(h1, "AN2009_098_eta1", fileStamp, det, true, 2); 
        plotEff(h1, "AN2009_098_ecalIso", fileStamp, det, true, 2, true, 4.2, 3.4);
        plotEff(h1, "AN2009_098_hcalIso", fileStamp, det, true, 2, true, 2.0, 1.3);
        plotEff(h1, "AN2009_098_tkIso", fileStamp, det, true, 2, true, 2.2, 1.1);


}

void plotResultsW(TString det, TString fileStamp)
{

        gROOT->ProcessLine(".L ~/tdrStyle.C");
        gROOT->ProcessLine("setTDRStyle()");

        // luminorm for 1pb-1
        // luminosity is already normalised to 1pb-1 in the looper
        HistogramUtilities h1("Results_" + fileStamp + ".root", 1.0);

        // W studies related
        //

        plotStack(h1, "weff_pt", "p_{T} (GeV)", fileStamp, det, 2);
        plotStack(h1, "weff_iso", "Isolation", fileStamp, det, 2);
        plotStack(h1, "weff_tcmet", "tcMET (GeV)", fileStamp, det, 2);
        plotStack(h1, "weff_jptpt", "Lead JPT p_{T} (GeV)", fileStamp, det, 2);
        plotStack(h1, "weff_tcmet_after_iso", "tcMET (GeV)", fileStamp, det, 2);
        plotStack(h1, "weff_jptpt_after_iso", "Lead JPT p_{T} (GeV)", fileStamp, det, 2);
	plotStack(h1, "weff_leastemjpt_after_iso", "Lowest EMFrac JPT", fileStamp, det, 2);

        plotStack(h1, "weff_d0corr_after_iso", "d0corr", fileStamp, det, 2);
        plotStack(h1, "weff_d0corr_after_iso_jpt", "d0corr", fileStamp, det, 2);


        plotStack(h1, "weff_leadjptphi_after_iso", "#Delta#phi{Lead JPT, electron}  (Degrees)", fileStamp, det, 2);
        plotStack(h1, "weff_tcmet_after_iso_jpt", "tcMET (GeV)", fileStamp, det, 2);
        plotStack(h1, "weff_leadjptphi_after_iso_jpt", "#Delta#phi{Lead JPT, electron}  (Degrees)", fileStamp, det, 2);
        plotStack(h1, "weff_leadjptphi_after_iso_jpt_tcmet", "#Delta#phi{Lead JPT, electron}  (Degrees)", fileStamp, det, 2);

        plotStack(h1, "weff_jptphimax_after_iso", "#Delta#phi_{Max}{JPT, electron}  (Degrees)", fileStamp, det, 2);
        plotStack(h1, "weff_tcmet_after_iso_jpt_conv", "tcMET (GeV)", fileStamp, det, 2);

        plotEff(h1, "weff_leadjptphi_after_iso", fileStamp, det, true, 2);
        plotEff(h1, "weff_jptpt_after_iso", fileStamp, det, true, 2);
        plotEff(h1, "weff_leadjptphi_after_iso_jpt", fileStamp, det, true, 2);
        plotEff(h1, "weff_leadjptphi_after_iso_jpt_tcmet", fileStamp, det, true, 2);

        plotEff(h1, "weff_jptphimax_after_iso", fileStamp, det, true, 2);
        plotEff(h1, "weff_jptpt_after_iso", fileStamp, det, true, 2);

	// distributions selected
	plotStack(h1, "weffs_sigmaIEtaIEta", "sigmaIEtaIEta", fileStamp, det);
        plotStack(h1, "weffbg_sigmaIEtaIEta", "sigmaIEtaIEta", fileStamp, det);

}

void plotResultsID(TString det, TString fileStamp)
{

	gROOT->ProcessLine(".L ~/tdrStyle.C");
	gROOT->ProcessLine("setTDRStyle()");
        gROOT->ProcessLine("gStyle->SetPalette(1)");

	// luminorm for 1pb-1
	// luminosity is already normalised to 1pb-1 in the looper
	HistogramUtilities h1("Results_" + fileStamp + ".root", 1.0);

	// electron id related
	//

        plotEff(h1, "dEtaIn", "IDStudy", det, true, 2, true, 0.007, 0.010);
        plotEff(h1, "dPhiIn", "IDStudy", det, true, 2, true, 0.020, 0.025);
        plotEff(h1, "dPhiInSigned", "IDStudy", det, true, 2, false);
        plotEff(h1, "hoe", "IDStudy", det, true, 2, true, 0.01, 0.01);
        plotEff(h1, "eopIn", "IDStudy", det, false, 4, true, 0.5, 0.5);
        plotEff(h1, "sigmaIEtaIEta", "IDStudy", det, true, 2, true, -1, 0.03);
        plotEff(h1, "sigmaIPhiIPhi", "IDStudy", det, true, 4, true);
        plotEff(h1, "E2x5Norm5x5","IDStudy", det, false, 1, false, 0.90, -1);
        plotEff(h1, "E1x5Norm5x5", "IDStudy", det, false, 4, false);
        plotEff(h1, "d0corr", "IDStudy", det, true, 4, true, 0.025, 0.035);

	// N-1 with respect to TasV1
	//
	// the distributions
        plotEff(h1, "dEtaInTasV1NM1", "IDStudy", det, true, 2, true, 0.007, 0.010);
        plotEff(h1, "dPhiInTasV1NM1", "IDStudy", det, true, 2, true, 0.020, 0.025);
        plotEff(h1, "hoeTasV1NM1", "IDStudy", det, true, 2, true, 0.01, 0.01);
        plotEff(h1, "sigmaIEtaIEtaTasV1NM1", "IDStudy", det, true, 2, true, -1, 0.03);
        plotEff(h1, "E2x5Norm5x5TasV1NM1","IDStudy", det, false, 1, false, 0.90, -1);
	// and the efficiency curves
        plotEffVar(h1, "dEtaInTasV1NM1_pt", det, "IDStudy", 4);
        plotEffVar(h1, "dPhiInTasV1NM1_pt", det, "IDStudy", 4);
        plotEffVar(h1, "hoeTasV1NM1_pt", det, "IDStudy", 4);
        plotEffVar(h1, "sigmaIEtaIEtaTasV1NM1_pt", det, "IDStudy", 4);
        plotEffVar(h1, "E2x5Norm5x5TasV1NM1_pt", det, "IDStudy", 4);
        plotEffVar(h1, "tasElectronV1_pt", det, "IDStudy", 4);

        plotEffVar(h1, "dEtaInTasV1NM1_eta", det, "IDStudy", 4);
        plotEffVar(h1, "dPhiInTasV1NM1_eta", det, "IDStudy", 4);
        plotEffVar(h1, "hoeTasV1NM1_eta", det, "IDStudy", 4);
        plotEffVar(h1, "sigmaIEtaIEtaTasV1NM1_eta", det, "IDStudy", 4);
        plotEffVar(h1, "E2x5Norm5x5TasV1NM1_eta", det, "IDStudy", 4);
        plotEffVar(h1, "tasElectronV1_eta", det, "IDStudy", 4);


	// the isolation part
	//
        plotEff(h1, "wwIsoAll", "IDStudy", det, true, 1, true);
        plotEff(h1, "tkIso03All", "IDStudy", det, true, 1, true);
        plotEff(h1, "ecalIso03All", "IDStudy", det, true, 1, true);
        plotEff(h1, "hcalIso03All", "IDStudy", det, true, 1, true);

        plotEff(h1, "tkIso03AllRe", "IDStudy", det, true, 1, true);
        plotEff(h1, "tkIso03AllReShVeto", "IDStudy", det, true, 1, true);
        plotEff(h1, "tkIso03AllReRel", "IDStudy", det, true, 1, true);
        plotEff(h1, "caloIso03All", "IDStudy", det, true, 1, true);
        plotEff(h1, "tkIso03AllReJura01", "IDStudy", det, true, 1, true);
        plotEff(h1, "tkIso03AllReJura02", "IDStudy", det, true, 1, true);
        plotEff(h1, "tkIso03AllReJura03", "IDStudy", det, true, 1, true);
	plotEff(h1, "tkIso03AllReJura01In015", "IDStudy", det, true, 1, true);

	// 2D stuff
	plot2DSB(h1, "tkIso03All2D", "p_{T} (GeV/c)", "tkIso03All", "IDStudy", det);
        plot2DSB(h1, "caloIso03All2D", "E_{T} (GeV)", "caloIso03All", "IDStudy", det);
        plot2DSB(h1, "tkIso03AllRedR2D", "dEta", "dPhi", "IDStudy", det);

	// N-1
        plotEff(h1, "tkIso03AllNM1", "IDStudy", det, true, 1, true);
        plotEff(h1, "tkIso03AllIDNM1", "IDStudy", det, true, 1, true);
        plotEff(h1, "tkIso03AllConvNM1", "IDStudy", det, true, 1, true);
        plotEff(h1, "tkIso03AllConvIDNM1", "IDStudy", det, true, 1, true);

        plotEff(h1, "ecalIso03AllNM1", "IDStudy", det, true, 1, true);
        plotEff(h1, "hcalIso03AllNM1", "IDStudy", det, true, 1, true);

	plotEff(h1, "tkIso03AllReJura01In015NM1", "IDStudy", det, true, 1, true);
        plotEff(h1, "tkIso03AllReJura01In015IDNM1", "IDStudy", det, true, 1, true);
        plotEff(h1, "tkIso03AllReJura01In015ConvNM1", "IDStudy", det, true, 1, true);
        plotEff(h1, "tkIso03AllReJura01In015ConvIDNM1", "IDStudy", det, true, 1, true);

        plotEff(h1, "tkIso03AllReShCutNM1", "IDStudy", det, true, 1, true);


}

