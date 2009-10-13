
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

void plotEff(HistogramUtilities &h1, TString name, TString saveName, TString det, bool ascending, int rebin, bool legendOnRight, float cutValEB, float cutValEE)
{

	std::cout << "[plotEff]" << std::endl;
	TH1F *h1_signal = h1.getHistogram(theSignal, name, "", det, rebin, "");
        TH1F *h1_background = h1.getHistogram(theBackground, name, "", det, rebin, "");
	TH1F *h1_both = h1.getHistogram(theBackground | theSignal, name, "", det, rebin, "");
	// latter bool is "ascending"
	bool normalise = true;
        TGraph *gr = (TGraph*)(eff_rej(*h1_signal, *h1_background, normalise, ascending).Clone());
        gr->SetTitle(name + ";Signal;Background");
	gr->SetMarkerStyle(23);

        TCanvas *c = new TCanvas();
        c->cd();
        gr->Draw("AP");
	gr->GetXaxis()->SetRangeUser(0, 1.1);
	gr->GetYaxis()->SetRangeUser(0, 1.1);
	Utilities::saveCanvas(c, "results/" + saveName + "_eff_" + name + "_" + det);

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
	std::cout << "[plotEff] Done" << std::endl;
}



void plotEffVar(HistogramUtilities &h1, TString name, TString saveName, Int_t rebin = 1)
{

        TCanvas *c = new TCanvas();
        TH1F *h1_eff = 0;
        TH1F *h1_total = 0;
        TH1F *h1_pass = 0;
	//TGraphAsymmErrors *gr_eff = new TGraphAsymmErrors();

        h1_total = h1.getHistogram(theSignal, name, "", "denom", rebin, "");
        h1_pass = h1.getHistogram(theSignal, name, "", "numer", rebin, "");
        //h1_total->Rebin(rebin);
        //h1_pass->Rebin(rebin);
	//gr_eff->BayesDivide(h1_pass, h1_total);
        h1_eff = (TH1F*)h1_pass->Clone();
        h1_eff->Reset();
        h1_eff->Divide(h1_pass, h1_total, 1.0, 1.0, "B");
        h1_eff->SetName(h1_pass->GetName());
        c->SetName(TString("c_s_") + h1_pass->GetName());
        c->cd();
        h1_eff->SetLineWidth(2);
        h1_eff->SetMarkerStyle(20);
	h1_eff->GetYaxis()->SetRangeUser(0, 1.1);
        h1_eff->Draw("E1");
	//gr_eff->Draw("AP");
        //gr_eff->GetYaxis()->SetRangeUser(0, 1.1);
        Utilities::saveCanvas(c, "results/" + saveName + "_effVar_s_" + name);


        h1_total = h1.getHistogram(theBackground, name, "", "denom", rebin, "");
        h1_pass = h1.getHistogram(theBackground, name, "", "numer", rebin, "");
	//h1_total->Rebin(rebin);
	//h1_pass->Rebin(rebin);
        //gr_eff->BayesDivide(h1_pass, h1_total);
        h1_eff->Reset();
        h1_eff->Divide(h1_pass, h1_total, 1.0, 1.0, "B");
        h1_eff->SetName(h1_pass->GetName());
        c->SetName(TString("c_bg_") + h1_pass->GetName());
        c->cd();
	h1_eff->SetLineWidth(2);
	h1_eff->SetMarkerStyle(20);
        h1_eff->GetYaxis()->SetRangeUser(0, 1.1);
        h1_eff->Draw("E1");
	//gr_eff->Draw("AP");
        //gr_eff->GetYaxis()->SetRangeUser(0, 1.1);
        Utilities::saveCanvas(c, "results/" + saveName + "_effVar_bg_" + name);

        delete c;
        delete h1_eff;
        delete h1_pass;
        delete h1_total;
	//delete gr_eff;
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
        Utilities::saveCanvas(c1, "results/" + saveName  + "_log_" + name + "_" + det);

	delete c1;
	delete st;
	delete lg_all;

}

void test()
{
        gROOT->ProcessLine(".L ~/tdrStyle.C");
        gROOT->ProcessLine("setTDRStyle()");
        HistogramUtilities h1("Results.root", 0.001);

        plotEffVar(h1, "dEtaIn_pt_ee", "", 4);
        plotEffVar(h1, "dEtaIn_eta_ee", "");
        plotEffVar(h1, "dEtaIn_phi_ee", "");
        plotEffVar(h1, "dEtaIn_pt_eb", "", 4);
        plotEffVar(h1, "dEtaIn_eta_eb", "");
        plotEffVar(h1, "dEtaIn_phi_eb", "");

        plotEffVar(h1, "dPhiIn_pt_ee", "", 4);
        plotEffVar(h1, "dPhiIn_eta_ee", "");
        plotEffVar(h1, "dPhiIn_phi_ee", "");
        plotEffVar(h1, "dPhiIn_pt_eb", "", 4);
        plotEffVar(h1, "dPhiIn_eta_eb", "");
        plotEffVar(h1, "dPhiIn_phi_eb", "");

        plotEffVar(h1, "hoe_pt_ee", "", 4);
        plotEffVar(h1, "hoe_eta_ee", "");
        plotEffVar(h1, "hoe_phi_ee", "");
        plotEffVar(h1, "hoe_pt_eb", "", 4);
        plotEffVar(h1, "hoe_eta_eb", "");
        plotEffVar(h1, "hoe_phi_eb", "");

        plotEffVar(h1, "sieie_pt_ee", "", 4);
        plotEffVar(h1, "sieie_eta_ee", "");
        plotEffVar(h1, "sieie_phi_ee", "");
        plotEffVar(h1, "sieie_pt_eb", "", 4);
        plotEffVar(h1, "sieie_eta_eb", "");
        plotEffVar(h1, "sieie_phi_eb", "");

        plotEffVar(h1, "robustTight_pt_eb", "");
        plotEffVar(h1, "robustTight_pt_ee", "");

        plotEffVar(h1, "classBasedTight_pt_eb", "");
        plotEffVar(h1, "classBasedTight_pt_ee", "");

        plotEffVar(h1, "eopInGT05_pt_ee", "", 4);
        plotEffVar(h1, "eopInGT05_pt_eb", "", 4);

        plotEffVar(h1, "eopInLT30_pt_ee", "", 4);
        plotEffVar(h1, "eopInLT30_pt_eb", "", 4);

}

void plotAllResults()
{
	plotResults("ee");
	plotResults("eb");
}

void plotAllResultsW()
{
        plotResultsW("ee", "iso10_jpt25_tcmet30");
        plotResultsW("eb", "iso10_jpt25_tcmet30");

        plotResultsW("ee", "iso10_jpt25_tcmet20");
        plotResultsW("eb", "iso10_jpt25_tcmet20");


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
        plotStack(h1, "weff_leadjptphi_after_iso", "#Delta#phi{Lead JPT, electron}  (Degrees)", fileStamp, det, 2);
        plotStack(h1, "weff_tcmet_after_iso_jpt", "tcMET (GeV)", fileStamp, det, 2);
        plotStack(h1, "weff_leadjptphi_after_iso_jpt", "#Delta#phi{Lead JPT, electron}  (Degrees)", fileStamp, det, 2);
        plotStack(h1, "weff_leadjptphi_after_iso_jpt_tcmet", "#Delta#phi{Lead JPT, electron}  (Degrees)", fileStamp, det, 2);

        plotStack(h1, "weff_jptphimax_after_iso", "#Delta#phi_{Max}{JPT, electron}  (Degrees)", fileStamp, det, 2);

        plotEff(h1, "weff_leadjptphi_after_iso", fileStamp, det, true, 2);
        plotEff(h1, "weff_jptpt_after_iso", fileStamp, det, true, 2);
        plotEff(h1, "weff_leadjptphi_after_iso_jpt", fileStamp, det, true, 2);
        plotEff(h1, "weff_leadjptphi_after_iso_jpt_tcmet", fileStamp, det, true, 2);

        plotEff(h1, "weff_jptphimax_after_iso", fileStamp, det, true, 2);

}

void plotResults(TString det)
{

	gROOT->ProcessLine(".L ~/tdrStyle.C");
	gROOT->ProcessLine("setTDRStyle()");

	// luminorm for 1pb-1
	// luminosity is already normalised to 1pb-1 in the looper
	HistogramUtilities h1("Results.root", 1.0);

	// electron id related
	//

        plotEff(h1, "dEtaIn", "", det, true, 2, true, 0.007, 0.010);
        plotEff(h1, "dPhiIn", "", det, true, 2, true, 0.020, 0.025);
        plotEff(h1, "dPhiInSigned", "", det, true, 2, true);
        plotEff(h1, "hoe", det, "", true, 2, true, 0.01, 0.01);
        plotEff(h1, "eopIn", det, "", false, 4, true, 0.5, 0.5);
        plotEff(h1, "sigmaIEtaIEta", "", det, true, 2, true, -1, 0.03);
        plotEff(h1, "sigmaIPhiIPhi", "", det, true, 4, true);
        plotEff(h1, "E2x5Norm5x5", "", det, false, 4, false, 0.90, -1);
        plotEff(h1, "E1x5Norm5x5", "", det, false, 4, false);
        plotEff(h1, "d0corr", det, "", true, 4, true, 0.025, 0.035);

	// the isolation part
        plotEff(h1, "wwIsoAll", "", det, true, 1, true);
        plotEff(h1, "tkIso03All", "", det, true, 1, true);
        plotEff(h1, "ecalIso03All", "", det, true, 1, true);
        plotEff(h1, "ecalIsoMod03All", "", det, true, 1, true);
        plotEff(h1, "hcalIso03All", "", det, true, 1, true);

}

