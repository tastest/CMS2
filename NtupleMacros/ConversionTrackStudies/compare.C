

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TROOT.h"

void saveCanvas(TCanvas *c1, const TString name, bool log)
{

	TString dir = "results/";
	c1->SetLogy(log);
	if (log) {
		c1->SaveAs(dir + name + "_log.png");
		c1->SaveAs(dir + name + "_log.root");
		c1->SaveAs(dir + name + "_log.C");
	}
	else {
		c1->SaveAs(dir + name + "_lin.png");
		c1->SaveAs(dir + name + "_lin.root");
		c1->SaveAs(dir + name + "_lin.C");
	}
}

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

void compareDataDataEcalSpecial(TFile *f_data, TFile *f_mc, TString saveName)
{

	// data

        TH1F *h1_pseudo_ecalIso03_data_eem = (TH1F*)f_data->Get("h1_pseudo_ecalIso03_eem");
        TH1F *h1_pseudo_ecalIso03_data_eep = (TH1F*)f_data->Get("h1_pseudo_ecalIso03_eep");
        formatHist(h1_pseudo_ecalIso03_data_eem, true);
        formatHist(h1_pseudo_ecalIso03_data_eep, true);
	h1_pseudo_ecalIso03_data_eem->SetMarkerColor(kBlack);
	h1_pseudo_ecalIso03_data_eem->SetLineColor(kBlack);

        TH1F *h1_pseudo_ecalIso03_data_ebm = (TH1F*)f_data->Get("h1_pseudo_ecalIso03_ebm");
        TH1F *h1_pseudo_ecalIso03_data_ebp = (TH1F*)f_data->Get("h1_pseudo_ecalIso03_ebp");
        formatHist(h1_pseudo_ecalIso03_data_ebm, true);
        formatHist(h1_pseudo_ecalIso03_data_ebp, true);
        h1_pseudo_ecalIso03_data_ebm->SetMarkerColor(kBlack);
        h1_pseudo_ecalIso03_data_ebm->SetLineColor(kBlack);

	// mc

        TH1F *h1_pseudo_ecalIso03_mc_eem = (TH1F*)f_mc->Get("h1_pseudo_ecalIso03_eem");
        TH1F *h1_pseudo_ecalIso03_mc_eep = (TH1F*)f_mc->Get("h1_pseudo_ecalIso03_eep");
        formatHist(h1_pseudo_ecalIso03_mc_eem, false);
        formatHist(h1_pseudo_ecalIso03_mc_eep, false);
        h1_pseudo_ecalIso03_mc_eem->SetMarkerColor(kBlack);
        h1_pseudo_ecalIso03_mc_eem->SetLineColor(kBlack);

        TH1F *h1_pseudo_ecalIso03_mc_ebm = (TH1F*)f_mc->Get("h1_pseudo_ecalIso03_ebm");
        TH1F *h1_pseudo_ecalIso03_mc_ebp = (TH1F*)f_mc->Get("h1_pseudo_ecalIso03_ebp");
        formatHist(h1_pseudo_ecalIso03_mc_ebm, false);
        formatHist(h1_pseudo_ecalIso03_mc_ebp, false);
        h1_pseudo_ecalIso03_mc_ebm->SetMarkerColor(kBlack);
        h1_pseudo_ecalIso03_mc_ebm->SetLineColor(kBlack);

	TH1F *h1_pseudo_ecalIso03_dataratio_ebpm = (TH1F*)h1_pseudo_ecalIso03_data_ebm->Clone();
	h1_pseudo_ecalIso03_dataratio_ebpm->Divide(h1_pseudo_ecalIso03_data_ebp);

        TH1F *h1_pseudo_ecalIso03_dataratio_eepm = (TH1F*)h1_pseudo_ecalIso03_data_eem->Clone();
	h1_pseudo_ecalIso03_dataratio_eepm->Divide(h1_pseudo_ecalIso03_data_eep);

        TH1F *h1_pseudo_ecalIso03_mcratio_ebpm = (TH1F*)h1_pseudo_ecalIso03_mc_ebm->Clone();
        h1_pseudo_ecalIso03_mcratio_ebpm->Divide(h1_pseudo_ecalIso03_mc_ebp);

        TH1F *h1_pseudo_ecalIso03_mcratio_eepm = (TH1F*)h1_pseudo_ecalIso03_mc_eem->Clone();
        h1_pseudo_ecalIso03_mcratio_eepm->Divide(h1_pseudo_ecalIso03_mc_eep);

        TCanvas *c1 = new TCanvas();
        c1->cd();

	h1_pseudo_ecalIso03_data_eem->Draw();
	h1_pseudo_ecalIso03_data_eep->Draw("SAME");
        saveCanvas(c1, saveName + "_pseudo_ecalIso03_data_eepm", false);
        saveCanvas(c1, saveName + "_pseudo_ecalIso03_data_eepm", true);

	h1_pseudo_ecalIso03_mcratio_ebpm->GetXaxis()->SetTitle("ecalIso03(EB-)/ecalIso03(EB+)");
        h1_pseudo_ecalIso03_mcratio_ebpm->Draw("HIST");
        h1_pseudo_ecalIso03_dataratio_ebpm->Draw("SAME E1");
	saveCanvas(c1, saveName + "_pseudo_ecalIso03_ratio_ebpm", false);

        h1_pseudo_ecalIso03_data_ebm->Draw();
        h1_pseudo_ecalIso03_data_ebp->Draw("SAME");
        saveCanvas(c1, saveName + "_pseudo_ecalIso03_data_ebpm", false);
        saveCanvas(c1, saveName + "_pseudo_ecalIso03_data_ebpm", true);

        h1_pseudo_ecalIso03_mcratio_eepm->GetXaxis()->SetTitle("ecalIso03(EE-)/ecalIso03(EE+)");
	h1_pseudo_ecalIso03_mcratio_eepm->Draw("HIST");
        h1_pseudo_ecalIso03_dataratio_eepm->Draw("SAME E1");
        saveCanvas(c1, saveName + "_pseudo_ecalIso03_ratio_eepm", false);	
	
	delete c1;
	delete h1_pseudo_ecalIso03_data_ebm;
	delete h1_pseudo_ecalIso03_data_ebp;
	delete h1_pseudo_ecalIso03_data_eem;
	delete h1_pseudo_ecalIso03_data_eep;

}

double histMax(TH1F *hist1, TH1F *hist2)
{

	if (hist1->GetMaximum(hist1->GetMaximum()) > hist2->GetMaximum(hist2->GetMaximum())) return hist1->GetMaximum(hist1->GetMaximum())*1.1;
	return hist2->GetMaximum(hist2->GetMaximum())*1.1;
}

void compareOneEcalSpecial(TFile *f_mc, TFile *f_data, TString saveName, TString det)
{

        TH1F *h1_pseudo_ecalIso03_mc = (TH1F*)f_mc->Get("h1_pseudo_ecalIso03_" + det);
        TH1F *h1_pseudo_ecalIso03_data = (TH1F*)f_data->Get("h1_pseudo_ecalIso03_" + det);
        formatHist(h1_pseudo_ecalIso03_mc, false);
        formatHist(h1_pseudo_ecalIso03_data, true);

        TCanvas *c1 = new TCanvas();
        c1->cd();

        h1_pseudo_ecalIso03_data->Draw();
        h1_pseudo_ecalIso03_mc->Scale(h1_pseudo_ecalIso03_data->Integral()/h1_pseudo_ecalIso03_mc->Integral());
        h1_pseudo_ecalIso03_mc->Draw("SAME HIST");
        h1_pseudo_ecalIso03_data->Draw("SAME");

        saveCanvas(c1, saveName + "_pseudo_ecalIso03_" + det, false);
        saveCanvas(c1, saveName + "_pseudo_ecalIso03_" + det, true);

	delete c1;
        delete h1_pseudo_ecalIso03_mc;
        delete h1_pseudo_ecalIso03_data;

}

void compareOne(TFile *f_mc, TFile *f_data, TString saveName, TString det)
{

	TH1F *h1_pseudo_ecalEta_mc = (TH1F*)f_mc->Get("h1_pseudo_ecalEta_" + det);
	TH1F *h1_pseudo_ecalEta_data = (TH1F*)f_data->Get("h1_pseudo_ecalEta_" + det);
	formatHist(h1_pseudo_ecalEta_mc, false);
	formatHist(h1_pseudo_ecalEta_data, true);

	TLegend *lg = new TLegend(0.6, 0.6, 0.9, 0.9);
        lg->SetFillColor(kWhite);
        lg->SetLineColor(kWhite);
        lg->SetShadowColor(kWhite);
        TString upperDet = det;
        upperDet.ToUpper();
        lg->AddEntry(h1_pseudo_ecalEta_mc, "MC (" + upperDet + ")", "f");
        lg->AddEntry(h1_pseudo_ecalEta_data, "Data (" + upperDet + ")", "lp");

	TH1F *h1_pseudo_ecalPhi_mc = (TH1F*)f_mc->Get("h1_pseudo_ecalPhi_" + det);
	TH1F *h1_pseudo_ecalPhi_data = (TH1F*)f_data->Get("h1_pseudo_ecalPhi_" + det);
	formatHist(h1_pseudo_ecalPhi_mc, false);
	formatHist(h1_pseudo_ecalPhi_data, true);

	TH1F *h1_pseudo_ecalIso03_qual_mc = (TH1F*)f_mc->Get("h1_pseudo_ecalIso03_qual_" + det);
	TH1F *h1_pseudo_ecalIso03_qual_data = (TH1F*)f_data->Get("h1_pseudo_ecalIso03_qual_" + det);
	formatHist(h1_pseudo_ecalIso03_qual_mc, false);
	formatHist(h1_pseudo_ecalIso03_qual_data, true);

        TH1F *h1_pseudo_ecalIso03_fakeSR_mc = (TH1F*)f_mc->Get("h1_pseudo_ecalIso03_fakeSR_" + det);
        TH1F *h1_pseudo_ecalIso03_fakeSR_data = (TH1F*)f_data->Get("h1_pseudo_ecalIso03_fakeSR_" + det);
        formatHist(h1_pseudo_ecalIso03_fakeSR_mc, false);
        formatHist(h1_pseudo_ecalIso03_fakeSR_data, true);

	TH1F *h1_pseudo_ecalIso03_mc = (TH1F*)f_mc->Get("h1_pseudo_ecalIso03_" + det);
	TH1F *h1_pseudo_ecalIso03_data = (TH1F*)f_data->Get("h1_pseudo_ecalIso03_" + det);
	formatHist(h1_pseudo_ecalIso03_mc, false);
	formatHist(h1_pseudo_ecalIso03_data, true);

        TH1F *h1_pseudo_ecalIso03_E015_mc = (TH1F*)f_mc->Get("h1_pseudo_ecalIso03_E015_" + det);
        TH1F *h1_pseudo_ecalIso03_E015_data = (TH1F*)f_data->Get("h1_pseudo_ecalIso03_E015_" + det);
        formatHist(h1_pseudo_ecalIso03_E015_mc, false);
        formatHist(h1_pseudo_ecalIso03_E015_data, true);

        TH1F *h1_pseudo_ecalIso03_E010_mc = (TH1F*)f_mc->Get("h1_pseudo_ecalIso03_E010_" + det);
        TH1F *h1_pseudo_ecalIso03_E010_data = (TH1F*)f_data->Get("h1_pseudo_ecalIso03_E010_" + det);
        formatHist(h1_pseudo_ecalIso03_E010_mc, false);
        formatHist(h1_pseudo_ecalIso03_E010_data, true);

        TH1F *h1_pseudo_ecalIso03_recHitEt_mc = (TH1F*)f_mc->Get("h1_pseudo_ecalIso03_recHitEt_" + det);
        TH1F *h1_pseudo_ecalIso03_recHitEt_data = (TH1F*)f_data->Get("h1_pseudo_ecalIso03_recHitEt_" + det);
        formatHist(h1_pseudo_ecalIso03_recHitEt_mc, false);
        formatHist(h1_pseudo_ecalIso03_recHitEt_data, true);

        TH1F *h1_pseudo_ecalIso03_recHitE_mc = (TH1F*)f_mc->Get("h1_pseudo_ecalIso03_recHitE_" + det);
        TH1F *h1_pseudo_ecalIso03_recHitE_data = (TH1F*)f_data->Get("h1_pseudo_ecalIso03_recHitE_" + det);
        formatHist(h1_pseudo_ecalIso03_recHitE_mc, false);
        formatHist(h1_pseudo_ecalIso03_recHitE_data, true);

        TH1F *h1_pseudo_ecalIso03_recHitN_mc = (TH1F*)f_mc->Get("h1_pseudo_ecalIso03_recHitN_" + det);
        TH1F *h1_pseudo_ecalIso03_recHitN_data = (TH1F*)f_data->Get("h1_pseudo_ecalIso03_recHitN_" + det);
        formatHist(h1_pseudo_ecalIso03_recHitN_mc, false);
        formatHist(h1_pseudo_ecalIso03_recHitN_data, true);

	TH1F *h1_pseudo_hcalD1Iso03_mc = (TH1F*)f_mc->Get("h1_pseudo_hcalD1Iso03_" + det);
	TH1F *h1_pseudo_hcalD1Iso03_data = (TH1F*)f_data->Get("h1_pseudo_hcalD1Iso03_" + det);
	formatHist(h1_pseudo_hcalD1Iso03_mc, false);
	formatHist(h1_pseudo_hcalD1Iso03_data, true);

	TH1F *h1_pseudo_hcalD2Iso03_mc = (TH1F*)f_mc->Get("h1_pseudo_hcalD2Iso03_" + det);
	TH1F *h1_pseudo_hcalD2Iso03_data = (TH1F*)f_data->Get("h1_pseudo_hcalD2Iso03_" + det);
	formatHist(h1_pseudo_hcalD2Iso03_mc, false);
	formatHist(h1_pseudo_hcalD2Iso03_data, true);

	TH1F *h1_pseudo_tkIso03_mc = (TH1F*)f_mc->Get("h1_pseudo_tkIso03_" + det);
	TH1F *h1_pseudo_tkIso03_data = (TH1F*)f_data->Get("h1_pseudo_tkIso03_" + det);
	formatHist(h1_pseudo_tkIso03_mc, false);
	formatHist(h1_pseudo_tkIso03_data, true);

	TH1F *h1_pseudo_dRClosestTower_mc = (TH1F*)f_mc->Get("h1_pseudo_dRClosestTower_" + det);
	TH1F *h1_pseudo_dRClosestTower_data = (TH1F*)f_data->Get("h1_pseudo_dRClosestTower_" + det);
	formatHist(h1_pseudo_dRClosestTower_mc, false);
	formatHist(h1_pseudo_dRClosestTower_data, true);

	TCanvas *c1 = new TCanvas();
	c1->cd();
	h1_pseudo_ecalEta_mc->GetYaxis()->SetRangeUser(0.0, h1_pseudo_ecalEta_mc->GetMaximum() * 1.5);
	h1_pseudo_ecalEta_mc->Scale(h1_pseudo_ecalEta_data->Integral()/h1_pseudo_ecalEta_mc->Integral());
	h1_pseudo_ecalEta_mc->Draw("HIST");
	h1_pseudo_ecalEta_data->Draw("SAME");
	saveCanvas(c1, saveName + "_pseudo_ecalEta_" + det, false);

	h1_pseudo_ecalPhi_mc->GetYaxis()->SetRangeUser(0.0, h1_pseudo_ecalPhi_mc->GetMaximum() * 1.5);
	h1_pseudo_ecalPhi_mc->Scale(h1_pseudo_ecalPhi_data->Integral()/h1_pseudo_ecalPhi_mc->Integral());
	h1_pseudo_ecalPhi_mc->Draw("HIST");
	h1_pseudo_ecalPhi_data->Draw("SAME");
	saveCanvas(c1, saveName + "_pseudo_ecalPhi_" + det, false);

	h1_pseudo_ecalIso03_qual_data->Draw();
	h1_pseudo_ecalIso03_qual_mc->Scale(h1_pseudo_ecalIso03_qual_data->Integral()/h1_pseudo_ecalIso03_qual_mc->Integral());
	h1_pseudo_ecalIso03_qual_mc->Draw("SAME HIST");
	h1_pseudo_ecalIso03_qual_data->Draw("SAME");
	saveCanvas(c1, saveName + "_pseudo_ecalIso03_qual_" + det, false);
	saveCanvas(c1, saveName + "_pseudo_ecalIso03_qual_" + det, true);

        h1_pseudo_ecalIso03_fakeSR_data->Draw();
        h1_pseudo_ecalIso03_fakeSR_mc->Scale(h1_pseudo_ecalIso03_fakeSR_data->Integral()/h1_pseudo_ecalIso03_fakeSR_mc->Integral());
        h1_pseudo_ecalIso03_fakeSR_mc->Draw("SAME HIST");
        h1_pseudo_ecalIso03_fakeSR_data->Draw("SAME");
        saveCanvas(c1, saveName + "_pseudo_ecalIso03_fakeSR_" + det, false);
        saveCanvas(c1, saveName + "_pseudo_ecalIso03_fakeSR_" + det, true);


	h1_pseudo_ecalIso03_data->Draw();
	h1_pseudo_ecalIso03_mc->Scale(h1_pseudo_ecalIso03_data->Integral()/h1_pseudo_ecalIso03_mc->Integral());
	h1_pseudo_ecalIso03_mc->Draw("SAME HIST");
	h1_pseudo_ecalIso03_data->Draw("SAME");
	saveCanvas(c1, saveName + "_pseudo_ecalIso03_" + det, false);
	saveCanvas(c1, saveName + "_pseudo_ecalIso03_" + det, true);

        h1_pseudo_ecalIso03_E015_data->Draw();
        h1_pseudo_ecalIso03_E015_mc->Scale(h1_pseudo_ecalIso03_E015_data->Integral()/h1_pseudo_ecalIso03_E015_mc->Integral());
        h1_pseudo_ecalIso03_E015_mc->Draw("SAME HIST");
        h1_pseudo_ecalIso03_E015_data->Draw("SAME");
        saveCanvas(c1, saveName + "_pseudo_ecalIso03_E015_" + det, false);
        saveCanvas(c1, saveName + "_pseudo_ecalIso03_E015_" + det, true);

        h1_pseudo_ecalIso03_E010_data->Draw();
        h1_pseudo_ecalIso03_E010_mc->Scale(h1_pseudo_ecalIso03_E010_data->Integral()/h1_pseudo_ecalIso03_E010_mc->Integral());
        h1_pseudo_ecalIso03_E010_mc->Draw("SAME HIST");
        h1_pseudo_ecalIso03_E010_data->Draw("SAME");
        saveCanvas(c1, saveName + "_pseudo_ecalIso03_E010_" + det, false);
        saveCanvas(c1, saveName + "_pseudo_ecalIso03_E010_" + det, true);

        h1_pseudo_ecalIso03_recHitEt_data->Draw();
        h1_pseudo_ecalIso03_recHitEt_mc->Scale(h1_pseudo_ecalIso03_recHitEt_data->Integral()/h1_pseudo_ecalIso03_recHitEt_mc->Integral());
        h1_pseudo_ecalIso03_recHitEt_mc->Draw("SAME HIST");
        h1_pseudo_ecalIso03_recHitEt_data->Draw("SAME");
        saveCanvas(c1, saveName + "_pseudo_ecalIso03_recHitEt_" + det, false);
        saveCanvas(c1, saveName + "_pseudo_ecalIso03_recHitEt_" + det, true);

        h1_pseudo_ecalIso03_recHitE_data->Draw();
        h1_pseudo_ecalIso03_recHitE_mc->Scale(h1_pseudo_ecalIso03_recHitE_data->Integral()/h1_pseudo_ecalIso03_recHitE_mc->Integral());
        h1_pseudo_ecalIso03_recHitE_mc->Draw("SAME HIST");
        h1_pseudo_ecalIso03_recHitE_data->Draw("SAME");
        saveCanvas(c1, saveName + "_pseudo_ecalIso03_recHitE_" + det, false);
        saveCanvas(c1, saveName + "_pseudo_ecalIso03_recHitE_" + det, true);


        h1_pseudo_ecalIso03_recHitN_data->Draw();
        h1_pseudo_ecalIso03_recHitN_mc->Scale(h1_pseudo_ecalIso03_recHitN_data->Integral()/h1_pseudo_ecalIso03_recHitN_mc->Integral());
        h1_pseudo_ecalIso03_recHitN_mc->Draw("SAME HIST");
        h1_pseudo_ecalIso03_recHitN_data->Draw("SAME");
        saveCanvas(c1, saveName + "_pseudo_ecalIso03_recHitN_" + det, false);
        saveCanvas(c1, saveName + "_pseudo_ecalIso03_recHitN_" + det, true);

	h1_pseudo_hcalD1Iso03_data->GetXaxis()->SetTitle("HCAL Iso (D1) (GeV)");
	h1_pseudo_hcalD1Iso03_data->GetYaxis()->SetTitle("Number of Cones");
	h1_pseudo_hcalD1Iso03_data->Draw();
	h1_pseudo_hcalD1Iso03_mc->Scale(h1_pseudo_hcalD1Iso03_data->Integral()/h1_pseudo_hcalD1Iso03_mc->Integral());
	h1_pseudo_hcalD1Iso03_mc->Draw("SAME HIST");
	h1_pseudo_hcalD1Iso03_data->Draw("SAME");
	lg->Draw();
	saveCanvas(c1, saveName + "_pseudo_hcalD1Iso03_" + det, false);
	saveCanvas(c1, saveName + "_pseudo_hcalD1Iso03_" + det, true);

        h1_pseudo_hcalD1Iso03_data->GetYaxis()->SetRangeUser(0.001, histMax(h1_pseudo_hcalD1Iso03_data, h1_pseudo_hcalD1Iso03_mc));
	h1_pseudo_hcalD1Iso03_data->GetXaxis()->SetRangeUser(0.1, 10);
        h1_pseudo_hcalD1Iso03_data->Draw();
        h1_pseudo_hcalD1Iso03_mc->Draw("SAME HIST");
        h1_pseudo_hcalD1Iso03_data->Draw("SAME");
       	lg->Draw();
	saveCanvas(c1, saveName + "_pseudo_hcalD1Iso03_zoom_" + det, false);

        h1_pseudo_hcalD2Iso03_data->GetXaxis()->SetTitle("HCAL Iso (D2) (GeV)");
        h1_pseudo_hcalD2Iso03_data->GetYaxis()->SetTitle("Number of Cones");
	h1_pseudo_hcalD2Iso03_data->Draw();
	h1_pseudo_hcalD2Iso03_mc->Scale(h1_pseudo_hcalD2Iso03_data->Integral()/h1_pseudo_hcalD2Iso03_mc->Integral());
	h1_pseudo_hcalD2Iso03_mc->Draw("SAME HIST");
	h1_pseudo_hcalD2Iso03_data->Draw("SAME");
	lg->Draw();
	saveCanvas(c1, saveName + "_pseudo_hcalD2Iso03_" + det, false);
	saveCanvas(c1, saveName + "_pseudo_hcalD2Iso03_" + det, true);

        h1_pseudo_hcalD2Iso03_data->GetXaxis()->SetRangeUser(0.1, 10);
        h1_pseudo_hcalD2Iso03_data->Draw();
        h1_pseudo_hcalD2Iso03_mc->Draw("SAME HIST");
        h1_pseudo_hcalD2Iso03_data->Draw("SAME");
	lg->Draw();
        saveCanvas(c1, saveName + "_pseudo_hcalD2Iso03_zoom_" + det, false);

        h1_pseudo_tkIso03_data->GetXaxis()->SetTitle("Track Iso (GeV)");
        h1_pseudo_tkIso03_data->GetYaxis()->SetTitle("Number of Cones");
	h1_pseudo_tkIso03_data->GetYaxis()->SetRangeUser(0.01, h1_pseudo_tkIso03_data->GetMaximum()*1.5);
	h1_pseudo_tkIso03_data->Draw();
	h1_pseudo_tkIso03_mc->Scale(h1_pseudo_tkIso03_data->Integral()/h1_pseudo_tkIso03_mc->Integral());
	h1_pseudo_tkIso03_mc->Draw("SAME HIST");
	h1_pseudo_tkIso03_data->Draw("SAME");
	lg->Draw();
	saveCanvas(c1, saveName + "_pseudo_tkIso03_" + det, false);
	saveCanvas(c1, saveName + "_pseudo_tkIso03_" + det, true);

        h1_pseudo_tkIso03_data->GetYaxis()->SetRangeUser(0.001, 20);
        h1_pseudo_tkIso03_data->GetXaxis()->SetRangeUser(0.1, 3.0);
        h1_pseudo_tkIso03_data->Draw();
        h1_pseudo_tkIso03_mc->Draw("SAME HIST");
        h1_pseudo_tkIso03_data->Draw("SAME");
        lg->Draw(); 
        saveCanvas(c1, saveName + "_pseudo_tkIso03_zoom_" + det, false);

	h1_pseudo_dRClosestTower_data->Draw();
	h1_pseudo_dRClosestTower_mc->Scale(h1_pseudo_dRClosestTower_data->Integral()/h1_pseudo_dRClosestTower_mc->Integral());
	h1_pseudo_dRClosestTower_mc->Draw("SAME HIST");
	h1_pseudo_dRClosestTower_data->Draw("SAME");
	saveCanvas(c1, saveName + "_pseudo_dRClosestTower_" + det, false);


	delete h1_pseudo_ecalEta_mc;
	delete h1_pseudo_ecalEta_data;

	delete h1_pseudo_ecalPhi_mc;
	delete h1_pseudo_ecalPhi_data;

	delete h1_pseudo_ecalIso03_qual_mc;
	delete h1_pseudo_ecalIso03_qual_data;

        delete h1_pseudo_ecalIso03_fakeSR_mc;
        delete h1_pseudo_ecalIso03_fakeSR_data;

	delete h1_pseudo_ecalIso03_mc;
	delete h1_pseudo_ecalIso03_data;

        delete h1_pseudo_ecalIso03_E015_mc;
        delete h1_pseudo_ecalIso03_E015_data;

        delete h1_pseudo_ecalIso03_E010_mc;
        delete h1_pseudo_ecalIso03_E010_data;

        delete h1_pseudo_ecalIso03_recHitEt_mc;
        delete h1_pseudo_ecalIso03_recHitEt_data;

        delete h1_pseudo_ecalIso03_recHitE_mc;
        delete h1_pseudo_ecalIso03_recHitE_data;

        delete h1_pseudo_ecalIso03_recHitN_mc;
        delete h1_pseudo_ecalIso03_recHitN_data;

	delete h1_pseudo_hcalD1Iso03_mc;
	delete h1_pseudo_hcalD1Iso03_data;

	delete h1_pseudo_hcalD2Iso03_mc;
	delete h1_pseudo_hcalD2Iso03_data;

	delete h1_pseudo_tkIso03_mc;
	delete h1_pseudo_tkIso03_data;

	delete h1_pseudo_dRClosestTower_mc;
	delete h1_pseudo_dRClosestTower_data;


}

void compare(TString ntuple_mc, TString ntuple_data, TString saveName)
{

	// set up style
	//
	gROOT->ProcessLine(".L tdrStyle.C");
	gROOT->ProcessLine("setTDRStyle()");
	//gStyle->SetOptStat(111111);

	// open file
	//
	TFile *f_mc = new TFile(ntuple_mc, "READ");
	TFile *f_data = new TFile(ntuple_data, "READ");
	gROOT->cd();

	compareOne(f_mc, f_data, saveName, "ee");
	compareOne(f_mc, f_data, saveName, "eb");

	compareOneEcalSpecial(f_mc, f_data, saveName, "eep");
        compareOneEcalSpecial(f_mc, f_data, saveName, "eem");
        compareOneEcalSpecial(f_mc, f_data, saveName, "ebp");
        compareOneEcalSpecial(f_mc, f_data, saveName, "ebm");

	// compare + and - ecal
	compareDataDataEcalSpecial(f_data, f_mc, saveName);

	delete f_mc;
	delete f_data;
}


