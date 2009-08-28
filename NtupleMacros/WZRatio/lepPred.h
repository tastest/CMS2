#include "TCanvas.h"
#include "TH1.h"
#include "TFile.h"
#include "TFile.h"
#include "Tools/HistogramUtilities.h"
#include "MetRes/plotMetRes.h"

void lepPred () 
{
     TFile *f_de = TFile::Open("Results_de.root");
     TFile *f_dm = TFile::Open("Results_dm.root");
     TFile *f_se = TFile::Open("Results_se.root");
     TFile *f_sm = TFile::Open("Results_sm.root");
     
     TH1F *h_ttbar_njets_de = (TH1F *)f_de->Get("ttbar_njets_all");
     TH1F *h_ttbar_njets_dm = (TH1F *)f_dm->Get("ttbar_njets_all");
     TH1F *h_ttbar_njets_se = (TH1F *)f_se->Get("ttbar_njets_all");
     TH1F *h_ttbar_njets_sm = (TH1F *)f_sm->Get("ttbar_njets_all");

     TCanvas *c_se = new TCanvas("c_se", "electron");
     TCanvas *c_sm = new TCanvas("c_sm", "muon");
     c_se->SetLogy();
     c_sm->SetLogy();

     TH1F *h_ttbar_njets_se_pred = (TH1F *)h_ttbar_njets_de->Clone();
     TH1F *h_ttbar_njets_sm_pred = (TH1F *)h_ttbar_njets_dm->Clone();
     h_ttbar_njets_se_pred->Clear();
     h_ttbar_njets_sm_pred->Clear();
     for (int i = 0; i <= h_ttbar_njets_de->GetNbinsX() + 1; ++i) {
	   if (i < 3) {
		 h_ttbar_njets_se_pred->SetBinContent(i, 0);
		 h_ttbar_njets_sm_pred->SetBinContent(i, 0);
	   }
	   h_ttbar_njets_se_pred->SetBinContent(i, h_ttbar_njets_de->GetBinContent(i - 2));
	   h_ttbar_njets_sm_pred->SetBinContent(i, h_ttbar_njets_dm->GetBinContent(i - 2));
     }
     
     h_ttbar_njets_se->SetLineColor(kBlack);
     h_ttbar_njets_sm->SetLineColor(kBlack);
     h_ttbar_njets_se_pred->SetLineColor(kRed);
     h_ttbar_njets_sm_pred->SetLineColor(kRed);
     c_se->cd();
     h_ttbar_njets_se->Divide(h_ttbar_njets_se_pred);
     h_ttbar_njets_se->Draw();
	 //h_ttbar_njets_se_pred->Draw("same");
     c_sm->cd();
     h_ttbar_njets_sm->Divide(h_ttbar_njets_sm_pred);
     h_ttbar_njets_sm->Draw();
	 //h_ttbar_njets_sm_pred->Draw("same");

	 sources_t theSources = sources_wzratio;
	 HistogramUtilities* hsinglee = new HistogramUtilities("Results_se.root");
	 HistogramUtilities* hsinglem = new HistogramUtilities("Results_sm.root");
	 HistogramUtilities* hsingle = new HistogramUtilities("Results_se.root", "Results_sm.root");
	 HistogramUtilities* hdoublee = new HistogramUtilities("Results_de.root");
	 HistogramUtilities* hdoublem = new HistogramUtilities("Results_dm.root");
	 HistogramUtilities* hdouble = new HistogramUtilities("Results_de.root", "Results_dm.root");
	 HistogramUtilities* hdoubleem = new HistogramUtilities("Results_em.root");
	 //HistogramUtilities* hdouble3 = new HistogramUtilities("Results_de.root", "Results_dm.root", "Results_em.root");

	 //Need this in order to not store clones, which could have same names in 2 files, thus breaking file_->Get()
	 TH1::AddDirectory(false); 

	 TLegend *lg_all = hsingle->getLegend(theSources, "njets", "", "all");
	 // TLegend *lg_all = hdoubleem->getLegend(theSources, "njets", "", "all");

	 make2fileStack(hsingle, lg_all, theSources, "njets", "", "all", true, "njets_all_e+m");
	 makeStack(hsinglee, lg_all, theSources, "njets", "", "all", true, "njets_all_e");
	 makeStack(hsinglem, lg_all, theSources, "njets", "", "all", true, "njets_all_m");
	 
	 make2fileStack(hdouble, lg_all, theSources, "njets", "", "all", true, "njets_all_ee+mm");
	 makeStack(hdoublee, lg_all, theSources, "njets", "", "all", true, "njets_all_ee");
	 makeStack(hdoublem, lg_all, theSources, "njets", "", "all", true, "njets_all_mm");

	 makeStack(hdoubleem, lg_all, theSources, "njets", "", "all", true, "njets_all_em");

	 //make3fileStack(hdouble3, lg_all, theSources, "njets", "", "all", true, "");

	 //THStack *s = hsingle->get2fileStack(theSources, "njets", "", "all");
	 //TCanvas *c_st = new TCanvas(s->GetName(), s->GetName());
	 //s->Draw();
	 //c_st->SaveAs((TString)c_st->GetName()+".png");
}
