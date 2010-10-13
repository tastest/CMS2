#include "TCanvas.h"
#include "TH1.h"
#include "TFile.h"
#include "TFile.h"

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
//      h_ttbar_njets_se_pred->Draw("same");
     c_sm->cd();
     h_ttbar_njets_sm->Divide(h_ttbar_njets_sm_pred);
     h_ttbar_njets_sm->Draw();
//      h_ttbar_njets_sm_pred->Draw("same");
}
     
     
