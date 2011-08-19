#include "Utilities.h"
#include "processDYEstResults.C"

void printR(TString nJets)
{

        gROOT->ProcessLine(".L ~/tdrStyle.C");
        gROOT->ProcessLine("setTDRStyle()");
        gStyle->SetOptTitle(1);

        TFile f("DYEstResults_ForWW_MET45_INCL.root");

	gROOT->cd();

        TH1F *h1_R_mm = getRHist(f, "dymm", nJets, "mm", true);
                h1_R_mm->SetMarkerStyle(24);
                h1_R_mm->SetLineColor(kBlue);
                h1_R_mm->SetLineWidth(2);
                h1_R_mm->GetXaxis()->SetTitle("MET Cut (GeV)");
        TH1F *h1_R_ee = getRHist(f, "dyee", nJets, "ee", true);
                h1_R_ee->SetMarkerStyle(26);
                h1_R_ee->SetLineColor(kRed);
                h1_R_ee->SetLineWidth(2);
        h1_R_mm->SetTitle("h1_R DY " + nJets);

        TCanvas *c0 = new TCanvas();
        c0->cd();
        h1_R_mm->GetXaxis()->SetRangeUser(20, 60);
        h1_R_mm->GetYaxis()->SetRangeUser(0, 1.5);
        h1_R_mm->Draw("HIST E1");
        h1_R_ee->Draw("Same HIST E1");

        TLegend *l1 = new TLegend(0.2, 0.7, 0.6, 0.9);
        l1->SetFillColor(kWhite);
        l1->SetLineColor(kWhite);
        l1->AddEntry(h1_R_mm, "dymm", "lp");
        l1->AddEntry(h1_R_ee, "dyee", "lp");
        l1->Draw();

        // look at ZZ
        TH1F *h1_R_ZZee = getRHist(f, "zz", nJets, "ee", false);
                h1_R_ZZee->SetMarkerStyle(28);
                h1_R_ZZee->SetLineColor(kRed);
                h1_R_ZZee->SetLineWidth(2);
        TH1F *h1_R_ZZmm = getRHist(f, "zz", nJets, "mm", false);
                h1_R_ZZmm->SetMarkerStyle(30);
                h1_R_ZZmm->SetLineColor(kBlue);
                h1_R_ZZmm->SetLineWidth(2);
        h1_R_ZZmm->SetTitle("h1_R ZZ " + nJets);

        TCanvas *c1 = new TCanvas();
        c1->cd();
        h1_R_ZZmm->GetXaxis()->SetRangeUser(20, 60);
        h1_R_ZZmm->GetYaxis()->SetRangeUser(0, 1.5);
        h1_R_ZZmm->Draw("HIST E1");
        h1_R_ZZee->Draw("HIST SAME E1");

        // look at WZ
        TH1F *h1_R_WZmm = getRHist(f, "wz", nJets, "mm", true);
                h1_R_WZmm->SetMarkerStyle(30);
                h1_R_WZmm->SetLineColor(kBlue);
                h1_R_WZmm->SetLineWidth(2);
        TH1F *h1_R_WZee = getRHist(f, "wz", nJets, "ee", true);
                h1_R_WZee->SetMarkerStyle(32);
                h1_R_WZee->SetLineColor(kRed);
                h1_R_WZee->SetLineWidth(2);
        h1_R_WZmm->SetTitle("h1_R WZ " + nJets);

        TCanvas *c2 = new TCanvas();
        c2->cd();
        h1_R_WZmm->GetXaxis()->SetRangeUser(20, 60);
        h1_R_WZmm->GetYaxis()->SetRangeUser(0, 1.5);
        h1_R_WZmm->Draw("HIST E1");
        h1_R_WZee->Draw("SAME HIST E1");

}

