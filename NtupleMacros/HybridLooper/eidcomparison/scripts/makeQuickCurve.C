
#include "TGraph.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TColor.h"
#include "TLegend.h"

#include "../../../Tools/histtools.cc"

void formatGraph(TGraph *g1, Int_t marker, Color_t color)
{
        g1->SetLineColor(color);
        g1->SetMarkerColor(color);
        g1->SetMarkerStyle(marker);
}

void formatHist(TH1F *h1, bool signal)
{
    if (signal) {
        h1->SetFillStyle(3002);
        h1->SetFillColor(kBlue);
        h1->SetLineColor(kBlue);
        h1->SetMarkerColor(kBlue);
    } else {
        h1->SetFillStyle(3002);
        h1->SetFillColor(kRed);
        h1->SetLineColor(kRed);
        h1->SetMarkerColor(kRed);
    }
}

void makeQuickCurve(TString extra_denominator)
{

    TFile f_signal("../baby_hww130.root", "READ");
    TFile f_bg("../baby_wjets.root", "READ");
    TTree *t_signal = (TTree*)f_signal.Get("tree");
    TTree *t_bg = (TTree*)f_bg.Get("tree");
    gROOT->cd();

    //
    // pf MVA
    //

    TH1F *h1_signal_pfmva = new TH1F("h1_signal_pfmva", "signal_pfmva", 500, -1.0, 1.0);
    TH1F *h1_bg_pfmva = new TH1F("h1_bg_pfmva", "bg_pfmva", 500, -1.0, 1.0);
    formatHist(h1_signal_pfmva, true);
    formatHist(h1_bg_pfmva, false);

    //
    // egamma LH
    //

    TH1F *h1_signal_lh = new TH1F("h1_signal_lh", "signal_lh", 500, 0.0, 1.0);
    TH1F *h1_bg_lh = new TH1F("h1_bg_lh", "bg_lh", 500, 0.0, 1.0);
    formatHist(h1_signal_lh, true);
    formatHist(h1_bg_lh, false);

    //
    // VBTF
    //
    TH1F *h1_signal_vbtf95 = new TH1F("h1_signal_vbtf95", "signal_vbtf95", 500, 0.0, 1.0);
    TH1F *h1_bg_vbtf95 = new TH1F("h1_bg_vbtf95", "bg_vbtf95", 500, 0.0, 1.0);
    formatHist(h1_signal_vbtf95, true);
    formatHist(h1_bg_vbtf95, false);

    TH1F *h1_signal_vbtf90 = new TH1F("h1_signal_vbtf90", "signal_vbtf90", 500, 0.0, 1.0);
    TH1F *h1_bg_vbtf90 = new TH1F("h1_bg_vbtf90", "bg_vbtf90", 500, 0.0, 1.0);
    formatHist(h1_signal_vbtf90, true);
    formatHist(h1_bg_vbtf90, false);

    TH1F *h1_signal_vbtf80 = new TH1F("h1_signal_vbtf80", "signal_vbtf80", 500, 0.0, 1.0);
    TH1F *h1_bg_vbtf80 = new TH1F("h1_bg_vbtf80", "bg_vbtf80", 500, 0.0, 1.0);
    formatHist(h1_signal_vbtf80, true);
    formatHist(h1_bg_vbtf80, false);

    TH1F *h1_signal_vbtf70 = new TH1F("h1_signal_vbtf70", "signal_vbtf70", 500, 0.0, 1.0);
    TH1F *h1_bg_vbtf70 = new TH1F("h1_bg_vbtf70", "bg_vbtf70", 500, 0.0, 1.0);
    formatHist(h1_signal_vbtf70, true);
    formatHist(h1_bg_vbtf70, false);

    //
    // set up denominators
    //

    TString denominator = "abs(reco_typelt) == 11 && " + extra_denominator;
    TString denominator_signal = "35/1000.0*weight*("+ denominator + " && abs(reco_mctypelt) == 11)";
    TString denominator_bg = "35/1000.0*weight*("+ denominator + " && abs(reco_mctypelt) != 11)";

    t_signal->Draw("reco_pfmvalt >> h1_signal_pfmva", denominator_signal);
    t_bg->Draw("reco_pfmvalt >> h1_bg_pfmva", denominator_bg);
    TGraph *effrej_pfmva = (TGraph*)eff_rej(*h1_signal_pfmva, *h1_bg_pfmva, false, false).Clone("effrej_pfmva");
    formatGraph(effrej_pfmva, 20, kBlack);

    t_signal->Draw("reco_lhlt >> h1_signal_lh", denominator_signal);
    t_bg->Draw("reco_lhlt >> h1_bg_lh", denominator_bg);
    TGraph *effrej_lh = (TGraph*)eff_rej(*h1_signal_lh, *h1_bg_lh, false, false).Clone("effrej_lh");
    formatGraph(effrej_lh, 20, kBlue);

    t_signal->Draw("reco_vbtf95lt >> h1_signal_vbtf95", denominator_signal);
    t_bg->Draw("reco_vbtf95lt >> h1_bg_vbtf95", denominator_bg);
    TGraph *effrej_vbtf95 = (TGraph*)eff_rej(*h1_signal_vbtf95, *h1_bg_vbtf95, false, false).Clone("effrej_vbtf95");
    formatGraph(effrej_vbtf95, 22, kRed);

    t_signal->Draw("reco_vbtf90lt >> h1_signal_vbtf90", denominator_signal);
    t_bg->Draw("reco_vbtf90lt >> h1_bg_vbtf90", denominator_bg);
    TGraph *effrej_vbtf90 = (TGraph*)eff_rej(*h1_signal_vbtf90, *h1_bg_vbtf90, false, false).Clone("effrej_vbtf90");
    formatGraph(effrej_vbtf90, 22, kRed + 1);

    t_signal->Draw("reco_vbtf80lt >> h1_signal_vbtf80", denominator_signal);
    t_bg->Draw("reco_vbtf80lt >> h1_bg_vbtf80", denominator_bg);
    TGraph *effrej_vbtf80 = (TGraph*)eff_rej(*h1_signal_vbtf80, *h1_bg_vbtf80, false, false).Clone("effrej_vbtf80");
    formatGraph(effrej_vbtf80, 22, kRed + 2);

    t_signal->Draw("reco_vbtf70lt >> h1_signal_vbtf70", denominator_signal);
    t_bg->Draw("reco_vbtf70lt >> h1_bg_vbtf70", denominator_bg);
    TGraph *effrej_vbtf70 = (TGraph*)eff_rej(*h1_signal_vbtf70, *h1_bg_vbtf70, false, false).Clone("effrej_vbtf70");
    formatGraph(effrej_vbtf70, 22, kRed + 3);

    TCanvas *c_pfmva = new TCanvas();
    c_pfmva->cd();
    h1_bg_pfmva->Draw();
    h1_signal_pfmva->Draw("SAME");

    TCanvas *c_lh = new TCanvas();
    c_lh->cd();
    h1_bg_lh->Draw();
    h1_signal_lh->Draw("SAME");

    TLegend *l1 = new TLegend(0.13, 0.47, 0.23, 0.87);
    l1->AddEntry(effrej_pfmva, "PF MVA", "lp");
    l1->AddEntry(effrej_lh, "e/#gamma LH", "lp");
    l1->AddEntry(effrej_vbtf95, "VBTF 95", "lp");
    l1->AddEntry(effrej_vbtf90, "VBTF 90", "lp");
    l1->AddEntry(effrej_vbtf80, "VBTF 80", "lp");
    l1->AddEntry(effrej_vbtf70, "VBTF 70", "lp");

    TCanvas *c_effrej = new TCanvas();
    c_effrej->cd();
    effrej_pfmva->Draw("AP");
    effrej_pfmva->GetYaxis()->SetTitle("Number of W+Jets");
    effrej_pfmva->GetXaxis()->SetTitle("Number of H->WW (M130)");
    effrej_pfmva->SetTitle(extra_denominator);
    effrej_lh->Draw("P");
    effrej_vbtf95->Draw("P");
    effrej_vbtf90->Draw("P");
    effrej_vbtf80->Draw("P");
    effrej_vbtf70->Draw("P");
    l1->Draw();

}

