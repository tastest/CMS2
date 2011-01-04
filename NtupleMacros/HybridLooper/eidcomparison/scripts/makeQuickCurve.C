
#include "TGraph.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TColor.h"

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

void makeQuickCurve()
{

    TFile f_signal("../baby_lm0.root", "READ");
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
    // set up denominators
    //

    TString denominator = "abs(reco_typelt) == 11 && reco_ptlt > 15.0 && ((reco_algolt&(1<<2))==(1<<2))";
    TString denominator_signal = "weight*("+ denominator + " && abs(reco_mctypelt) == 11)";
    TString denominator_bg = "weight*("+ denominator + " && abs(reco_mctypelt) != 11)";

    t_signal->Draw("reco_pfmvalt >> h1_signal_pfmva", denominator_signal);
    t_bg->Draw("reco_pfmvalt >> h1_bg_pfmva", denominator_bg);
    TGraph *effrej_pfmva = (TGraph*)eff_rej(*h1_signal_pfmva, *h1_bg_pfmva, false, false).Clone("effrej_pfmva");
    formatGraph(effrej_pfmva, 20, kBlack);

    t_signal->Draw("reco_lhlt >> h1_signal_lh", denominator_signal);
    t_bg->Draw("reco_lhlt >> h1_bg_lh", denominator_bg);
    TGraph *effrej_lh = (TGraph*)eff_rej(*h1_signal_lh, *h1_bg_lh, false, false).Clone("effrej_lh");
    formatGraph(effrej_lh, 20, kBlue);

    TCanvas *c_pfmva = new TCanvas();
    c_pfmva->cd();
    h1_bg_pfmva->Draw();
    h1_signal_pfmva->Draw("SAME");

    TCanvas *c_lh = new TCanvas();
    c_lh->cd();
    h1_bg_lh->Draw();
    h1_signal_lh->Draw("SAME");

    TCanvas *c_effrej = new TCanvas();
    c_effrej->cd();
    effrej_pfmva->Draw("AP");
    effrej_lh->Draw("P");

}

