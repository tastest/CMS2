#include "PlotUtilities.h"

#include "TROOT.h"
#include "TFile.h"
#include "TList.h"
#include "TObject.h"
#include "TObjArray.h"
#include "TH2.h"
#include "TKey.h"
#include "TRegexp.h"
#include <TIterator.h>
#include "TClass.h"
#include "TLegend.h"

#include <iostream>
#include <cmath>

H cumulate (const H &in, bool increasing)
{
    H h_out(in.GetName() + TString("tmp"), in.GetTitle(), in.GetNbinsX(),
            in.GetBinLowEdge(1), in.GetBinLowEdge(in.GetNbinsX() + 1));
    h_out.Sumw2();
    h_out.SetFillColor(in.GetFillColor());
    h_out.SetFillStyle(in.GetFillStyle());
    h_out.SetLineStyle(in.GetLineStyle());
    h_out.SetLineColor(in.GetLineColor());
    h_out.SetLineWidth(in.GetLineWidth());
    double sum = 0;
    double err2 = 0;
    if (increasing) {
        for (int j = 0; j <= in.GetNbinsX() + 1; ++j) {
            sum += in.GetBinContent(j);
            err2 += in.GetBinError(j)*in.GetBinError(j);
            h_out.SetBinContent(j, sum);
            h_out.SetBinError(j, sqrt(err2));
        }
    } else {
        for (int j = in.GetNbinsX() + 1; j >= 0; --j) {
            sum += in.GetBinContent(j);
            err2 += in.GetBinError(j)*in.GetBinError(j);
            h_out.SetBinContent(j, sum);
            h_out.SetBinError(j, sqrt(err2));
        }
    }
    return h_out;
}

TGraph eff_rej (const H &signal, H &background, bool normalize, bool increasing)
{
    H sig = *(TH1F*)signal.Clone("h_tmp_s");
    if (normalize)
        sig.Scale(1 / sig.Integral(0, sig.GetNbinsX() + 1));
    H bg = *(TH1F*)background.Clone("h_tmp_bg");
    if (normalize)
        bg.Scale(1 / bg.Integral(0, bg.GetNbinsX() + 1));
    H sig_cum = cumulate(sig, increasing);
    H bg_cum = cumulate(bg, increasing);
    TGraph ret(signal.GetNbinsX());
    for (int i = 1; i <= signal.GetNbinsX(); ++i) {
        const double x = sig_cum.GetBinCenter(i);
        const double s = sig_cum.GetBinContent(i);
        const double b = bg_cum.GetBinContent(bg_cum.FindBin(x));
        ret.SetPoint(i - 1, b, s); // gotta love offsets
        //        printf("point %d: %f s, %f b\n", i, s, b);
    }
    return ret;
}

void deleteHistos() {
    // Delete all existing histograms in memory
    TObject* obj;
    TList* list = gDirectory->GetList() ;
    TIterator* iter = list->MakeIterator();
    while ( (obj = (iter->Next())) ) {
        if (obj->IsA()->InheritsFrom(TH1::Class()) ||
                obj->IsA()->InheritsFrom(TH2::Class()) ) {
                delete obj;
            }
    }
}

void saveHist(const char* filename, const char* pat)
{

    TList* list = gDirectory->GetList() ;
    TIterator* iter = list->MakeIterator();
    TRegexp re(pat,kTRUE) ;

    TFile outf(filename,"RECREATE") ;
    printf("[PlotUtilities::saveHist] Saving histograms to %s\n", filename);

    int counter = 0;
    while(TObject* obj=iter->Next()) {

        // don't save TH1Keys objects
        // this is a bug fudge
        if (TString(obj->GetName()).Contains("histokeys")) continue;

        // save other stuff
        if (TString(obj->GetName()).Index(re)>=0) {
            obj->Write() ;
            ++counter;
        }

    }

    printf("[PlotUtilities::saveHist] Saved %i histograms\n", counter);
    outf.Close() ;
    delete iter ;
}

TCanvas *ComparePlots(TFile *f, const char *hist1, const char *hist2, const char *label1, const char *label2, unsigned int rebin, bool norm, bool log, unsigned int opt)
{

    // get hists
    TH1F *h1 = (TH1F*)f->Get(hist1)->Clone(hist1);
    TH1F *h2 = (TH1F*)f->Get(hist2)->Clone(hist2);

    // overflow
    h1->SetBinContent(h1->GetNbinsX(), h1->GetBinContent(h1->GetNbinsX()) + h1->GetBinContent(h1->GetNbinsX()+1));
    h1->SetBinContent(h1->GetNbinsX()+1, 0.0);
    h2->SetBinContent(h2->GetNbinsX(), h2->GetBinContent(h2->GetNbinsX()) + h2->GetBinContent(h2->GetNbinsX()+1));
    h2->SetBinContent(h2->GetNbinsX()+1, 0.0);

    // rebin
    h1->Rebin(rebin);
    h2->Rebin(rebin);

    // normalize
    if (norm) h1->Scale(h2->Integral(0, h2->GetNbinsX() + 1) / h1->Integral(0, h1->GetNbinsX() + 1));

    // format
    if (opt == 1) {
        h1->SetLineWidth(2);
        h1->SetLineColor(kRed);
        h2->SetLineColor(kBlue);
        h2->SetFillColor(kBlue);
        h2->SetFillStyle(0);
        h2->SetLineWidth(2);
        h2->SetMarkerStyle(20);
        h2->SetMarkerColor(kBlue);
    } 
    if (opt == 2) {
        h1->SetLineWidth(2);
        h1->SetLineColor(kRed);
        h1->SetFillColor(kRed);
        h1->SetFillStyle(3004);
        h2->SetFillStyle(0);
        h2->SetLineWidth(2);
        h2->SetMarkerStyle(20);
        h2->SetMarkerColor(kBlack);
        h2->SetLineColor(kBlack);
    }

    // legend
    TLegend *l1 = new TLegend(0.70, 0.85, 0.98, 0.98);
    l1->SetLineColor(kWhite);
    l1->SetFillColor(kWhite);
    l1->SetShadowColor(kWhite);
    l1->AddEntry(h1, label1, "f");
    l1->AddEntry(h2, label2, "f");

    // draw
    TCanvas *c1 = new TCanvas();
    if (log) c1->SetLogy();
    c1->cd();
    h1->Draw("HIST");
    h2->Draw("SAME E1");
    double yMax = TMath::Max(h1->GetMaximum(), h2->GetMaximum());
    h1->GetYaxis()->SetRangeUser(0.0, TMath::Max(yMax + 2*sqrt(yMax), yMax + 0.3*yMax));
    if (log) h1->GetYaxis()->SetRangeUser(0.1, yMax * 1000);
    l1->Draw();

    return c1;

}

TGraph GetROC(TFile *f, const char *hist1, const char *hist2, bool increasing)
{

    // get hists
    TH1F *sig = (TH1F*)f->Get(hist1)->Clone(hist1);
    TH1F *bkg = (TH1F*)f->Get(hist2)->Clone(hist2);
    
    // overflow
    sig->SetBinContent(sig->GetNbinsX(), sig->GetBinContent(sig->GetNbinsX()) + sig->GetBinContent(sig->GetNbinsX()+1));
    sig->SetBinContent(sig->GetNbinsX()+1, 0.0);
    bkg->SetBinContent(bkg->GetNbinsX(), bkg->GetBinContent(bkg->GetNbinsX()) + bkg->GetBinContent(bkg->GetNbinsX()+1));
    bkg->SetBinContent(bkg->GetNbinsX()+1, 0.0);

    TGraph gr_roc = eff_rej(*sig, *bkg, true, increasing);
    return gr_roc;
}


TGraph GetEff(TFile *f, const char *hist1, bool increasing)
{

    // get hists
    TH1F *sig = (TH1F*)f->Get(hist1)->Clone(hist1);

    // overflow
    sig->SetBinContent(sig->GetNbinsX(), sig->GetBinContent(sig->GetNbinsX()) + sig->GetBinContent(sig->GetNbinsX()+1));
    sig->SetBinContent(sig->GetNbinsX()+1, 0.0);

    TGraph ret(sig->GetNbinsX());
    for (int i = 1; i <= sig->GetNbinsX(); ++i) {
        float val = sig->GetBinLowEdge(i);
        float ntotal = sig->Integral(1, sig->GetNbinsX());
        float npass = sig->Integral(i, sig->GetNbinsX());
        if (increasing) npass = ntotal - npass;
        ret.SetPoint(i - 1, val, npass/ntotal);
    }

    return ret;
}


