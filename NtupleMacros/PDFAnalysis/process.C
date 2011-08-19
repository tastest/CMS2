#include "TFile.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TLegend.h"

#include <string>
#include <math.h>
#include <iostream>

enum PDFFamily {
    PDF_CTEQ,
    PDF_CTEQ_AS,
    PDF_MSTW,
};

void FormatLegend(TLegend *lg) {
    lg->SetLineColor(kWhite);
    lg->SetShadowColor(kWhite);
    lg->SetFillColor(kWhite);
}

//
// CT10 alpha_S variation analysis
//

void AnalyseCT10as(TFile *f, std::string sampleName, TH1F *h1, TH1F *h1_up, TH1F *h1_down)
{

    TH1F *histArr[10];
    for (unsigned int i = 0; i < 10; ++i) {
        histArr[i] = (TH1F*)f->Get(Form("%s_CT10as_%i", sampleName.c_str(), i));
    }

    for (Int_t bin = 1; bin <= h1->GetNbinsX(); ++bin)
    {

        const unsigned int max_set = 3;
        Double_t X0 = histArr[5]->GetBinContent(bin);
        Double_t binCenter = histArr[5]->GetBinCenter(bin);
        if (X0 == 0) continue;
        Double_t Xi_up = histArr[10-max_set]->GetBinContent(bin);
        Double_t Xi_down = histArr[max_set]->GetBinContent(bin);
        Double_t plus_max = sqrt(pow(TMath::Max(TMath::Max(Xi_up - X0, Xi_down - X0), 0.0), 2))/1.645;
        Double_t minus_max = sqrt(pow(TMath::Max(TMath::Max(X0 - Xi_up, X0 - Xi_down), 0.0), 2))/1.645;
        h1->SetBinContent(bin, X0);
        h1_down->SetBinContent(bin, -1*minus_max/X0);
        h1_up->SetBinContent(bin, plus_max/X0);
    }

    for (unsigned int i = 0; i < 10; ++i) {
        delete histArr[i];
    }
}

//
// MSTW/CT10 intrinsic uncertainty analysis
//

void AnalyseCTMSTW(TFile *f, std::string sampleName, std::string pdfName, const unsigned int &nsets, TH1F *h1, TH1F *h1_up, TH1F *h1_down)
{

    TH1F *histArr[nsets];
    for (unsigned int i = 0; i < nsets; ++i) {
        histArr[i] = (TH1F*)f->Get(Form("%s_%s_%i", sampleName.c_str(), pdfName.c_str(), i));
    }

    for (Int_t bin = 1; bin <= h1->GetNbinsX(); ++bin)
    {

        Double_t X0 = histArr[0]->GetBinContent(bin);
        Double_t binCenter = histArr[0]->GetBinCenter(bin);
        Double_t plus_max = 0.0;
        Double_t minus_max = 0.0;

        if (X0 == 0) continue;

        for (unsigned int subset = 0; subset < ((nsets - 1)/2); ++subset) 
        {
            Double_t Xi_up = histArr[(subset*2) + 1]->GetBinContent(bin);
            Double_t Xi_down = histArr[(subset*2) + 2]->GetBinContent(bin);
            plus_max += pow(TMath::Max(TMath::Max(Xi_up - X0, Xi_down - X0), 0.0), 2);
            minus_max += pow(TMath::Max(TMath::Max(X0 - Xi_up, X0 - Xi_down), 0.0), 2);
        }
        plus_max = sqrt(plus_max);
        minus_max = sqrt(minus_max);
        h1->SetBinContent(bin, X0);
        h1_down->SetBinContent(bin, -1*minus_max/X0);
        h1_up->SetBinContent(bin, plus_max/X0);
    }

    for (unsigned int i = 0; i < nsets; ++i) {
        delete histArr[i];
    }
}

//
// NNPDF intrinsic uncertainty analysis
//

void AnalyseNNPDF(TFile *f, std::string sampleName, TH1F *h1, TH1F *h1_up, TH1F *h1_down)
{

    //
    // first get the average value of the set of replicas with central value alpha_S
    //

    TH1F *h1_temp = (TH1F*)f->Get(Form("%s_NNPDF20_100_0", sampleName.c_str()));
    for (Int_t bin = 1; bin <= h1->GetNbinsX(); ++bin) {
        h1->SetBinContent(bin, h1_temp->GetBinContent(bin));
    }

    //
    // now read in the sets with larger alpha_s
    // and also smaller
    //

    std::vector<TH1F*> histArr;
    unsigned int nsets_up = 0;
    unsigned int nsets_down = 0;
    for (unsigned int i = 1; i < 101; ++i)
        histArr.push_back((TH1F*)f->Get(Form("%s_NNPDF20_100_%i", sampleName.c_str(), i)));
    for (unsigned int i = 0; i < 72; ++i) {
        histArr.push_back((TH1F*)f->Get(Form("%s_NNPDF20_as_0120_100_%i", sampleName.c_str(), i)));
        histArr.push_back((TH1F*)f->Get(Form("%s_NNPDF20_as_0118_100_%i", sampleName.c_str(), i)));
    }
    for (unsigned int i = 0; i < 27; ++i) {
        histArr.push_back((TH1F*)f->Get(Form("%s_NNPDF20_as_0121_100_%i", sampleName.c_str(), i)));
        histArr.push_back((TH1F*)f->Get(Form("%s_NNPDF20_as_0117_100_%i", sampleName.c_str(), i)));
    }
    for (unsigned int i = 0; i < 5; ++i) {
        histArr.push_back((TH1F*)f->Get(Form("%s_NNPDF20_as_0122_100_%i", sampleName.c_str(), i)));
        histArr.push_back((TH1F*)f->Get(Form("%s_NNPDF20_as_0116_100_%i", sampleName.c_str(), i)));
    }

    //
    // make the bands
    //

    for (unsigned int bin = 0; bin < h1->GetNbinsX() + 1; ++bin) {

        float X0 = h1->GetBinContent(bin);
        if (X0 == 0) continue;

        float delta_up = 0.0;
        float delta_down = 0.0;
        for (unsigned int i = 0; i < histArr.size(); ++i) {
            if (histArr[i]->GetBinContent(bin) > h1->GetBinContent(bin)) {
                delta_up += pow(histArr[i]->GetBinContent(bin) - h1->GetBinContent(bin), 2);
                ++nsets_up;
            }
            else {
                delta_down += pow(h1->GetBinContent(bin) - histArr[i]->GetBinContent(bin), 2);
                ++nsets_down;
            }
        }

        delta_up = sqrt( (1.0/(float(nsets_up) - 1.0)) * delta_up);
        delta_down = sqrt( (1.0/(float(nsets_down) - 1.0)) * delta_down);
        h1_up->SetBinContent(bin, delta_up/X0);
        h1_down->SetBinContent(bin, -1*delta_down/X0);
    }

    //
    // clean up
    //

    delete h1_temp;
    for (unsigned int i = 0; i < histArr.size(); ++i) delete histArr[i];

    //
    // now compare to the sets with smaller and larger alpha_S
    //

}

//std::string sampleName = "nn_hww115_ww";
void process(std::string sampleName)
{
    gROOT->ProcessLine(".L ~/tdrStyle.C");
    gROOT->ProcessLine("setTDRStyle()");

    TFile *f = new TFile("histos_mc.root", "READ");
    gROOT->cd();

    Int_t nbins = 40;
    Float_t min = -1.0;
    Float_t max = 1.0;

    //
    // CTEQ
    //

    // central
    TH1F *h1_CT10   = new TH1F(Form("%s_h1_CT10", sampleName.c_str()),  "centre", nbins, min, max);
    TH1F *h1_CT10_up   = new TH1F(Form("%s_h1_CT10_up", sampleName.c_str()),  "centre up", nbins, min, max);
    TH1F *h1_CT10_down  = new TH1F(Form("%s_h1_CT10_down", sampleName.c_str()), "centre down", nbins, min, max);
    AnalyseCTMSTW(f, sampleName, "CT10", 53, h1_CT10, h1_CT10_up, h1_CT10_down);

    // alpha_s
    TH1F *h1_CT10_alpha_S   = new TH1F(Form("%s_h1_CT10_alpha_S", sampleName.c_str()),  "alpha_S", nbins, min, max);
    TH1F *h1_CT10_alpha_S_up   = new TH1F(Form("%s_h1_CT10_alpha_S_up", sampleName.c_str()),  "alpha_S up", nbins, min, max);
    TH1F *h1_CT10_alpha_S_down  = new TH1F(Form("%s_h1_CT10_alpha_S_down", sampleName.c_str()), "alpha_S down", nbins, min, max);
    AnalyseCT10as(f, sampleName, h1_CT10_alpha_S, h1_CT10_alpha_S_up, h1_CT10_alpha_S_down);
    h1_CT10_alpha_S_up->SetLineColor(kRed);
    h1_CT10_alpha_S_down->SetLineColor(kRed);
    h1_CT10_alpha_S_down->SetLineStyle(kDashed);

    TH1F *h1_CT10_env_up   = new TH1F(Form("%s_h1_CT10_env_up", sampleName.c_str()),  "CT10 up", nbins, min, max);    
    TH1F *h1_CT10_env_down   = new TH1F(Form("%s_h1_CT10_env_down", sampleName.c_str()),  "CT10 down", nbins, min, max);
    h1_CT10_env_up->SetLineColor(kRed);
    h1_CT10_env_down->SetLineColor(kRed);
    for (unsigned int bin = 0; bin < h1_CT10->GetNbinsX() + 1; ++bin) {
        float up = sqrt(pow(h1_CT10_up->GetBinContent(bin), 2) 
                            + pow(h1_CT10_alpha_S_up->GetBinContent(bin), 2));
        float down = -1*sqrt(pow(h1_CT10_down->GetBinContent(bin), 2) 
                            + pow(h1_CT10_alpha_S_down->GetBinContent(bin), 2));
        h1_CT10_env_up->SetBinContent(bin, up);
        h1_CT10_env_down->SetBinContent(bin, down);
    }

    TLegend *lg_cteq = new TLegend(0.2, 0.2, 0.65, 0.4);
    lg_cteq->AddEntry(h1_CT10_alpha_S_up, "CT10 alpha_S +68%", "l");
    lg_cteq->AddEntry(h1_CT10_up, "CT10", "l");
    lg_cteq->AddEntry(h1_CT10_alpha_S_down, "CT10 alpha_S -68%", "l");
    FormatLegend(lg_cteq);

    TCanvas *c_CT10as = new TCanvas();
    c_CT10as->cd();
    h1_CT10_alpha_S_up->Draw();
    h1_CT10_up->Draw("SAME");
    h1_CT10_down->Draw("SAME");
    h1_CT10_alpha_S_down->Draw("SAME");
    h1_CT10_alpha_S_up->GetYaxis()->SetRangeUser(-0.2, 0.2);
    lg_cteq->Draw();
    c_CT10as->SaveAs(Form("results/%s_CT10_alphaS.png", sampleName.c_str()));

    TLegend *lg_cteqenv = new TLegend(0.2, 0.2, 0.65, 0.4);
    lg_cteqenv->AddEntry(h1_CT10_env_up, "CT10 PDF+alpha_S", "l");
    lg_cteqenv->AddEntry(h1_CT10_up, "CT10 PDF", "l");
    FormatLegend(lg_cteqenv);

    TCanvas *c_CT10env = new TCanvas();
    c_CT10env->cd();
    h1_CT10_env_up->Draw();
    h1_CT10_env_down->Draw("SAME");
    h1_CT10_up->Draw("SAME");
    h1_CT10_down->Draw("SAME");
    h1_CT10_env_up->GetYaxis()->SetRangeUser(-0.2, 0.2);
    lg_cteqenv->Draw();
    c_CT10env->SaveAs(Form("results/%s_CT10_envelope.png", sampleName.c_str()));
    
    //
    // MSTW
    //

    TH1F *h1_MSTW   = new TH1F(Form("%s_h1_MSTW", sampleName.c_str()),  "centre", nbins, min, max);
    TH1F *h1_MSTW_up   = new TH1F(Form("%s_h1_MSTW_up", sampleName.c_str()),  "centre up", nbins, min, max);
    TH1F *h1_MSTW_down  = new TH1F(Form("%s_h1_MSTW_down", sampleName.c_str()), "centre down", nbins, min, max);
    AnalyseCTMSTW(f, sampleName, "MSTW2008nlo68cl", 41, h1_MSTW, h1_MSTW_up, h1_MSTW_down);

    TH1F *h1_MSTW_asup   = new TH1F(Form("%s_h1_MSTW_asup", sampleName.c_str()),  "asup", nbins, min, max);
    TH1F *h1_MSTW_asup_up   = new TH1F(Form("%s_h1_MSTW_asup_up", sampleName.c_str()),  "asup up", nbins, min, max);
    TH1F *h1_MSTW_asup_down  = new TH1F(Form("%s_h1_MSTW_asup_down", sampleName.c_str()), "asup down", nbins, min, max);
    AnalyseCTMSTW(f, sampleName, "MSTW2008nlo68cl_asmz+68cl", 41, h1_MSTW_asup, h1_MSTW_asup_up, h1_MSTW_asup_down);
    h1_MSTW_asup_up->SetLineColor(kRed);
    h1_MSTW_asup_down->SetLineColor(kRed);

    TH1F *h1_MSTW_asdown   = new TH1F(Form("%s_h1_MSTW_asdown", sampleName.c_str()),  "asdown", nbins, min, max);
    TH1F *h1_MSTW_asdown_up   = new TH1F(Form("%s_h1_MSTW_asdown_up", sampleName.c_str()),  "asdown up", nbins, min, max);
    TH1F *h1_MSTW_asdown_down  = new TH1F(Form("%s_h1_MSTW_asdown_down", sampleName.c_str()), "asdown down", nbins, min, max);
    AnalyseCTMSTW(f, sampleName, "MSTW2008nlo68cl_asmz-68cl", 41, h1_MSTW_asdown, h1_MSTW_asdown_up, h1_MSTW_asdown_down);
    h1_MSTW_asdown_up->SetLineColor(kRed);
    h1_MSTW_asdown_up->SetLineStyle(kDashed);
    h1_MSTW_asdown_down->SetLineColor(kRed);
    h1_MSTW_asdown_down->SetLineStyle(kDashed);

    TH1F *h1_MSTW_asup05   = new TH1F(Form("%s_h1_MSTW_asup05", sampleName.c_str()),  "asup05", nbins, min, max);
    TH1F *h1_MSTW_asup05_up   = new TH1F(Form("%s_h1_MSTW_asup05_up", sampleName.c_str()),  "asup05 up", nbins, min, max);
    TH1F *h1_MSTW_asup05_down  = new TH1F(Form("%s_h1_MSTW_asup05_down", sampleName.c_str()), "asup05 down", nbins, min, max);
    AnalyseCTMSTW(f, sampleName, "MSTW2008nlo68cl_asmz+68clhalf", 41, h1_MSTW_asup05, h1_MSTW_asup05_up, h1_MSTW_asup05_down);
    h1_MSTW_asup05_up->SetLineColor(kBlue);
    h1_MSTW_asup05_down->SetLineColor(kBlue);

    TH1F *h1_MSTW_asdown05   = new TH1F(Form("%s_h1_MSTW_asdown05", sampleName.c_str()),  "asdown05", nbins, min, max);
    TH1F *h1_MSTW_asdown05_up   = new TH1F(Form("%s_h1_MSTW_asdown05_up", sampleName.c_str()),  "asdown05 up", nbins, min, max);
    TH1F *h1_MSTW_asdown05_down  = new TH1F(Form("%s_h1_MSTW_asdown05_down", sampleName.c_str()), "asdown05 down", nbins, min, max);
    AnalyseCTMSTW(f, sampleName, "MSTW2008nlo68cl_asmz-68clhalf", 41, h1_MSTW_asdown05, h1_MSTW_asdown05_up, h1_MSTW_asdown05_down);
    h1_MSTW_asdown05_up->SetLineColor(kBlue);
    h1_MSTW_asdown05_up->SetLineStyle(kDashed);
    h1_MSTW_asdown05_down->SetLineColor(kBlue);
    h1_MSTW_asdown05_down->SetLineStyle(kDashed);

    // now combine to get the PDF + alpha_S uncertainty
    // according to the envelope described in AN2011/055
    TH1F *h1_MSTW_env_up   = new TH1F(Form("%s_h1_MSTW_env_up", sampleName.c_str()),  "MSTW up", nbins, min, max);    
    TH1F *h1_MSTW_env_down   = new TH1F(Form("%s_h1_MSTW_env_down", sampleName.c_str()),  "MSTW down", nbins, min, max);
    h1_MSTW_env_up->SetLineColor(kRed);
    h1_MSTW_env_down->SetLineColor(kRed);
    for (unsigned int bin = 0; bin < h1_MSTW->GetNbinsX() + 1; ++bin) {
        float delta_A_up = h1_MSTW_up->GetBinContent(bin);
        float delta_A_asup = h1_MSTW_asup_up->GetBinContent(bin);
        float delta_A_asup05 = h1_MSTW_asup05_up->GetBinContent(bin);
        float delta_A_down = h1_MSTW_down->GetBinContent(bin);
        float delta_A_asdown = h1_MSTW_asdown_down->GetBinContent(bin);
        float delta_A_asdown05 = h1_MSTW_asdown05_down->GetBinContent(bin);
        h1_MSTW_env_up->SetBinContent(bin, TMath::Max(TMath::Max(delta_A_up, delta_A_asup), delta_A_asup05));
        h1_MSTW_env_down->SetBinContent(bin, TMath::Min(TMath::Min(delta_A_down, delta_A_asdown), delta_A_asdown05));
    }

    TLegend *lg_mstw = new TLegend(0.2, 0.2, 0.65, 0.4);
    lg_mstw->AddEntry(h1_MSTW_asup_up, "MSTW2008nlo68cl_asmz+68cl", "l");
    lg_mstw->AddEntry(h1_MSTW_asup05_up, "MSTW2008nlo68cl_asmz+68clhalf", "l");
    lg_mstw->AddEntry(h1_MSTW, "MSTW2008nlo68cl_asmz", "l");
    lg_mstw->AddEntry(h1_MSTW_asdown05_up, "MSTW2008nlo68cl_asmz-68clhalf", "l");
    lg_mstw->AddEntry(h1_MSTW_asdown_up, "MSTW2008nlo68cl_asmz-68cl", "l");
    FormatLegend(lg_mstw);

    TCanvas *c_MSTWas = new TCanvas();
    c_MSTWas->cd();
    h1_MSTW_asup_up->Draw();
    h1_MSTW_up->Draw("SAME");
    h1_MSTW_down->Draw("SAME");
    h1_MSTW_asup_down->Draw("SAME");
    h1_MSTW_asdown_up->Draw("SAME");
    h1_MSTW_asdown_down->Draw("SAME");
    h1_MSTW_asup05_up->Draw("SAME");
    h1_MSTW_asup05_down->Draw("SAME");
    h1_MSTW_asdown05_up->Draw("SAME");
    h1_MSTW_asdown05_down->Draw("SAME");
    h1_MSTW_asup_up->GetYaxis()->SetRangeUser(-0.2, 0.2);
    lg_mstw->Draw();
    c_MSTWas->SaveAs(Form("results/%s_MSTW_alphaS.png", sampleName.c_str()));

    TLegend *lg_mstwenv = new TLegend(0.2, 0.2, 0.65, 0.4);
    lg_mstwenv->AddEntry(h1_MSTW_up, "MSTW PDF", "l");
    lg_mstwenv->AddEntry(h1_MSTW_env_up, "MSTW PDF+alpha_S", "l");
    FormatLegend(lg_mstwenv);

    TCanvas *c_MSTWenv = new TCanvas();
    c_MSTWenv->cd();
    h1_MSTW_env_up->Draw();
    h1_MSTW_env_down->Draw("SAME");
    h1_MSTW_up->Draw("SAME");
    h1_MSTW_down->Draw("SAME");
    h1_MSTW_env_up->GetYaxis()->SetRangeUser(-0.2, 0.2);
    lg_mstwenv->Draw();
    c_MSTWenv->SaveAs(Form("results/%s_MSTW_envelope.png", sampleName.c_str()));

    //
    // NNPDF
    // 

    TH1F *h1_NNPDF   = new TH1F(Form("%s_h1_NNPDF", sampleName.c_str()),  "centre", nbins, min, max);
    TH1F *h1_NNPDF_up   = new TH1F(Form("%s_h1_NNPDF_up", sampleName.c_str()),  "centre up", nbins, min, max);
    TH1F *h1_NNPDF_down  = new TH1F(Form("%s_h1_NNPDF_down", sampleName.c_str()), "centre down", nbins, min, max);
    AnalyseNNPDF(f, sampleName, h1_NNPDF, h1_NNPDF_up, h1_NNPDF_down);

    TLegend *lg_nnpdf = new TLegend(0.2, 0.2, 0.65, 0.4);
    lg_nnpdf->AddEntry(h1_NNPDF_up, "NNPDF 2.0 PDF+alpha_S", "l");
    FormatLegend(lg_nnpdf);

    TCanvas *c_NNPDF = new TCanvas();
    c_NNPDF->cd();
    h1_NNPDF_up->Draw();
    h1_NNPDF_down->Draw("SAME");
    h1_NNPDF_up->GetYaxis()->SetRangeUser(-0.2, 0.2);
    lg_nnpdf->Draw("SAME");
    c_NNPDF->SaveAs(Form("results/%s_NNPDF_envelope.png", sampleName.c_str()));

    //
    //
    //

    // ! GRAND TOTAL !

    //
    // compare overall shapes
    //

    TH1F *h1_CTEQ6M = (TH1F*)f->Get(Form("%s_cteq6mE_0", sampleName.c_str()))->Clone("h1_CTEQ6M");
    TH1F *h1_NNPDF_diff = (TH1F*)h1_NNPDF->Clone("h1_NNPDF_diff");
    TH1F *h1_CT10_diff = (TH1F*)h1_CT10->Clone("h1_CT10_diff");
    TH1F *h1_MSTW_diff = (TH1F*)h1_MSTW->Clone("h1_MSTW_diff");
    h1_NNPDF_diff->SetLineColor(kBlack);
    h1_CT10_diff->SetLineColor(kRed);
    h1_MSTW_diff->SetLineColor(kBlue);

    TH1F *h1_NNPDF_diff_env_up = (TH1F*)h1_NNPDF->Clone("h1_NNPDF_diff_env");
    TH1F *h1_CT10_diff_env_up = (TH1F*)h1_CT10->Clone("h1_CT10_diff_env");
    TH1F *h1_MSTW_diff_env_up = (TH1F*)h1_MSTW->Clone("h1_MSTW_diff_env");
    h1_NNPDF_diff_env_up->SetLineColor(kBlack);
    h1_CT10_diff_env_up->SetLineColor(kRed);
    h1_MSTW_diff_env_up->SetLineColor(kBlue);
    TH1F *h1_NNPDF_diff_env_down = (TH1F*)h1_NNPDF->Clone("h1_NNPDF_diff_env");
    TH1F *h1_CT10_diff_env_down = (TH1F*)h1_CT10->Clone("h1_CT10_diff_env");
    TH1F *h1_MSTW_diff_env_down = (TH1F*)h1_MSTW->Clone("h1_MSTW_diff_env");
    h1_NNPDF_diff_env_down->SetLineColor(kBlack);
    h1_CT10_diff_env_down->SetLineColor(kRed);
    h1_MSTW_diff_env_down->SetLineColor(kBlue);
    h1_NNPDF_diff_env_down->SetLineStyle(kDashed);
    h1_CT10_diff_env_down->SetLineStyle(kDashed);
    h1_MSTW_diff_env_down->SetLineStyle(kDashed);

    TH1F *h1_env_up   = new TH1F(Form("%s_h1_env_up", sampleName.c_str()),  "up", nbins, min, max);
    TH1F *h1_env_down   = new TH1F(Form("%s_h1_env_down", sampleName.c_str()),  "down", nbins, min, max);
    h1_env_up->SetLineWidth(2);
    h1_env_down->SetLineWidth(2);

    for (unsigned int bin = 0; bin < h1_CTEQ6M->GetNbinsX() + 1; ++bin) {
        float A0 = h1_CTEQ6M->GetBinContent(bin);
        h1_NNPDF_diff->SetBinContent(bin, h1_NNPDF_diff->GetBinContent(bin) - A0);
        h1_CT10_diff->SetBinContent(bin, h1_CT10_diff->GetBinContent(bin) - A0);
        h1_MSTW_diff->SetBinContent(bin, h1_MSTW_diff->GetBinContent(bin) - A0);

        h1_NNPDF_diff_env_up->SetBinContent(bin, h1_NNPDF_diff->GetBinContent(bin) + h1_NNPDF_up->GetBinContent(bin));
        h1_CT10_diff_env_up->SetBinContent(bin, h1_CT10_diff->GetBinContent(bin) + h1_CT10_up->GetBinContent(bin));
        h1_MSTW_diff_env_up->SetBinContent(bin, h1_MSTW_diff->GetBinContent(bin) + h1_MSTW_up->GetBinContent(bin));

        h1_NNPDF_diff_env_down->SetBinContent(bin, h1_NNPDF_diff->GetBinContent(bin) + h1_NNPDF_down->GetBinContent(bin));
        h1_CT10_diff_env_down->SetBinContent(bin, h1_CT10_diff->GetBinContent(bin) + h1_CT10_down->GetBinContent(bin));
        h1_MSTW_diff_env_down->SetBinContent(bin, h1_MSTW_diff->GetBinContent(bin) + h1_MSTW_down->GetBinContent(bin));

        h1_env_up->SetBinContent(bin, TMath::Max(TMath::Max(h1_NNPDF_diff_env_up->GetBinContent(bin), h1_CT10_diff_env_up->GetBinContent(bin)), h1_MSTW_diff_env_up->GetBinContent(bin)));
        h1_env_down->SetBinContent(bin, TMath::Min(TMath::Min(h1_NNPDF_diff_env_down->GetBinContent(bin), h1_CT10_diff_env_down->GetBinContent(bin)), h1_MSTW_diff_env_down->GetBinContent(bin)));

    }

    TLegend *lg_env = new TLegend(0.2, 0.2, 0.65, 0.4);
    lg_env->AddEntry(h1_env_up, "Overall Uncertainty", "l");
    lg_env->AddEntry(h1_NNPDF_diff_env_up, "NNPDF", "l");
    lg_env->AddEntry(h1_CT10_diff_env_up, "CT10", "l");
    lg_env->AddEntry(h1_MSTW_diff_env_up, "MSTW", "l");
    FormatLegend(lg_env);

    TCanvas *c_central = new TCanvas();
    c_central->cd();
    h1_env_up->Draw();
    h1_env_down->Draw("SAME");
    h1_NNPDF_diff_env_up->Draw("SAME");
    h1_CT10_diff_env_up->Draw("SAME");
    h1_MSTW_diff_env_up->Draw("SAME");
    h1_NNPDF_diff_env_down->Draw("SAME");
    h1_CT10_diff_env_down->Draw("SAME");
    h1_MSTW_diff_env_down->Draw("SAME");
    h1_env_up->GetYaxis()->SetRangeUser(-0.2, 0.2);
    lg_env->Draw();
    c_central->SaveAs(Form("results/%s_envelope.png", sampleName.c_str()));

    f->Close();
}

