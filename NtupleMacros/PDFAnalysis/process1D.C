#include "TLine.h"
#include "TFile.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"

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

    for (Int_t binx = 1; binx <= h1->GetNbinsX(); ++binx)
    {

        const unsigned int max_set = 3;
        Double_t X0 = histArr[5]->GetBinContent(binx);
        //Double_t binCenter = histArr[5]->GetBinCenter(binx);
        if (X0 == 0) continue;
        Double_t Xi_up = histArr[10-max_set]->GetBinContent(binx);
        Double_t Xi_down = histArr[max_set]->GetBinContent(binx);
        Double_t plus_max = sqrt(pow(TMath::Max(TMath::Max(Xi_up - X0, Xi_down - X0), 0.0), 2))/1.645;
        Double_t minus_max = sqrt(pow(TMath::Max(TMath::Max(X0 - Xi_up, X0 - Xi_down), 0.0), 2))/1.645;
        h1->SetBinContent(binx, X0);
        h1_down->SetBinContent(binx, -1*minus_max);
        h1_up->SetBinContent(binx, plus_max);
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

    for (Int_t binx = 1; binx <= h1->GetNbinsX(); ++binx)
    {

        Double_t X0 = histArr[0]->GetBinContent(binx);
        //Double_t binCenter = histArr[0]->GetBinCenter(binx);
        Double_t plus_max = 0.0;
        Double_t minus_max = 0.0;

        if (X0 == 0) continue;

        for (unsigned int subset = 0; subset < ((nsets - 1)/2); ++subset) 
        {
            Double_t Xi_up = histArr[(subset*2) + 1]->GetBinContent(binx);
            Double_t Xi_down = histArr[(subset*2) + 2]->GetBinContent(binx);
            plus_max += pow(TMath::Max(TMath::Max(Xi_up - X0, Xi_down - X0), 0.0), 2);
            minus_max += pow(TMath::Max(TMath::Max(X0 - Xi_up, X0 - Xi_down), 0.0), 2);
        }
        plus_max = sqrt(plus_max);
        minus_max = sqrt(minus_max);
        h1->SetBinContent(binx, X0);
        h1_down->SetBinContent(binx, -1*minus_max);
        h1_up->SetBinContent(binx, plus_max);
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
    for (Int_t binx= 1; binx <= h1->GetNbinsX(); ++binx) {
        h1->SetBinContent(binx, h1_temp->GetBinContent(binx));
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

    for (Int_t binx = 1; binx <= h1->GetNbinsX(); ++binx)
    {

        float X0 = h1->GetBinContent(binx);
        if (X0 == 0) continue;

        float delta_up = 0.0;
        float delta_down = 0.0;
        for (unsigned int i = 0; i < histArr.size(); ++i) {
            if (histArr[i]->GetBinContent(binx) > h1->GetBinContent(binx)) {
                delta_up += pow(histArr[i]->GetBinContent(binx) - h1->GetBinContent(binx), 2);
                ++nsets_up;
            }
            else {
                delta_down += pow(h1->GetBinContent(binx) - histArr[i]->GetBinContent(binx), 2);
                ++nsets_down;
            }
        }

        delta_up = sqrt( (1.0/(float(nsets_up) - 1.0)) * delta_up);
        delta_down = sqrt( (1.0/(float(nsets_down) - 1.0)) * delta_down);
        h1_up->SetBinContent(binx, delta_up);
        h1_down->SetBinContent(binx, -1*delta_down);
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

// sampleName e.g. qqWW_DF_0j...
// genPdf e.g. cteq6ll_0 or CT10_5...
void process1D(std::string fileName, std::string sampleName, std::string genPdf)
{

    gSystem->mkdir("results");
    gROOT->ProcessLine(".L ~/tdrStyle.C");
    gROOT->ProcessLine("setTDRStyle()");
    gStyle->SetPalette(1);
    gStyle->SetPaintTextFormat("4.3f");

    TFile *f = new TFile(Form("%s.root", fileName.c_str()), "READ");
    gROOT->cd();

    TH1F *h1_temp = (TH1F*)f->Get(Form("%s_%s", sampleName.c_str(), genPdf.c_str()))->Clone("temp");
    Int_t nbinsx    = h1_temp->GetXaxis()->GetNbins();
    Float_t xmin    = h1_temp->GetXaxis()->GetXmin();
    Float_t xmax    = h1_temp->GetXaxis()->GetXmax();
    delete h1_temp;

    //
    // CTEQ
    //

    // central
    TH1F *h1_CT10   = new TH1F(Form("%s_h1_CT10", sampleName.c_str()),  "centre", nbinsx, xmin, xmax);
    TH1F *h1_CT10_up   = new TH1F(Form("%s_h1_CT10_up", sampleName.c_str()),  "centre up", nbinsx, xmin, xmax);
    TH1F *h1_CT10_down  = new TH1F(Form("%s_h1_CT10_down", sampleName.c_str()), "centre down", nbinsx, xmin, xmax);
    AnalyseCTMSTW(f, sampleName, "CT10", 53, h1_CT10, h1_CT10_up, h1_CT10_down);

    // alpha_s
    TH1F *h1_CT10_alpha_S   = new TH1F(Form("%s_h1_CT10_alpha_S", sampleName.c_str()),  "alpha_S", nbinsx, xmin, xmax);
    TH1F *h1_CT10_alpha_S_up   = new TH1F(Form("%s_h1_CT10_alpha_S_up", sampleName.c_str()),  "alpha_S up", nbinsx, xmin, xmax);
    TH1F *h1_CT10_alpha_S_down  = new TH1F(Form("%s_h1_CT10_alpha_S_down", sampleName.c_str()), "alpha_S down", nbinsx, xmin, xmax);

    AnalyseCT10as(f, sampleName, h1_CT10_alpha_S, h1_CT10_alpha_S_up, h1_CT10_alpha_S_down);
    h1_CT10_alpha_S_up->SetLineColor(kRed);
    h1_CT10_alpha_S_down->SetLineColor(kRed);
    h1_CT10_alpha_S_down->SetLineStyle(kDashed);


    TH1F *h1_CT10_env_up   = new TH1F(Form("%s_h1_CT10_env_up", sampleName.c_str()),  "CT10 up", nbinsx, xmin, xmax);    
    TH1F *h1_CT10_env_down   = new TH1F(Form("%s_h1_CT10_env_down", sampleName.c_str()),  "CT10 down", nbinsx, xmin, xmax);
    h1_CT10_env_up->SetLineColor(kRed);
    h1_CT10_env_down->SetLineColor(kRed);
    for (int binx = 0; binx < h1_CT10->GetNbinsX() + 1; ++binx) {
        float up = sqrt(pow(h1_CT10_up->GetBinContent(binx), 2) 
                + pow(h1_CT10_alpha_S_up->GetBinContent(binx), 2));
        float down = -1*sqrt(pow(h1_CT10_down->GetBinContent(binx), 2) 
                + pow(h1_CT10_alpha_S_down->GetBinContent(binx), 2));
        h1_CT10_env_up->SetBinContent(binx, up);
        h1_CT10_env_down->SetBinContent(binx, down);
    }

    //
    // MSTW
    //

    TH1F *h1_MSTW   = new TH1F(Form("%s_h1_MSTW", sampleName.c_str()),  "centre", nbinsx, xmin, xmax);
    TH1F *h1_MSTW_up   = new TH1F(Form("%s_h1_MSTW_up", sampleName.c_str()),  "centre up", nbinsx, xmin, xmax);
    TH1F *h1_MSTW_down  = new TH1F(Form("%s_h1_MSTW_down", sampleName.c_str()), "centre down", nbinsx, xmin, xmax);
    AnalyseCTMSTW(f, sampleName, "MSTW2008nlo68cl", 41, h1_MSTW, h1_MSTW_up, h1_MSTW_down);

    TH1F *h1_MSTW_asup   = new TH1F(Form("%s_h1_MSTW_asup", sampleName.c_str()),  "asup", nbinsx, xmin, xmax);
    TH1F *h1_MSTW_asup_up   = new TH1F(Form("%s_h1_MSTW_asup_up", sampleName.c_str()),  "asup up", nbinsx, xmin, xmax);
    TH1F *h1_MSTW_asup_down  = new TH1F(Form("%s_h1_MSTW_asup_down", sampleName.c_str()), "asup down", nbinsx, xmin, xmax);
    AnalyseCTMSTW(f, sampleName, "MSTW2008nlo68cl_asmz+68cl", 41, h1_MSTW_asup, h1_MSTW_asup_up, h1_MSTW_asup_down);
    h1_MSTW_asup_up->SetLineColor(kRed);
    h1_MSTW_asup_down->SetLineColor(kRed);

    TH1F *h1_MSTW_asdown   = new TH1F(Form("%s_h1_MSTW_asdown", sampleName.c_str()),  "asdown", nbinsx, xmin, xmax);
    TH1F *h1_MSTW_asdown_up   = new TH1F(Form("%s_h1_MSTW_asdown_up", sampleName.c_str()),  "asdown up", nbinsx, xmin, xmax);
    TH1F *h1_MSTW_asdown_down  = new TH1F(Form("%s_h1_MSTW_asdown_down", sampleName.c_str()), "asdown down", nbinsx, xmin, xmax);
    AnalyseCTMSTW(f, sampleName, "MSTW2008nlo68cl_asmz-68cl", 41, h1_MSTW_asdown, h1_MSTW_asdown_up, h1_MSTW_asdown_down);
    h1_MSTW_asdown_up->SetLineColor(kRed);
    h1_MSTW_asdown_up->SetLineStyle(kDashed);
    h1_MSTW_asdown_down->SetLineColor(kRed);
    h1_MSTW_asdown_down->SetLineStyle(kDashed);

    TH1F *h1_MSTW_asup05   = new TH1F(Form("%s_h1_MSTW_asup05", sampleName.c_str()),  "asup05", nbinsx, xmin, xmax);
    TH1F *h1_MSTW_asup05_up   = new TH1F(Form("%s_h1_MSTW_asup05_up", sampleName.c_str()),  "asup05 up", nbinsx, xmin, xmax);
    TH1F *h1_MSTW_asup05_down  = new TH1F(Form("%s_h1_MSTW_asup05_down", sampleName.c_str()), "asup05 down", nbinsx, xmin, xmax);
    AnalyseCTMSTW(f, sampleName, "MSTW2008nlo68cl_asmz+68clhalf", 41, h1_MSTW_asup05, h1_MSTW_asup05_up, h1_MSTW_asup05_down);
    h1_MSTW_asup05_up->SetLineColor(kBlue);
    h1_MSTW_asup05_down->SetLineColor(kBlue);

    TH1F *h1_MSTW_asdown05   = new TH1F(Form("%s_h1_MSTW_asdown05", sampleName.c_str()),  "asdown05", nbinsx, xmin, xmax);
    TH1F *h1_MSTW_asdown05_up   = new TH1F(Form("%s_h1_MSTW_asdown05_up", sampleName.c_str()),  "asdown05 up", nbinsx, xmin, xmax);
    TH1F *h1_MSTW_asdown05_down  = new TH1F(Form("%s_h1_MSTW_asdown05_down", sampleName.c_str()), "asdown05 down", nbinsx, xmin, xmax);
    AnalyseCTMSTW(f, sampleName, "MSTW2008nlo68cl_asmz-68clhalf", 41, h1_MSTW_asdown05, h1_MSTW_asdown05_up, h1_MSTW_asdown05_down);
    h1_MSTW_asdown05_up->SetLineColor(kBlue);
    h1_MSTW_asdown05_up->SetLineStyle(kDashed);
    h1_MSTW_asdown05_down->SetLineColor(kBlue);
    h1_MSTW_asdown05_down->SetLineStyle(kDashed);

    // now combine to get the PDF + alpha_S uncertainty
    // according to the envelope described in AN2011/055
    TH1F *h1_MSTW_env_up   = new TH1F(Form("%s_h1_MSTW_env_up", sampleName.c_str()),  "MSTW up", nbinsx, xmin, xmax);    
    TH1F *h1_MSTW_env_down   = new TH1F(Form("%s_h1_MSTW_env_down", sampleName.c_str()),  "MSTW down", nbinsx, xmin, xmax);
    h1_MSTW_env_up->SetLineColor(kRed);
    h1_MSTW_env_down->SetLineColor(kRed);
    for (int binx = 0; binx < h1_MSTW->GetNbinsX() + 1; ++binx) {
        float delta_A_up = h1_MSTW_up->GetBinContent(binx);
        float delta_A_asup = h1_MSTW_asup_up->GetBinContent(binx);
        float delta_A_asup05 = h1_MSTW_asup05_up->GetBinContent(binx);
        float delta_A_down = h1_MSTW_down->GetBinContent(binx);
        float delta_A_asdown = h1_MSTW_asdown_down->GetBinContent(binx);
        float delta_A_asdown05 = h1_MSTW_asdown05_down->GetBinContent(binx);
        h1_MSTW_env_up->SetBinContent(binx, TMath::Max(TMath::Max(delta_A_up, delta_A_asup), delta_A_asup05));
        h1_MSTW_env_down->SetBinContent(binx, TMath::Min(TMath::Min(delta_A_down, delta_A_asdown), delta_A_asdown05));
    }

    //
    // NNPDF
    // 

    TH1F *h1_NNPDF   = new TH1F(Form("%s_h1_NNPDF", sampleName.c_str()),  "centre", nbinsx, xmin, xmax);
    TH1F *h1_NNPDF_up   = new TH1F(Form("%s_h1_NNPDF_up", sampleName.c_str()),  "centre up", nbinsx, xmin, xmax);
    TH1F *h1_NNPDF_down  = new TH1F(Form("%s_h1_NNPDF_down", sampleName.c_str()), "centre down", nbinsx, xmin, xmax);
    AnalyseNNPDF(f, sampleName, h1_NNPDF, h1_NNPDF_up, h1_NNPDF_down);

    //
    // ! GRAND TOTAL !
    //

    TH1F *h1_GenPDF = (TH1F*)f->Get(Form("%s_%s", sampleName.c_str(), genPdf.c_str()))->Clone("GenPDF");
    TH1F *h1_env_up   = new TH1F(Form("%s_h1_env_up", sampleName.c_str()),  "up", nbinsx, xmin, xmax);
    TH1F *h1_env_down   = new TH1F(Form("%s_h1_env_down", sampleName.c_str()),  "down", nbinsx, xmin, xmax);
    h1_env_up->SetFillColor(kGray);
    h1_env_up->SetLineColor(kWhite);
    h1_env_down->SetFillColor(10);
    h1_env_down->SetLineColor(kWhite);

    for (int binx = 0; binx < h1_GenPDF->GetNbinsX() + 1; ++binx) {

        // get up and down variations including both
        // - intrinsic pdf variation from central value
        // - variation of central value from generated set
        Float_t up_NNPDF    = h1_NNPDF->GetBinContent(binx) + h1_NNPDF_up->GetBinContent(binx);
        Float_t up_CT10     = h1_CT10->GetBinContent(binx) + h1_CT10_env_up->GetBinContent(binx);
        Float_t up_MSTW     = h1_MSTW->GetBinContent(binx) + h1_MSTW_env_up->GetBinContent(binx);
        Float_t down_NNPDF  = h1_NNPDF->GetBinContent(binx) + h1_NNPDF_down->GetBinContent(binx);
        Float_t down_CT10   = h1_CT10->GetBinContent(binx) + h1_CT10_env_down->GetBinContent(binx);
        Float_t down_MSTW   = h1_MSTW->GetBinContent(binx) + h1_MSTW_env_down->GetBinContent(binx);

        // get the envelope of envelopes...
        h1_env_up->SetBinContent(binx,    TMath::Max(TMath::Max(up_NNPDF,     up_CT10),   up_MSTW));
        h1_env_down->SetBinContent(binx,  TMath::Min(TMath::Min(down_NNPDF,   down_CT10), down_MSTW));
    }

    //
    // ... scale variations to central shape
    //

    // get the mid point of the up and down envelopes
    // as the central alternate pdf shape
    TH1F* h1_env_midpoint = (TH1F*)h1_env_up->Clone("h1_env_midpoint");
    h1_env_midpoint->SetTitle("midpoint");
    h1_env_midpoint->Add(h1_env_down, -1.0);
    h1_env_midpoint->Scale(0.5);
    h1_env_midpoint->Add(h1_env_down, +1.0);

    // now scale the central alternate pdf shape
    // to the generated shape in the sample
    Float_t scaleFactor = h1_GenPDF->Integral(0, h1_GenPDF->GetXaxis()->GetNbins())
        / h1_env_midpoint->Integral(0, h1_env_midpoint->GetXaxis()->GetNbins());
    h1_env_midpoint->Scale(scaleFactor);
    h1_env_up->Scale(scaleFactor);
    h1_env_down->Scale(scaleFactor);
    h1_CT10->Scale(scaleFactor);
    h1_CT10_env_up->Scale(scaleFactor);
    h1_CT10_env_down->Scale(scaleFactor);
    h1_MSTW->Scale(scaleFactor);
    h1_MSTW_env_up->Scale(scaleFactor);
    h1_MSTW_env_down->Scale(scaleFactor);
    h1_NNPDF->Scale(scaleFactor);
    h1_NNPDF_up->Scale(scaleFactor);
    h1_NNPDF_down->Scale(scaleFactor);

    // compute the ratios
    TH1F *h1_ratio = (TH1F*)h1_env_midpoint->Clone("h1_ratio");
    h1_ratio->Divide(h1_GenPDF);
    TH1F *h1_ratioUp = (TH1F*)h1_env_up->Clone("h1_ratioUp");
    h1_ratioUp->Divide(h1_GenPDF);
    TH1F *h1_ratioDown = (TH1F*)h1_env_down->Clone("h1_ratioDown");
    h1_ratioDown->Divide(h1_GenPDF);

    TH1F *h1_ratioCT10 = (TH1F*)h1_CT10->Clone("h1_ratioCT10");
    h1_ratioCT10->Divide(h1_GenPDF);
    TH1F *h1_ratioCT10Up = (TH1F*)h1_CT10_env_up->Clone("h1_ratioCT10Up");
    h1_ratioCT10Up->Divide(h1_GenPDF);
    TH1F *h1_ratioCT10Down = (TH1F*)h1_CT10_env_down->Clone("h1_ratioCT10Down");
    h1_ratioCT10Down->Divide(h1_GenPDF);

    TH1F *h1_ratioMSTW = (TH1F*)h1_MSTW->Clone("h1_ratioMSTW");
    h1_ratioMSTW->Divide(h1_GenPDF);
    TH1F *h1_ratioMSTWUp = (TH1F*)h1_MSTW_env_up->Clone("h1_ratioMSTWUp");
    h1_ratioMSTWUp->Divide(h1_GenPDF);
    TH1F *h1_ratioMSTWDown = (TH1F*)h1_MSTW_env_down->Clone("h1_ratioMSTWDown");
    h1_ratioMSTWDown->Divide(h1_GenPDF);

    TH1F *h1_ratioNNPDF = (TH1F*)h1_NNPDF->Clone("h1_ratioNNPDF");
    h1_ratioNNPDF->Divide(h1_GenPDF);
    TH1F *h1_ratioNNPDFUp = (TH1F*)h1_NNPDF_up->Clone("h1_ratioNNPDFUp");
    h1_ratioNNPDFUp->Divide(h1_GenPDF);
    TH1F *h1_ratioNNPDFDown = (TH1F*)h1_NNPDF_down->Clone("h1_ratioNNPDFDown");
    h1_ratioNNPDFDown->Divide(h1_GenPDF);

    //
    // save it 
    //

    TFile f_out(Form("results/PDFUncertainty_%s_%s.root",  fileName.c_str(), sampleName.c_str()), "RECREATE");
    f_out.cd();
    h1_ratio->Write();
    h1_ratioUp->Write();
    h1_ratioDown->Write();
    f_out.Close();

    //
    //
    //

    // sigh, this is getting tiring...
    // I need some tea
    TCanvas *c_darjeeling = new TCanvas("c_darjeeling", "c_darjeeling", 1000, 600);
    c_darjeeling->SetRightMargin(0.15);
    c_darjeeling->cd();

    TLine *line = new TLine(h1_ratio->GetXaxis()->GetXmin(), 1.0, h1_ratio->GetXaxis()->GetXmax(), 1.0);
    line->SetLineStyle(kDashed);
    line->SetLineColor(kBlack);

    // individual CT10 
    h1_ratioCT10Up->Draw("HIST");
    h1_ratioCT10Up->SetLineStyle(kDashed);
    h1_ratioCT10Up->SetLineColor(kRed);
    h1_ratioCT10Up->GetYaxis()->SetRangeUser(-0.2, 1.2);
    h1_ratioCT10Up->GetYaxis()->SetTitle("Ratio {Alternate/Default}");
    h1_ratioCT10Down->Draw("HIST SAME");
    h1_ratioCT10Down->SetLineStyle(kDashed);
    h1_ratioCT10Down->SetLineColor(kRed);
    h1_ratioCT10->Draw("HIST SAME");
    h1_ratioCT10->SetLineWidth(2);
    h1_ratioCT10->SetLineColor(kRed);
    h1_ratioCT10->SetFillStyle(0);
    line->Draw();
    c_darjeeling->RedrawAxis();
    c_darjeeling->SaveAs(Form("results/%s_%s_alternateShapeRatioCT10.png", fileName.c_str(), sampleName.c_str()));

    h1_ratioCT10Down->Draw("HIST");
    c_darjeeling->SaveAs("results/test.png");

    // individual MSTW
    h1_ratioMSTWUp->Draw("HIST");
    h1_ratioMSTWUp->SetLineStyle(kDashed);
    h1_ratioMSTWUp->SetLineColor(kBlue);
    h1_ratioMSTWUp->GetYaxis()->SetRangeUser(-0.2, 1.2);
    h1_ratioMSTWUp->GetYaxis()->SetTitle("Ratio {Alternate/Default}");
    h1_ratioMSTWDown->Draw("HIST SAME");
    h1_ratioMSTWDown->SetLineStyle(kDashed);
    h1_ratioMSTWDown->SetLineColor(kBlue);
    h1_ratioMSTW->Draw("HIST SAME");
    h1_ratioMSTW->SetLineWidth(2);
    h1_ratioMSTW->SetLineColor(kBlue);
    h1_ratioMSTW->SetFillStyle(0);
    line->Draw();
    c_darjeeling->RedrawAxis();
    c_darjeeling->SaveAs(Form("results/%s_%s_alternateShapeRatioMSTW.png", fileName.c_str(), sampleName.c_str()));

    // individual NNPDF
    h1_ratioNNPDFUp->Draw("HIST");
    h1_ratioNNPDFUp->SetLineStyle(kDashed);
    h1_ratioNNPDFUp->SetLineColor(kBlack);
    h1_ratioNNPDFUp->GetYaxis()->SetRangeUser(-0.2, 1.2);
    h1_ratioNNPDFUp->GetYaxis()->SetTitle("Ratio {Alternate/Default}");
    h1_ratioNNPDFDown->Draw("HIST SAME");
    h1_ratioNNPDFDown->SetLineStyle(kDashed);
    h1_ratioNNPDFDown->SetLineColor(kBlack);
    h1_ratioNNPDF->Draw("HIST SAME");
    h1_ratioNNPDF->SetLineWidth(2);
    h1_ratioNNPDF->SetLineColor(kBlack);
    h1_ratioNNPDF->SetFillStyle(0);
    line->Draw();
    c_darjeeling->RedrawAxis();
    c_darjeeling->SaveAs(Form("results/%s_%s_alternateShapeRatioNNPDF.png", fileName.c_str(), sampleName.c_str()));

    // and its ratio
    Int_t palette[7];
    palette[0] = kWhite;
    for (unsigned int i=0;i<7;++i){
        palette[i] = 18-i;
    }
    gStyle->SetPalette(7,palette);
    gStyle->SetPaintTextFormat("4.2f");

    // ratio of central and alternate...
    c_darjeeling->SetLogy(0);
    TLegend *l1 = new TLegend(0.2, 0.6, 0.6, 0.85);
    l1->SetLineColor(kWhite);
    l1->SetShadowColor(kWhite);
    l1->SetFillColor(kWhite);
    l1->AddEntry(h1_ratio, "PDF variation midpoint", "l");
    l1->AddEntry(h1_ratioUp, "PDF variation +/- 1#sigma", "f");

    h1_ratioUp->Draw("HIST");
    h1_ratioUp->GetYaxis()->SetRangeUser(0.7, 1.5);
    h1_ratioUp->GetYaxis()->SetTitle("Ratio {Alternate/Default}");
    h1_ratioDown->Draw("HIST SAME");
    h1_ratio->Draw("HIST SAME");
    h1_ratio->SetLineWidth(2);
    h1_ratio->SetLineColor(kBlack);
    h1_ratio->SetFillStyle(0);
    line->Draw();
    l1->Draw();
    c_darjeeling->RedrawAxis();
    c_darjeeling->SaveAs(Form("results/%s_%s_alternateShapeRatio.png", fileName.c_str(), sampleName.c_str()));


    l1->AddEntry(h1_ratioCT10, "CT10 variation midpoint", "l");
    l1->AddEntry(h1_ratioMSTW, "MSTW variation midpoint", "l");
    l1->AddEntry(h1_ratioNNPDF, "NNPDF variation midpoint", "l");
    h1_ratioUp->Draw("HIST");
    h1_ratioDown->Draw("HIST SAME");
    h1_ratio->Draw("HIST SAME");
    h1_ratioCT10->Draw("HIST SAME");
    h1_ratioCT10->SetLineWidth(2);
    h1_ratioCT10->SetLineColor(kRed);
    h1_ratioCT10->SetLineStyle(kDashed);
    h1_ratioCT10->SetFillStyle(0);
    h1_ratioMSTW->Draw("HIST SAME");
    h1_ratioMSTW->SetLineWidth(2);
    h1_ratioMSTW->SetLineColor(kBlue);
    h1_ratioMSTW->SetLineStyle(kDashed+1);
    h1_ratioMSTW->SetFillStyle(0);
    h1_ratioNNPDF->Draw("HIST SAME");
    h1_ratioNNPDF->SetLineWidth(2);
    h1_ratioNNPDF->SetLineColor(kBlack);
    h1_ratioNNPDF->SetLineStyle(kDashed+2);
    h1_ratioNNPDF->SetFillStyle(0);
    line->Draw();
    l1->Draw();
    c_darjeeling->SaveAs(Form("results/%s_%s_alternateShapeRatioDetailed.png", fileName.c_str(), sampleName.c_str()));

    f->Close();
}

