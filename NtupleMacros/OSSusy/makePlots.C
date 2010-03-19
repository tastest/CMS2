#include "TROOT.h"
#include "TH1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <iostream>
#include <iomanip>
#include <math.h>
#include <string>
#include "TLegendEntry.h"
//#include "histtools.h"

namespace hist {
  void add(const char* outHistName, const char* patORpfx0);
  TLegend* legend(TCanvas* canvas, Option_t* option = "lpf", Bool_t addColor = kFALSE, Int_t token = -1, Float_t xmin = 0.75, Float_t ymin = 0.75, Float_t xmax = 0.99, Float_t ymax = 0.99);
}

////void loadHist(const char* filename, const char* pfx = 0, const char* pat = "*", Bool_t doAdd = kFALSE);
void loadHist(const char* filename, const char* directory = 0, const char* pfx = 0, const char* pat = "*", Bool_t doAdd = kFALSE);

void makePlots(string filename, string met, string jet) {

  gROOT->SetStyle("Plain");
  gROOT->ProcessLine("gStyle->SetOptStat(0)");

  const char samples[19][6] = {"ttdil", "DYee", "wjets", "tW", "ww", "wz", "zz", "LM0", "LM1", "LM2", "LM3", "LM4", "LM5", "LM6", "LM7", "LM8", "LM9", "LM10", "LM11"};

  const unsigned int nSMSamples = 7;
  const unsigned int nLMSamples = 12;
  
  //TFile* file = new TFile(filename.c_str());

  loadHist(filename.c_str(), 0, 0, "*dilPt*", kTRUE);
  loadHist(filename.c_str(), 0, 0, "*tcmet*", kTRUE);
  loadHist(filename.c_str(), 0, 0, "*metmuonjes*", kTRUE);

  std::cout << "adding histograms..." << std::endl;

  hist::add("sm_dilPt_2j_all", "^[^L].*_hdilPt_2j_all$");
  //hist::add("sm_dilPt_2j_all", "^LM11_hdilPt_2j_all$");
  hist::add("sm_tcmet_2j_all", "^[^L].*_htcmet_2j_all$");
  hist::add("sm_metmuonjes_2j_all", "^[^L].*_hmetmuonjes_2j_all$");

  TH1F* sm_dilPt_jpts = (TH1F*)gDirectory->Get("sm_dilPt_2j_all");
  
  if(sm_dilPt_jpts == 0) return;

  TH1F* sm_tcmet_jpts = (TH1F*)gDirectory->Get("sm_tcmet_2j_all");

  if(sm_tcmet_jpts == 0) return;

  TH1F* sm_metmuonjes_calo = (TH1F*)gDirectory->Get("sm_metmuonjes_2j_all");

  if(sm_metmuonjes_calo == 0) return;

  std::cout << "making predictions..." << std::endl;

  int bin50  = sm_dilPt_jpts->FindBin(50);
  int bin100 = sm_dilPt_jpts->FindBin(100);
  int bin175 = sm_dilPt_jpts->FindBin(175);
  double scale = sm_dilPt_jpts->Integral(bin50, 101) / sm_dilPt_jpts->Integral(0, 101);
  sm_dilPt_jpts->Scale( 1 / scale );
  double dilpt_prediction100 = sm_dilPt_jpts->Integral(bin100, 101);
  double dilpt_prediction175 = sm_dilPt_jpts->Integral(bin175, 101);

  std::cout << "the scale factor is: " << scale << std::endl;

  TCanvas* c1 = new TCanvas;
  c1->cd();
  sm_dilPt_jpts->SetMarkerColor(kBlue);
  sm_dilPt_jpts->SetLineColor(kBlue);
  sm_dilPt_jpts->SetMarkerStyle(22);
  sm_dilPt_jpts->Draw();

  TH1F* craphisto = new TH1F();
  craphisto->SetFillColor(0);
  craphisto->SetMarkerColor(0);
  craphisto->SetLineColor(0);
  craphisto->Draw("sames");
  
  if( met == "tcmet" ) {
    sm_tcmet_jpts->SetMarkerColor(kRed);
    sm_tcmet_jpts->SetLineColor(kRed);
    sm_tcmet_jpts->SetMarkerStyle(24);
    sm_tcmet_jpts->Draw("sames");
  }

  else if( met == "metmuonjes" ) {
    sm_metmuonjes_calo->SetMarkerColor(kRed);
    sm_metmuonjes_calo->SetLineColor(kRed);
    sm_metmuonjes_calo->SetMarkerStyle(24);
    sm_metmuonjes_calo->Draw("sames");
  }

  sm_dilPt_jpts->SetXTitle("p_{T} (GeV)");
  sm_dilPt_jpts->SetTitle("");
  
  TLegend* legend = hist::legend(c1, "p");
  TList* list = legend->GetListOfPrimitives();
  TLegendEntry* leg1 = (TLegendEntry*)list->At(0);

  //if( met == "tcmet" ) {
    leg1->SetLabel("Predicted SM MET distribution");
    TLegendEntry* leg3 = (TLegendEntry*)list->At(1);
    leg3->SetLabel("from dilepton p_{T} distribution");
    TLegendEntry* leg2 = (TLegendEntry*)list->At(2);
    leg2->SetLabel("Observed SM MET distribution");
  //}
/*
  if( met == "metmuonjes" ) {
    leg1->SetLabel("Predicted SM MET distribution");
    TLegendEntry* leg3 = (TLegendEntry*)list->At(1);
    leg3->SetLabel("from dilepton p_{T} distribution");
    TLegendEntry* leg2 = (TLegendEntry*)list->At(2);
    leg2->SetLabel("Observed SM MET distribution");
  }
  */

  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  legend->Draw();

  std::cout << dilpt_prediction100 << " events predicted from dilepton pt distribution for " << met.c_str() << " > 100 GeV" << std::endl;
  std::cout << dilpt_prediction175 << " events predicted from dilepton pt distribution for " << met.c_str() << " > 175 GeV" << std::endl;

  double met_observed100 = -1;
  double met_observed175 = -1;

  if(met == "tcmet") {
    std::cout << "working on tcmet/jpts..." << std::endl;
    int metbin100 = sm_tcmet_jpts->FindBin(100);
    int metbin175 = sm_tcmet_jpts->FindBin(175);
    met_observed100 = sm_tcmet_jpts->Integral(metbin100, 101);
    met_observed175 = sm_tcmet_jpts->Integral(metbin175, 101);
  }
  else if(met == "metmuonjes") {
    std::cout << "working on metmuonsjes/calo..." << std::endl;
    int metbin100 = sm_metmuonjes_calo->FindBin(100);
    int metbin175 = sm_metmuonjes_calo->FindBin(175);
    met_observed100 = sm_metmuonjes_calo->Integral(metbin100, 101);
    met_observed175 = sm_metmuonjes_calo->Integral(metbin175, 101);
  }

  std::cout << met_observed100 << " events observed in " << met.c_str() << " distribution above 100 GeV for SM." << std::endl;  
  std::cout << met_observed175 << " events observed in " << met.c_str() << " distribution above 175 GeV for SM." << std::endl;  

  /*
  TCanvas* canvas = new TCanvas();
  TLegend* legend = new TLegend(0.60,0.62,0.90,0.80);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  */
  
  std::cout << "\n\n" << std::endl;  
  std::cout << setw(15) << "Sample" << setw(25) << "Observed (100)" << setw(25) << "Expected (100)" << setw(25) << "Observed (175)" << setw(25) << "Expected (175)" << std::endl;

  for( unsigned int i = nSMSamples; i < 19; ++i ) {

    TH1F* hmet   = (TH1F*)gDirectory->Get(Form("%s%s%s%s", samples[i], "_h", met.c_str(), "_2j_all"));

    if (hmet   == 0) return;

    int lmbin100 = hmet->FindBin(100);
    int lmbin175 = hmet->FindBin(175);
    double lm_observed100 = hmet->Integral(lmbin100, 101);
    double lm_observed175 = hmet->Integral(lmbin175, 101);

    std::cout << setw(15) << samples[i] << setw(25) << setprecision(3) << lm_observed100+met_observed100 << setw(25) << setprecision(3) << dilpt_prediction100 << setw(25) << setprecision(3) << lm_observed175+met_observed175 << setw(25) << setprecision(3) << dilpt_prediction175 << std::endl;

  }

  std::cout << "\n\n" << std::endl;  
  std::cout << setw(15) << "Sample" << setw(25) << "Observed (100)" << setw(25) << "Observed (175)" << std::endl;

  for( unsigned int j = 0; j < nSMSamples; ++j ) {

    TH1F* hmet   = (TH1F*)gDirectory->Get(Form("%s%s%s%s", samples[j], "_h", met.c_str(), "_2j_all"));

    if (hmet   == 0) return;

    int smbin100 = hmet->FindBin(100);
    int smbin175 = hmet->FindBin(175);
    double sm_observed100 = hmet->Integral(smbin100, 101);
    double sm_observed175 = hmet->Integral(smbin175, 101);

    std::cout << setw(15) << samples[j] << setw(25) << setprecision(3) << sm_observed100 << setw(25) << setprecision(3) << sm_observed175 << std::endl;


  }

    /*
  hdilpt->SetLineColor(1);
  hdilpt->SetLineWidth(2);
    */
//    TH1F* hmet = (TH1F*)(file->Get(Form("%s%s%s%s", samples[i], "_h", met.c_str(), "_2j_all")));
//
//    if (hmet == 0) return;
    /*
  hmet->SetLineColor(2);
  hmet->SetLineWidth(2);
  hmet->SetFillColor(0);

  canvas->cd();

  hdilpt->Draw();
  hmet->Draw("sames");

  hdilpt->SetXTitle("|p_{T}| (GeV)");
  hdilpt->SetTitle(Form("%s%s", "Overlay of dilepton pT and ", met.c_str()));

  legend->AddEntry(hdilpt, "Dilepton pt", "l");
  legend->AddEntry(hmet, Form("%s", met.c_str(), "l"));
  legend->Draw();
    */


//    canvas->SaveAs(Form("%s%s%s", met.c_str(), jet.c_str(), "_top_prediction_dist.root"));
//    canvas->SaveAs(Form("%s%s%s", met.c_str(), jet.c_str(), "_top_prediction_dist.pdf"));

}
