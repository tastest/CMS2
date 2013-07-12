#include "TH1D.h"
#include "TString.h"
#include "TCanvas.h"
#include <iostream>
#include "TFile.h"
#include "TDirectory.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TLatex.h"
#include <fstream>
#include "tdrStyle.C"
#include "CommonFunctions.C"
#include <vector>

using namespace std;

TH1D* histo1;
TH1D* histo2;
TH1D* histo3;
TH1D* datahisto;




void makejsbiasplots(TString histname = "ttMasspull",  int drawlog =0, double rangelow = -1, double rangehigh = 3, double rangeylow = 0, double rangeyhigh = 0, int rebin = 1,TString FName1 = "results/hist_usePtGt2020_hypDisamb_usepfMET_usepfJets_useOS_vetoHypMassLt12_requireBTag_sortJetCandidatesbyPt_generalLeptonVeto_createBabyNtuples_applylepIDCuts_applylepIsoCuts_vetoZmass_veto2Jets_vetoMET.root"){
  setTDRStyle();

  bool is2D = rangeyhigh == 0 ? 0 : 1; 
  std::cout << "Opening " << FName1.Data() << "\n";
  TFile *f_1         = TFile::Open(FName1.Data());  
  histo1 = (TH1D*)f_1->Get(Form("ttdil_h%s_allj_all", histname.Data()));
  std::cout << "hist " << histo1->GetName() << " with entries " << histo1->GetEntries() << std::endl;
  
  gStyle->SetOptStat(1001001100);
  if(drawlog && !is2D) gStyle->SetOptLogy(1);
  else gStyle->SetOptLogy(0);

  if(drawlog && is2D) gStyle->SetOptLogz(1);
  else gStyle->SetOptLogz(0);

  if(!is2D) histo1->Rebin(rebin);
  histo1->GetXaxis()->SetRangeUser(rangelow,rangehigh);
  if(is2D) histo1->GetYaxis()->SetRangeUser(rangeylow,rangeyhigh);
  histo1->SetLineColor(kBlue);
  histo1-> SetFillColor(0);

  TCanvas *c1 = new TCanvas();
  c1->cd();

  if(!is2D) histo1->Draw("hist");
  else histo1->Draw("COLZ");

  /*
  TLegend *leg = new TLegend(0.74,0.76,0.92,0.92);
  leg->SetFillStyle(0);
  leg->SetLineColor(0);
  leg->SetTextSize(0.035);
  leg->SetFillStyle(0);
  leg->AddEntry(histo1, "ttdil","l");

  leg->Draw("same");
  */

/*
  TPaveText *pt1 = new TPaveText(0.20, 0.76, 0.40, 0.91, "brNDC");
  pt1->SetName("pt1name");
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(0);
  
  TText *blah;
  blah = pt1->AddText("CMS Preliminary, 5.0 fb^{-1} at  #sqrt{s}=7 TeV");
  blah->SetTextSize(0.032);
  blah->SetTextAlign(11);  

  pt1->Draw();
*/

  if(!drawlog) c1->Print(Form("%s.pdf", histname.Data()));
  else c1->Print(Form("%s_log.pdf", histname.Data()));

  f_1->Close();
  
}
