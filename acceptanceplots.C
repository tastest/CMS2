#include "TH1D.h"
#include "TH2D.h"
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

TH1D* hnumerator;
TH1D* hdenominator;
TH1D* hacceptance;

TH2D* hnumerator2d_mtt;
TH2D* hdenominator2d_mtt;
TH2D* hnumerator2drebinned_mtt;
TH2D* hdenominator2drebinned_mtt;
TH2D* hacceptance2drebinned_mtt;

TH2D* hnumerator2d_ttpt;
TH2D* hdenominator2d_ttpt;
TH2D* hnumerator2drebinned_ttpt;
TH2D* hdenominator2drebinned_ttpt;
TH2D* hacceptance2drebinned_ttpt;

TH2D* hnumerator2d_ttrapidity2;
TH2D* hdenominator2d_ttrapidity2;
TH2D* hnumerator2drebinned_ttrapidity2;
TH2D* hdenominator2drebinned_ttrapidity2;
TH2D* hacceptance2drebinned_ttrapidity2;

void acceptanceplots(TString histname = "lepAzimAsym", bool drawnorm = false, TString FName1 = "results/hist_usePtGt2020_hypDisamb_usepfMET_usepfJets_useOS_vetoHypMassLt12_requireBTag_sortJetCandidatesbyPt_generalLeptonVeto_createBabyNtuples_applylepIDCuts_applylepIsoCuts_vetoZmass_veto2Jets_vetoMET.root", TString FName2 = "results/hist_noCuts.root"){
  setTDRStyle();

  std::cout << "Opening " << FName1.Data() << "\n";
  TFile *f_1         = TFile::Open(FName1.Data());  
  hnumerator = (TH1D*)f_1->Get(Form("ttdil_h%sGen_allj_all", histname.Data())); 
  hnumerator2d_mtt = (TH2D*)f_1->Get(Form("ttdil_h%sGen2d_allj_all", histname.Data())); 
  hnumerator2d_ttpt = (TH2D*)f_1->Get(Form("ttdil_h%sttpTGen2d_allj_all", histname.Data())); 
  hnumerator2d_ttrapidity2 = (TH2D*)f_1->Get(Form("ttdil_h%sttRapidity2Gen2d_allj_all", histname.Data())); 

  std::cout << "Opening " << FName2.Data() << "\n";  
  TFile *f_2         = TFile::Open(FName2.Data());  
  hdenominator = (TH1D*)f_2->Get(Form("ttdil_h%sGen_allj_all", histname.Data()));
  hdenominator2d_mtt = (TH2D*)f_2->Get(Form("ttdil_h%sGen2d_allj_all", histname.Data())); 
  hdenominator2d_ttpt = (TH2D*)f_2->Get(Form("ttdil_h%sttpTGen2d_allj_all", histname.Data())); 
  hdenominator2d_ttrapidity2 = (TH2D*)f_2->Get(Form("ttdil_h%sttRapidity2Gen2d_allj_all", histname.Data())); 

  std::cout << "Opened " << Form("ttdil_h%sGen_allj_all", histname.Data()) << " and "<< Form("ttdil_h%sGen2d_allj_all", histname.Data()) <<"\n";

  Double_t pi = 3.141592653589793;
  Double_t bins_lepChargeAsym[] =  { -4., -0.8, -0.4, 0., 0.4, 0.8, 4.}; 
  Double_t bins_lepAzimAsym[] = {-1., -0.8, -0.4, 0., 0.4, 0.8, 1.}; 
  Double_t bins_lepAzimAsym2[] = {0., 4.*pi/20., 7.*pi/20., 10.*pi/20., 13.*pi/20., 16.*pi/20., pi}; 
  Double_t bins_topCosTheta[] = {-1., -0.7, -0.4, 0., 0.4, 0.7, 1.}; 
  Double_t bins_pseudorapiditydiff[] =  { -4., -1.0, -0.5, 0., 0.5, 1.0, 4.}; 
  Double_t bins_rapiditydiff[] =  { -4., -0.8, -0.3, 0., 0.3, 0.8, 4.}; 
  Double_t bins_rapiditydiffMarco[] =  { -4., -0.7, -0.3, 0., 0.3, 0.7, 4.}; 
  Double_t bins_lepCosTheta[] = {-1., -0.6, -0.3, 0., 0.3, 0.6, 1.}; 
  Double_t bins_topSpinCorr[] = {-1., -0.5, -0.2, 0., 0.2, 0.5, 1.}; 

  Double_t bins_lepChargeAsym_for2D[] =  { -4., 0., 4.};
  Double_t bins_lepAzimAsym_for2D[] = {-1., 0., 1.};
  Double_t bins_lepAzimAsym2_for2D[] = {0., pi/2., pi};
  Double_t bins_topCosTheta_for2D[] = {-1., 0., 1.};
  Double_t bins_pseudorapiditydiff_for2D[] =  { -4., 0., 4.};
  Double_t bins_rapiditydiff_for2D[] =  { -4., 0., 4.};
  Double_t bins_rapiditydiffMarco_for2D[] =  { -4., 0., 4.};
  Double_t bins_lepCosTheta_for2D[] = {-1., 0., 1.};
  Double_t bins_topSpinCorr_for2D[] = {-1., 0., 1.};

  Double_t binsmtt[] = {0., 410., 510., 1200.}; 
  Double_t binsttpt[] = {0., 24., 52., 300}; 
  Double_t binsttrapidity2[] = {0., 0.3, 0.7, 3.0}; 
  Double_t bins[7];
  Double_t binsfor2D[3];


  if(histname == "lepChargeAsym") memcpy(bins,bins_lepChargeAsym,7*8);
  if(histname == "lepAzimAsym") memcpy(bins,bins_lepAzimAsym,7*8);
  if(histname == "lepAzimAsym2") memcpy(bins,bins_lepAzimAsym2,7*8);
  if(histname == "topCosTheta") memcpy(bins,bins_topCosTheta,7*8);
  if(histname == "pseudorapiditydiff") memcpy(bins,bins_pseudorapiditydiff,7*8);
  if(histname == "rapiditydiff") memcpy(bins,bins_rapiditydiff,7*8);
  if(histname == "rapiditydiffMarco") memcpy(bins,bins_rapiditydiffMarco,7*8);
  if(histname == "lepCosTheta" || histname == "lepPlusCosTheta" || histname == "lepMinusCosTheta") memcpy(bins,bins_lepCosTheta,7*8);
  if(histname == "topSpinCorr") memcpy(bins,bins_topSpinCorr,7*8);

  if(histname == "lepChargeAsym") memcpy(binsfor2D,bins_lepChargeAsym_for2D,3*8);
  if(histname == "lepAzimAsym") memcpy(binsfor2D,bins_lepAzimAsym_for2D,3*8);
  if(histname == "lepAzimAsym2") memcpy(binsfor2D,bins_lepAzimAsym2_for2D,3*8);
  if(histname == "topCosTheta") memcpy(binsfor2D,bins_topCosTheta_for2D,3*8);
  if(histname == "pseudorapiditydiff") memcpy(binsfor2D,bins_pseudorapiditydiff_for2D,3*8);
  if(histname == "rapiditydiff") memcpy(binsfor2D,bins_rapiditydiff_for2D,3*8);
  if(histname == "rapiditydiffMarco") memcpy(binsfor2D,bins_rapiditydiffMarco_for2D,3*8);
  if(histname == "lepCosTheta" || histname == "lepPlusCosTheta" || histname == "lepMinusCosTheta") memcpy(binsfor2D,bins_lepCosTheta_for2D,3*8);
  if(histname == "topSpinCorr") memcpy(binsfor2D,bins_topSpinCorr_for2D,3*8);

  hnumerator = (TH1D*) hnumerator->Rebin(6,Form("numerator_%s", histname.Data()),bins);
  hdenominator = (TH1D*) hdenominator->Rebin(6,Form("denominator_%s", histname.Data()),bins);

  hnumerator2drebinned_mtt = new TH2D(Form("numerator_%s_mtt", histname.Data()),Form("numerator_%s_mtt", histname.Data()),2,binsfor2D,3, binsmtt);
  TAxis *xaxis = hnumerator2d_mtt->GetXaxis();
  TAxis *yaxis = hnumerator2d_mtt->GetYaxis();
  for (int j=1;j<=yaxis->GetNbins();j++) {
    for (int i=1;i<=xaxis->GetNbins();i++) {
      hnumerator2drebinned_mtt->Fill(xaxis->GetBinCenter(i),yaxis->GetBinCenter(j), hnumerator2d_mtt->GetBinContent(i,j));
    }
  }
  hdenominator2drebinned_mtt = new TH2D(Form("denominator_%s_mtt", histname.Data()),Form("denominator_%s_mtt", histname.Data()),2,binsfor2D,3, binsmtt);
  TAxis *xaxisd = hdenominator2d_mtt->GetXaxis();
  TAxis *yaxisd = hdenominator2d_mtt->GetYaxis();
  for (int j=1;j<=yaxisd->GetNbins();j++) {
    for (int i=1;i<=xaxisd->GetNbins();i++) {
      hdenominator2drebinned_mtt->Fill(xaxisd->GetBinCenter(i),yaxisd->GetBinCenter(j), hdenominator2d_mtt->GetBinContent(i,j));
    }
  }

  hnumerator2drebinned_ttpt = new TH2D(Form("numerator_%s_ttpt", histname.Data()),Form("numerator_%s_ttpt", histname.Data()),2,binsfor2D,3, binsttpt);
  xaxis = hnumerator2d_ttpt->GetXaxis();
  yaxis = hnumerator2d_ttpt->GetYaxis();
  for (int j=1;j<=yaxis->GetNbins();j++) {
    for (int i=1;i<=xaxis->GetNbins();i++) {
      hnumerator2drebinned_ttpt->Fill(xaxis->GetBinCenter(i),yaxis->GetBinCenter(j), hnumerator2d_ttpt->GetBinContent(i,j));
    }
  }
  hdenominator2drebinned_ttpt = new TH2D(Form("denominator_%s_ttpt", histname.Data()),Form("denominator_%s_ttpt", histname.Data()),2,binsfor2D,3, binsttpt);
  xaxisd = hdenominator2d_ttpt->GetXaxis();
  yaxisd = hdenominator2d_ttpt->GetYaxis();
  for (int j=1;j<=yaxisd->GetNbins();j++) {
    for (int i=1;i<=xaxisd->GetNbins();i++) {
      hdenominator2drebinned_ttpt->Fill(xaxisd->GetBinCenter(i),yaxisd->GetBinCenter(j), hdenominator2d_ttpt->GetBinContent(i,j));
    }
  }

  hnumerator2drebinned_ttrapidity2 = new TH2D(Form("numerator_%s_ttrapidity2", histname.Data()),Form("numerator_%s_ttrapidity2", histname.Data()),2,binsfor2D,3, binsttrapidity2);
  xaxis = hnumerator2d_ttrapidity2->GetXaxis();
  yaxis = hnumerator2d_ttrapidity2->GetYaxis();
  for (int j=1;j<=yaxis->GetNbins();j++) {
    for (int i=1;i<=xaxis->GetNbins();i++) {
      hnumerator2drebinned_ttrapidity2->Fill(xaxis->GetBinCenter(i),yaxis->GetBinCenter(j), hnumerator2d_ttrapidity2->GetBinContent(i,j));
    }
  }
  hdenominator2drebinned_ttrapidity2 = new TH2D(Form("denominator_%s_ttrapidity2", histname.Data()),Form("denominator_%s_ttrapidity2", histname.Data()),2,binsfor2D,3, binsttrapidity2);
  xaxisd = hdenominator2d_ttrapidity2->GetXaxis();
  yaxisd = hdenominator2d_ttrapidity2->GetYaxis();
  for (int j=1;j<=yaxisd->GetNbins();j++) {
    for (int i=1;i<=xaxisd->GetNbins();i++) {
      hdenominator2drebinned_ttrapidity2->Fill(xaxisd->GetBinCenter(i),yaxisd->GetBinCenter(j), hdenominator2d_ttrapidity2->GetBinContent(i,j));
    }
  }


  TString accepthistname = "accept_";
  accepthistname += histname;

  TFile *output = new TFile(Form("%s.root", accepthistname.Data()), "RECREATE");  

  hacceptance =  (TH1D*) hnumerator->Clone(accepthistname.Data());
  hacceptance->SetTitle(accepthistname.Data());
  hacceptance->Reset();
  hacceptance->Divide(hnumerator,hdenominator,1., 1.);

  hacceptance2drebinned_mtt =  (TH2D*) hnumerator2drebinned_mtt->Clone( Form("%s_mtt", accepthistname.Data()) );
  hacceptance2drebinned_mtt->Reset();
  hacceptance2drebinned_mtt->SetTitle(Form("%s_mtt", accepthistname.Data()));
  hacceptance2drebinned_mtt->Divide(hnumerator2drebinned_mtt,hdenominator2drebinned_mtt,1., 1.);

  hacceptance2drebinned_ttpt =  (TH2D*) hnumerator2drebinned_ttpt->Clone( Form("%s_ttpt", accepthistname.Data()) );
  hacceptance2drebinned_ttpt->Reset();
  hacceptance2drebinned_ttpt->SetTitle(Form("%s_ttpt", accepthistname.Data()));
  hacceptance2drebinned_ttpt->Divide(hnumerator2drebinned_ttpt,hdenominator2drebinned_ttpt,1., 1.);

  hacceptance2drebinned_ttrapidity2 =  (TH2D*) hnumerator2drebinned_ttrapidity2->Clone( Form("%s_ttrapidity2", accepthistname.Data()) );
  hacceptance2drebinned_ttrapidity2->Reset();
  hacceptance2drebinned_ttrapidity2->SetTitle(Form("%s_ttrapidity2", accepthistname.Data()));
  hacceptance2drebinned_ttrapidity2->Divide(hnumerator2drebinned_ttrapidity2,hdenominator2drebinned_ttrapidity2,1., 1.);

  hnumerator->SetLineColor(kBlue);
  hnumerator-> SetFillColor(0);
  hnumerator->SetMarkerColor(kBlue);
  hdenominator->SetLineColor(kRed);
  hdenominator->SetMarkerColor(kRed);
  hdenominator-> SetFillColor(0);
  hacceptance->SetLineColor(kBlack);
  hacceptance->SetMarkerColor(kBlack);
  hacceptance-> SetFillColor(0);

  gStyle->SetPaintTextFormat("6.4f");

  TCanvas *c1 = new TCanvas();
  c1->cd();

  hacceptance->SetMaximum(1.25*hacceptance->GetMaximum());
  if(hacceptance->GetMinimum() <0.15 *hacceptance->GetMaximum() ) hacceptance->SetMinimum(0.);  
  if(hacceptance->GetMinimum() > 0.) hacceptance->SetMinimum(0.75*hacceptance->GetMinimum() );  

  hacceptance->GetYaxis()->SetTitle("Acceptance");
  hacceptance->GetYaxis()->SetTitleOffset(1.6);

  if(histname.Contains("lepPlusCosTheta") ) {
    hacceptance->GetXaxis()->SetTitle("cos(#theta^{+}_{l})");
  }
  if(histname.Contains("lepMinusCosTheta") ) {
    hacceptance->GetXaxis()->SetTitle("cos(#theta^{-}_{l})");
  }
  if(histname.Contains("lepChargeAsym") ) {
    hacceptance->GetXaxis()->SetTitle(" |#eta_{l^{+}}| - |#eta_{l^{-}}| ");
  }

  if(!drawnorm){
    hacceptance->Draw("hist TEXT00E");
  }
  else{
    hacceptance->DrawNormalized("hist");
    hnumerator->DrawNormalized("histsame");
    hdenominator->DrawNormalized("histsame");
  }

  TLegend *leg = new TLegend(0.74,0.86,0.90,0.92);
  leg->SetFillStyle(0);
  leg->SetLineColor(0);
  leg->SetTextSize(0.032);
  leg->SetFillStyle(0);
  leg->AddEntry(hacceptance, "acceptance","l");
  if(drawnorm) leg->AddEntry(hnumerator, "numerator","l");
  if(drawnorm) leg->AddEntry(hdenominator, "denominator","l");

  TPaveText *pt1 = new TPaveText(0.18, 0.86, 0.40, 0.92, "brNDC");
  pt1->SetName("pt1name");
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(0);

  TText *blah;
  blah = pt1->AddText("CMS Preliminary, 5.0 fb^{-1} at  #sqrt{s}=7 TeV");
  blah->SetTextSize(0.032);
  blah->SetTextAlign(11);  

  pt1->Draw();

  c1->Print(Form("%s.pdf", accepthistname.Data()));

  hacceptance->Write();
  hdenominator->Write();
  hnumerator->Write();

  hdenominator2d_mtt->Write();
  hnumerator2d_mtt->Write();
  hacceptance2drebinned_mtt->Write();
  hdenominator2drebinned_mtt->Write();
  hnumerator2drebinned_mtt->Write();

  hdenominator2d_ttpt->Write();
  hnumerator2d_ttpt->Write();
  hacceptance2drebinned_ttpt->Write();
  hdenominator2drebinned_ttpt->Write();
  hnumerator2drebinned_ttpt->Write();

  hdenominator2d_ttrapidity2->Write();
  hnumerator2d_ttrapidity2->Write();
  hacceptance2drebinned_ttrapidity2->Write();
  hdenominator2drebinned_ttrapidity2->Write();
  hnumerator2drebinned_ttrapidity2->Write();

  f_1->Close();
  f_2->Close();
  output->Close();

}
