#ifndef COMMONFUNCTIONS_C
#define COMMONFUNCTIONS_C

#include "TFile.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TArrow.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "THStack.h"
#include <iostream>
#include <stdlib.h>


float GetEntries(const TH1F *h, int lowBin, int highBin) {
  
  if(lowBin < 0 || highBin > h->GetNbinsX()+1) {
    cout << "Bins out of range. Setting the lowBin to the underflow and the highBin to the overflow" << endl;
    lowBin = 0;
    highBin = h->GetNbinsX();
  }
  
  float nentries = 0;
  for(int i = lowBin; i < highBin+1; i++) 
    nentries = nentries + h->GetBinContent(i);
  
  return nentries;
}

float GetTotalError(const TH1F *h, int lowBin, int highBin) {
  
  if(lowBin < 0 || highBin > h->GetNbinsX()+1) {
    cout << "Bins out of range. Setting the lowBin to the underflow and the highBin to the overflow" << endl;
    lowBin = 0;
    highBin = h->GetNbinsX();
  }
  
  float err2 = 0;
  for(int i = lowBin; i < highBin+1; i++) 
    err2 = err2 + pow(h->GetBinError(i),2);
  
  return sqrt(err2);
  
}


std::string formatFloat(double x, const char* formatS) {
  std::string xS = Form(Form("%s", formatS),x);
  double xB = atof(xS.c_str());
  if (x>0 && xB==0){
    xS = Form(" %6.1g",x);
  }
  return xS;
}


void setStyle(TH1F *& hist, int rebin, bool scale, Color_t color, Style_t marker, TString x_title, TString y_title, float xMin, float xMax)
{
  hist->Rebin(rebin);
  if(scale)
    hist->Scale(1.0/rebin);
  hist->SetLineColor(color);
  hist->SetLineWidth(3);
  hist->SetMarkerColor(color);
  hist->SetMarkerStyle(marker);
  hist->SetMarkerSize(1.0);
  hist->GetXaxis()->SetTitle(x_title);
  hist->GetYaxis()->SetTitle(y_title);
  hist->SetTitleOffset(1.3,"Y");
  hist->GetXaxis()->SetRangeUser(xMin,xMax);
  hist->SetStats(0);
  hist->SetFillColor(0);
  //hist->SetBinContent(1, hist->GetBinContent(0)+hist->GetBinContent(1));
}

void setStyle(TH1F *& hist, int rebin, Color_t color, Style_t marker, TString x_title, TString y_title, float xMin, float xMax)
{
  hist->Rebin(rebin);
  hist->Scale(1.0/rebin);
  hist->SetLineColor(color);
  hist->SetLineWidth(3);
  hist->SetMarkerColor(color);
  hist->SetMarkerStyle(marker);
  hist->SetMarkerSize(1.0);
  hist->GetXaxis()->SetTitle(x_title);
  hist->GetYaxis()->SetTitle(y_title);
  hist->SetTitle("");
  hist->SetTitleOffset(1.8,"Y");
  hist->GetXaxis()->SetRangeUser(xMin,xMax);
  hist->SetStats(0);
  hist->SetFillColor(0);
}


TArrow *getArrow(const int nHist, TH1F **hist, float cutVal)
{
  float max (0);
  for(int i=0; i<nHist;i++) 
    max = max > hist[i]->GetMaximum() ? max : hist[i]->GetMaximum();
    
  TArrow *arr_cut = new TArrow(cutVal, max, cutVal, 0, 0.08, "|>");
  arr_cut->SetAngle(25.0);
  arr_cut->SetLineWidth(3.0);
  arr_cut->SetLineColor(kBlack);
  arr_cut->SetFillColor(kBlack);
  arr_cut->SetLineStyle(kDashed);
  arr_cut->SetFillStyle(3001);
  return arr_cut;
}

TLegend *getLegend(const int nHist, TH1F**hist, TString* & label, double x1, double y1, double x2, double y2) {
  TLegend *leg = new TLegend(x1, y1, x2, y2, "", "brNDC");
  for(int i = 0;i<nHist;i++) 
    leg->AddEntry(hist[i], label[i], "lp");
  leg->SetFillColor(0);
  return leg;
}

int getBin(TH1F* h, float x) 
{
  int bin_sel(-1);
  for(int i=0;i<h->GetNbinsX();i++) {
    if(x>h->GetBinLowEdge(i) && x <= h->GetBinLowEdge(i) + h->GetBinWidth(i)) 
      bin_sel = i;
  }
  return bin_sel;
}




void drawPlot(TH1F** hist, const int nHist, TString plotName, bool logy, std::vector<TString> legend)
{
  TCanvas *c1 = new TCanvas("c1","c1");
  c1->cd();  
  gPad->SetGridy();
  //gPad->SetGridx();
    // save and draw the plots
  TLegend *leg = new TLegend(0.45, 0.15, 0.85, 0.35);
  leg->SetFillColor(0);

  for(int i=0;i<nHist;i++) {
    leg->AddEntry(hist[i], legend[i], "l");
  }

  if(logy) gPad->SetLogy();  
  for(int i=0;i<nHist;i++) {
    if (i==0) hist[i]->Draw("h");
    else
      hist[i]->Draw("sameh");
  }
  leg->Draw("SAME"); 
  c1->SaveAs(TString("epsfiles/"+plotName+".eps"));
  c1->SaveAs(TString("pngfiles/"+plotName+".png"));

  delete leg;  
  delete c1;
}


void getEff(TH1F*  hist, TH1F* hist_eff) {
  for (int i=2;i<hist->GetNbinsX()+2;i++) {
    double denom = hist->Integral(0, hist->GetNbinsX()+1);
    Double_t num = hist->Integral(0, i-1);
    Double_t eff = num/denom;
    hist_eff->SetBinContent(i-1,eff);
    hist_eff->SetBinError( i-1, sqrt(eff*(1-eff)/hist->GetEntries()));
  }
}

void fixRangeY(TH1* r,TH1* s){
  double ymin = (r->GetBinContent(r->GetMinimumBin()) < s->GetBinContent(s->GetMinimumBin())) ? 
    r->GetBinContent(r->GetMinimumBin()) : s->GetBinContent(s->GetMinimumBin());
  double ymax = (r->GetBinContent(r->GetMaximumBin()) > s->GetBinContent(s->GetMaximumBin())) ?
    r->GetBinContent(r->GetMaximumBin()) : s->GetBinContent(s->GetMaximumBin());
  r->GetYaxis()->SetRangeUser(ymin*0.9,ymax*1.1);
  s->GetYaxis()->SetRangeUser(ymin*0.9,ymax*1.1);
}


void fixRangeY(TH1 **h, const int nHist) {
  double ymin = h[0]->GetBinContent(h[0]->GetMinimumBin());
  double ymax =  h[0]->GetBinContent(h[0]->GetMaximumBin());
  for(int i=1;i<nHist;i++) {
    ymin = ymin < h[i]->GetBinContent(h[i]->GetMinimumBin()) ? ymin : h[i]->GetBinContent(h[i]->GetMinimumBin());
    ymax = ymax > h[i]->GetBinContent(h[i]->GetMaximumBin()) ? ymax : h[i]->GetBinContent(h[i]->GetMaximumBin());
  }
  for (int i=0;i<nHist;i++) 
    h[i]->GetYaxis()->SetRangeUser(ymin*0.9,ymax*1.1);
}

#endif
