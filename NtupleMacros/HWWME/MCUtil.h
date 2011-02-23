//
// == This Class is designed to access the event based MC truth 
//    information and filling the lepton selection efficiency, 
//    system boost (kT), and fakerate histograms
//    Yanyan Gao (ygao@fnal.gov) 2011/2/20
// ==

#ifndef MCUTIL_H
#define MCUTIL_H

#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TAxis.h"

TFile *mcutilFile_;

TH2F  *els_numer_mc;
TH2F  *els_denom_mc;
TH2F  *els_eff_mc;

TH1F  *els_numer_mc_eta;
TH1F  *els_denom_mc_eta;
TH1F  *els_eff_mc_eta;

TH1F  *els_numer_mc_pt;
TH1F  *els_denom_mc_pt;
TH1F  *els_eff_mc_pt;

TH2F  *mus_numer_mc;
TH2F  *mus_denom_mc;
TH2F  *mus_eff_mc;

TH1F  *mus_numer_mc_eta;
TH1F  *mus_denom_mc_eta;
TH1F  *mus_eff_mc_eta;

TH1F  *mus_numer_mc_pt;
TH1F  *mus_denom_mc_pt;
TH1F  *mus_eff_mc_pt;

TH1F *kx;
TH1F *ky;
TH1F *kt;

void InitMCUtilHist(std::string process) {
  cout << "InitMCUtilHist()"<<endl; 

  if(process != "")
    mcutilFile_ = TFile::Open(string(process + "_MCUtil.root").c_str(),"RECREATE");
  else 
    mcutilFile_ = TFile::Open("MCUtil.root","RECREATE");

  els_numer_mc = new TH2F("els_numer_mc", "els_numer_mc", 20, -2.5, 2.5, 20, 10, 100);
  els_denom_mc = new TH2F("els_denom_mc", "els_denom_mc", 20, -2.5, 2.5, 20, 10, 100);
  els_eff_mc = new TH2F("els_eff_mc", "els_eff_mc", 20, -2.5, 2.5, 20, 10, 100);

  els_numer_mc_eta = new TH1F("els_numer_mc_eta", "els_numer_mc_eta", 20, -2.5, 2.5);
  els_denom_mc_eta = new TH1F("els_denom_mc_eta", "els_denom_mc_eta", 20, -2.5, 2.5);
  els_eff_mc_eta = new TH1F("els_eff_mc_eta", "els_eff_mc_eta", 20, -2.5, 2.5);

  els_numer_mc_pt = new TH1F("els_numer_mc_pt", "els_numer_mc_pt", 20, 10, 100);
  els_denom_mc_pt = new TH1F("els_denom_mc_pt", "els_denom_mc_pt", 20, 10, 100);
  els_eff_mc_pt = new TH1F("els_eff_mc_pt", "els_eff_mc_pt", 20, 10, 100);

  mus_numer_mc = new TH2F("mus_numer_mc", "mus_numer_mc", 20, -2.5, 2.5, 20, 10, 100);
  mus_denom_mc = new TH2F("mus_denom_mc", "mus_denom_mc", 20, -2.5, 2.5, 20, 10, 100);
  mus_eff_mc = new TH2F("mus_eff_mc", "mus_eff_mc", 20, -2.5, 2.5, 20, 10, 100);

  mus_numer_mc_eta = new TH1F("mus_numer_mc_eta", "mus_numer_mc_eta", 20, -2.5, 2.5);
  mus_denom_mc_eta = new TH1F("mus_denom_mc_eta", "mus_denom_mc_eta", 20, -2.5, 2.5);
  mus_eff_mc_eta = new TH1F("mus_eff_mc_eta", "mus_eff_mc_eta", 20, -2.5, 2.5);

  mus_numer_mc_pt = new TH1F("mus_numer_mc_pt", "mus_numer_mc_pt", 20, 10, 100);
  mus_denom_mc_pt = new TH1F("mus_denom_mc_pt", "mus_denom_mc_pt", 20, 10, 100);
  mus_eff_mc_pt = new TH1F("mus_eff_mc_pt", "mus_eff_mc_pt", 20, 10, 100);

  kx = new TH1F("kx", "kx", 50, -50, 50);
  ky = new TH1F("ky", "ky", 50, -50, 50);
  kt = new TH1F("kt", "kt", 100, 0, 100);
}
 
void FillEffHist(TString process, double weight);
void FillKtHist(TString process, double weight);

// Utility Functions
unsigned int getVVType();
bool isIdentified( TString process);

void getEff(double & numer, double & denom, double & eff, double & efferr ) 
{
  if (denom == 0.0) return;
  eff = numer/denom;
  efferr = sqrt(eff*(1-eff)/denom);
}

// 
// This function assumes the two 2D histograms are filled with event counts rather 
// than weighted event counts
// 

void Fill2DEffHist(TH2F* & hist_numer, TH2F* & hist_denom, TH2F* & hist_eff) {
  // cout << "Fill2DEffHist()" << endl;
  if(!hist_numer | !hist_denom ) return;
  // check the hist_numer and hist_denom have the same bin structures
  if(   hist_numer->GetNbinsX() != hist_denom->GetNbinsX() 
	|| hist_numer->GetNbinsY() != hist_denom->GetNbinsY() ) return;
  if(   hist_numer->GetXaxis()->GetXmin() != hist_denom->GetXaxis()->GetXmin() 
	|| hist_numer->GetYaxis()->GetXmin() != hist_denom->GetYaxis()->GetXmin()) return;
  if(   hist_numer->GetXaxis()->GetXmax() != hist_denom->GetXaxis()->GetXmax() 
	|| hist_numer->GetYaxis()->GetXmax() != hist_denom->GetYaxis()->GetXmax()) return;
  
  // include the under/over-flow bins
  for(int iX=0;iX<hist_numer->GetNbinsX()+2;iX++) {
    for(int iY=0;iY<hist_numer->GetNbinsY()+2;iY++) {
      double numer = hist_numer->GetBinContent(iX, iY);
      double denom = hist_denom->GetBinContent(iX, iY);
      double eff(0.), efferr (0.); 
      if(denom!=0) 
	getEff(numer, denom, eff, efferr);
      hist_eff->SetBinContent(iX, iY, eff);
      hist_eff->SetBinError(iX, iY, efferr);
    }
  }
}


// 
// This function assumes the two 1D histograms are filled with event counts rather 
// than weighted event counts
// 

void Fill1DEffHist(TH1F* & hist_numer, TH1F* & hist_denom, TH1F* & hist_eff) {
  // cout << "Fill1DEffHist()" << endl;
  if(!hist_numer | !hist_denom ) return;
  // check the hist_numer and hist_denom have the same bin structures
  if(    hist_numer->GetNbinsX() != hist_denom->GetNbinsX() 
	 || hist_numer->GetXaxis()->GetXmin() != hist_denom->GetXaxis()->GetXmin()
	 || hist_numer->GetXaxis()->GetXmax() != hist_denom->GetXaxis()->GetXmax() ) return;
 
  for(int iX=0;iX<hist_numer->GetNbinsX()+2;iX++) {
    double numer = hist_numer->GetBinContent(iX);
    double denom = hist_denom->GetBinContent(iX);
    double eff(0.), efferr (0.); 
    if(denom!=0) 
      getEff(numer, denom, eff, efferr);

    hist_eff->SetBinContent(iX, eff);
    hist_eff->SetBinError(iX, efferr);
  }
}


void saveMCUtilOutput()
{
    Fill2DEffHist(els_numer_mc, els_denom_mc, els_eff_mc);
    Fill1DEffHist(els_numer_mc_eta, els_denom_mc_eta, els_eff_mc_eta);
    Fill1DEffHist(els_numer_mc_pt, els_denom_mc_pt, els_eff_mc_pt);
    
    Fill2DEffHist(mus_numer_mc, mus_denom_mc, mus_eff_mc);
    Fill1DEffHist(mus_numer_mc_eta, mus_denom_mc_eta, mus_eff_mc_eta);
    Fill1DEffHist(mus_numer_mc_pt, mus_denom_mc_pt, mus_eff_mc_pt);

    mcutilFile_->cd();

    els_denom_mc->Write();
    els_numer_mc->Write();
    els_eff_mc->Write();

    els_denom_mc_eta->Write();
    els_numer_mc_eta->Write();
    els_eff_mc_eta->Write();

    els_denom_mc_pt->Write();
    els_numer_mc_pt->Write();
    els_eff_mc_pt->Write();

    mus_denom_mc->Write();
    mus_numer_mc->Write();
    mus_eff_mc->Write();

    mus_denom_mc_eta->Write();
    mus_numer_mc_eta->Write();
    mus_eff_mc_eta->Write();

    mus_denom_mc_pt->Write();
    mus_numer_mc_pt->Write();
    mus_eff_mc_pt->Write();

    kx->Write();
    ky->Write();
    kt->Write();
    
    mcutilFile_->Close();
  
}


#endif
