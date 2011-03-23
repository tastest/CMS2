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

TH2F  *els_numer_mc_;
TH2F  *els_denom_mc_;
TH2F  *els_eff_mc_;

TH1F  *els_numer_mc_eta_;
TH1F  *els_denom_mc_eta_;
TH1F  *els_eff_mc_eta_;

TH1F  *els_numer_mc_pt_;
TH1F  *els_denom_mc_pt_;
TH1F  *els_eff_mc_pt_;

TH2F  *mus_numer_mc_;
TH2F  *mus_denom_mc_;
TH2F  *mus_eff_mc_;

TH1F  *mus_numer_mc_eta_;
TH1F  *mus_denom_mc_eta_;
TH1F  *mus_eff_mc_eta_;

TH1F  *mus_numer_mc_pt_;
TH1F  *mus_denom_mc_pt_;
TH1F  *mus_eff_mc_pt_;

TH1F *kx_;
TH1F *ky_;

// These histograms are used for the yield counting in ME
TH1F *dilmass_ee_;
TH1F *dilmass_mm_;
TH1F *dilmass_em_;

void InitMCUtilHist(const char* process, TFile *utilFile_) {
  cout << "InitMCUtilHist(): " << " process = " << process <<endl; 

  utilFile_->cd();
  els_numer_mc_ = new TH2F(Form("%s_heleNumer", process), "Number of electron numerators", 20, -2.5, 2.5, 20, 0, 100);
  els_denom_mc_ = new TH2F(Form("%s_heleDenom",process),  "Number of electron denominators", 20, -2.5, 2.5, 20, 0, 100);
  els_eff_mc_ = new TH2F(Form("%s_heleEff",process),  "Electron Efficiency", 20, -2.5, 2.5, 20, 0, 100);

  els_numer_mc_eta_ = new TH1F(Form("%s_heleNumerEta",process), "Number of electron numerators", 20, -2.5, 2.5);
  els_denom_mc_eta_ = new TH1F(Form("%s_heleDenomEta",process), "Number of electron denominators", 20, -2.5, 2.5);
  els_eff_mc_eta_ = new TH1F(Form("%s_heleEffEta",process), "Electron Efficiency", 20, -2.5, 2.5);

  els_numer_mc_pt_ = new TH1F(Form("%s_heleNumerPt",process), "Number of electron numerators", 20, 0, 100);
  els_denom_mc_pt_ = new TH1F(Form("%s_heleDenomPt",process), "Number of electron denominators", 20, 0, 100);
  els_eff_mc_pt_ = new TH1F(Form("%s_heleEffPt",process), "Electron Efficiency", 20, 0, 100);

  mus_numer_mc_ = new TH2F(Form("%s_hmuNumer",process), "Number of muon numerators", 20, -2.5, 2.5, 20, 0, 100);
  mus_denom_mc_ = new TH2F(Form("%s_hmuDenom",process), "Number of muon denominators", 20, -2.5, 2.5, 20, 0, 100);
  mus_eff_mc_ = new TH2F(Form("%s_hmuEff",process), "Muon Efficiency", 20, -2.5, 2.5, 20, 0, 100);

  mus_numer_mc_eta_ = new TH1F(Form("%s_hmuNumerEta",process), "Number of muon numerators", 20, -2.5, 2.5);
  mus_denom_mc_eta_ = new TH1F(Form("%s_hmuDenomEta",process), "Number of muon denominators", 20, -2.5, 2.5);
  mus_eff_mc_eta_ = new TH1F(Form("%s_hmuEffEta",process),  "Muon Efficiency", 20, -2.5, 2.5);

  mus_numer_mc_pt_ = new TH1F(Form("%s_hmuNumerPt",process), "Number of muon numerators", 20, 0, 100);
  mus_denom_mc_pt_ = new TH1F(Form("%s_hmuDenomPt",process), "Number of muon denominators", 20, 0, 100);
  mus_eff_mc_pt_ = new TH1F(Form("%s_hmuEffPt",process), "Muon Efficiency", 20, 0, 100);

  kx_ = new TH1F(Form("%s_kx",process), "System Boost in X", 50, -50, 50);
  ky_ = new TH1F(Form("%s_ky",process), "System Boost in Y", 50, -50, 50);

  dilmass_ee_ = new TH1F(Form("%s_dilmass_ee",process), "Di-lepton mass in ee", 50, 0, 250);
  dilmass_ee_ ->Sumw2();
  dilmass_em_ = new TH1F(Form("%s_dilmass_em",process), "Di-lepton mass in em", 50, 0, 250);
  dilmass_em_ ->Sumw2();
  dilmass_mm_ = new TH1F(Form("%s_dilmass_mm",process), "Di-lepton mass in mm", 50, 0, 250);
  dilmass_mm_ ->Sumw2();


}
 
void FillEffHist(TString process, double weight);
void FillKtHist(TString process, double weight);

// Utility Functions
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

void Fill2DEffHist(TH2F* hist_numer, TH2F* hist_denom, TH2F* hist_eff) {
  cout << "Fill2DEffHist()" << endl;
  if(!hist_numer || !hist_denom ) return;
  
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

void Fill1DEffHist(TH1F* hist_numer, TH1F* hist_denom, TH1F* hist_eff) {
  // cout << "Fill1DEffHist()" << endl;
  if(!hist_numer || !hist_denom ) return;
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


void saveMCUtilOutput(TFile *utilFile_)
{
  cout << "saveMCUtilHist()" << endl;
  utilFile_->cd();
  
  Fill2DEffHist(els_numer_mc_, els_denom_mc_, els_eff_mc_);
  Fill1DEffHist(els_numer_mc_eta_, els_denom_mc_eta_, els_eff_mc_eta_);
  Fill1DEffHist(els_numer_mc_pt_, els_denom_mc_pt_, els_eff_mc_pt_);
  
  Fill2DEffHist(mus_numer_mc_, mus_denom_mc_, mus_eff_mc_);
  Fill1DEffHist(mus_numer_mc_eta_, mus_denom_mc_eta_, mus_eff_mc_eta_);
  Fill1DEffHist(mus_numer_mc_pt_, mus_denom_mc_pt_, mus_eff_mc_pt_);
  
  els_denom_mc_->Write();
  els_numer_mc_->Write();
  els_eff_mc_->Write();
  
  els_denom_mc_eta_->Write();
  els_numer_mc_eta_->Write();
  els_eff_mc_eta_->Write();
  
  els_denom_mc_pt_->Write();
  els_numer_mc_pt_->Write();
  els_eff_mc_pt_->Write();
  
  mus_denom_mc_->Write();
  mus_numer_mc_->Write();
  mus_eff_mc_->Write();
  
  mus_denom_mc_eta_->Write();
  mus_numer_mc_eta_->Write();
  mus_eff_mc_eta_->Write();
  
  mus_denom_mc_pt_->Write();
  mus_numer_mc_pt_->Write();
  mus_eff_mc_pt_->Write();
  
  kx_->Write();
  ky_->Write();

  dilmass_ee_->Write();
  dilmass_mm_->Write();
  dilmass_em_->Write();
  
}


#endif
