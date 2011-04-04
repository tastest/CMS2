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


// Histograms to calculate real lepton efficiencies

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

// System boost

TH1F *kx_;
TH1F *ky_;

// Histograms to calculate the probablitiy of a parton -> lepton

TH2F *parton_;
TH2F *els_fake_;
TH2F *els_genfr_;
TH2F *mus_fake_;
TH2F *mus_genfr_;

TH1F *parton_eta_;
TH1F *els_fake_eta_;
TH1F *els_genfr_eta_;
TH1F *mus_fake_eta_;
TH1F *mus_genfr_eta_;

TH1F *parton_pt_;
TH1F *els_fake_pt_;
TH1F *els_genfr_pt_;
TH1F *mus_fake_pt_;
TH1F *mus_genfr_pt_;

// electron - letpon FO response 

TH2F *els_fo_parton_;
TH2F *mus_fo_parton_;

// Histograms to calculate the probablitiy of a lepton FO -> lepton

TH2F *els_fo_;
TH2F *els_good_;
TH2F *els_fr_;

TH2F *mus_fo_;
TH2F *mus_good_;
TH2F *mus_fr_;


// Not clean area, diagonistic histograms

TH1F *els_fo_parton_dEta_;
TH1F *els_fo_parton_dPhi_;
TH1F *els_fo_MCType_;
TH1F *els_fo_partonID_;

TH1F *mus_fo_parton_dEta_;
TH1F *mus_fo_parton_dPhi_;
TH1F *mus_fo_MCType_;
TH1F *mus_fo_partonID_;

void InitMCUtilHist(const char* process, TFile *utilFile_) {
  cout << "InitMCUtilHist(): " << " process = " << process <<endl; 

  utilFile_->cd();  

  if( TString(process) == "ww") {
  
    els_numer_mc_ = new TH2F(Form("%s_heleNumer", process), "electron numerators", 20, -2.5, 2.5, 20, 0, 100);
    els_denom_mc_ = new TH2F(Form("%s_heleDenom",process),  "electron denominators", 20, -2.5, 2.5, 20, 0, 100);
    els_eff_mc_ = new TH2F(Form("%s_heleEff",process),  "Electron Efficiency", 20, -2.5, 2.5, 20, 0, 100);
    
    els_numer_mc_eta_ = new TH1F(Form("%s_heleNumerEta",process), "electron numerators", 20, -2.5, 2.5);
    els_denom_mc_eta_ = new TH1F(Form("%s_heleDenomEta",process), "electron denominators", 20, -2.5, 2.5);
    els_eff_mc_eta_ = new TH1F(Form("%s_heleEffEta",process), "Electron Efficiency", 20, -2.5, 2.5);
    
    els_numer_mc_pt_ = new TH1F(Form("%s_heleNumerPt",process), "electron numerators", 20, 0, 100);
    els_denom_mc_pt_ = new TH1F(Form("%s_heleDenomPt",process), "electron denominators", 20, 0, 100);
    els_eff_mc_pt_ = new TH1F(Form("%s_heleEffPt",process), "Electron Efficiency", 20, 0, 100);
    
    mus_numer_mc_ = new TH2F(Form("%s_hmuNumer",process), "muon numerators", 20, -2.5, 2.5, 20, 0, 100);
    mus_denom_mc_ = new TH2F(Form("%s_hmuDenom",process), "muon denominators", 20, -2.5, 2.5, 20, 0, 100);
    mus_eff_mc_ = new TH2F(Form("%s_hmuEff",process), "Muon Efficiency", 20, -2.5, 2.5, 20, 0, 100);
    
    mus_numer_mc_eta_ = new TH1F(Form("%s_hmuNumerEta",process), "muon numerators", 20, -2.5, 2.5);
    mus_denom_mc_eta_ = new TH1F(Form("%s_hmuDenomEta",process), "muon denominators", 20, -2.5, 2.5);
    mus_eff_mc_eta_ = new TH1F(Form("%s_hmuEffEta",process),  "Muon Efficiency", 20, -2.5, 2.5);
    
    mus_numer_mc_pt_ = new TH1F(Form("%s_hmuNumerPt",process), "muon numerators", 20, 0, 100);
    mus_denom_mc_pt_ = new TH1F(Form("%s_hmuDenomPt",process), "muon denominators", 20, 0, 100);
    mus_eff_mc_pt_ = new TH1F(Form("%s_hmuEffPt",process), "Muon Efficiency", 20, 0, 100);
    
  }
  
  kx_ = new TH1F(Form("%s_kx",process), "System Boost in X", 50, -50, 50);
  ky_ = new TH1F(Form("%s_ky",process), "System Boost in Y", 50, -50, 50);

  if (TString(process) == "wjets") {
    
    const Double_t ptbins_fakerate[10] = {10.,15.,20.,25.,30.,35., 40., 50., 75., 100.};
    const int nptbins = 9;
    const Double_t etabins_fakerate[6] = {0.0, 0.5, 1.0, 1.479, 2.0, 2.5};
    const int netabins = 5;
    
    // parton to lepton FO

    parton_ = new TH2F(Form("%s_hparton",process), "Generator Parton", netabins,etabins_fakerate,nptbins,ptbins_fakerate);
    els_fake_ = new TH2F(Form("%s_heleFake",process), "FO Electrons",  netabins,etabins_fakerate,nptbins,ptbins_fakerate);
    mus_fake_ = new TH2F(Form("%s_hmuFake",process), "FO Muons", netabins,etabins_fakerate,nptbins,ptbins_fakerate);

    els_genfr_ = new TH2F(Form("%s_heleGenFR",process), "Parton to Electron Generator FR",  netabins,etabins_fakerate,nptbins,ptbins_fakerate);
    els_fo_parton_ = new TH2F(Form("%s_heleFOResponse",process), "matched parton pt / Electron FO pT", nptbins,ptbins_fakerate, 10, 0.5, 4.5);
    mus_genfr_ = new TH2F(Form("%s_hmuGenFR",process), "Parton to Muon Generator FR", netabins,etabins_fakerate,nptbins,ptbins_fakerate);
    mus_fo_parton_ = new TH2F(Form("%s_hmuFOResponse",process), "matched parton pt / Muon FO pT", nptbins,ptbins_fakerate, 10, 0.5, 4.5);

    // lepton FO to lepton 

    els_fo_ = new TH2F(Form("%s_heleFO",process), "Electron FO", netabins,etabins_fakerate,nptbins,ptbins_fakerate);
    els_good_ = new TH2F(Form("%s_heleGood",process), "Electron Passing Analysis Cuts", netabins,etabins_fakerate,nptbins,ptbins_fakerate);
    els_fr_ = new TH2F(Form("%s_heleFR",process), "Electron FO to Good Electron FR",  netabins,etabins_fakerate,nptbins,ptbins_fakerate);

    mus_fo_ = new TH2F(Form("%s_hmuFO",process), "Muon FO", netabins,etabins_fakerate,nptbins,ptbins_fakerate);
    mus_good_ = new TH2F(Form("%s_hmuGood",process), "Muon Passing Analysis Cuts", netabins,etabins_fakerate,nptbins,ptbins_fakerate);
    mus_fr_ = new TH2F(Form("%s_hmuFR",process), "Muon FO to Good Muon FR",  netabins,etabins_fakerate,nptbins,ptbins_fakerate);
	
    // 1-d parton - lepton FO histograms
    
    parton_eta_ = new TH1F(Form("%s_hpartonEta",process), "Generator Parton Eta", netabins, etabins_fakerate);
    els_fake_eta_ = new TH1F(Form("%s_heleFakeEta",process), "FO electrons", netabins, etabins_fakerate);
    els_genfr_eta_ = new TH1F(Form("%s_heleGenFREta",process), "Parton to Electron Generator FR", netabins, etabins_fakerate);
    mus_fake_eta_ = new TH1F(Form("%s_hmuFakeEta",process), "FO muons", netabins, etabins_fakerate);
    mus_genfr_eta_ = new TH1F(Form("%s_hmuGenFREta",process), "Parton to Muon Generator FR", netabins, etabins_fakerate);
    
    
    parton_pt_ = new TH1F(Form("%s_hpartonPt",process), "Generator Parton Pt", nptbins,ptbins_fakerate);
    els_fake_pt_ = new TH1F(Form("%s_heleFakePt",process), "FO electrons",nptbins,ptbins_fakerate);
    els_genfr_pt_ = new TH1F(Form("%s_heleGenFRPt",process), "Parton to Electron Generator FR", nptbins,ptbins_fakerate);
    mus_fake_pt_ = new TH1F(Form("%s_hmuFakePt",process), "FO muons", nptbins,ptbins_fakerate);
    mus_genfr_pt_ = new TH1F(Form("%s_hmuGenFRPt",process), "Parton to Muon Generator FR", nptbins,ptbins_fakerate);
     
    // Diagonistic histograms
    els_fo_parton_dEta_ = new TH1F(Form("%s_hdEtaEleFOParton",process), "#Delta#eta(Lepton FO, Parton)", 20,0,0.2);
    els_fo_parton_dPhi_ = new TH1F(Form("%s_hdPhiEleFOParton",process), "#Delta#eta(Lepton FO, Parton)", 20,0,0.2);
    els_fo_MCType_ = new TH1F(Form("%s_heleFakeMCType",process), "Electron FO MC Type", 7,-1,6);
    els_fo_partonID_ = new TH1F(Form("%s_helePartonID",process), "Electron FO Matched Parton ID",35, -10, 25);

    mus_fo_parton_dEta_ = new TH1F(Form("%s_hdEtaMuFOParton",process), "#Delta#eta(Lepton FO, Parton)", 20,0,0.2);
    mus_fo_parton_dPhi_ = new TH1F(Form("%s_hdPhiMuFOParton",process), "#Delta#eta(Lepton FO, Parton)", 20,0,0.2);
    mus_fo_MCType_ = new TH1F(Form("%s_hmuFakeMCType",process), "Muctron FO MC Type", 7,-1,6);
    mus_fo_partonID_ = new TH1F(Form("%s_hmuPartonID",process), "Muctron FO Matched Parton ID",35, -10, 25);

    
  }
  
}
 
void fillEffHist(const char* process, double weight);
void fillKtHist(const char* process, double weight);
void fillFOHist();
void findClosestEleFO(LorentzVector v_parton, double& minDR, int& idx_minDR);
void findClosestMuFO(LorentzVector v_parton, double& minDR, int& idx_minDR);

// Utility Functions
bool isIdentified(const char* process);

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

void fill2DEffHist(TH2F* hist_numer, TH2F* hist_denom, TH2F* hist_eff) {
  cout << "fill2DEffHist()" << endl;
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

void fill1DEffHist(TH1F* hist_numer, TH1F* hist_denom, TH1F* hist_eff) {
  // cout << "fill1DEffHist()" << endl;
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


void saveMCUtilOutput(const char* process, TFile *utilFile_)
{
  cout << "saveMCUtilHist()" << endl;
  utilFile_->cd();
  
  if(TString(process) == "ww") {
    fill2DEffHist(els_numer_mc_, els_denom_mc_, els_eff_mc_);
    fill1DEffHist(els_numer_mc_eta_, els_denom_mc_eta_, els_eff_mc_eta_);
    fill1DEffHist(els_numer_mc_pt_, els_denom_mc_pt_, els_eff_mc_pt_);
    
    fill2DEffHist(mus_numer_mc_, mus_denom_mc_, mus_eff_mc_);
    fill1DEffHist(mus_numer_mc_eta_, mus_denom_mc_eta_, mus_eff_mc_eta_);
    fill1DEffHist(mus_numer_mc_pt_, mus_denom_mc_pt_, mus_eff_mc_pt_);
    
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
  }
  
  kx_->Write();
  ky_->Write();

  // Fake rate related histograms filled only for wjets MC
  if(TString(process) == "wjets") {
    fill2DEffHist(els_fake_, parton_, els_genfr_);
    fill2DEffHist(mus_fake_, parton_, mus_genfr_);
    fill1DEffHist(els_fake_eta_, parton_eta_, els_genfr_eta_);
    fill1DEffHist(mus_fake_eta_, parton_eta_, mus_genfr_eta_);
    fill1DEffHist(els_fake_pt_, parton_pt_, els_genfr_pt_);
    fill1DEffHist(mus_fake_pt_, parton_pt_, mus_genfr_pt_);
   
    fill2DEffHist(els_good_, els_fo_, els_fr_);
    fill2DEffHist(mus_good_, mus_fo_, mus_fr_);
 
    parton_->Write();
    els_fake_->Write();
    els_genfr_->Write();
    mus_fake_->Write();
    mus_genfr_->Write();

    parton_eta_->Write();
    els_fake_eta_->Write();
    els_genfr_eta_->Write();
    mus_fake_eta_->Write();
    mus_genfr_eta_->Write();
    
    parton_pt_->Write();
    els_fake_pt_->Write();
    els_genfr_pt_->Write();
    mus_fake_pt_->Write();
    mus_genfr_pt_->Write();
    
    els_fo_->Write();
    els_good_->Write();
    els_fr_->Write();

    mus_fo_->Write();
    mus_good_->Write();
    mus_fr_->Write();
    
    els_fo_parton_->Write();
    mus_fo_parton_->Write();
    
    els_fo_parton_dEta_->Write();
    els_fo_parton_dPhi_->Write();
    els_fo_MCType_->Write();
    els_fo_partonID_->Write();
    
    mus_fo_parton_dEta_->Write();
    mus_fo_parton_dPhi_->Write();
    mus_fo_MCType_->Write();
    mus_fo_partonID_->Write();
  }
  
 

}


#endif
