#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include <iostream>
#include <iomanip>
#include <string>
#include "wwtypes.cc"

void getEstimate(const TH1F* h, double& estimate, double& error, const TH1F* qcdShape){
  estimate = 0;
  error = 0;
  double sideband = 0;
  double sideband_err = 0;
  for ( Int_t x=1; x<=h->GetNbinsX(); ++x ){
    if ( h->GetBinLowEdge(x) >= 0.2 && h->GetBinLowEdge(x) <= 1.0 ){
      sideband += h->GetBinContent(x);
      sideband_err += h->GetBinError(x)*h->GetBinError(x);
    }
  } 
  estimate = sideband*0.1/0.8;
  error = sqrt(sideband_err)*0.1/0.8;
}

void printLine( const char*name,
		TH1F* DYee, TH1F* DYmm, TH1F* DYtt, TH1F* tt, TH1F* wjets, TH1F* wz, TH1F* zz, TH1F* ww, TH1F* tw,
		TH1F* qcdShape)
{
  string pm = "+/-";
  // string pm = "&plusmn;";
  
  cout << "|" << Form(" %12s ",name) << "|";
  double integral(0);
  double error(0);
  if (DYee) {
    getEstimate(DYee, integral, error, qcdShape);
    cout << Form(" %5.2f%s%4.2f ", integral, pm.c_str(), error);
  } 
  else
    cout << "    skipped   ";
  cout << "|";
  if (DYmm) {
    getEstimate(DYmm, integral, error, qcdShape);
    cout << Form(" %5.2f%s%4.2f ", integral, pm.c_str(), error);
  }
  else
    cout << "    skipped   ";
  cout << "|";
  if (DYtt) {
    getEstimate(DYtt, integral, error, qcdShape);
    cout << Form(" %5.2f%s%4.2f ", integral, pm.c_str(), error);
  }
  else
    cout << "    skipped   ";
  cout << "|";
  if (tt) {
    getEstimate(tt, integral, error, qcdShape);
    cout << Form(" %5.2f%s%4.2f ", integral, pm.c_str(), error);
  }
  else
    cout << "    skipped   ";
  cout << "|";
  if (wjets) {
    getEstimate(wjets, integral, error, qcdShape);
    cout << Form(" %5.2f%s%4.2f ", integral, pm.c_str(), error);
  }
  else
    cout << "    skipped   ";
  cout << "|";
  if (wz) {
    getEstimate(wz, integral, error, qcdShape);
    cout << Form(" %5.2f%s%4.2f ", integral, pm.c_str(), error);
  }
  else
    cout << "    skipped   ";
  cout << "|";
  if (zz) {
    getEstimate(zz, integral, error, qcdShape);
    cout << Form(" %5.2f%s%4.2f ", integral, pm.c_str(), error);
  }
  else
    cout << "    skipped   ";
  cout << "|";
  if (ww) {
    getEstimate(ww, integral, error, qcdShape);
    cout << Form(" %5.2f%s%4.2f ", integral, pm.c_str(), error);
  }
  else
    cout << "    skipped   ";
  cout << "|";
  if (tw) {
    getEstimate(tw, integral, error, qcdShape);
    cout << Form(" %5.2f%s%4.2f ", integral, pm.c_str(), error);
  }
  else
    cout << "    skipped   ";
  cout << "|" <<endl;
  
}
void estimateJetBkg_sideband(const char* file = "processed_data.root")
{
  using namespace std;
  
  //hist file:
  TFile *ftt = TFile::Open(file);
  assert(ftt);

  //do the emu row:
  TH1F *DYee_mm  = dynamic_cast<TH1F*>(ftt->Get("dyee_hMuRelIso_mm"));
  TH1F *DYee_em  = dynamic_cast<TH1F*>(ftt->Get("dyee_hMuRelIso_em"));
  TH1F *DYmm_mm  = dynamic_cast<TH1F*>(ftt->Get("dymm_hMuRelIso_mm"));
  TH1F *DYmm_em  = dynamic_cast<TH1F*>(ftt->Get("dymm_hMuRelIso_em"));
  TH1F *DYtt_mm  = dynamic_cast<TH1F*>(ftt->Get("dytt_hMuRelIso_mm"));
  TH1F *DYtt_em  = dynamic_cast<TH1F*>(ftt->Get("dytt_hMuRelIso_em"));
  TH1F *tt_mm    = dynamic_cast<TH1F*>(ftt->Get("ttbar_hMuRelIso_mm"));
  TH1F *tt_em    = dynamic_cast<TH1F*>(ftt->Get("ttbar_hMuRelIso_em"));
  TH1F *wjets_mm = dynamic_cast<TH1F*>(ftt->Get("wjets_hMuRelIso_mm"));
  TH1F *wjets_em = dynamic_cast<TH1F*>(ftt->Get("wjets_hMuRelIso_em"));
  TH1F *wz_mm    = dynamic_cast<TH1F*>(ftt->Get("wz_hMuRelIso_mm"));
  TH1F *wz_em    = dynamic_cast<TH1F*>(ftt->Get("wz_hMuRelIso_em"));
  TH1F *zz_mm    = dynamic_cast<TH1F*>(ftt->Get("zz_hMuRelIso_mm"));
  TH1F *zz_em    = dynamic_cast<TH1F*>(ftt->Get("zz_hMuRelIso_em"));
  TH1F *ww_mm    = dynamic_cast<TH1F*>(ftt->Get("ww_hMuRelIso_mm"));
  TH1F *ww_em    = dynamic_cast<TH1F*>(ftt->Get("ww_hMuRelIso_em"));
  TH1F *tw_mm    = dynamic_cast<TH1F*>(ftt->Get("tw_hMuRelIso_mm"));
  TH1F *tw_em    = dynamic_cast<TH1F*>(ftt->Get("tw_hMuRelIso_em"));
  
  cout << "Jet induced muon fake background estimation ([0.2,1.0]*0.1/0.8):" <<endl;
  cout << "\n" << Form("| %12s | %12s | %12s | %12s | %12s | %12s | %12s | %12s | %12s | %12s |",
		       "", "*DY ee*","*DY mumu*","*DY tautau*","*ttbar*","*Wjets*","*WZ*","*ZZ*","*WW*","*TW*")
       << endl;
  printLine("mm - flat", DYee_mm, DYmm_mm, DYtt_mm, tt_mm, wjets_mm, wz_mm, zz_mm, ww_mm, tw_mm, 0);
  printLine("em - flat", DYee_em, DYmm_em, DYtt_em, tt_em, wjets_em, wz_em, zz_em, ww_em, tw_em, 0);
  cout <<endl;
  /*
  //
  // Apply fake rates
  //
  TH1F *qcd_el_fake  = dynamic_cast<TH1F*>(ftt->Get("qcd_hFakableRateSingleElectron"));
  TH1F *qcd_el_final = dynamic_cast<TH1F*>(ftt->Get("qcd_hFinalRateSingleElectron"));
  if ( qcd_el_fake && qcd_el_final ){
    TH1F *qcd_el_fake_not_final = (TH1F*)qcd_el_fake->Clone("qcd_el_fake_not_final");
  
    // nominal fake rate
    TH1F *fakeRate = (TH1F*)qcd_el_final->Clone("qcd_el_fakeRate");
    fakeRate->Divide(qcd_el_fake);

    // fake_rate/(1-fake_rate)
    TH1F *qcdShape = (TH1F*)qcd_el_final->Clone("qcd_el_qcdShape");
    qcdShape->Divide(qcd_el_fake_not_final);

    cout <<"Jet induced electron fake background estimation:\n" <<endl;    
    printLine("ee bkg", DYee_mm, DYmm_mm, DYtt_mm, tt_mm, wjets_mm, wz_mm, zz_mm, ww_mm, tw_mm, qcdShape);
    printLine("em bkg", DYee_em, DYmm_em, DYtt_em, tt_em, wjets_em, wz_em, zz_em, ww_em, tw_em, qcdShape);
    cout <<endl;
  } else {
    cout << "Fake rates are not available to estimate jet induced background. Please re-run with QCD samples" << endl;
  }
  */
}
