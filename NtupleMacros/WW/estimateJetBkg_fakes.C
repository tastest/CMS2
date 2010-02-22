#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include <iostream>
#include <iomanip>
#include <string>
#include "wwtypes.cc"

void getIntegral(const TH2F* h, double& integral, double& error, const TH2F* fakeRate2){
  integral = 0;
  error = 0;
  if ( fakeRate2 == 0 ){
    for ( Int_t x=0; x<h->GetNbinsX(); ++x )
      for ( Int_t y=0; y<h->GetNbinsY(); ++y ){
	integral += h->GetBinContent(x+1,y+1);
	error += h->GetBinError(x+1,y+1)*h->GetBinError(x+1,y+1);
      }
  } else {
    for ( Int_t x=0; x<h->GetNbinsX(); ++x )
      for ( Int_t y=0; y<h->GetNbinsY(); ++y ){
	integral += h->GetBinContent(x+1,y+1)*fakeRate2->GetBinContent(x+1,y+1);
	error += h->GetBinError(x+1,y+1)*h->GetBinError(x+1,y+1)*fakeRate2->GetBinContent(x+1,y+1)*fakeRate2->GetBinContent(x+1,y+1)+
	  h->GetBinContent(x+1,y+1)*h->GetBinContent(x+1,y+1)*fakeRate2->GetBinError(x+1,y+1)*fakeRate2->GetBinError(x+1,y+1);
      }
  }
  error = sqrt(error);
}

void printLine( const char*name,
		TH2F* DYee, TH2F* DYmm, TH2F* DYtt, TH2F* tt, TH2F* wjets, TH2F* wz, TH2F* zz, TH2F* ww, TH2F* tw,
		TH2F* fakeRate2)
{
  string pm = "+/-";
  // string pm = "&plusmn;";
  
  cout << "|" << Form(" %12s ",name) << "|";
  double integral(0);
  double error(0);
  if (DYee) {
    getIntegral(DYee, integral, error, fakeRate2);
    cout << Form(" %5.2f%s%4.2f ", integral, pm.c_str(), error);
  } 
  else
    cout << "    skipped   ";
  cout << "|";
  if (DYmm) {
    getIntegral(DYmm, integral, error, fakeRate2);
    cout << Form(" %5.2f%s%4.2f ", integral, pm.c_str(), error);
  }
  else
    cout << "    skipped   ";
  cout << "|";
  if (DYtt) {
    getIntegral(DYtt, integral, error, fakeRate2);
    cout << Form(" %5.2f%s%4.2f ", integral, pm.c_str(), error);
  }
  else
    cout << "    skipped   ";
  cout << "|";
  if (tt) {
    getIntegral(tt, integral, error, fakeRate2);
    cout << Form(" %5.2f%s%4.2f ", integral, pm.c_str(), error);
  }
  else
    cout << "    skipped   ";
  cout << "|";
  if (wjets) {
    getIntegral(wjets, integral, error, fakeRate2);
    cout << Form(" %5.2f%s%4.2f ", integral, pm.c_str(), error);
  }
  else
    cout << "    skipped   ";
  cout << "|";
  if (wz) {
    getIntegral(wz, integral, error, fakeRate2);
    cout << Form(" %5.2f%s%4.2f ", integral, pm.c_str(), error);
  }
  else
    cout << "    skipped   ";
  cout << "|";
  if (zz) {
    getIntegral(zz, integral, error, fakeRate2);
    cout << Form(" %5.2f%s%4.2f ", integral, pm.c_str(), error);
  }
  else
    cout << "    skipped   ";
  cout << "|";
  if (ww) {
    getIntegral(ww, integral, error, fakeRate2);
    cout << Form(" %5.2f%s%4.2f ", integral, pm.c_str(), error);
  }
  else
    cout << "    skipped   ";
  cout << "|";
  if (tw) {
    getIntegral(tw, integral, error, fakeRate2);
    cout << Form(" %5.2f%s%4.2f ", integral, pm.c_str(), error);
  }
  else
    cout << "    skipped   ";
  cout << "|" <<endl;
  
}
void estimateJetBkg_fakes(const char* file = "processed_data.root")
{
  using namespace std;
  
  //hist file:
  TFile *ftt = TFile::Open(file);
  assert(ftt);

  //do the emu row:
  TH2F *DYee_ee  = dynamic_cast<TH2F*>(ftt->Get("dyee_helFRfakable_ee"));
  TH2F *DYee_em  = dynamic_cast<TH2F*>(ftt->Get("dyee_helFRfakable_em"));
  TH2F *DYmm_ee  = dynamic_cast<TH2F*>(ftt->Get("dymm_helFRfakable_ee"));
  TH2F *DYmm_em  = dynamic_cast<TH2F*>(ftt->Get("dymm_helFRfakable_em"));
  TH2F *DYtt_ee  = dynamic_cast<TH2F*>(ftt->Get("dytt_helFRfakable_ee"));
  TH2F *DYtt_em  = dynamic_cast<TH2F*>(ftt->Get("dytt_helFRfakable_em"));
  TH2F *tt_ee    = dynamic_cast<TH2F*>(ftt->Get("ttbar_helFRfakable_ee"));
  TH2F *tt_em    = dynamic_cast<TH2F*>(ftt->Get("ttbar_helFRfakable_em"));
  TH2F *wjets_ee = dynamic_cast<TH2F*>(ftt->Get("wjets_helFRfakable_ee"));
  TH2F *wjets_em = dynamic_cast<TH2F*>(ftt->Get("wjets_helFRfakable_em"));
  TH2F *wz_ee    = dynamic_cast<TH2F*>(ftt->Get("wz_helFRfakable_ee"));
  TH2F *wz_em    = dynamic_cast<TH2F*>(ftt->Get("wz_helFRfakable_em"));
  TH2F *zz_ee    = dynamic_cast<TH2F*>(ftt->Get("zz_helFRfakable_ee"));
  TH2F *zz_em    = dynamic_cast<TH2F*>(ftt->Get("zz_helFRfakable_em"));
  TH2F *ww_ee    = dynamic_cast<TH2F*>(ftt->Get("ww_helFRfakable_ee"));
  TH2F *ww_em    = dynamic_cast<TH2F*>(ftt->Get("ww_helFRfakable_em"));
  TH2F *tw_ee    = dynamic_cast<TH2F*>(ftt->Get("tw_helFRfakable_ee"));
  TH2F *tw_em    = dynamic_cast<TH2F*>(ftt->Get("tw_helFRfakable_em"));
  
  cout <<"Electron fake count:" <<endl;
  cout << "\n" << Form("| %12s | %12s | %12s | %12s | %12s | %12s | %12s | %12s | %12s | %12s |",
		       "", "*DY ee*","*DY mumu*","*DY tautau*","*ttbar*","*Wjets*","*WZ*","*ZZ*","*WW*","*TW*")
       << endl;
  printLine("ee fakes", DYee_ee, DYmm_ee, DYtt_ee, tt_ee, wjets_ee, wz_ee, zz_ee, ww_ee, tw_ee, 0);
  printLine("em fakes", DYee_em, DYmm_em, DYtt_em, tt_em, wjets_em, wz_em, zz_em, ww_em, tw_em, 0);
  cout <<endl;

  //
  // Apply fake rates
  //
  TH2F *qcd_el_fake  = dynamic_cast<TH2F*>(ftt->Get("qcd_hFakableRateSingleElectron"));
  TH2F *qcd_el_final = dynamic_cast<TH2F*>(ftt->Get("qcd_hFinalRateSingleElectron"));
  if ( qcd_el_fake && qcd_el_final ){
    TH2F *qcd_el_fake_not_final = (TH2F*)qcd_el_fake->Clone("qcd_el_fake_not_final");
  
    // nominal fake rate
    TH2F *fakeRate = (TH2F*)qcd_el_final->Clone("qcd_el_fakeRate");
    fakeRate->Divide(qcd_el_fake);

    // fake_rate/(1-fake_rate)
    TH2F *fakeRate2 = (TH2F*)qcd_el_final->Clone("qcd_el_fakeRate2");
    fakeRate2->Divide(qcd_el_fake_not_final);

    cout <<"Jet induced electron fake background estimation:\n" <<endl;    
    printLine("ee bkg", DYee_ee, DYmm_ee, DYtt_ee, tt_ee, wjets_ee, wz_ee, zz_ee, ww_ee, tw_ee, fakeRate2);
    printLine("em bkg", DYee_em, DYmm_em, DYtt_em, tt_em, wjets_em, wz_em, zz_em, ww_em, tw_em, fakeRate2);
    cout <<endl;
  } else {
    cout << "Fake rates are not available to estimate jet induced background. Please re-run with QCD samples" << endl;
  }
  
}
