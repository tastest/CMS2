#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TSystem.h"
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
  
  // 
  // fakerate version 
  // "WW" or "fakerate"
  //
  string jetbkgversion = "WW";
  if (gSystem->Getenv("JETBKGVERSION")){
    jetbkgversion = gSystem->Getenv("JETBKGVERSION");
    cout << "JetBkgVersion: " << jetbkgversion << endl;
  }

  //hist file:
  TFile *ftt = TFile::Open(file);
  assert(ftt);
  
  TFile *el_fakeRateFile = TFile::Open("$CMS2_LOCATION/NtupleMacros/data/el_FR_3X.root", "read");
  assert(el_fakeRateFile);
    
  //do the emu row:
  string prefix = "";
  if(jetbkgversion == "fakerate") prefix = "_fakerate";
  TH2F *DYee_ee  = dynamic_cast<TH2F*>(ftt->Get(TString("dyee_helFRfakable"+prefix+"_ee")));
  TH2F *DYee_em  = dynamic_cast<TH2F*>(ftt->Get(TString("dyee_helFRfakable"+prefix+"_em")));
  TH2F *DYmm_ee  = dynamic_cast<TH2F*>(ftt->Get(TString("dymm_helFRfakable"+prefix+"_ee")));
  TH2F *DYmm_em  = dynamic_cast<TH2F*>(ftt->Get(TString("dymm_helFRfakable"+prefix+"_em")));
  TH2F *DYtt_ee  = dynamic_cast<TH2F*>(ftt->Get(TString("dytt_helFRfakable"+prefix+"_ee")));
  TH2F *DYtt_em  = dynamic_cast<TH2F*>(ftt->Get(TString("dytt_helFRfakable"+prefix+"_em")));
  TH2F *tt_ee    = dynamic_cast<TH2F*>(ftt->Get(TString("ttbar_helFRfakable"+prefix+"_ee")));
  TH2F *tt_em    = dynamic_cast<TH2F*>(ftt->Get(TString("ttbar_helFRfakable"+prefix+"_em")));
  TH2F *wjets_ee = dynamic_cast<TH2F*>(ftt->Get(TString("wjets_helFRfakable"+prefix+"_ee")));
  TH2F *wjets_em = dynamic_cast<TH2F*>(ftt->Get(TString("wjets_helFRfakable"+prefix+"_em")));
  TH2F *wz_ee    = dynamic_cast<TH2F*>(ftt->Get(TString("wz_helFRfakable"+prefix+"_ee")));
  TH2F *wz_em    = dynamic_cast<TH2F*>(ftt->Get(TString("wz_helFRfakable"+prefix+"_em")));
  TH2F *zz_ee    = dynamic_cast<TH2F*>(ftt->Get(TString("zz_helFRfakable"+prefix+"_ee")));
  TH2F *zz_em    = dynamic_cast<TH2F*>(ftt->Get(TString("zz_helFRfakable"+prefix+"_em")));
  TH2F *ww_ee    = dynamic_cast<TH2F*>(ftt->Get(TString("ww_helFRfakable"+prefix+"_ee")));
  TH2F *ww_em    = dynamic_cast<TH2F*>(ftt->Get(TString("ww_helFRfakable"+prefix+"_em")));
  TH2F *tw_ee    = dynamic_cast<TH2F*>(ftt->Get(TString("tw_helFRfakable"+prefix+"_ee")));
  TH2F *tw_em    = dynamic_cast<TH2F*>(ftt->Get(TString("tw_helFRfakable"+prefix+"_em")));
  
  cout <<"Electron fake count ("<<jetbkgversion<<") :" <<endl;
  cout << "\n" << Form("| %12s | %12s | %12s | %12s | %12s | %12s | %12s | %12s | %12s | %12s |",
		       "", "*DY ee*","*DY mumu*","*DY tautau*","*ttbar*","*Wjets*","*WZ*","*ZZ*","*WW*","*TW*")
       << endl;
  printLine("ee fakes", DYee_ee, DYmm_ee, DYtt_ee, tt_ee, wjets_ee, wz_ee, zz_ee, ww_ee, tw_ee, 0);
  printLine("em fakes", DYee_em, DYmm_em, DYtt_em, tt_em, wjets_em, wz_em, zz_em, ww_em, tw_em, 0);
  cout <<endl;

  //
  // Apply fake rates 
  //

  TH2F *fakeRate2;
  // If it is from WW code, calcualte from scratch the fake rate
  if(jetbkgversion=="WW") {
    TH2F *qcd_el_fake  = dynamic_cast<TH2F*>(ftt->Get("qcd_hFakableRateSingleElectron"));
    TH2F *qcd_el_final = dynamic_cast<TH2F*>(ftt->Get("qcd_hFinalRateSingleElectron"));
    if ( qcd_el_fake && qcd_el_final ){
      TH2F *qcd_el_fake_not_final = (TH2F*)qcd_el_fake->Clone("qcd_el_fake_not_final");
      
      // nominal fake rate
      //TH2F *fakeRate = (TH2F*)qcd_el_final->Clone("qcd_el_fakeRate");
      //fakeRate->Divide(qcd_el_fake);

      // fake_rate/(1-fake_rate)
      //fakeRate2 = (TH2F*)qcd_el_final->Clone("qcd_el_fakeRate2");
      //fakeRate2->Divide(qcd_el_fake_not_final);
   
      fakeRate2 = (TH2F*)qcd_el_final->Clone("fakeRate2");
      fakeRate2->Divide(qcd_el_final, qcd_el_fake, 1., 1., "B");
    
    } else {
    cout << "Fake rates are not available to estimate jet induced background. Please re-run with QCD samples" << endl;
    }
  }
  
  if(jetbkgversion=="fakerate") {
    fakeRate2     = dynamic_cast<TH2F *>( el_fakeRateFile->Get("QCD30_el_v2_cand02_FR_etavspt") );
    if(!fakeRate2) 
      cout << "Fake rates are not available from NtupleMacros/data to estimate jet induced background. Please check" << endl;
  }
  
  cout <<"Jet induced electron fake background estimation ("<<jetbkgversion<<") :\n" <<endl;    
  printLine("ee bkg", DYee_ee, DYmm_ee, DYtt_ee, tt_ee, wjets_ee, wz_ee, zz_ee, ww_ee, tw_ee, fakeRate2);
  printLine("em bkg", DYee_em, DYmm_em, DYtt_em, tt_em, wjets_em, wz_em, zz_em, ww_em, tw_em, fakeRate2);
  cout <<endl;
}
