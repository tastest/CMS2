// 
// It is assumed that fake inputs are calculated as number of events
// when one of the lepton fails nominal selection, but satisfy the fakable
// and the othe lepton satisfies the nominal selection
//
// The fake rate ratio is assumed to be a straight ratio of number of
// leptons that passed the nominal selection devided by number of those
// that passed the fakable object requirements.
// 

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TSystem.h"
#include <iostream>
#include <iomanip>
#include <string>
#include "wwtypes.cc"

bool CheckBinLimits(const TArrayD* h1Array, const TArrayD* h2Array)
{
  Int_t fN = h1Array->fN;
  if ( fN != 0 ) {
    if ( h2Array->fN != fN ) {
      printf("DifferentBinLimits\n");
      return false;
    }
    else {
      for ( int i = 0; i < fN; ++i ) {
	if ( fabs(h1Array->GetAt(i)-h2Array->GetAt(i)) > 1E-3 ) {
	  printf("DifferentBinLimits\n");
	  return false;
	}
      }
    }
  }

  return true;
}

bool CheckConsistency(const TH1* h1, const TH1* h2)
{
  // Check histogram compatibility   Int_t nbinsx = h1->GetNbinsX();
  Int_t nbinsx = h1->GetNbinsX();
  Int_t nbinsy = h1->GetNbinsY();
  Int_t nbinsz = h1->GetNbinsZ();

  // Check whether the histograms have the same number of bins.
  if (nbinsx != h2->GetNbinsX() || nbinsy != h2->GetNbinsY() || nbinsz != h2->GetNbinsZ()) {
    printf("DifferentNumberOfBins\n");
    return false;
  }
  // Check that the axis limits of the histograms are the same
  if (h1->fXaxis.GetXmin() != h2->fXaxis.GetXmin() ||
      h1->fXaxis.GetXmax() != h2->fXaxis.GetXmax() ||
      h1->fYaxis.GetXmin() != h2->fYaxis.GetXmin() ||
      h1->fYaxis.GetXmax() != h2->fYaxis.GetXmax() ||
      h1->fZaxis.GetXmin() != h2->fZaxis.GetXmin() ||
      h1->fZaxis.GetXmax() != h2->fZaxis.GetXmax()) {
    printf("DifferentAxisLimits\n");
    return false;
  }

  bool ret = true;

  ret &= CheckBinLimits(h1->GetXaxis()->GetXbins(), h2->GetXaxis()->GetXbins());
  ret &= CheckBinLimits(h1->GetYaxis()->GetXbins(), h2->GetYaxis()->GetXbins());
  ret &= CheckBinLimits(h1->GetZaxis()->GetXbins(), h2->GetZaxis()->GetXbins());

  return ret;
}

//___________________________________________________________________________

void getIntegral(const TH2F* h, double& integral, double& error, const TH2F* fakeRate){
  integral = 0;
  error = 0;

  // first check if histogram have the same binning
  if ( fakeRate && ! CheckConsistency(h,fakeRate) ){
    printf("Error: histograms are not consistent. Abort");
    return;
  }

  if ( fakeRate == 0 ){
    for ( Int_t x=0; x<h->GetNbinsX(); ++x )
      for ( Int_t y=0; y<h->GetNbinsY(); ++y ){
	integral += h->GetBinContent(x+1,y+1);
	error += h->GetBinError(x+1,y+1)*h->GetBinError(x+1,y+1);
      }
  } else {
    for ( Int_t x=0; x<h->GetNbinsX(); ++x )
      for ( Int_t y=0; y<h->GetNbinsY(); ++y ){
	double bin    = h->GetBinContent(x+1,y+1);
	double binerr = h->GetBinError(x+1,y+1);
	double fr     = fakeRate->GetBinContent(x+1,y+1);
	double frerr  = fakeRate->GetBinError(x+1,y+1);
	integral += bin*fr/(1-fr);
	error += pow(binerr*fr/(1-fr),2) + pow(bin*frerr/(1-fr)/(1-fr), 2);
      }
  }
  error = sqrt(error);
}

void printLine( const char*name,
		TH2F* DYee, TH2F* DYmm, TH2F* DYtt, TH2F* tt, TH2F* wjets, TH2F* wz, TH2F* zz, TH2F* ww, TH2F* tw, TH2F* data,
		TH2F* fakeRate)
{
  string pm = "+/-";
  // string pm = "&plusmn;";
  
  cout << "|" << Form(" %12s ",name) << "|";
  double integral(0);
  double error(0);
  if (DYee) {
    getIntegral(DYee, integral, error, fakeRate);
    cout << Form(" %5.2f%s%4.2f ", integral, pm.c_str(), error);
  } 
  else
    cout << "    skipped   ";
  cout << "|";
  if (DYmm) {
    getIntegral(DYmm, integral, error, fakeRate);
    cout << Form(" %5.2f%s%4.2f ", integral, pm.c_str(), error);
  }
  else
    cout << "    skipped   ";
  cout << "|";
  if (DYtt) {
    getIntegral(DYtt, integral, error, fakeRate);
    cout << Form(" %5.2f%s%4.2f ", integral, pm.c_str(), error);
  }
  else
    cout << "    skipped   ";
  cout << "|";
  if (tt) {
    getIntegral(tt, integral, error, fakeRate);
    cout << Form(" %5.2f%s%4.2f ", integral, pm.c_str(), error);
  }
  else
    cout << "    skipped   ";
  cout << "|";
  if (wjets) {
    getIntegral(wjets, integral, error, fakeRate);
    cout << Form(" %5.2f%s%4.2f ", integral, pm.c_str(), error);
  }
  else
    cout << "    skipped   ";
  cout << "|";
  if (wz) {
    getIntegral(wz, integral, error, fakeRate);
    cout << Form(" %5.2f%s%4.2f ", integral, pm.c_str(), error);
  }
  else
    cout << "    skipped   ";
  cout << "|";
  if (zz) {
    getIntegral(zz, integral, error, fakeRate);
    cout << Form(" %5.2f%s%4.2f ", integral, pm.c_str(), error);
  }
  else
    cout << "    skipped   ";
  cout << "|";
  if (ww) {
    getIntegral(ww, integral, error, fakeRate);
    cout << Form(" %5.2f%s%4.2f ", integral, pm.c_str(), error);
  }
  else
    cout << "    skipped   ";
  cout << "|";
  if (tw) {
    getIntegral(tw, integral, error, fakeRate);
    cout << Form(" %5.2f%s%4.2f ", integral, pm.c_str(), error);
  }
  else
    cout << "    skipped   ";
  cout << "|";
  if (data) {
    getIntegral(data, integral, error, fakeRate);
    cout << Form(" %5.2f%s%4.2f ", integral, pm.c_str(), error);
  }
  else
    cout << "    skipped   ";
  cout << "|" <<endl;
  
}

void estimateJetBkg_fakes(const char* file = "processed_data.root", 
			  const char* el_fakerate_file = 0,
			  const char* el_fakerate_name = 0,
			  const char* mu_fakerate_file = 0,
			  const char* mu_fakerate_name = 0)
{
  using namespace std;
  
  // 
  // fakerate version 
  // "WW" or "fakerate"
  //
  bool externalElFakeRate = (el_fakerate_file != 0) && (el_fakerate_name != 0);
  bool externalMuFakeRate = (mu_fakerate_file != 0) && (mu_fakerate_name != 0);

  //hist file:
  TFile *ftt = TFile::Open(file);
  assert(ftt);
  
  string prefix = "";
  if(externalElFakeRate) prefix = "_fakerate";
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
  TH2F *data_ee  = dynamic_cast<TH2F*>(ftt->Get(TString("data_helFRfakable"+prefix+"_ee")));
  TH2F *data_em  = dynamic_cast<TH2F*>(ftt->Get(TString("data_helFRfakable"+prefix+"_em")));

  if( externalElFakeRate )
    cout <<"Electron fake count for external fakerates:" <<endl;
  else
    cout <<"Electron fake count for ad-hoc fakerates:" <<endl;

  cout << "\n" << Form("| %12s | %12s | %12s | %12s | %12s | %12s | %12s | %12s | %12s | %12s | %12s |",
		       "", "*DY ee*","*DY mumu*","*DY tautau*","*ttbar*","*Wjets*","*WZ*","*ZZ*","*WW*","*TW*","*Data*")
       << endl;
  printLine("ee fakes", DYee_ee, DYmm_ee, DYtt_ee, tt_ee, wjets_ee, wz_ee, zz_ee, ww_ee, tw_ee, data_ee, 0);
  printLine("em fakes", DYee_em, DYmm_em, DYtt_em, tt_em, wjets_em, wz_em, zz_em, ww_em, tw_em, data_em, 0);
  cout <<endl;

  //
  // Apply fake rates 
  //

  TH2F *fakeRate(0);
  // If it is from WW code, calcualte from scratch the fake rate
  if( ! externalElFakeRate ){
    TH2F *qcd_el_fake  = dynamic_cast<TH2F*>(ftt->Get("qcd_hFakableRateSingleElectron"));
    TH2F *qcd_el_final = dynamic_cast<TH2F*>(ftt->Get("qcd_hFinalRateSingleElectron"));
    if ( qcd_el_fake && qcd_el_final ){
      //fakeRate = (TH2F*)qcd_el_final->Clone("fakeRate");
      //fakeRate->Divide(qcd_el_fake);
      fakeRate = (TH2F*)qcd_el_final->Clone("fakeRate");
      fakeRate->Divide(qcd_el_final, qcd_el_fake, 1., 1., "B");
    } else {
      cout << "Fake rates are not available to estimate jet induced background. Please re-run with QCD samples" << endl;
    }
    cout <<"Jet induced electron fake background estimation (ad-hoc) :\n" <<endl;    
  } else {
    TFile *el_fakeRateFile = TFile::Open(el_fakerate_file);
    // TFile *el_fakeRateFile = TFile::Open("../data/el_FR_3X.root");
    assert(el_fakeRateFile);
    // fakeRate = dynamic_cast<TH2F*>( el_fakeRateFile->Get("QCD30_el_v2_cand02_FR_etavspt") );
    fakeRate = dynamic_cast<TH2F*>( el_fakeRateFile->Get(el_fakerate_name) );
    if(!fakeRate) 
      cout << "Fake rates are not available from NtupleMacros/data to estimate jet induced background. Please check" << endl;
    cout <<"Jet induced electron fake background estimation (external fake rate) :\n" <<endl;    
  }
  printLine("ee bkg", DYee_ee, DYmm_ee, DYtt_ee, tt_ee, wjets_ee, wz_ee, zz_ee, ww_ee, tw_ee, data_ee, fakeRate);
  printLine("em bkg", DYee_em, DYmm_em, DYtt_em, tt_em, wjets_em, wz_em, zz_em, ww_em, tw_em, data_em, fakeRate);
  cout <<endl;

  //
  // Now the muons
  //

  prefix = "";
  if(externalMuFakeRate) prefix = "_fakerate";
  TH2F *DYee_mm  = dynamic_cast<TH2F*>(ftt->Get(TString("dyee_hmuFRfakable"+prefix+"_mm")));
  TH2F *DYee_me  = dynamic_cast<TH2F*>(ftt->Get(TString("dyee_hmuFRfakable"+prefix+"_em")));
  TH2F *DYmm_mm  = dynamic_cast<TH2F*>(ftt->Get(TString("dymm_hmuFRfakable"+prefix+"_mm")));
  TH2F *DYmm_me  = dynamic_cast<TH2F*>(ftt->Get(TString("dymm_hmuFRfakable"+prefix+"_em")));
  TH2F *DYtt_mm  = dynamic_cast<TH2F*>(ftt->Get(TString("dytt_hmuFRfakable"+prefix+"_mm")));
  TH2F *DYtt_me  = dynamic_cast<TH2F*>(ftt->Get(TString("dytt_hmuFRfakable"+prefix+"_em")));
  TH2F *tt_mm    = dynamic_cast<TH2F*>(ftt->Get(TString("ttbar_hmuFRfakable"+prefix+"_mm")));
  TH2F *tt_me    = dynamic_cast<TH2F*>(ftt->Get(TString("ttbar_hmuFRfakable"+prefix+"_em")));
  TH2F *wjets_mm = dynamic_cast<TH2F*>(ftt->Get(TString("wjets_hmuFRfakable"+prefix+"_mm")));
  TH2F *wjets_me = dynamic_cast<TH2F*>(ftt->Get(TString("wjets_hmuFRfakable"+prefix+"_em")));
  TH2F *wz_mm    = dynamic_cast<TH2F*>(ftt->Get(TString("wz_hmuFRfakable"+prefix+"_mm")));
  TH2F *wz_me    = dynamic_cast<TH2F*>(ftt->Get(TString("wz_hmuFRfakable"+prefix+"_em")));
  TH2F *zz_mm    = dynamic_cast<TH2F*>(ftt->Get(TString("zz_hmuFRfakable"+prefix+"_mm")));
  TH2F *zz_me    = dynamic_cast<TH2F*>(ftt->Get(TString("zz_hmuFRfakable"+prefix+"_em")));
  TH2F *ww_mm    = dynamic_cast<TH2F*>(ftt->Get(TString("ww_hmuFRfakable"+prefix+"_mm")));
  TH2F *ww_me    = dynamic_cast<TH2F*>(ftt->Get(TString("ww_hmuFRfakable"+prefix+"_em")));
  TH2F *tw_mm    = dynamic_cast<TH2F*>(ftt->Get(TString("tw_hmuFRfakable"+prefix+"_mm")));
  TH2F *tw_me    = dynamic_cast<TH2F*>(ftt->Get(TString("tw_hmuFRfakable"+prefix+"_em")));
  TH2F *data_mm  = dynamic_cast<TH2F*>(ftt->Get(TString("data_hmuFRfakable"+prefix+"_mm")));
  TH2F *data_me  = dynamic_cast<TH2F*>(ftt->Get(TString("data_hmuFRfakable"+prefix+"_em")));

  if( externalMuFakeRate )
    cout <<"Muon fake count for external fakerates:" <<endl;
  else
    cout <<"Muon fake count for ad-hoc fakerates:" <<endl;

  cout << "\n" << Form("| %12s | %12s | %12s | %12s | %12s | %12s | %12s | %12s | %12s | %12s | %12s |",
		       "", "*DY ee*","*DY mumu*","*DY tautau*","*ttbar*","*Wjets*","*WZ*","*ZZ*","*WW*","*TW*","*Data*")
       << endl;
  printLine("mm fakes", DYee_mm, DYmm_mm, DYtt_mm, tt_mm, wjets_mm, wz_mm, zz_mm, ww_mm, tw_mm, data_mm, 0);
  printLine("em fakes", DYee_me, DYmm_me, DYtt_me, tt_me, wjets_me, wz_me, zz_me, ww_me, tw_me, data_me, 0);
  cout <<endl;

  //
  // Apply fake rates 
  //

  fakeRate = 0;
  // If it is from WW code, calculate from scratch the fake rate
  if( ! externalMuFakeRate ){
    TH2F *qcd_mu_fake  = dynamic_cast<TH2F*>(ftt->Get("qcd_hFakableRateSingleMuon"));
    TH2F *qcd_mu_final = dynamic_cast<TH2F*>(ftt->Get("qcd_hFinalRateSingleMuon"));
    if ( qcd_mu_fake && qcd_mu_final ){
      //fakeRate = (TH2F*)qcd_mu_final->Clone("fakeRate");
      //fakeRate->Divide(qcd_mu_fake);
      fakeRate = (TH2F*)qcd_mu_final->Clone("fakeRate");
      fakeRate->Divide(qcd_mu_final, qcd_mu_fake, 1., 1., "B");
    } else {
      cout << "Fake rates are not available to estimate jet induced background. Please re-run with QCD samples" << endl;
    }
    cout <<"Jet induced muon fake background estimation (ad-hoc) :\n" <<endl;    
  } else {
    TFile *mu_fakeRateFile = TFile::Open(mu_fakerate_file);
    assert(mu_fakeRateFile);
    fakeRate = dynamic_cast<TH2F*>( mu_fakeRateFile->Get(mu_fakerate_name) );
    if(!fakeRate) 
      cout << "Fake rates are not available from NtupleMacros/data to estimate jet induced background. Please check" << endl;
    cout <<"Jet induced muon fake background estimation (external fake rate) :\n" <<endl;    
  }
  printLine("mm bkg", DYee_mm, DYmm_mm, DYtt_mm, tt_mm, wjets_mm, wz_mm, zz_mm, ww_mm, tw_mm, data_mm, fakeRate);
  printLine("em bkg", DYee_me, DYmm_me, DYtt_me, tt_me, wjets_me, wz_me, zz_me, ww_me, tw_me, data_me, fakeRate);
  cout <<endl;
}
