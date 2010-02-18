#include "TSystem.h"
#include "TFile.h"
#include "TH2.h"
#include "fakerates.h"
#include "CMS2.h"
#include "electronSelections.h"
#include "muonSelections.h"

/*muons*/
class TH2F &fakeRateMuon (enum fakeRateVersion);
class TH2F &fakeRateErrorMuon (enum fakeRateVersion);
bool   isFakeDenominatorMuon_v1 (int);
double muFakeProb_v1 (int);
double muFakeProbErr_v1 (int);
static TFile *mu_fakeRateFile = 0;
static TH2F  *mu_fakeRate = 0;
static TH2F  *mu_fakeRate_err = 0;
bool isFakeableMuon (int i_mu, enum fakeRateVersion version){
  if(version == mu_v1) return  isFakeDenominatorMuon_v1(i_mu);
  else {
    std::cout<<"isFakeableMuon: invalid fakeRateVersion given. Check it!"<<std::endl;
    return false;
  }
}
double muFakeProb (int i_mu, enum fakeRateVersion version){
  if(version == mu_v1) return muFakeProb_v1(i_mu);
  else {
    std::cout<<"muFakeProb: invalid muon fakeRateVersion given. Check it!"<<std::endl;
    return -999.;
  }
}
double muFakeProbErr (int i_mu, enum fakeRateVersion version){
  if(version == mu_v1) return muFakeProbErr_v1(i_mu);
  else {
    std::cout<<"muFakeProbErr: invalid muon fakeRateVersion given. Check it!"<<std::endl;
    return -999.;
  }
}

double muFakeProb_v1 (int i_mu){
  float prob = 0.0;
  float prob_error = 0.0;
  TH2F *theFakeRate = &fakeRateMuon(mu_v1);
  TH2F *theFakeRateErr = &fakeRateErrorMuon(mu_v1);
  // cut definition
  float pt = cms2.mus_p4()[i_mu].Pt();
  float eta = fabs(cms2.mus_p4()[i_mu].Eta());
  float upperEdge = theFakeRate->GetYaxis()->GetBinLowEdge(theFakeRate->GetYaxis()->GetNbins()) + theFakeRate->GetYaxis()->GetBinWidth(theFakeRate->GetYaxis()->GetNbins()) - 0.001;

  if ( pt > upperEdge )
    pt = upperEdge;
  prob = theFakeRate->GetBinContent(theFakeRate->FindBin(eta,pt));
  prob_error = theFakeRateErr->GetBinContent(theFakeRateErr->FindBin(eta,pt));

  if (prob>1.0 || prob<0.0) {
    std::cout<<"ERROR FROM MU FAKE RATE!!! prob = " << prob << std::endl;
  }
  if (prob==0.0){
    std::cout<<"ERROR FROM MU FAKE RATE!!! prob = " << prob
             <<" for Et = " <<cms2.mus_p4()[i_mu].Pt()
             <<" and Eta = " <<cms2.mus_p4()[i_mu].Eta()
             << std::endl;
  }
  return prob;
}

double muFakeProbErr_v1 (int i_mu){
  float prob = 0.0;
  float prob_error = 0.0;
  TH2F *theFakeRate = &fakeRateMuon(mu_v1);
  TH2F *theFakeRateErr = &fakeRateErrorMuon(mu_v1);
  // cut definition
  float pt = cms2.mus_p4()[i_mu].Pt();
  float eta = fabs(cms2.mus_p4()[i_mu].Eta());
  float upperEdge = theFakeRate->GetYaxis()->GetBinLowEdge(theFakeRate->GetYaxis()->GetNbins()) + theFakeRate->GetYaxis()->GetBinWidth(theFakeRate->GetYaxis()->GetNbins()) - 0.001;

  if ( pt > upperEdge )
    pt = upperEdge;
  prob = theFakeRate->GetBinContent(theFakeRate->FindBin(eta,pt));
  prob_error = theFakeRateErr->GetBinContent(theFakeRateErr->FindBin(eta,pt));

  if (prob>1.0 || prob<0.0) {
    std::cout<<"ERROR FROM MU FAKE RATE!!! prob = " << prob << std::endl;
  }
  if (prob==0.0){
    std::cout<<"ERROR FROM MU FAKE RATE!!! prob = " << prob
             <<" for Et = " <<cms2.mus_p4()[i_mu].Pt()
             <<" and Eta = " <<cms2.mus_p4()[i_mu].Eta()
             << std::endl;
  }
  return prob_error;
}
bool isFakeDenominatorMuon_v1 (int index) {
    if ( cms2.mus_p4()[index].pt() < 10.) {
      std::cout << "muonID ERROR: requested muon is too low pt,  Abort." << std::endl;
      return false;
    }
    if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.5)  return false; // eta cut
    //if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) return false; //glb fit chisq
    if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 50) return false;
    if (((cms2.mus_type().at(index)) & (1<<1)) == 0)    return false; // global muon
    if (((cms2.mus_type().at(index)) & (1<<2)) == 0)    return false; // tracker muon
    if (cms2.mus_validHits().at(index) < 11)            return false; // # of tracker hits
    //if (cms2.mus_iso_ecalvetoDep().at(index) > 4)       return false; // ECalE < 4 
    //if (cms2.mus_iso_hcalvetoDep().at(index) > 6)       return false; // HCalE < 6 
    if (cms2.mus_iso_ecalvetoDep().at(index) > 10)       return false;
    if (cms2.mus_iso_hcalvetoDep().at(index) > 12)       return false;
    if (cms2.mus_gfit_validSTAHits().at(index) == 0)    return false; // Glb fit must have hits in mu chambers
    //if (TMath::Abs(cms2.mus_d0corr().at(index)) > 0.02) return false; // d0 from beamspot
    //if (muonIsoValue(index) > 0.1)                      return false; // Isolation cut
    if (TMath::Abs(cms2.mus_d0corr().at(index)) > 0.2) return false; // d0 from beamspot
    if (muonIsoValue(index) > 0.4)                      return false; // Isolation cut
    return true;

}
TH2F &fakeRateMuon (enum fakeRateVersion version){
  if ( mu_fakeRateFile == 0 ) {
    mu_fakeRateFile = TFile::Open("$CMS2_LOCATION/NtupleMacros/data/mu_FR_3X.root", "read");
    if ( mu_fakeRateFile == 0 ) {
      std::cout << "$CMS2_LOCATION/NtupleMacros/data/QCDFRplots-v1.root could not be found!!" << std::endl;
      std::cout << "Please make sure that $CMS2_LOCATION points to your CMS2 directory and that" << std::endl;
      std::cout << "$CMS2_LOCATION/NtupleMacros/data/QCDFRplots-v1.root exists!" << std::endl;
      gSystem->Exit(1);
    }
    mu_fakeRate       = dynamic_cast<TH2F *>( mu_fakeRateFile->Get("QCD30_mu_FR_etavspt") );
    mu_fakeRate_err   = dynamic_cast<TH2F *>( mu_fakeRateFile->Get("QCD30_mu_FRErr_etavspt") );
  }
  return *mu_fakeRate;
}
TH2F &fakeRateErrorMuon (enum fakeRateVersion version){
  if ( mu_fakeRateFile == 0 ) {
    mu_fakeRateFile = TFile::Open("$CMS2_LOCATION/NtupleMacros/data/mu_FR_3X.root", "read");
    if ( mu_fakeRateFile == 0 ) {
      std::cout << "$CMS2_LOCATION/NtupleMacros/data/QCDFRplots-v1.root could not be found!!" << std::endl;
      std::cout << "Please make sure that $CMS2_LOCATION points to your CMS2 directory and that" << std::endl;
      std::cout << "$CMS2_LOCATION/NtupleMacros/data/QCDFRplots-v1.root exists!" << std::endl;
      gSystem->Exit(1);
    }
    mu_fakeRate_err = dynamic_cast<TH2F *>(mu_fakeRateFile->Get("QCD30_mu_FRErr_etavspt"));
  } 
  return *mu_fakeRate_err;
}



/* electrons */
class TH2F &fakeRateEl (enum fakeRateVersion);
class TH2F &fakeRateErrorEl (enum fakeRateVersion);
bool   isFakeDenominatorElectron_v1 (int);
double elFakeProb_v1 (int);
double elFakeProbErr_v1 (int);
bool   isFakeDenominatorElectron_v2 (int);
double elFakeProb_v2 (int);
double elFakeProbErr_v2 (int);
bool   isFakeDenominatorElectron_v3 (int);
double elFakeProb_v3 (int);
double elFakeProbErr_v3 (int);
static TFile *el_fakeRateFile = 0;
static TH2F  *el_fakeRate_v1 = 0;
static TH2F  *el_fakeRate_v2 = 0;
static TH2F  *el_fakeRate_v3 = 0;
static TFile *el_fakeRateErrorFile = 0;
static TH2F  *el_fakeRateErr_v1 = 0;
static TH2F  *el_fakeRateErr_v2 = 0;
static TH2F  *el_fakeRateErr_v3 = 0;

bool isFakeableElectron (int i_el, enum fakeRateVersion version)
{
  if(version == el_v1) return isFakeDenominatorElectron_v1(i_el);
  if(version == el_v2) return isFakeDenominatorElectron_v2(i_el);
  if(version == el_v3) return isFakeDenominatorElectron_v3(i_el);
  else {
    std::cout<<"isFakeable: invalid fakeRateVersion given. Check it!"<<std::endl;
    return false;
  }
}
double elFakeProb (int i_el, enum fakeRateVersion version)
{
  if(version == el_v1) return elFakeProb_v1(i_el);
  if(version == el_v2) return elFakeProb_v2(i_el);
  if(version == el_v3) return elFakeProb_v3(i_el);
  else {
    std::cout<<"elFakeProb: invalid fakeRateVersion given. Check it!"<<std::endl;
    return -999.;
  }
}

double elFakeProbErr (int i_el, enum fakeRateVersion version)
{
  if(version == el_v1) return elFakeProbErr_v1(i_el);
  if(version == el_v2) return elFakeProbErr_v2(i_el);
  if(version == el_v3) return elFakeProbErr_v3(i_el);
  else {
    std::cout<<"elFakeProbErr: invalid fakeRateVersion given. Check it!"<<std::endl;
    return -999.;
  }
}

// --------------------------------
// electrons v1
// --------------------------------
bool isFakeDenominatorElectron_v1 (int index) {
 	if (!cms2.els_type()[index] & (1<<ISECALDRIVEN)) return false;
	if (fabs(cms2.els_p4()[index].eta()) > 2.5) return false;
  if (!electronId_noMuon(index)) return false;
  //if (!electronId_cand01(index)) return false;
  //if (!electronImpact_cand01(index)) return false;
  //if (electronIsolation_relsusy_cand1(index, true) > 0.10) return false;
  if (electronIsolation_relsusy_cand1(index, true) > 0.40) return false;
  if (isFromConversionPartnerTrack(index)) return false;
  return true;
}
double elFakeProb_v1 (int i_el){
  float prob = 0.0;
  float prob_error = 0.0;
  TH2F *theFakeRate = &fakeRateEl(el_v1);
  TH2F *theFakeRateErr = &fakeRateErrorEl(el_v1);
  // cut definition
  float pt = cms2.els_p4()[i_el].Pt();
  float eta = fabs(cms2.els_p4()[i_el].Eta());
  float upperEdge = theFakeRate->GetYaxis()->GetBinLowEdge(theFakeRate->GetYaxis()->GetNbins()) + theFakeRate->GetYaxis()->GetBinWidth(theFakeRate->GetYaxis()->GetNbins()) - 0.001;

  if ( pt > upperEdge ) pt = upperEdge;
     
  prob       = theFakeRate->GetBinContent(theFakeRate->FindBin(eta,pt));
  prob_error = theFakeRateErr->GetBinContent(theFakeRateErr->FindBin(eta,pt));
     
  if (prob>1.0 || prob<0.0) {
	  std::cout<<"ERROR FROM EL FAKE RATE!!! prob = " << prob << std::endl;
  }
  if (prob==0.0){
	  std::cout<<"ERROR FROM EL FAKE RATE!!! prob = " << prob
             <<" for Et = " <<cms2.els_p4()[i_el].Pt()
             <<" and Eta = " <<cms2.els_p4()[i_el].Eta()
             << std::endl;
  }
  return prob;
}

double elFakeProbErr_v1 (int i_el){
  float prob = 0.0;
  float prob_error = 0.0;
  TH2F *theFakeRate = &fakeRateEl(el_v1);
  TH2F *theFakeRateErr = &fakeRateErrorEl(el_v1);
  // cut definition
  float pt = cms2.els_p4()[i_el].Pt();
  float eta = fabs(cms2.els_p4()[i_el].Eta());
  float upperEdge = theFakeRate->GetYaxis()->GetBinLowEdge(theFakeRate->GetYaxis()->GetNbins()) + theFakeRate->GetYaxis()->GetBinWidth(theFakeRate->GetYaxis()->GetNbins()) - 0.001;

  if ( pt > upperEdge ) pt = upperEdge;
     
  prob       = theFakeRate->GetBinContent(theFakeRate->FindBin(eta,pt));
  prob_error = theFakeRateErr->GetBinContent(theFakeRateErr->FindBin(eta,pt));
     
  if (prob>1.0 || prob<0.0) {
	  std::cout<<"ERROR FROM EL FAKE RATE!!! prob = " << prob << std::endl;
  }
  if (prob==0.0){
	  std::cout<<"ERROR FROM EL FAKE RATE!!! prob = " << prob
             <<" for Et = " <<cms2.els_p4()[i_el].Pt()
             <<" and Eta = " <<cms2.els_p4()[i_el].Eta()
             << std::endl;
  }
  return prob_error;
}

// --------------------------------
// electrons v2
// --------------------------------
bool isFakeDenominatorElectron_v2 (int index) {
 	if (!cms2.els_type()[index] & (1<<ISECALDRIVEN)) return false;
	if (fabs(cms2.els_p4()[index].eta()) > 2.5) return false;
  if (!electronId_noMuon(index)) return false;
  //if (!electronId_cand01(index)) return false;
  //if (!electronImpact_cand01(index)) return false;
  if (electronIsolation_relsusy_cand1(index, true) > 0.10) return false;
  if (isFromConversionPartnerTrack(index)) return false;
  return true;
}
double elFakeProb_v2 (int i_el){
  float prob = 0.0;
  float prob_error = 0.0;
  TH2F *theFakeRate = &fakeRateEl(el_v2);
  TH2F *theFakeRateErr = &fakeRateErrorEl(el_v2);
  // cut definition
  float pt = cms2.els_p4()[i_el].Pt();
  float eta = fabs(cms2.els_p4()[i_el].Eta());
  float upperEdge = theFakeRate->GetYaxis()->GetBinLowEdge(theFakeRate->GetYaxis()->GetNbins()) + theFakeRate->GetYaxis()->GetBinWidth(theFakeRate->GetYaxis()->GetNbins()) - 0.001;

  if ( pt > upperEdge ) pt = upperEdge;
     
  prob       = theFakeRate->GetBinContent(theFakeRate->FindBin(eta,pt));
  prob_error = theFakeRateErr->GetBinContent(theFakeRateErr->FindBin(eta,pt));
     
  if (prob>1.0 || prob<0.0) {
	  std::cout<<"ERROR FROM EL FAKE RATE!!! prob = " << prob << std::endl;
  }
  if (prob==0.0){
	  std::cout<<"ERROR FROM EL FAKE RATE!!! prob = " << prob
             <<" for Et = " <<cms2.els_p4()[i_el].Pt()
             <<" and Eta = " <<cms2.els_p4()[i_el].Eta()
             << std::endl;
  }
  return prob;
}

double elFakeProbErr_v2 (int i_el){
  float prob = 0.0;
  float prob_error = 0.0;
  TH2F *theFakeRate = &fakeRateEl(el_v2);
  TH2F *theFakeRateErr = &fakeRateErrorEl(el_v2);
  // cut definition
  float pt = cms2.els_p4()[i_el].Pt();
  float eta = fabs(cms2.els_p4()[i_el].Eta());
  float upperEdge = theFakeRate->GetYaxis()->GetBinLowEdge(theFakeRate->GetYaxis()->GetNbins()) + theFakeRate->GetYaxis()->GetBinWidth(theFakeRate->GetYaxis()->GetNbins()) - 0.001;

  if ( pt > upperEdge ) pt = upperEdge;
     
  prob       = theFakeRate->GetBinContent(theFakeRate->FindBin(eta,pt));
  prob_error = theFakeRateErr->GetBinContent(theFakeRateErr->FindBin(eta,pt));
     
  if (prob>1.0 || prob<0.0) {
	  std::cout<<"ERROR FROM EL FAKE RATE!!! prob = " << prob << std::endl;
  }
  if (prob==0.0){
	  std::cout<<"ERROR FROM EL FAKE RATE!!! prob = " << prob
             <<" for Et = " <<cms2.els_p4()[i_el].Pt()
             <<" and Eta = " <<cms2.els_p4()[i_el].Eta()
             << std::endl;
  }
  return prob_error;
}
// --------------------------------
// electrons v3
// --------------------------------
bool isFakeDenominatorElectron_v3 (int index) {
 	if (!cms2.els_type()[index] & (1<<ISECALDRIVEN)) return false;
  if (fabs(cms2.els_p4()[index].eta()) > 2.5) return false;
  if (!electronId_noMuon(index)) return false;
  if (!electronId_cand01(index)) return false;
  //if (!electronImpact_cand01(index)) return false;
  //if (electronIsolation_relsusy_cand1(index, true) > 0.10) return false;
  if (electronIsolation_relsusy_cand1(index, true) > 0.40) return false;
  if (isFromConversionPartnerTrack(index)) return false;
  return true;
}
double elFakeProb_v3 (int i_el){
  float prob = 0.0;
  float prob_error = 0.0;
  TH2F *theFakeRate = &fakeRateEl(el_v3);
  TH2F *theFakeRateErr = &fakeRateErrorEl(el_v3);
  // cut definition
  float pt = cms2.els_p4()[i_el].Pt();
  float eta = fabs(cms2.els_p4()[i_el].Eta());
  float upperEdge = theFakeRate->GetYaxis()->GetBinLowEdge(theFakeRate->GetYaxis()->GetNbins()) + theFakeRate->GetYaxis()->GetBinWidth(theFakeRate->GetYaxis()->GetNbins()) - 0.001;

  if ( pt > upperEdge ) pt = upperEdge;
     
  prob       = theFakeRate->GetBinContent(theFakeRate->FindBin(eta,pt));
  prob_error = theFakeRateErr->GetBinContent(theFakeRateErr->FindBin(eta,pt));
     
  if (prob>1.0 || prob<0.0) {
	  std::cout<<"ERROR FROM EL FAKE RATE!!! prob = " << prob << std::endl;
  }
  if (prob==0.0){
	  std::cout<<"ERROR FROM EL FAKE RATE!!! prob = " << prob
             <<" for Et = " <<cms2.els_p4()[i_el].Pt()
             <<" and Eta = " <<cms2.els_p4()[i_el].Eta()
             << std::endl;
  }
  return prob;
}


double elFakeProbErr_v3 (int i_el){
  float prob = 0.0;
  float prob_error = 0.0;
  TH2F *theFakeRate = &fakeRateEl(el_v3);
  TH2F *theFakeRateErr = &fakeRateErrorEl(el_v3);
  // cut definition
  float pt = cms2.els_p4()[i_el].Pt();
  float eta = fabs(cms2.els_p4()[i_el].Eta());
  float upperEdge = theFakeRate->GetYaxis()->GetBinLowEdge(theFakeRate->GetYaxis()->GetNbins()) + theFakeRate->GetYaxis()->GetBinWidth(theFakeRate->GetYaxis()->GetNbins()) - 0.001;

  if ( pt > upperEdge ) pt = upperEdge;
     
  prob       = theFakeRate->GetBinContent(theFakeRate->FindBin(eta,pt));
  prob_error = theFakeRateErr->GetBinContent(theFakeRateErr->FindBin(eta,pt));
     
  if (prob>1.0 || prob<0.0) {
	  std::cout<<"ERROR FROM EL FAKE RATE!!! prob = " << prob << std::endl;
  }
  if (prob==0.0){
	  std::cout<<"ERROR FROM EL FAKE RATE!!! prob = " << prob
             <<" for Et = " <<cms2.els_p4()[i_el].Pt()
             <<" and Eta = " <<cms2.els_p4()[i_el].Eta()
             << std::endl;
  }
  return prob_error;
}
// --------------------------------
/* new root file access for 3X */
// --------------------------------
TH2F &fakeRateEl (enum fakeRateVersion version){
  if ( el_fakeRateFile == 0 ) {
    el_fakeRateFile = TFile::Open("$CMS2_LOCATION/NtupleMacros/data/el_FR_3X.root", "read");
    if ( el_fakeRateFile == 0 ) {
      std::cout << "$CMS2_LOCATION/NtupleMacros/data/el_FR_3X.root could not be found!!" << std::endl;
      std::cout << "Please make sure that $CMS2_LOCATION points to your CMS2 directory and that" << std::endl;
      std::cout << "$CMS2_LOCATION/NtupleMacros/data/el_FR_3X.root exists!" << std::endl;
      gSystem->Exit(1);
    }
    el_fakeRate_v1   = dynamic_cast<TH2F *>( el_fakeRateFile->Get("QCD30_el_IDn_ISO_04_FRptvseta") );
    el_fakeRate_v2   = dynamic_cast<TH2F *>( el_fakeRateFile->Get("QCD30_el_IDn_ISO_01_FRptvseta") );
    el_fakeRate_v3   = dynamic_cast<TH2F *>( el_fakeRateFile->Get("QCD30_el_IDy_ISO_04_FRptvseta") );
 }
 if( version == el_v1 ){ 
      return *el_fakeRate_v1;
 } else if( version == el_v2 ){
      return *el_fakeRate_v2;
 } else if( version == el_v3 ){ 
      return *el_fakeRate_v3;
 } 
 cout << "ERROR: unknown electron version" << endl;	
 gSystem->Exit(1);
 return *(TH2F*)0;
}
TH2F &fakeRateErrorEl (enum fakeRateVersion version){
 if ( el_fakeRateErrorFile == 0 ) {
    el_fakeRateErrorFile = TFile::Open("$CMS2_LOCATION/NtupleMacros/data/el_FR_3X.root", "read");
    if ( el_fakeRateErrorFile == 0 ) {
      std::cout << "$CMS2_LOCATION/NtupleMacros/data/el_FR_3X.root could not be found!!" << std::endl;
      std::cout << "Please make sure that $CMS2_LOCATION points to your CMS2 directory and that" << std::endl;
      std::cout << "$CMS2_LOCATION/NtupleMacros/data/el_FR_3X.root exists!" << std::endl;
      gSystem->Exit(1);
    }
    el_fakeRateErr_v1   = dynamic_cast<TH2F *>( el_fakeRateFile->Get("QCD30_el_IDn_ISO_04_FRErrptvseta") );
    el_fakeRateErr_v2   = dynamic_cast<TH2F *>( el_fakeRateFile->Get("QCD30_el_IDn_ISO_01_FRErrptvseta") );
    el_fakeRateErr_v3   = dynamic_cast<TH2F *>( el_fakeRateFile->Get("QCD30_el_IDy_ISO_04_FRErrptvseta") );
  } 
  if( version == el_v1 ){
    return *el_fakeRateErr_v1;
  } else if( version == el_v2 ){
    return *el_fakeRateErr_v2;
  } else if( version == el_v3 ){
    return *el_fakeRateErr_v3;
  } 
  cout << "ERROR: unknown electron version" << endl;
  gSystem->Exit(1);
  return *(TH2F *)0;
} 




// -----------------------------------------
// MC helper functions for fakerate tests:
// -----------------------------------------
int elFakeMCCategory(int i_el) {
  int category = -1;
  if(
     (abs(cms2.els_mc_id()[i_el])       == 11  && 
      abs(cms2.els_mc_motherid()[i_el]) == 22) ||
     (abs(cms2.els_mc_id()[i_el])       == 22) ||
     (abs(cms2.els_mc_id()[i_el])        > 100 && 
      abs(cms2.els_mc_id()[i_el])        < 200)||
     (abs(cms2.els_mc_id()[i_el])       == 11  && 
      abs(cms2.els_mc_motherid()[i_el]) == 111)
     ) {
    // electrons from gamma (conversion)
    category = 1;
  }
  else if(
          (abs(cms2.els_mc_id()[i_el]) > 200     && 
           abs(cms2.els_mc_id()[i_el]) < 400  )  ||
          (abs(cms2.els_mc_id()[i_el]) > 2000    && 
           abs(cms2.els_mc_id()[i_el]) < 4000 )  ||
          (abs(cms2.els_mc_id()[i_el]) == 11 && 
           abs(cms2.els_mc_motherid()[i_el]) > 200     && 
           abs(cms2.els_mc_motherid()[i_el]) < 400  )  || 
          (abs(cms2.els_mc_id()[i_el]) == 11 && 
           abs(cms2.els_mc_motherid()[i_el]) > 2000    && 
           abs(cms2.els_mc_motherid()[i_el]) < 4000 )  
          ) {
    // electron candidate or its mother is a light hadron
    category = 2;
  }
  else if( ( abs(cms2.els_mc_id()[i_el]) == 11 
             && abs(cms2.els_mc_motherid()[i_el]) >=400
             && abs(cms2.els_mc_motherid()[i_el]) <=600 )  || 
           ( abs(cms2.els_mc_id()[i_el]) == 11 
             && abs(cms2.els_mc_motherid()[i_el]) >=4000
             && abs(cms2.els_mc_motherid()[i_el]) <=6000 )
           ) {
    // heavy hadrons
    category = 3;
  }
  else {
    // the rest
    category = 4;
  }
  return category;
}
int muFakeMCCategory(int i_mu) {
  int category = -1;
  
  if( // punchthrough / sailthrough
     (abs(cms2.mus_mc_id()[i_mu]) != 13 )
     ) {
    category = 1;
  }
  else if( 
          abs(cms2.mus_mc_id()[i_mu]) == 13 && 
          abs(cms2.mus_mc_motherid()[i_mu]) < 400 
          ) {
    // light hadrons
    category = 2;
  }
  else if(
          ( abs(cms2.mus_mc_id()[i_mu]) == 13          &&
            abs(cms2.mus_mc_motherid()[i_mu]) >=400   &&
            abs(cms2.mus_mc_motherid()[i_mu]) <=600 ) || 
          ( abs(cms2.mus_mc_id()[i_mu]) == 13          &&
            abs(cms2.mus_mc_motherid()[i_mu]) >=4000   &&
            abs(cms2.mus_mc_motherid()[i_mu]) <=6000 )
          ) {
    // heavy hadrons
    category = 3;
  }
  else {
    // the rest
    category = 4;
  }
  return category;
}
