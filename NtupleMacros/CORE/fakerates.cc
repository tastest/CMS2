// ROOT includes
#include "TSystem.h"
#include "TFile.h"
#include "TH2.h"

// TAS includes
#include "fakerates.h"
#include "CMS2.h"
#include "electronSelections.h"
#include "muonSelections.h"


/*muons*/
class TH2F &fakeRateMuon (enum fakeRateVersion);
static TFile *mu_fakeRateFile = 0;
static TH2F  *mu_fakeRate = 0;
static TH2F  *mu_ttbar_fakeRate = 0;
bool isFakeableMuon (int index, enum fakeRateVersion version){

if( version == mu_v1 ){

  // return denominator for current version
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
else if ( version == mu_ttbar ){
    if ( cms2.mus_p4()[index].pt() < 10.) {
      std::cout << "muonID ERROR: requested muon is too low pt,  Abort." << std::endl;
      return false;
    }
    if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.5)  return false; // eta cut
    //if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) return false; //glb fit chisq
    if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 50) return false; //glb fit chisq
    if (((cms2.mus_type().at(index)) & (1<<1)) == 0)    return false; // global muon
    if (((cms2.mus_type().at(index)) & (1<<2)) == 0)    return false; // tracker muon
    if (cms2.mus_validHits().at(index) < 11)            return false; // # of tracker hits
    //if (TMath::Abs(cms2.mus_d0corr().at(index)) > 0.02) return false; // d0 from beamspot
    if (TMath::Abs(cms2.mus_d0corr().at(index)) > 0.2) return false; // d0 from beamspot
    return true;
} 
else {
  cout << "ERROR: unhandled muon version" << endl;
  return false;
}

}
double muFakeProb (int i_mu, enum fakeRateVersion version){

  // initializations
  float prob = 0.0;
  float pt = cms2.mus_p4()[i_mu].Pt();
  float eta = fabs(cms2.mus_p4()[i_mu].Eta());

  // Get FR(eta,pt) from 2D hist
  TH2F *theFakeRate = &fakeRateMuon(version);
  float upperEdge = theFakeRate->GetYaxis()->GetBinLowEdge(theFakeRate->GetYaxis()->GetNbins()) 
                    + theFakeRate->GetYaxis()->GetBinWidth(theFakeRate->GetYaxis()->GetNbins()) - 0.001;
  if ( pt > upperEdge ) pt = upperEdge;
  prob = theFakeRate->GetBinContent(theFakeRate->FindBin(eta,pt));

  // sanity check
  if ( prob > 1.0 || prob < 0.0 ) std::cout << "ERROR FROM MU FAKE RATE!!! prob = " << prob << std::endl;
  if (prob==0.0) std::cout << "ERROR FROM MU FAKE RATE!!! prob = " << prob <<" for Et = " <<cms2.mus_p4()[i_mu].Pt()
                           << " and Eta = " <<cms2.mus_p4()[i_mu].Eta() << std::endl;

  //
  return prob;
}
double muFakeProbErr (int i_mu, enum fakeRateVersion version){

  // initializations
  float prob = 0.0;
  float prob_error = 0.0;
  float pt = cms2.mus_p4()[i_mu].Pt();
  float eta = fabs(cms2.mus_p4()[i_mu].Eta());

  // Get FR(eta,pt) & error(eta,pt) from 2D hist 
  TH2F *theFakeRate = &fakeRateMuon(version);
  float upperEdge = theFakeRate->GetYaxis()->GetBinLowEdge(theFakeRate->GetYaxis()->GetNbins()) + theFakeRate->GetYaxis()->GetBinWidth(theFakeRate->GetYaxis()->GetNbins()) - 0.001;
  if ( pt > upperEdge ) pt = upperEdge;
  prob        = theFakeRate->GetBinContent( theFakeRate->FindBin(eta,pt) );
  prob_error  = theFakeRate->GetBinError( theFakeRate->FindBin(eta,pt) );

  // sanity check
  if ( prob > 1.0 || prob < 0.0 ) std::cout << "ERROR FROM MU FAKE RATE!!! prob = " << prob << std::endl;
  if (prob==0.0) std::cout << "ERROR FROM MU FAKE RATE!!! prob = " << prob <<" for Et = " <<cms2.mus_p4()[i_mu].Pt()
                           << " and Eta = " <<cms2.mus_p4()[i_mu].Eta() << std::endl;
  //
  return prob_error;
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
    mu_ttbar_fakeRate = dynamic_cast<TH2F *>( mu_fakeRateFile->Get("QCD30_mu_ttbar_FR_etavspt") );
  }
  if(version == mu_v1){
    return *mu_fakeRate;
  } else if (version == mu_ttbar){
    return *mu_ttbar_fakeRate;
  } 
  cout << "ERROR: unknown muon version" << endl;
  gSystem->Exit(1); 
  return *(TH2F*)0; 
}

/* electrons */
class TH2F &fakeRateEl (enum fakeRateVersion);
double elFakeProb_test (int, enum fakeRateVersion);
double elFakeProbErr_test (int, enum fakeRateVersion);

// file containing FR & errors for all supported denominators
static TFile *el_fakeRateFile = 0;

// histograms for supported denominators
static TH2F  *el_fakeRate_ttbar = 0;

static TH2F  *el_fakeRate_v1_cand01 = 0;
static TH2F  *el_fakeRate_v1_cand02 = 0;
static TH2F  *el_fakeRate_v1_cand02flip = 0;
static TH2F  *el_fakeRate_v2_cand01 = 0;
static TH2F  *el_fakeRate_v2_cand02 = 0;
static TH2F  *el_fakeRate_v2_cand02flip = 0;
static TH2F  *el_fakeRate_v3_cand01 = 0;
static TH2F  *el_fakeRate_v3_cand02 = 0;
static TH2F  *el_fakeRate_v3_cand02flip = 0;

// denominator selection
bool isFakeableElectron (int index, enum fakeRateVersion version){

// remove ID, ISO, IP
  if(version == el_v1_cand01){
    if (!cms2.els_type()[index] & (1<<ISECALDRIVEN)) return false;
    if (fabs(cms2.els_p4()[index].eta()) > 2.5) return false;
    if (!electronId_noMuon(index)) return false;
    //if (!electronId_cand01(index)) return false;
    //if (!electronImpact_cand01(index)) return false;
    //if (electronIsolation_relsusy_cand1(index, true) > 0.10) return false;
    if (electronIsolation_relsusy_cand1(index, true) > 0.40) return false;
    if (isFromConversionPartnerTrack(index)) return false;
    return true;
  } else if(version == el_v1_cand02){
    if (!cms2.els_type()[index] & (1<<ISECALDRIVEN)) return false;
    if (fabs(cms2.els_p4()[index].eta()) > 2.5) return false;
    if (!electronId_noMuon(index)) return false;
    //  if (!electronId_cand02(index)) return false;
    // extra on request of FKW
    //  if (!electronId_extra(index)) return false;
    //  if (!electronImpact_cand01(index)) return false;
    //  if (electronIsolation_relsusy_cand1(index, true) > 0.10) return false;
    if (electronIsolation_relsusy_cand1(index, true) > 0.40) return false;
    if (isFromConversionPartnerTrack(index)) return false;
    // on request of FKW we have also here:
    //  if (isChargeFlip(index)) return false;
    return true;
  } else if(version == el_v1_cand02flip){
    if (!cms2.els_type()[index] & (1<<ISECALDRIVEN)) return false;
    if (fabs(cms2.els_p4()[index].eta()) > 2.5) return false;
    if (!electronId_noMuon(index)) return false;
    //  if (!electronId_cand02(index)) return false;
    // extra on request of FKW
    //  if (!electronId_extra(index)) return false;
    //  if (!electronImpact_cand01(index)) return false;
    //  if (electronIsolation_relsusy_cand1(index, true) > 0.10) return false;
    if (electronIsolation_relsusy_cand1(index, true) > 0.40) return false;
    if (isFromConversionPartnerTrack(index)) return false;
    // on request of FKW we have also here:
    if (isChargeFlip(index)) return false;
    return true;
  } 
// remove ID, IP
    else if(version == el_v2_cand01){
    if (!cms2.els_type()[index] & (1<<ISECALDRIVEN)) return false;
    if (fabs(cms2.els_p4()[index].eta()) > 2.5) return false;
    if (!electronId_noMuon(index)) return false;
    //if (!electronId_cand01(index)) return false;
    //if (!electronImpact_cand01(index)) return false;
    if (electronIsolation_relsusy_cand1(index, true) > 0.10) return false;
    if (isFromConversionPartnerTrack(index)) return false;
    return true;
  } else if(version == el_v2_cand02){
    if (!cms2.els_type()[index] & (1<<ISECALDRIVEN)) return false;
    if (fabs(cms2.els_p4()[index].eta()) > 2.5) return false;
    if (!electronId_noMuon(index)) return false;
    //  if (!electronId_cand02(index)) return false;
    // extra on request of FKW
    //  if (!electronId_extra(index)) return false;
    //  if (!electronImpact_cand01(index)) return false;
    if (electronIsolation_relsusy_cand1(index, true) > 0.10) return false;
    if (isFromConversionPartnerTrack(index)) return false;
    // on request of FKW we have also here:
    // if (isChargeFlip(index)) return false;
    return true;
  } else if(version == el_v2_cand02flip){
    if (!cms2.els_type()[index] & (1<<ISECALDRIVEN)) return false;
    if (fabs(cms2.els_p4()[index].eta()) > 2.5) return false;
    if (!electronId_noMuon(index)) return false;
    //  if (!electronId_cand02(index)) return false;
    // extra on request of FKW
    //  if (!electronId_extra(index)) return false;
    //  if (!electronImpact_cand01(index)) return false;
    if (electronIsolation_relsusy_cand1(index, true) > 0.10) return false;
    if (isFromConversionPartnerTrack(index)) return false;
    // on request of FKW we have also here:
    if (isChargeFlip(index)) return false;
    return true;
  } 
// remove ISO, IP
    else if(version == el_v3_cand01){
    if (!cms2.els_type()[index] & (1<<ISECALDRIVEN)) return false;
    if (fabs(cms2.els_p4()[index].eta()) > 2.5) return false;
    if (!electronId_noMuon(index)) return false;
    if (!electronId_cand01(index)) return false;
    //if (!electronImpact_cand01(index)) return false;
    //if (electronIsolation_relsusy_cand1(index, true) > 0.10) return false;
    if (electronIsolation_relsusy_cand1(index, true) > 0.40) return false;
    if (isFromConversionPartnerTrack(index)) return false;
    return true;
  } else if(version == el_v3_cand02){
    if (!cms2.els_type()[index] & (1<<ISECALDRIVEN)) return false;
    if (fabs(cms2.els_p4()[index].eta()) > 2.5) return false;
    if (!electronId_noMuon(index)) return false;
    if (!electronId_cand02(index)) return false;
    // extra on request of FKW
    if (!electronId_extra(index)) return false;
    //  if (!electronImpact_cand01(index)) return false;
    //  if (electronIsolation_relsusy_cand1(index, true) > 0.10) return false;
    if (electronIsolation_relsusy_cand1(index, true) > 0.40) return false;
    if (isFromConversionPartnerTrack(index)) return false;
    // on request of FKW we have also here:
    //  if (isChargeFlip(index)) return false;
    return true;
  } else if(version == el_v3_cand02flip){
    if (!cms2.els_type()[index] & (1<<ISECALDRIVEN)) return false;
    if (fabs(cms2.els_p4()[index].eta()) > 2.5) return false;
    if (!electronId_noMuon(index)) return false;
    if (!electronId_cand02(index)) return false;
    // extra on request of FKW
    if (!electronId_extra(index)) return false;
    //  if (!electronImpact_cand01(index)) return false;
    //  if (electronIsolation_relsusy_cand1(index, true) > 0.10) return false;
    if (electronIsolation_relsusy_cand1(index, true) > 0.40) return false;
    if (isFromConversionPartnerTrack(index)) return false;
    // on request of FKW we have also here:
    if (isChargeFlip(index)) return false;
    return true;
  } 
  else if( version == el_ttbar ){
    if (!cms2.els_type()[index] & (1<<ISECALDRIVEN)) return false;
    if (fabs(cms2.els_p4()[index].eta()) > 2.5) return false;
    if (!electronId_noMuon(index)) return false;
    //if (!electronId_cand01(index)) return false;
    //if (!electronImpact_cand01(index)) return false; 
    if (electronIsolation_relsusy_cand1(index, true) > 0.10) return false;
    if (isFromConversionPartnerTrack(index)) return false;
    return true;
  }
  else {
    std::cout<<"isFakeable: invalid fakeRateVersion given. Check it!"<<std::endl;
    return false;
  }
}

double elFakeProb (int i_el, enum fakeRateVersion version){

  // initialization
  float prob = 0.0;
  float pt = cms2.els_p4().at(i_el).Pt();
  float eta = fabs(cms2.els_p4().at(i_el).Eta());

  // get the FR & from 2D hist
  TH2F *theFakeRate = &fakeRateEl( version );
  float upperEdge =  theFakeRate->GetYaxis()->GetBinLowEdge(theFakeRate->GetYaxis()->GetNbins()) 
                    + theFakeRate->GetYaxis()->GetBinWidth(theFakeRate->GetYaxis()->GetNbins()) - 0.001;
  if ( pt > upperEdge ) pt = upperEdge;
  prob       = theFakeRate->GetBinContent( theFakeRate->FindBin(eta,pt) );

  // sanity checks on FR and error
  if ( prob > 1.0 || prob < 0.0) std::cout << "ERROR FROM EL FAKE RATE!!! prob = " << prob << std::endl;
  if (prob==0.0) std::cout << "ERROR FROM EL FAKE RATE!!! prob = " << prob << " for Et = " << cms2.els_p4()[i_el].Pt() 
                           << " and Eta = " <<cms2.els_p4()[i_el].Eta() << std::endl;
  // return FR
  return prob;
}

double elFakeProbErr (int i_el, enum fakeRateVersion version){

  // initialization
  float prob = 0.0;
  float prob_error = 0.0;
  float pt = cms2.els_p4().at(i_el).Pt();
  float eta = fabs(cms2.els_p4().at(i_el).Eta());
  
  // get the FR & error from 2D hist
  TH2F *theFakeRate = &fakeRateEl( version );
  float upperEdge = theFakeRate->GetYaxis()->GetBinLowEdge(theFakeRate->GetYaxis()->GetNbins()) 
                    + theFakeRate->GetYaxis()->GetBinWidth(theFakeRate->GetYaxis()->GetNbins()) - 0.001;
  if ( pt > upperEdge ) pt = upperEdge;
  prob       = theFakeRate->GetBinContent(theFakeRate->FindBin(eta,pt));
  prob_error = theFakeRate->GetBinError(theFakeRate->FindBin(eta,pt));
 
  // sanity checks on FR and error
  if ( prob > 1.0 || prob < 0.0) std::cout << "ERROR FROM EL FAKE RATE!!! prob = " << prob << std::endl;
  if (prob==0.0) std::cout << "ERROR FROM EL FAKE RATE!!! prob = " << prob << " for Et = " << cms2.els_p4()[i_el].Pt() 
                           << " and Eta = " <<cms2.els_p4()[i_el].Eta() << std::endl;
  
  //
  return prob_error;
}

// read fake rate & errors from ROOT file in cvs NtupleMacros/data
TH2F &fakeRateEl (enum fakeRateVersion version){
  if ( el_fakeRateFile == 0 ) {
    el_fakeRateFile = TFile::Open("$CMS2_LOCATION/NtupleMacros/data/el_FR_3X.root", "read");
    if ( el_fakeRateFile == 0 ) {
      std::cout << "$CMS2_LOCATION/NtupleMacros/data/el_FR_3X.root could not be found!!" << std::endl;
      std::cout << "Please make sure that $CMS2_LOCATION points to your CMS2 directory and that" << std::endl;
      std::cout << "$CMS2_LOCATION/NtupleMacros/data/el_FR_3X.root exists!" << std::endl;
      gSystem->Exit(1);
    }
    el_fakeRate_v1_cand01     = dynamic_cast<TH2F *>( el_fakeRateFile->Get("QCD30_el_v1_cand01_FR_etavspt") );
    el_fakeRate_v1_cand02     = dynamic_cast<TH2F *>( el_fakeRateFile->Get("QCD30_el_v1_cand02_FR_etavspt") );
    el_fakeRate_v1_cand02flip = dynamic_cast<TH2F *>( el_fakeRateFile->Get("QCD30_el_v1_cand02flip_FR_etavspt") );
    el_fakeRate_v2_cand01     = dynamic_cast<TH2F *>( el_fakeRateFile->Get("QCD30_el_v2_cand01_FR_etavspt") );
    el_fakeRate_v2_cand02     = dynamic_cast<TH2F *>( el_fakeRateFile->Get("QCD30_el_v2_cand02_FR_etavspt") );
    el_fakeRate_v2_cand02flip = dynamic_cast<TH2F *>( el_fakeRateFile->Get("QCD30_el_v2_cand02flip_FR_etavspt") );
    el_fakeRate_v3_cand01     = dynamic_cast<TH2F *>( el_fakeRateFile->Get("QCD30_el_v3_cand01_FR_etavspt") );
    el_fakeRate_v3_cand02     = dynamic_cast<TH2F *>( el_fakeRateFile->Get("QCD30_el_v3_cand02_FR_etavspt") );
    el_fakeRate_v3_cand02flip = dynamic_cast<TH2F *>( el_fakeRateFile->Get("QCD30_el_v3_cand02flip_FR_etavspt") );
    el_fakeRate_ttbar = dynamic_cast<TH2F *>( el_fakeRateFile->Get("QCD30_el_ttbar_FR_etavspt") );
  }
  if( version == el_v1_cand01 ){ 
    return *el_fakeRate_v1_cand01;
  } else if( version == el_v1_cand02 ){ 
    return *el_fakeRate_v1_cand02;
  } else if( version == el_v1_cand02flip ){ 
    return *el_fakeRate_v1_cand02flip;
  } else if( version == el_v2_cand01 ){
    return *el_fakeRate_v2_cand01;
  } else if( version == el_v2_cand02 ){
    return *el_fakeRate_v2_cand02;
  } else if( version == el_v2_cand02flip ){
    return *el_fakeRate_v2_cand02flip;
  } else if( version == el_v3_cand01 ){ 
    return *el_fakeRate_v3_cand01;
  } else if( version == el_v3_cand02 ){ 
    return *el_fakeRate_v3_cand02;
  } else if( version == el_v3_cand02flip ){ 
    return *el_fakeRate_v3_cand02flip;
  } else if( version == el_ttbar ){ 
    return *el_fakeRate_ttbar;
  }
  cout << "ERROR: unknown electron version" << endl;	
  gSystem->Exit(1);
  return *(TH2F*)0;
}

void PrintTH2F( TH2F* theFakeRate ){
  for(int j=0; j <= 1+theFakeRate->GetNbinsX(); j++){
    cout << j << "\t";
  }
  cout << endl;
  for(int i=0; i <= 1+theFakeRate->GetNbinsY(); i++){
    cout << i << "\t";
    for(int j=0; j <= 1+theFakeRate->GetNbinsX(); j++){
      cout << Form("%.07f", theFakeRate->GetBinContent(j,1+theFakeRate->GetNbinsY()-i) ) << "\t";
    }
    cout << endl;
  }
  return;
} 
  
//
void PrintErrTH2F( TH2F* theFakeRate ){
  for(int j=0; j <= 1+theFakeRate->GetNbinsX(); j++){
    cout << j << "\t";
  }
  cout << endl;
  for(int i=0; i <= 1+theFakeRate->GetNbinsY(); i++){
    cout << i << "\t";
    for(int j=0; j <= 1+theFakeRate->GetNbinsX(); j++){
      cout << Form("%.07f", theFakeRate->GetBinError(j,1+theFakeRate->GetNbinsY()-i) ) << "\t";
    }
    cout << endl; 
  }
  return;
} 

