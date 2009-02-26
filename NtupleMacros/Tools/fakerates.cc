#include "TFile.h"
#include "TSystem.h"
#include "TH2.h"
#include "CORE/CMS2.h"
#include "CORE/selections.h"
#include "Tools/fakerates.h"

static TFile *el_fakeRateFile_v2_2 = 0;
static TH2F  *el_fakeRate_v2_2 = 0;

static TFile *el_fakeRateFile_v5 = 0;
static TH2F  *el_fakeRate_v5 = 0;
static TH2F  *el_fakeRate_err_v5 = 0;

static TFile *el_fakeRateFile_v6 = 0;
static TH2F  *el_fakeRate_v6 = 0;
static TH2F  *el_fakeRate_err_v6 = 0;

static TFile *el_fakeRateFile_v7 = 0;
static TH2F  *el_fakeRate_v7 = 0;
static TH2F  *el_fakeRate_err_v7 = 0;

double elFakeProb_v2_2 (int i_el, int add_error_times)
{
     float prob = 0.0;
     float prob_error = 0.0;
     TH2F *theFakeRate = &fakeRate();
     // cut definition
     float pt = cms2.els_p4()[i_el].Pt();
     float upperEdge = theFakeRate->GetYaxis()->GetBinLowEdge(theFakeRate->GetYaxis()->GetNbins()) + theFakeRate->GetYaxis()->GetBinWidth(theFakeRate->GetYaxis()->GetNbins()) - 0.001;
     if ( pt > upperEdge )
       pt = upperEdge;
     prob = theFakeRate->GetBinContent(theFakeRate->FindBin(cms2.els_p4()[i_el].Eta(),pt));
     prob_error =
	  theFakeRate->GetBinError(theFakeRate->FindBin(cms2.els_p4()[i_el].Eta(),pt));
     
     if (prob>1.0 || prob<0.0) {
// 	  std::cout<<"ERROR FROM FAKE RATE!!! prob = " << prob << std::endl;
     }
     if (prob==0.0){
// 	  std::cout<<"ERROR FROM FAKE RATE!!! prob = " << prob
// 		   <<" for Et = " <<cms2.els_p4()[i_el].Pt()
// 		   <<" and Eta = " <<cms2.els_p4()[i_el].Eta()
// 		   << std::endl;
     }
     return prob+add_error_times*prob_error;
}

// Is the i-th electron in the electron block a fakeable object?
// For now: return true
bool isFakeable_v2_2 (int i_el) {
     //
     // returns true if input fulfills certain cuts
     //
     
     // cut definition
     float et_cut        = 0.;
     float pt_cut        = 15.;
     float eta_cut       = 2.5;
     //   float tkIso_cut     = 10.; // was 50
     float iso_ratio_cut = 0.92; //
     float eOverP_cut    = 999999.99;
     float hOverE_cut    = 0.2;
     
     float iso_ratio = 0.0;
     
     if( (cms2.els_p4()[i_el].Pt()+cms2.els_tkIso()[i_el]) > 0.0 ) 
	  iso_ratio = cms2.els_p4()[i_el].Pt()/(cms2.els_p4()[i_el].Pt()+cms2.els_tkIso()[i_el]);
     else iso_ratio = 0.0; // reject events with 0 momentum - do we have thses at all?
	  
     bool result = true;
     
     if (cms2.els_closestMuon().at(i_el) != -1)		result = false;
     if ( cms2.els_eSC()[i_el]      < et_cut )                 result = false;
     if ( cms2.els_p4()[i_el].Pt()  < pt_cut )                 result = false;
     if ( TMath::Abs(cms2.els_p4()[i_el].Eta()) > eta_cut )      result = false;
     //   // previous iso requirement, use this OR the one below!
     //   if ( cms2.els_tkIso()[i_el]    > tkIso_cut )              result = false;
     //new isolation requirement
     if ( iso_ratio               < iso_ratio_cut )          result = false;
     if ( cms2.els_eOverPIn()[i_el] > eOverP_cut )             result = false;
     if ( cms2.els_hOverE()[i_el]   > hOverE_cut )             result = false;
     
     return result;
}

bool isNumeratorElectron_v2_2 (int index, int type) { // 0=loose, 1=tight, for pass4: 1=loose, 2=tight
     //
     // returns true if input fulfills certain cuts
     //
     
     // cut definition
     float et_cut        = 0.;
     float pt_cut        = 15;
     float eta_cut       = 2.5;
     //   float tkIso_cut     = 5.;
     //isolation requirement
     float iso_ratio_cut = 0.92; // use alternatively to iso cut above
     float eOverP_cut    = 999999.99;
     //float eOverP_cut   = 3.;
     float hOverE_cut    = 0.2;
     float d0_cut        = 0.025;

     float iso_ratio =
	  cms2.els_p4()[index].Pt()/(cms2.els_p4()[index].Pt()+cms2.els_tkIso()[index]);
     
     unsigned int njets_cut        = 1;  // require at least N jets
     //  float HLT_jet_approx = 30.0; // require leading jet to be larger than N GeV v2_3
     //  float HLT_jet_approx = 60.0; // require leading jet to be larger than N GeV v2_4
     //  float HLT_jet_approx = 90.0; // require leading jet to be larger than N GeV v2_5
//      float HLT_jet_approx = 120.0; // require leading jet to be larger than N GeV v2_6
     //  float HLT_dijet_approx = 15.0; // require leading dijet to be larger than N=(et1+et2)/2 GeV
     
     bool result = true;
     
     if (cms2.els_closestMuon().at(index) != -1)		result = false;
     if (cms2.evt_njets() < njets_cut)                              result = false;
//      if ( jets_p4->at(0).Pt() < HLT_jet_approx )             result = false;
     //  if ( (jets_p4->at(0).Pt()+jets_p4->at(1).Pt())/2. < HLT_jet_approx )     result = false; // hmm! did not require a jet to be preset in the first case..
     
     if ( cms2.els_eSC()[index]      < et_cut )                 result = false;
     if ( cms2.els_p4()[index].Pt()  < pt_cut )                 result = false;
     if ( TMath::Abs(cms2.els_p4()[index].Eta()) > eta_cut )      result = false;
     //   // previous iso requirement, use this OR the one below!
     //   if ( els_tkIso->at(index)    > tkIso_cut )              result = false;
     //new isolation requirement
     if ( iso_ratio               < iso_ratio_cut )          result = false;
     if ( cms2.els_eOverPIn()[index] > eOverP_cut )             result = false;
     if ( cms2.els_hOverE()[index]   > hOverE_cut )             result = false;
     // add additional cleaning cuts (from FKW) 080324
     if ( TMath::Abs(cms2.els_d0corr()[index])  > d0_cut )            result = false;
     
     /*
       bool IdCuts = cut_verysimple(els_dEtaIn->at(index),
       els_dPhiIn->at(index),
       els_hOverE->at(index),
       els_eSeedOverPOut->at(index),
       els_sigmaEtaEta->at(index),
       els_p4->at(index).Eta());
     */
     
     /*
       bool IdCuts = electron_selection(index, type);
       if (!IdCuts) result = false;
     */
     
     // _pass4 has 3 types - 0=robust, 1=loose, 2=tight
     // - need to adjust in all places here
     bool IdCuts = cms2.els_tightId()[index];
     if (!IdCuts) result = false;
     
     return result;
}

double elFakeProb_v5 (int i_el, int add_error_times)
{
     float prob = 0.0;
     float prob_error = 0.0;
     TH2F *theFakeRate = &fakeRate();
     TH2F *theFakeRateErr = &fakeRateError();
     // cut definition
     float pt = cms2.els_p4()[i_el].Pt();
     float upperEdge = theFakeRate->GetYaxis()->GetBinLowEdge(theFakeRate->GetYaxis()->GetNbins()) + theFakeRate->GetYaxis()->GetBinWidth(theFakeRate->GetYaxis()->GetNbins()) - 0.001;
     if ( pt > upperEdge )
       pt = upperEdge;
     prob = theFakeRate->GetBinContent(theFakeRate->FindBin(cms2.els_p4()[i_el].Eta(),pt));
     prob_error =
	  theFakeRateErr->GetBinContent(theFakeRateErr->FindBin(cms2.els_p4()[i_el].Eta(),pt));
     
     if (prob>1.0 || prob<0.0) {
// 	  std::cout<<"ERROR FROM FAKE RATE!!! prob = " << prob << std::endl;
     }
     if (prob==0.0){
// 	  std::cout<<"ERROR FROM FAKE RATE!!! prob = " << prob
// 		   <<" for Et = " <<cms2.els_p4()[i_el].Pt()
// 		   <<" and Eta = " <<cms2.els_p4()[i_el].Eta()
// 		   << std::endl;
     }
     return prob+add_error_times*prob_error;
}

bool isFakeDenominatorElectron_v5 (int index) 
{
  //
  // returns true if input fulfills certain cuts
  //

  // cut definition
  float pt_cut        		= 15.;
  float eta_cut       		= 2.5;
  float hOverE_cut    		= 0.2;
  bool  use_calo_iso            = false;

  bool result = true;

  if (cms2.els_closestMuon().at(index) != -1)		result = false;
  if ( cms2.els_p4()[index].Pt()  < pt_cut )            result = false;
  if ( TMath::Abs(cms2.els_p4()[index].Eta()) > eta_cut ) result = false;
  if ( !passElectronIsolation(index,use_calo_iso) )          	result = false;
  //  if ( !passElectronIsolationLoose(index,true) )          	result = false; //v5_2
  if ( !passElectronIsolationLoose2(index,true) )          	result = false; //v5_4
  if ( cms2.els_hOverE()[index]   > hOverE_cut )        result = false;

  return result;

}

bool isFakeNumeratorElectron_v5 (int index, int type) 
{ 
  //
  // 1=loose, 2=tight
  //
  // returns true if input fulfills certain cuts
  //
  
  // cut definition
  float pt_cut        		= 15;
  float eta_cut       		= 2.5;
  bool  use_calo_iso            = true;

  bool result = true;

  if (cms2.els_closestMuon().at(index) != -1)		result = false;
  if ( cms2.els_p4()[index].Pt()  < pt_cut )                 result = false;
  if ( TMath::Abs(cms2.els_p4()[index].Eta()) > eta_cut )      result = false;
  if ( !passElectronIsolation(index,use_calo_iso) )          	result = false;
  if ( type == 1 ) {
    // loose
    if ( !goodLooseElectronWithoutIsolation(index) )   result = false;
  } else if ( type == 2 ) {
    // tight
    if ( !goodElectronWithoutIsolation(index) )   result = false;
  } else {
    cout << "WARNING: wrong electron type detected, please select loose (1) or tight (2)" << endl;
  }

  return result;
  
}

double elFakeProb_v6 (int i_el, int add_error_times)
{
     float prob = 0.0;
     float prob_error = 0.0;
     TH2F *theFakeRate = &fakeRate();
     TH2F *theFakeRateErr = &fakeRateError();
     // cut definition
     float pt = cms2.els_p4()[i_el].Pt();
     float upperEdge = theFakeRate->GetYaxis()->GetBinLowEdge(theFakeRate->GetYaxis()->GetNbins()) + theFakeRate->GetYaxis()->GetBinWidth(theFakeRate->GetYaxis()->GetNbins()) - 0.001;
     if ( pt > upperEdge )
       pt = upperEdge;
     prob = theFakeRate->GetBinContent(theFakeRate->FindBin(cms2.els_p4()[i_el].Eta(),pt));
     prob_error =
	  theFakeRateErr->GetBinContent(theFakeRateErr->FindBin(cms2.els_p4()[i_el].Eta(),pt));
     
     if (prob>1.0 || prob<0.0) {
// 	  std::cout<<"ERROR FROM FAKE RATE!!! prob = " << prob << std::endl;
     }
     if (prob==0.0){
// 	  std::cout<<"ERROR FROM FAKE RATE!!! prob = " << prob
// 		   <<" for Et = " <<cms2.els_p4()[i_el].Pt()
// 		   <<" and Eta = " <<cms2.els_p4()[i_el].Eta()
// 		   << std::endl;
     }
     return prob+add_error_times*prob_error;
}

bool isFakeDenominatorElectron_v6 (int index) 
{
  //
  // returns true if input fulfills certain cuts
  //

  // cut definition
  float pt_cut        		= 15.;
  float eta_cut       		= 2.5;
  float hOverE_cut    		= 0.2;
  bool  use_calo_iso            = false;

  bool result = true;

  if (cms2.els_closestMuon().at(index) != -1)		result = false;
  if ( cms2.els_p4()[index].Pt()  < pt_cut )            result = false;
  if ( TMath::Abs(cms2.els_p4()[index].Eta()) > eta_cut ) result = false;
  if ( !passElectronIsolation(index,use_calo_iso) )          	result = false;
  //  if ( !passElectronIsolationLoose(index,true) )          	result = false; //v5_2
  if ( !passElectronIsolationLoose2(index,true) )          	result = false; //v5_4
  if ( cms2.els_hOverE()[index]   > hOverE_cut )        result = false;

  return result;

}

bool isFakeNumeratorElectron_v6 (int index, int type) 
{ 
  //
  // 1=loose, 2=tight
  //
  // returns true if input fulfills certain cuts
  //
  
  // cut definition
  float pt_cut        		= 15;
  float eta_cut       		= 2.5;
  bool  use_calo_iso            = true;

  bool result = true;

  if (cms2.els_closestMuon().at(index) != -1)		result = false;
  if ( cms2.els_p4()[index].Pt()  < pt_cut )                 result = false;
  if ( TMath::Abs(cms2.els_p4()[index].Eta()) > eta_cut )      result = false;
  if ( !passElectronIsolation(index,use_calo_iso) )          	result = false;
  if ( type == 1 ) {
    // loose
    if ( !goodLooseElectronWithoutIsolation(index) )   result = false;
  } else if ( type == 2 ) {
    // tight
    if ( !goodElectronWithoutIsolation(index) )   result = false;
  } else {
    cout << "WARNING: wrong electron type detected, please select loose (1) or tight (2)" << endl;
  }

  return result;
  
}
//asdasd
double elFakeProb_v7 (int i_el, int add_error_times)
{
     float prob = 0.0;
     float prob_error = 0.0;
     TH2F *theFakeRate = &fakeRate();
     TH2F *theFakeRateErr = &fakeRateError();
     // cut definition
     float pt = cms2.els_p4()[i_el].Pt();
     float upperEdge = theFakeRate->GetYaxis()->GetBinLowEdge(theFakeRate->GetYaxis()->GetNbins()) + theFakeRate->GetYaxis()->GetBinWidth(theFakeRate->GetYaxis()->GetNbins()) - 0.001;
     if ( pt > upperEdge )
       pt = upperEdge;
     prob = theFakeRate->GetBinContent(theFakeRate->FindBin(cms2.els_p4()[i_el].Eta(),pt));
     prob_error =
	  theFakeRateErr->GetBinContent(theFakeRateErr->FindBin(cms2.els_p4()[i_el].Eta(),pt));
     
     if (prob>1.0 || prob<0.0) {
	  std::cout<<"ERROR FROM FAKE RATE!!! prob = " << prob << std::endl;
     }
     if (prob==0.0){
	  std::cout<<"ERROR FROM FAKE RATE!!! prob = " << prob
		   <<" for Et = " <<cms2.els_p4()[i_el].Pt()
		   <<" and Eta = " <<cms2.els_p4()[i_el].Eta()
		   << std::endl;
     }
     return prob+add_error_times*prob_error;
}

bool isFakeDenominatorElectron_v7 (int index) 
{
  //
  // returns true if input fulfills certain cuts
  //

  // cut definition
  float pt_cut        		= 20.;
  float eta_cut       		= 2.5;
  float hOverE_cut    		= 0.2;
  bool  use_calo_iso            = false;

  bool result = true;

  if (cms2.els_closestMuon().at(index) != -1)		result = false;
  if ( cms2.els_p4()[index].Pt()  < pt_cut )            result = false;
  if ( TMath::Abs(cms2.els_p4()[index].Eta()) > eta_cut ) result = false;
  if ( !passElectronIsolation(index,use_calo_iso) )          	result = false;
  //  if ( !passElectronIsolationLoose(index,true) )          	result = false; //v5_2
  if ( !passElectronIsolationLoose2(index,true) )          	result = false; //v5_4
  if ( cms2.els_hOverE()[index]   > hOverE_cut )        result = false;

  return result;

}

bool isFakeNumeratorElectron_v7 (int index, int type) 
{ 
  //
  // 1=loose, 2=tight
  //
  // returns true if input fulfills certain cuts
  //
  
  // cut definition
  float pt_cut        		= 20;
  float eta_cut       		= 2.5;
  bool  use_calo_iso            = true;

  bool result = true;

  if (cms2.els_closestMuon().at(index) != -1)		result = false;
  if ( cms2.els_p4()[index].Pt()  < pt_cut )                 result = false;
  if ( TMath::Abs(cms2.els_p4()[index].Eta()) > eta_cut )      result = false;
  if ( !passElectronIsolation(index,use_calo_iso) )          	result = false;
  if ( type == 1 ) {
    // loose
    if ( !goodLooseElectronWithoutIsolation(index) )   result = false;
  } else if ( type == 2 ) {
    // tight
    if ( !goodElectronWithoutIsolation(index) )   result = false;
  } else {
    cout << "WARNING: wrong electron type detected, please select loose (1) or tight (2)" << endl;
  }

  return result;
  
}
bool isFakeDenominatorMuon_v1 (int index) 
{
  //
  // returns true if input fulfills certain cuts
  //

  // cut definition
  float pt_cut        		= 20.;
  float eta_cut       		= 2.5;

  bool result = true;

  // need: 

  // - global muon
  if ( (cms2.mus_type().at(index)&0x2)==0 )                                result = false;

  // - pt > 20 GeV
  if ( cms2.mus_p4()[index].Pt()  < pt_cut )                               result = false;

  // - d0corr < 0.1 cm -> loosened to 0.2!!
  if (TMath::Abs(cms2.mus_d0corr().at(index))   > 0.2)                     result = false;

  // - chi2/ndf < 20 (?)
  if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) > 20.) result = false;

  // - isoSumPt < 15 -> needs revision - would go for same def as for electrons?
  // - changed iso to mu_rel_iso > 0.75
  if ( !passMuonIsolationLoose(index) )          	                   result = false; 

  //  if (cms2.mus_validHits().at(index) < 11)                             result = false;

  if ( TMath::Abs(cms2.mus_p4()[index].Eta()) > eta_cut )                  result = false;

  return result;

}

bool isFakeNumeratorMuon_v1 (int index, int type) 
{ 
  //
  // returns true if input fulfills certain cuts
  //
  
  // cut definition
  float pt_cut        		= 20;
  float eta_cut       		= 2.5;

  bool result = true;

  // Jakes original definition
  // - global muon
  // - pt > 20 GeV
  // - d0corr < 0.1 cm
  // - chi2/ndf < 10
  // - isoSumPt < 3 -> needs revision - would go for same def as for electrons?

  // now using WW analysis selection

  if ( cms2.mus_p4()[index].Pt()  < pt_cut )              result = false;
  if ( TMath::Abs(cms2.mus_p4()[index].Eta()) > eta_cut ) result = false;
  if ( !passMuonIsolation(index) )          	          result = false;
  if ( !goodMuonWithoutIsolation(index) )                 result = false;

  return result;
  
}

#define USE_V7

bool isFakeable (int i_el)
{
#ifdef USE_V5
  return isFakeDenominatorElectron_v5(i_el);
#endif
#ifdef USE_V6
  return isFakeDenominatorElectron_v6(i_el);
#endif
#ifdef USE_V7
  return isFakeDenominatorElectron_v7(i_el);
#else
  return isFakeable_v2_2(i_el);
#endif
}

bool isFakeableMuon (int i_mu)
{
  return  isFakeDenominatorMuon_v1(i_mu);
}

double elFakeProb (int i_el, int add_error_times)
{
#ifdef USE_V5
     return elFakeProb_v5(i_el, add_error_times);
#endif
#ifdef USE_V6
     return elFakeProb_v6(i_el, add_error_times);
#endif
#ifdef USE_V7
     return elFakeProb_v7(i_el, add_error_times);
#else
     return elFakeProb_v2_2(i_el, add_error_times);
#endif
}

double muFakeProb (int i_mu, int add_error_times)
{
  //     return muFakeProb_v1(i_mu, add_error_times);
     return -999.99;
}

bool isNumeratorElectron (int index, int type)
{
#ifdef USE_V5
     return isFakeNumeratorElectron_v5(index, 2);
#endif
#ifdef USE_V6
     return isFakeNumeratorElectron_v6(index, 2);
#endif
#ifdef USE_V7
     return isFakeNumeratorElectron_v7(index, 2);
#else
     return isNumeratorElectron_v2_2(index, type);
#endif
}

bool isNumeratorMuon (int index, int type)
{
     return isFakeNumeratorMuon_v1(index, 2);
}

TH2F &fakeRate ()
{
#ifdef USE_V5
     if ( el_fakeRateFile_v5 == 0 ) {
	  el_fakeRateFile_v5 = TFile::Open("$CMS2_LOCATION/NtupleMacros/data/fakeRates-v5_5.root", "read"); 
	  if ( el_fakeRateFile_v5 == 0 ) {
	       std::cout << "$CMS2_LOCATION/NtupleMacros/data/fakeRates-v5_5.root could not be found!!" << std::endl;
	       std::cout << "Please make sure that $CMS2_LOCATION points to your CMS2 directory and that" << std::endl;
	       std::cout << "$CMS2_LOCATION/NtupleMacros/data/fakeRates-v5_5.root exists!" << std::endl;
	       gSystem->Exit(1);
	  }
	  el_fakeRate_v5 = dynamic_cast<TH2F *>(el_fakeRateFile_v5->Get("fakeRateTemplate_wo_leading_elt_fakeRatesFull"));
	  el_fakeRate_err_v5 = dynamic_cast<TH2F *>(el_fakeRateFile_v5->Get("fakeRateTemplateError_wo_leading_elt_fakeRatesFull"));
     }
     return *el_fakeRate_v5;
#endif
#ifdef USE_V6
     if ( el_fakeRateFile_v6 == 0 ) {
	  el_fakeRateFile_v6 = TFile::Open("$CMS2_LOCATION/NtupleMacros/data/fakeRates-v6_6.root", "read"); 
	  if ( el_fakeRateFile_v6 == 0 ) {
	       std::cout << "$CMS2_LOCATION/NtupleMacros/data/fakeRates-v6_6.root could not be found!!" << std::endl;
	       std::cout << "Please make sure that $CMS2_LOCATION points to your CMS2 directory and that" << std::endl;
	       std::cout << "$CMS2_LOCATION/NtupleMacros/data/fakeRates-v6_6.root exists!" << std::endl;
	       gSystem->Exit(1);
	  }
 	  el_fakeRate_v6 = dynamic_cast<TH2F *>(el_fakeRateFile_v6->Get("fakeRateTemplate_wo_leading_elt_EleFakes"));
 	  el_fakeRate_err_v6 = dynamic_cast<TH2F *>(el_fakeRateFile_v6->Get("fakeRateTemplateError_wo_leading_elt_EleFakes"));
     }
     return *el_fakeRate_v6;
#endif
#ifdef USE_V7
     if ( el_fakeRateFile_v7 == 0 ) {
	  el_fakeRateFile_v7 = TFile::Open("$CMS2_LOCATION/NtupleMacros/data/fakeRates-v7_2.root", "read"); 
	  if ( el_fakeRateFile_v7 == 0 ) {
	       std::cout << "$CMS2_LOCATION/NtupleMacros/data/fakeRates-v7_2.root could not be found!!" << std::endl;
	       std::cout << "Please make sure that $CMS2_LOCATION points to your CMS2 directory and that" << std::endl;
	       std::cout << "$CMS2_LOCATION/NtupleMacros/data/fakeRates-v7_2.root exists!" << std::endl;
	       gSystem->Exit(1);
	  }
	  el_fakeRate_v7 = dynamic_cast<TH2F *>(el_fakeRateFile_v7->Get("fakeRateTemplate_elt_EleFakes"));
	  el_fakeRate_err_v7 = dynamic_cast<TH2F *>(el_fakeRateFile_v7->Get("fakeRateTemplateError_elt_EleFakes"));
     }
     return *el_fakeRate_v7;
#else
     if ( el_fakeRateFile_v2_2 == 0 ) {
	  el_fakeRateFile_v2_2 = TFile::Open("$CMS2_LOCATION/NtupleMacros/data/fakeRates-v2_2_allpt.root", "read");
	  if ( el_fakeRateFile_v2_2 == 0 ) {
	       std::cout << "$CMS2_LOCATION/NtupleMacros/data/fakeRates-v2_2_allpt.root could not be found!!" << std::endl;
	       std::cout << "Please make sure that $CMS2_LOCATION points to your CMS2 directory and that" << std::endl;
	       std::cout << "$CMS2_LOCATION/NtupleMacros/data/fakeRates-v2_2_allpt.root exists!" << std::endl;
	       gSystem->Exit(1);
	  }
	  el_fakeRate_v2_2 = dynamic_cast<TH2F *>(el_fakeRateFile_v2_2->Get("fakeRate_wo_leading_elt_qcd"));
     }
     return *el_fakeRate_v2_2;
#endif
}

TH2F &fakeRateError ()
{
#ifdef USE_V5
     if ( el_fakeRateFile_v5 == 0 ) {
	  el_fakeRateFile_v5 = TFile::Open("$CMS2_LOCATION/NtupleMacros/data/fakeRates-v5_5.root", "read"); 
	  if ( el_fakeRateFile_v5 == 0 ) {
	       std::cout << "$CMS2_LOCATION/NtupleMacros/data/fakeRates-v5_5.root could not be found!!" << std::endl;
	       std::cout << "Please make sure that $CMS2_LOCATION points to your CMS2 directory and that" << std::endl;
	       std::cout << "$CMS2_LOCATION/NtupleMacros/data/fakeRates-v5_5.root exists!" << std::endl;
	       gSystem->Exit(1);
	  }
	  el_fakeRate_v5 = dynamic_cast<TH2F *>(el_fakeRateFile_v5->Get("fakeRateTemplate_wo_leading_elt_fakeRatesFull"));
	  el_fakeRate_err_v5 = dynamic_cast<TH2F *>(el_fakeRateFile_v5->Get("fakeRateTemplateError_wo_leading_elt_fakeRatesFull"));
     }
     return *el_fakeRate_err_v5;
#endif
#ifdef USE_V6
     if ( el_fakeRateFile_v6 == 0 ) {
	  el_fakeRateFile_v6 = TFile::Open("$CMS2_LOCATION/NtupleMacros/data/fakeRates-v6_6.root", "read"); 
	  if ( el_fakeRateFile_v6 == 0 ) {
	       std::cout << "$CMS2_LOCATION/NtupleMacros/data/fakeRates-v6_6.root could not be found!!" << std::endl;
	       std::cout << "Please make sure that $CMS2_LOCATION points to your CMS2 directory and that" << std::endl;
	       std::cout << "$CMS2_LOCATION/NtupleMacros/data/fakeRates-v6_6.root exists!" << std::endl;
	       gSystem->Exit(1);
	  }
 	  el_fakeRate_v6 = dynamic_cast<TH2F *>(el_fakeRateFile_v6->Get("fakeRateTemplate_wo_leading_elt_EleFakes"));
 	  el_fakeRate_err_v6 = dynamic_cast<TH2F *>(el_fakeRateFile_v6->Get("fakeRateTemplateError_wo_leading_elt_EleFakes"));
     }
     return *el_fakeRate_err_v6;
#endif
#ifdef USE_V7
     if ( el_fakeRateFile_v7 == 0 ) {
	  el_fakeRateFile_v7 = TFile::Open("$CMS2_LOCATION/NtupleMacros/data/fakeRates-v7_2.root", "read"); 
	  if ( el_fakeRateFile_v7 == 0 ) {
	       std::cout << "$CMS2_LOCATION/NtupleMacros/data/fakeRates-v7_2.root could not be found!!" << std::endl;
	       std::cout << "Please make sure that $CMS2_LOCATION points to your CMS2 directory and that" << std::endl;
	       std::cout << "$CMS2_LOCATION/NtupleMacros/data/fakeRates-v7_2.root exists!" << std::endl;
	       gSystem->Exit(1);
	  }
	  el_fakeRate_v7 = dynamic_cast<TH2F *>(el_fakeRateFile_v7->Get("fakeRateTemplate_elt_EleFakes"));
	  el_fakeRate_err_v7 = dynamic_cast<TH2F *>(el_fakeRateFile_v7->Get("fakeRateTemplateError_elt_EleFakes"));
     }
     return *el_fakeRate_err_v7;
#else
     assert("use the bin errors in the fake rate histo instead of error histogram" && 0);
#endif
}

