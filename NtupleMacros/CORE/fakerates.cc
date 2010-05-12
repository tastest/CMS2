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
char* mu_filename;
static TH2F  *mu_fakeRate = 0;
static TH2F  *mu_ttbar_fakeRate = 0;

bool isFakeableMuon (int index, enum fakeRateVersion version) {

    if ( version == mu_v1 ) {
        return muonId(index, muonSelectionFO_mu_v1);
    } 
    else if ( version == mu_ttbar ) {
        return muonId(index, muonSelectionFO_mu_ttbar);
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

void SetMuonFile(char * filename){
    mu_filename = filename;
    cout << endl << "*Using muon fake rates from: " << mu_filename << endl;
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
        //mu_fakeRateFile = TFile::Open("$CMS2_LOCATION/NtupleMacros/data/mu_FR_3X.root", "read");
        mu_fakeRateFile = TFile::Open(mu_filename, "read");
        if ( mu_fakeRateFile == 0 ) {
            std::cout << "$CMS2_LOCATION/NtupleMacros/data/QCDFRplots-v1.root could not be found!!" << std::endl;
            std::cout << "Please make sure that $CMS2_LOCATION points to your CMS2 directory and that" << std::endl;
            std::cout << "$CMS2_LOCATION/NtupleMacros/data/QCDFRplots-v1.root exists!" << std::endl;
            gSystem->Exit(1);
        }
        mu_fakeRate       = dynamic_cast<TH2F *>( mu_fakeRateFile->Get("mu_FR_etavspt") );
        mu_ttbar_fakeRate = dynamic_cast<TH2F *>( mu_fakeRateFile->Get("mu_ttbar_FR_etavspt") );
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
char* el_filename;

// histograms for supported denominators
static TH2F  *el_fakeRate_ttbar_v1 = 0;
static TH2F  *el_fakeRate_ttbar_v2 = 0;
static TH2F  *el_fakeRate_ttbar_v3 = 0;

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

    // "V1"
    // remove ID, ISO, IP
    if(version == el_v1_cand01) {
        return pass_electronSelection(index, electronSelectionFO_el_v1_cand01);
    } 
    else if(version == el_v1_cand02) {
        return pass_electronSelection(index, electronSelectionFO_el_v1_cand02);
    } 
    else if(version == el_v1_cand02flip){
        return pass_electronSelection(index, electronSelectionFO_el_v1_cand02flip);
    } 

    // "V2"
    // remove ID, IP
    else if(version == el_v2_cand01) {
        return pass_electronSelection(index, electronSelectionFO_el_v2_cand01);    
    }
    else if(version == el_v2_cand02) {
        return pass_electronSelection(index, electronSelectionFO_el_v2_cand02);
    } 
    else if(version == el_v2_cand02flip) {
        return pass_electronSelection(index, electronSelectionFO_el_v2_cand02flip);
    }

    // "V3" 
    // remove ISO, IP
    else if(version == el_v3_cand01) { 
        return pass_electronSelection(index, electronSelectionFO_el_v3_cand01);
    } 
    else if(version == el_v3_cand02) {
        return pass_electronSelection(index, electronSelectionFO_el_v3_cand02);
    } 
    else if(version == el_v3_cand02flip) {
        return pass_electronSelection(index, electronSelectionFO_el_v3_cand02flip);
    } 

    // "TTBAR V1,2,3"
    else if( version == el_ttbar_v1 ) {
        return pass_electronSelection(index, electronSelectionFO_el_ttbar_v1);
    }
    else if( version == el_ttbar_v2 ) {
        return pass_electronSelection(index, electronSelectionFO_el_ttbar_v2);
    }
    else if( version == el_ttbar_v3 ) {
        return pass_electronSelection(index, electronSelectionFO_el_ttbar_v3);
    }

    // ERROR!
    else {
        std::cout<<"isFakeable: invalid fakeRateVersion given. Check it!"<<std::endl;
        return false;
    }
}

double elFakeProb (int i_el, enum fakeRateVersion version ){

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

void SetElectronFile(char * filename){
    el_filename = filename;
    cout << endl << "*Using electron fake rates from: " << el_filename << endl;
}

// read fake rate & errors from ROOT file in cvs NtupleMacros/data
TH2F &fakeRateEl (enum fakeRateVersion version ){
    if ( el_fakeRateFile == 0 ) {
        //el_fakeRateFile = TFile::Open("$CMS2_LOCATION/NtupleMacros/data/el_FR_3X.root", "read");
        //el_fakeRateFile = TFile::Open("/home/users/dbarge/FakeRates/CMSSW_3_5_6__CMS2_V03-03-23/el_FR_3X.root", "read");
        el_fakeRateFile = TFile::Open( el_filename, "read");
        if ( el_fakeRateFile == 0 ) {
            std::cout << "$CMS2_LOCATION/NtupleMacros/data/el_FR_3X.root could not be found!!" << std::endl;
            std::cout << "Please make sure that $CMS2_LOCATION points to your CMS2 directory and that" << std::endl;
            std::cout << "$CMS2_LOCATION/NtupleMacros/data/el_FR_3X.root exists!" << std::endl;
            gSystem->Exit(1);
        }
        el_fakeRate_v1_cand01     = dynamic_cast<TH2F *>( el_fakeRateFile->Get( "el_v1_cand01_Cleaned_FR_etavspt" ) );
        el_fakeRate_v1_cand02     = dynamic_cast<TH2F *>( el_fakeRateFile->Get( "el_v1_cand02_Cleaned_FR_etavspt" ) );
        el_fakeRate_v1_cand02flip = dynamic_cast<TH2F *>( el_fakeRateFile->Get( "el_v1_cand02flip_Cleaned_FR_etavspt" ) );
        el_fakeRate_v2_cand01     = dynamic_cast<TH2F *>( el_fakeRateFile->Get( "el_v2_cand01_Cleaned_FR_etavspt" ) );
        el_fakeRate_v2_cand02     = dynamic_cast<TH2F *>( el_fakeRateFile->Get( "el_v2_cand02_Cleaned_FR_etavspt" ) );
        el_fakeRate_v2_cand02flip = dynamic_cast<TH2F *>( el_fakeRateFile->Get( "el_v2_cand02flip_Cleaned_FR_etavspt" ) );
        el_fakeRate_v3_cand01     = dynamic_cast<TH2F *>( el_fakeRateFile->Get( "el_v3_cand01_Cleaned_FR_etavspt" ) );
        el_fakeRate_v3_cand02     = dynamic_cast<TH2F *>( el_fakeRateFile->Get( "el_v3_cand02_Cleaned_FR_etavspt" ) );
        el_fakeRate_v3_cand02flip = dynamic_cast<TH2F *>( el_fakeRateFile->Get( "el_v3_cand02flip_Cleaned_FR_etavspt" ) );
        el_fakeRate_ttbar_v1 = dynamic_cast<TH2F *>( el_fakeRateFile->Get( "el_ttbar_v1_Cleaned_FR_etavspt" ) );
        el_fakeRate_ttbar_v2 = dynamic_cast<TH2F *>( el_fakeRateFile->Get( "el_ttbar_v2_Cleaned_FR_etavspt" ) );
        el_fakeRate_ttbar_v3 = dynamic_cast<TH2F *>( el_fakeRateFile->Get( "el_ttbar_v3_Cleaned_FR_etavspt" ) );
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
    } 
    else if( version == el_ttbar_v1 ){ 
        return *el_fakeRate_ttbar_v1;
    } 
    else if( version == el_ttbar_v2 ){ 
        return *el_fakeRate_ttbar_v2;
    } 
    else if( version == el_ttbar_v3 ){ 
        return *el_fakeRate_ttbar_v3;
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

