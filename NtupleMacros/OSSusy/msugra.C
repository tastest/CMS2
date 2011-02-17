#include <algorithm>
#include <iostream>
#include <map>
#include <vector>
#include <sstream>
#include "TChain.h"
#include "TChainElement.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TProfile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TRandom3.h"
#include "TLegend.h"
#include "TLatex.h"

#include "cl95cms_landsberg.c"

using namespace std;

//-------------------------------------------
//parameters for CL95 function
//-------------------------------------------

const Double_t ilum            = 34.0;  // lumi
const Double_t slum            = 0.;    // lumi uncertainty (=0 b/c uncertainty is included in sig acceptance)
const Double_t eff             = 1.;    // sig efficiency
const Double_t bck             = 1.40;  // expected background
//const Double_t bck             = 0.00;  // expected background
const Double_t sbck            = 0.77;  // background error
const int      n               = 1;     // observed yield
const int      nuissanceModel  = 1;     // nuissance model (0 - Gaussian, 1 - lognormal, 2 - gamma)

//-------------------------------------------
// uncertainties
//-------------------------------------------

const float lumierr = 0.11;  // 11% lumi error
const float leperr  = 0.05;  // 5% error from lepton efficiency
const float pdferr  = 0.13;  // PDF uncertainty


//-------------------------------------------
//ABCD yields
//-------------------------------------------

const float yieldA  =  12.;
const float yieldB  =  37.;
const float yieldC  =   4.;


//calculate the expected UL? this is time-consuming so turn off if not needed
const bool calculateExpectedUL = true;

//rebin the TH2 yield histos (NOT RECOMMENDED)
const int   rebin   = 1;     // rebin yield histos





//const float nev     = 4.7;   // UL assuming 20% acceptance uncertainty
//const float kfact   = 1.4;   // approx k-factor
//const float fudge   = 1.;    // pb/fb conversion



TH1F* getCurve               ( TH2I *hist , char* name );
TGraphErrors* getCurve_TGraph( TH2I *hist , char* name );

void msugra( char* filename ){

  TFile *outfile = new TFile("exclusion/exclusion.root","RECREATE");

 //-----------------------------------
 // Here we load the yield histograms
 //-----------------------------------
 TFile *f = TFile::Open( filename );
 char* prefix = "LMscan10_";
 //char* prefix = "LMscan3_";

 TH2F* hyield     = (TH2F*) f->Get(Form("%slmgridyield",prefix));
 TH2F* hyield_k   = (TH2F*) f->Get(Form("%slmgridyield_k",prefix));
 TH2F* hyield_kup = (TH2F*) f->Get(Form("%slmgridyield_kup",prefix));
 TH2F* hyield_kdn = (TH2F*) f->Get(Form("%slmgridyield_kdn",prefix));
 TH2F* hyield_jup = (TH2F*) f->Get(Form("%slmgridyield_jup",prefix));
 TH2F* hyield_jdn = (TH2F*) f->Get(Form("%slmgridyield_jdn",prefix));

 TH2F* hyield_A = (TH2F*) f->Get(Form("%slmgridyield_A",prefix));
 TH2F* hyield_B = (TH2F*) f->Get(Form("%slmgridyield_B",prefix));
 TH2F* hyield_C = (TH2F*) f->Get(Form("%slmgridyield_C",prefix));
 TH2F* hyield_D = (TH2F*) f->Get(Form("%slmgridyield_D",prefix));

 TH2F* hUL_NLO     = (TH2F*) hyield_k->Clone("hUL_NLO");
 TH2F* hUL_NLO_exp = (TH2F*) hyield_k->Clone("hUL_NLO_exp");
 TH2F* hUL_LO      = (TH2F*) hyield_k->Clone("hUL_LO");

 TH1F* htotuncertainty = new TH1F("htotuncertainty","",200,0,2);

 hUL_NLO->Reset();
 hUL_LO->Reset();
 hUL_NLO_exp->Reset();

 if( rebin > 1 ){
   
   hyield->RebinX(rebin);
   hyield->RebinY(rebin);
   hyield->Scale( 1./(rebin*rebin) );

   hyield_k->RebinX(rebin);
   hyield_k->RebinY(rebin);
   hyield_k->Scale( 1./(rebin*rebin) );

   hUL_NLO->RebinX(rebin);
   hUL_NLO->RebinY(rebin);
   hUL_NLO->Scale( 1./(rebin*rebin) );

   hUL_NLO_exp->RebinX(rebin);
   hUL_NLO_exp->RebinY(rebin);
   hUL_NLO_exp->Scale( 1./(rebin*rebin) );

   hUL_LO->RebinX(rebin);
   hUL_LO->RebinY(rebin);
   hUL_LO->Scale( 1./(rebin*rebin) );

   hyield_kup->RebinX(rebin);
   hyield_kup->RebinY(rebin);
   hyield_kup->Scale( 1./(rebin*rebin) );

   hyield_kdn->RebinX(rebin);
   hyield_kdn->RebinY(rebin);
   hyield_kdn->Scale( 1./(rebin*rebin) );

   hyield_jup->RebinX(rebin);
   hyield_jup->RebinY(rebin);
   hyield_jup->Scale( 1./(rebin*rebin) );

   hyield_jdn->RebinX(rebin);
   hyield_jdn->RebinY(rebin);
   hyield_jdn->Scale( 1./(rebin*rebin) );

   hyield_A->RebinX(rebin);
   hyield_A->RebinY(rebin);
   hyield_A->Scale( 1./(rebin*rebin) );

   hyield_B->RebinX(rebin);
   hyield_B->RebinY(rebin);
   hyield_B->Scale( 1./(rebin*rebin) );

   hyield_C->RebinX(rebin);
   hyield_C->RebinY(rebin);
   hyield_C->Scale( 1./(rebin*rebin) );

   hyield_D->RebinX(rebin);
   hyield_D->RebinY(rebin);
   hyield_D->Scale( 1./(rebin*rebin) );

 }


 const unsigned int nm0bins  = hyield->GetXaxis()->GetNbins();
 const unsigned int nm12bins = hyield->GetYaxis()->GetNbins();

 //-----------------------------------------------
 //calculate acceptance error and UL at each scan point
 //-----------------------------------------------
 float accerr_NLO[nm0bins][nm12bins];
 Double_t ul_NLO[nm0bins][nm12bins];
 Double_t ul_NLO_SC[nm0bins][nm12bins];
 Double_t ul_NLO_nobkg[nm0bins][nm12bins];
 Double_t ul_NLO_exp[nm0bins][nm12bins];

 float accerr_LO[nm0bins][nm12bins];
 Double_t ul_LO[nm0bins][nm12bins];

 for( unsigned int m0bin = 1 ; m0bin <= nm0bins ; ++m0bin ){
   for( unsigned int m12bin = 1 ; m12bin <= nm12bins ; ++m12bin ){
     accerr_NLO[m0bin-1][m12bin-1] = 0.;
     accerr_LO[m0bin-1][m12bin-1]  = 0.;
     ul_NLO[m0bin-1][m12bin-1]     = 9999.;
     ul_NLO_SC[m0bin-1][m12bin-1]  = 9999.;
     ul_NLO_nobkg[m0bin-1][m12bin-1]  = 9999.;
     ul_NLO_exp[m0bin-1][m12bin-1] = 9999.;
     ul_LO[m0bin-1][m12bin-1]      = 9999.;
   }
 }

 //for( unsigned int m0bin = 1 ; m0bin <= 5 ; ++m0bin ){
 //for( unsigned int m12bin = 1 ; m12bin <= 5 ; ++m12bin ){
     
 for( unsigned int m0bin = 1  ; m0bin  <= nm0bins  ; ++m0bin ){
   for( unsigned int m12bin = 1 ; m12bin <= nm12bins ; ++m12bin ){

 //for( unsigned int m0bin = 18 ; m0bin <=18  ; ++m0bin ){
 //for( unsigned int m12bin = 12 ; m12bin <= 12 ; ++m12bin ){

 //for( unsigned int m0bin  = 5 ; m0bin  <= 7 ; ++m0bin ){
 //for( unsigned int m12bin = 5 ; m12bin <= 7 ; ++m12bin ){
     
     //get yields
     float yield     = hyield->GetBinContent( m0bin , m12bin );
     float yield_k   = hyield_k->GetBinContent( m0bin , m12bin );
     float yield_kup = hyield_kup->GetBinContent( m0bin , m12bin );
     float yield_kdn = hyield_kdn->GetBinContent( m0bin , m12bin );
     float yield_jup = hyield_jup->GetBinContent( m0bin , m12bin );
     float yield_jdn = hyield_jdn->GetBinContent( m0bin , m12bin );


//      if( fabs( hyield->GetXaxis()->GetBinCenter(m0bin) - 60 ) < 0.1 && fabs ( hyield->GetYaxis()->GetBinCenter(m12bin) -250 ) < 0.1 ){
//        cout << "Found LM1" << endl;
//      }
//      else{ continue; }


     cout << endl << endl << "------------------------------------------------------------------------" << endl;
     cout << endl << "m0 " << m0bin-1 << " m12 " << m12bin-1 << endl;
     cout << hyield->GetXaxis()->GetBinCenter(m0bin) << " " << hyield->GetYaxis()->GetBinCenter(m12bin) << endl;

     //skip bins with 0 yield
     if( fabs( yield_k ) < 1.e-10 ){
       cout << "zero yield, skipping!" << endl;
       continue; 
     }

//      if( hyield->GetXaxis()->GetBinCenter(m0bin) > 500 ){
//        cout << "m0 > 500 GeV, skipping" << endl;
//        continue;
//      }


     //-------------------------------------------------
     // this can save a lot of time, turned off here
     //-------------------------------------------------

     /*
     //a point with LO yield > 10 is definitely excluded
     if( yield > 10. ){
       ul_NLO[m0bin-1][m12bin-1]        = -1;
       ul_NLO_SC[m0bin-1][m12bin-1]     = -1;
       ul_NLO_nobkg[m0bin-1][m12bin-1]  = -1;
       ul_NLO_exp[m0bin-1][m12bin-1]    = -1;
       ul_LO[m0bin-1][m12bin-1]         = -1;
       hUL_NLO->SetBinContent     ( m0bin , m12bin , -1 );
       hUL_NLO_exp->SetBinContent ( m0bin , m12bin , -1 );
       hUL_LO->SetBinContent      ( m0bin , m12bin , -1 );
       cout << "yield " << yield << " point is excluded, skipping" << endl;
       continue;
     }

     //a point with NLO yield < 3 is definitely not excluded
     if( yield_k < 3. ){
       ul_NLO[m0bin-1][m12bin-1]        = 100;
       ul_NLO_SC[m0bin-1][m12bin-1]     = 100;
       ul_NLO_nobkg[m0bin-1][m12bin-1]  = 100;
       ul_NLO_exp[m0bin-1][m12bin-1]    = 100;
       ul_LO[m0bin-1][m12bin-1]         = 100;
       cout << "yield " << yield << " point is NOT excluded, skipping" << endl;
       hUL_NLO->SetBinContent     ( m0bin , m12bin , 100 );
       hUL_NLO_exp->SetBinContent ( m0bin , m12bin , 100 );
       hUL_LO->SetBinContent      ( m0bin , m12bin , 100 );
       continue;
     }
     */

     //uncertainty from k-factor
     float kerr_up   = fabs( ( yield_kup - yield_k ) / yield_k );
     float kerr_dn   = fabs( ( yield_kdn - yield_k ) / yield_k );
     float kerr      = TMath::Max( kerr_up , kerr_dn );
     
     //uncertainty from JES
     float jerr_up   = fabs( ( yield_jup - yield_k ) / yield_k );
     float jerr_dn   = fabs( ( yield_jdn - yield_k ) / yield_k );
     float jerr      = TMath::Max( jerr_up , jerr_dn );
     
     //add up NLO uncertainties (including k-factor uncertainty)
     float err2_NLO = 0.;
     err2_NLO += kerr * kerr;                                              //k-factor
     err2_NLO += jerr * jerr;                                              //JES
     err2_NLO += lumierr * lumierr;                                        //lumi
     err2_NLO += leperr * leperr;                                          //lep efficiency
     err2_NLO += pdferr * pdferr;                                          //PDF uncertainty
     accerr_NLO[m0bin-1][m12bin-1] = err2_NLO > 0 ? sqrt( err2_NLO ) : 0.; //total

     htotuncertainty->Fill( accerr_NLO[m0bin-1][m12bin-1] );

     
     //calculate observed NLO UL (including k-factor uncertainty)
     cout << "NLO UL: CL95( " << ilum << " , " << slum << " , " << eff << " , " << accerr_NLO[m0bin-1][m12bin-1] << " , " 
	  << bck << " , " << sbck << " , " << n << " , " << "false" << " , " << nuissanceModel << " )" << endl;
     ul_NLO[m0bin-1][m12bin-1] = ilum * CL95( ilum, slum, eff, accerr_NLO[m0bin-1][m12bin-1], bck, sbck, n, false, nuissanceModel );

     //calculate observed NLO UL (assume 0 bkg)
     cout << "NLO UL (no bkg): CL95( " << ilum << " , " << slum << " , " << eff << " , " << accerr_NLO[m0bin-1][m12bin-1] << " , " 
	  << 0.01 << " , " << 0.01 << " , " << n << " , " << "false" << " , " << nuissanceModel << " )" << endl;
     ul_NLO_nobkg[m0bin-1][m12bin-1] = ilum * CL95( ilum, slum, eff, accerr_NLO[m0bin-1][m12bin-1], 0.01, 0.01, n, false, nuissanceModel );

     /*
       //----------------------------------------------------------------
       //this is for signal contamination study using ABCD method
       //----------------------------------------------------------------


     //calculate observed NLO UL (with signal contamination)
     float A = yieldA - hyield_A->GetBinContent( m0bin , m12bin );
     float B = yieldB - hyield_B->GetBinContent( m0bin , m12bin );
     float C = yieldC - hyield_C->GetBinContent( m0bin , m12bin );
     float D = hyield_D->GetBinContent( m0bin , m12bin );

     cout << "Asusy " << hyield_A->GetBinContent( m0bin , m12bin ) << endl;
     cout << "Bsusy " << hyield_B->GetBinContent( m0bin , m12bin ) << endl;
     cout << "Csusy " << hyield_C->GetBinContent( m0bin , m12bin ) << endl;
     cout << "Dsusy " << hyield_D->GetBinContent( m0bin , m12bin ) << endl;

     cout << "A     " << A << endl;
     cout << "B     " << B << endl;
     cout << "C     " << C << endl;
     cout << "D     " << D << endl;

     if( fabs( D - yield_k ) > 1.e-10 ){
       cout << "yield " << yield_k << " region D " << D << endl;
     }

     float pred    = 0.;
     float prederr = 0.;

     if( A <= 0. || B <= 0. || C <= 0. ){
       pred    = 0.01;
       prederr = 0.01;
     }else{
       pred    = A * C / B;
       float prederr_stat = pred * sqrt( 1./A + 1./B + 1./C ); 
       float prederr_syst = pred * 0.2;
       prederr = sqrt( pow( prederr_stat , 2) + pow( prederr_syst , 2) );
     }

     cout << "pred  " << pred << " +/- " << prederr << endl;

     cout << "NLO UL (signal contamination): CL95( " << ilum << " , " << slum << " , " << eff << " , " << accerr_NLO[m0bin-1][m12bin-1] << " , " 
	  << pred << " , " << prederr << " , " << n << " , " << "false" << " , " << nuissanceModel << " )" << endl;
     ul_NLO_SC[m0bin-1][m12bin-1] = ilum * CL95( ilum, slum, eff, accerr_NLO[m0bin-1][m12bin-1], pred, prederr, n, false, nuissanceModel );

     cout << "UL with sig cont " << ul_NLO_SC[m0bin-1][m12bin-1] << endl;

     */

     //calculate expected NLO UL (including k-factor uncertainty)
     if( calculateExpectedUL ){
       if( accerr_NLO[m0bin-1][m12bin-1] > 0.5 ){
	 cout << "Large acceptance error! " << accerr_NLO[m0bin-1][m12bin-1] << endl;
	 cout << "Setting UL to 100!" << endl;
	 ul_NLO_exp[m0bin-1][m12bin-1] = 100.;
       }else{
	 cout << "Expected UL: CLA( " << ilum << " , " << slum << " , " << eff << " , " << accerr_NLO[m0bin-1][m12bin-1] 
	      << " , " << bck << " , " << sbck << " , " << nuissanceModel << " )" << endl;
	 ul_NLO_exp[m0bin-1][m12bin-1] = ilum * CLA(ilum, slum, eff, accerr_NLO[m0bin-1][m12bin-1], bck, sbck, nuissanceModel );
       }
     }
     else{
       ul_NLO_exp[m0bin-1][m12bin-1] = 100.;
     }
     

     //add up LO uncertainties (NOT including k-factor uncertainty)
     float err2_LO = 0.;
     err2_LO += jerr * jerr;                                               //JES
     err2_LO += lumierr * lumierr;                                         //lumi
     err2_LO += leperr * leperr;                                           //lep efficiency
     err2_LO += pdferr * pdferr;                                           //PDF uncertainty
     accerr_LO[m0bin-1][m12bin-1] = err2_LO > 0 ? sqrt( err2_LO ) : 0.;    //total

     /*     
     //calculate observed LO UL (NOT including k-factor uncertainty)
     cout << "LO UL CL95( " << ilum << " , " << slum << " , " << eff << " , " << accerr_LO[m0bin-1][m12bin-1] << " , " 
	  << bck << " , " << sbck << " , " << n << " , " << "false" << " , " << nuissanceModel << " )" << endl;
     ul_LO[m0bin-1][m12bin-1] = ilum * CL95( ilum, slum, eff, accerr_LO[m0bin-1][m12bin-1], bck, sbck, n, false, nuissanceModel );
     */

     //printout errors and UL
     cout << "yield            " << yield   << endl;
     cout << "yield * K        " << yield_k << endl;
     cout << "yield * Kup      " << yield_kup << endl;
     cout << "yield * Kdn      " << yield_kdn << endl;
     cout << "yield * JESup    " << yield_jup << endl;
     cout << "yield * JESdn    " << yield_jdn << endl;
     cout << "K error          " << kerr << endl;
     cout << "JES error        " << jerr << endl;
     cout << "total error NLO  " << accerr_NLO[m0bin-1][m12bin-1] << endl;
     cout << "NLO UL           " << ul_NLO[m0bin-1][m12bin-1] << endl;
     cout << "NLO UL sig cont  " << ul_NLO_SC[m0bin-1][m12bin-1] << endl;
     cout << "NLO UL no bkg    " << ul_NLO_nobkg[m0bin-1][m12bin-1] << endl;
     cout << "NLO UL exp       " << ul_NLO_exp[m0bin-1][m12bin-1] << endl;
     cout << "total error LO   " << accerr_LO[m0bin-1][m12bin-1] << endl;
     cout << "LO UL            " << ul_LO[m0bin-1][m12bin-1] << endl;
     
     hUL_NLO->SetBinContent( m0bin , m12bin , ul_NLO[m0bin-1][m12bin-1] );
     hUL_NLO_exp->SetBinContent( m0bin , m12bin , ul_NLO_exp[m0bin-1][m12bin-1] );
     hUL_LO->SetBinContent( m0bin , m12bin , ul_LO[m0bin-1][m12bin-1] );
   }
 }
 
 float xmin = hyield->GetXaxis()->GetXmin();
 float xmax = hyield->GetXaxis()->GetXmax();
 float ymin = hyield->GetYaxis()->GetXmin();
 float ymax = hyield->GetYaxis()->GetXmax();
 int nx     = hyield->GetXaxis()->GetNbins();
 int ny     = hyield->GetYaxis()->GetNbins();

 TH2I* hexcl_NLO_obs       = new TH2I("hexcl_NLO_obs",       "Observed NLO Exclusion",nx,xmin,xmax,ny,ymin,ymax);
 TH2I* hexcl_NLO_obs_SC    = new TH2I("hexcl_NLO_obs_SC",    "Observed NLO Exclusion (Sig Cont)",nx,xmin,xmax,ny,ymin,ymax);
 TH2I* hexcl_NLO_obs_nobkg = new TH2I("hexcl_NLO_obs_nobkg", "Observed NLO Exclusion (No Bkg)",nx,xmin,xmax,ny,ymin,ymax);
 TH2I* hexcl_LO_obs        = new TH2I("hexcl_LO_obs",        "Observed LO Exclusion", nx,xmin,xmax,ny,ymin,ymax);
 TH2I* hexcl_NLO_exp       = new TH2I("hexcl_NLO_exp",       "Expected NLO Exclusion",nx,xmin,xmax,ny,ymin,ymax);

 for( unsigned int m0bin  = 1 ; m0bin  <= nm0bins  ; ++m0bin  ){
   for( unsigned int m12bin = 1 ; m12bin <= nm12bins ; ++m12bin ){

     //for( unsigned int m0bin = 18 ; m0bin <= 18 ; ++m0bin ){
     //  for( unsigned int m12bin = 12 ; m12bin <= 12 ; ++m12bin ){
     
     //for( unsigned int m0bin = 5 ; m0bin <= 7 ; ++m0bin ){
     //for( unsigned int m12bin = 5 ; m12bin <= 7 ; ++m12bin ){
     
     hexcl_NLO_obs->SetBinContent( m0bin , m12bin , 0 );
     hexcl_NLO_obs_SC->SetBinContent( m0bin , m12bin , 0 );
     hexcl_NLO_obs_nobkg->SetBinContent( m0bin , m12bin , 0 );
     hexcl_LO_obs->SetBinContent( m0bin , m12bin , 0 );
     hexcl_NLO_exp->SetBinContent( m0bin , m12bin , 0 );

     //NLO observed
     int excluded_NLO_obs = 0;
     if( hyield_k->GetBinContent( m0bin , m12bin ) > ul_NLO[m0bin-1][m12bin-1] ) excluded_NLO_obs = 1;
     hexcl_NLO_obs->SetBinContent( m0bin , m12bin , excluded_NLO_obs );

     //NLO observed (sig cont)
     int excluded_NLO_obs_SC = 0;
     if( hyield_k->GetBinContent( m0bin , m12bin ) > ul_NLO_SC[m0bin-1][m12bin-1] ) excluded_NLO_obs_SC = 1;
     hexcl_NLO_obs_SC->SetBinContent( m0bin , m12bin , excluded_NLO_obs_SC );

     //NLO observed (nobkg)
     int excluded_NLO_obs_nobkg = 0;
     if( hyield_k->GetBinContent( m0bin , m12bin ) > ul_NLO_nobkg[m0bin-1][m12bin-1] ) excluded_NLO_obs_nobkg = 1;
     hexcl_NLO_obs_nobkg->SetBinContent( m0bin , m12bin , excluded_NLO_obs_nobkg );

     //NLO expected
     int excluded_NLO_exp = 0;
     if( hyield_k->GetBinContent( m0bin , m12bin ) > ul_NLO_exp[m0bin-1][m12bin-1] ) excluded_NLO_exp = 1;
     hexcl_NLO_exp->SetBinContent( m0bin , m12bin , excluded_NLO_exp );

     //LO observed
     int excluded_LO_obs = 0;
     if( hyield->GetBinContent( m0bin , m12bin ) > ul_LO[m0bin-1][m12bin-1] ) excluded_LO_obs = 1;
     hexcl_LO_obs->SetBinContent( m0bin , m12bin , excluded_LO_obs );     

     if( hyield_k->GetBinContent( m0bin , m12bin ) > 0. ){
     
       cout << endl << "m0 " << m0bin-1 << " m12 " << m12bin-1 << endl;
       cout << "NLO Yield     " << hyield_k->GetBinContent( m0bin , m12bin ) << endl;
       
       cout << "NLO UL        " << ul_NLO[m0bin-1][m12bin-1] << endl;
       cout << "Excluded?     " << excluded_NLO_obs << endl;
       
       cout << "NLO UL SC     " << ul_NLO_SC[m0bin-1][m12bin-1] << endl;
       cout << "Excluded?     " << excluded_NLO_obs_SC << endl;

       cout << "NLO UL no bkg " << ul_NLO_nobkg[m0bin-1][m12bin-1] << endl;
       cout << "Excluded?     " << excluded_NLO_obs_nobkg << endl;

       cout << "NLO UL exp    " << ul_NLO_exp[m0bin-1][m12bin-1] << endl;
       cout << "Excluded?     " << excluded_NLO_exp << endl;
       
       cout << "LO UL         " << ul_LO[m0bin-1][m12bin-1] << endl;
       cout << "Excluded?     " << excluded_LO_obs << endl;
       
     }
   }
 }

 TH1F*         limit_NLO_obs      = getCurve(        hexcl_NLO_obs , "limit_NLO_obs");
 TGraphErrors* limitgraph_NLO_obs = getCurve_TGraph( hexcl_NLO_obs , "limitgraph_NLO_obs");

 TH1F*         limit_NLO_obs_SC      = getCurve(        hexcl_NLO_obs_SC , "limit_NLO_obs_SC");
 TGraphErrors* limitgraph_NLO_obs_SC = getCurve_TGraph( hexcl_NLO_obs_SC , "limitgraph_NLO_obs_SC");

 TH1F*         limit_NLO_obs_nobkg      = getCurve(        hexcl_NLO_obs_nobkg , "limit_NLO_obs_nobkg");
 TGraphErrors* limitgraph_NLO_obs_nobkg = getCurve_TGraph( hexcl_NLO_obs_nobkg , "limitgraph_NLO_obs_nobkg");

 TH1F*         limit_NLO_exp      = getCurve(        hexcl_NLO_exp , "limit_NLO_exp");
 TGraphErrors* limitgraph_NLO_exp = getCurve_TGraph( hexcl_NLO_exp , "limitgraph_NLO_exp");

 TH1F*         limit_LO_obs       = getCurve(        hexcl_LO_obs , "limit_LO_obs");
 TGraphErrors* limitgraph_LO_obs  = getCurve_TGraph( hexcl_LO_obs , "limitgraph_LO_obs");




 //excluded points
 TCanvas *c1 = new TCanvas("c1","",1200,400);
 c1->Divide(3,1);

 TLatex *t = new TLatex();
 t->SetNDC();

 c1->cd(1);
 hexcl_NLO_obs->Draw("colz");
 t->DrawLatex(0.5,0.8,"observed NLO");

 c1->cd(2);
 hexcl_LO_obs->Draw("colz");
 t->DrawLatex(0.5,0.8,"observed LO");

 c1->cd(3);
 hexcl_NLO_exp->Draw("colz");
 t->DrawLatex(0.5,0.8,"expected NLO");
 
 c1->Print("exclusion/exclusion.png");

 //exclusion TH1
 TCanvas *c2 = new TCanvas("c2","",800,600);
 c2->cd();

 limit_NLO_obs->SetMinimum(90.);
 limit_NLO_obs->SetTitle("Exclusion Curve in mSUGRA Space");
 limit_NLO_obs->GetXaxis()->SetTitle("m_{0} (GeV)");
 limit_NLO_obs->GetYaxis()->SetTitle("m_{1/2} (GeV)");
 //limit_NLO_obs->SetLineWidth(2);
 //limit_NLO_obs->SetLineColor(2);
 limit_NLO_obs->GetXaxis()->SetRangeUser(0,500);
 limit_NLO_obs->GetYaxis()->SetRangeUser(100,400);
 //limit_NLO_obs->Draw("c");
 limit_NLO_obs->Draw();
 limit_LO_obs->SetLineColor(2);
 limit_NLO_exp->SetLineColor(4);
 limit_LO_obs->Draw("same");
 limit_NLO_exp->Draw("same");

 TLegend *leg = new TLegend(0.6,0.6,0.8,0.8);
 leg->AddEntry(limit_NLO_obs, "NLO obs","l");
 leg->AddEntry(limit_LO_obs,  "LO obs","l");
 leg->AddEntry(limit_NLO_exp, "NLO exp","l");
 leg->SetBorderSize(1);
 leg->SetFillColor(0);
 leg->Draw();

 c2->Print("exclusion/limit_hist.png");

 TCanvas *c3 = new TCanvas("c3","",800,600);
 c3->cd();

 limitgraph_LO_obs->SetMarkerColor(2);
 limitgraph_NLO_exp->SetMarkerColor(4);
 limitgraph_NLO_obs->Draw("AP");
 limitgraph_LO_obs->Draw("sameP");
 limitgraph_NLO_exp->Draw("sameP");
 leg->Draw();

 c3->Print("exclusion/limit_graph.png");

 outfile->cd();

 hexcl_NLO_obs->Write();
 limit_NLO_obs->Write();
 limitgraph_NLO_obs->Write();

 hexcl_NLO_obs_SC->Write();
 limit_NLO_obs_SC->Write();
 limitgraph_NLO_obs_SC->Write();

 hexcl_NLO_obs_nobkg->Write();
 limit_NLO_obs_nobkg->Write();
 limitgraph_NLO_obs_nobkg->Write();

 hexcl_LO_obs->Write();
 limit_LO_obs->Write();
 limitgraph_LO_obs->Write();

 hexcl_NLO_exp->Write();
 limit_NLO_exp->Write();
 limitgraph_NLO_exp->Write();

 hUL_NLO->Write();
 hUL_NLO_exp->Write();
 hUL_LO->Write();
 hyield_k->Write();

 htotuncertainty->Write();
 outfile->Close();
 
}

TH1F* getCurve( TH2I *hist , char* name){

 //--------------------------------------------
 // Get the parameters of the yield histogram
 //--------------------------------------------
 float xmin = hist->GetXaxis()->GetXmin();
 float xmax = hist->GetXaxis()->GetXmax();
 float ymin = hist->GetYaxis()->GetXmin();
 float ymax = hist->GetYaxis()->GetXmax();
 int nx     = hist->GetXaxis()->GetNbins();
 int ny     = hist->GetYaxis()->GetNbins();

 cout << "xmin = " << xmin <<endl;
 cout << "xmax = " << xmax <<endl;
 cout << "ymin = " << ymin <<endl;
 cout << "ymax = " << ymax <<endl;
 cout << "nx  = " << nx << endl;
 cout << "ny  = " << ny << endl;

 //---------------------------------------------------------------
 // We book a 1D histogram to keep the results in
 //---------------------------------------------------------------
 TH1F* limit = new TH1F(name,name,nx,xmin,xmax);
 //---------------------------------------------------------------
 // Use Sanjay's method, scan from the top and hope for the best
 //---------------------------------------------------------------
 float ybinsize = (ymax-ymin)/ny;
 float xbinsize = (xmax-xmin)/nx;
 for (int ix=1; ix<=nx; ix++) {
   float x = xmin + ix*xbinsize - 0.5*xbinsize;
   bool foundOne = false;

   for (int iy=ny; iy>0; iy--) {
     float this_ = hist->GetBinContent(ix,iy);
     //this_ = kfact*fudge*this_;
     //if (this_ > nev) {
     if (this_ > 0.5) {
       float yupperedge = ymin + iy*ybinsize;
       limit->Fill(x,yupperedge);
       // cout << yupperedge << " " << iy << endl;
       foundOne=true;
       break;
     }
   } //close iy loop
   //    if (!foundOne) limit->Fill(x,ymin);

 }   //close ix loop
 
 return limit;
}

TGraphErrors* getCurve_TGraph( TH2I *hist , char* name ){

 //--------------------------------------------
 // Get the parameters of the yield histogram
 //--------------------------------------------
 float xmin = hist->GetXaxis()->GetXmin();
 float xmax = hist->GetXaxis()->GetXmax();
 float ymin = hist->GetYaxis()->GetXmin();
 float ymax = hist->GetYaxis()->GetXmax();
 int nx     = hist->GetXaxis()->GetNbins();
 int ny     = hist->GetYaxis()->GetNbins();

 cout << "xmin = " << xmin <<endl;
 cout << "xmax = " << xmax <<endl;
 cout << "ymin = " << ymin <<endl;
 cout << "ymax = " << ymax <<endl;
 cout << "nx  = " << nx << endl;
 cout << "ny  = " << ny << endl;

 const unsigned int npoints = nx;
 float xpoint[npoints];
 float ypoint[npoints];
 float ex[npoints];
 float ey[npoints];

 for( unsigned int i = 0 ; i < npoints ; ++i ){
   xpoint[i] = 0.;
   ypoint[i] = 0.;
   ex[i]     = 0.;
   ey[i]     = 0.;
 }


 //---------------------------------------------------------------
 // Use Sanjay's method, scan from the top and hope for the best
 //---------------------------------------------------------------
 float ybinsize = (ymax-ymin)/ny;
 float xbinsize = (xmax-xmin)/nx;
 for (int ix=1; ix<=nx; ix++) {
   float x = xmin + ix*xbinsize - 0.5*xbinsize;
   bool foundOne = false;

   for (int iy=ny; iy>0; iy--) {
     float this_ = hist->GetBinContent(ix,iy);
     //this_ = kfact*fudge*this_;
     //if (this_ > nev) {
     if (this_ > 0.5) {
       float yupperedge = ymin + iy*ybinsize;
       // limit->Fill(x,yupperedge);
       xpoint[ix-1] = x;
       ypoint[ix-1] = yupperedge;
       // cout << yupperedge << " " << iy << endl;
       foundOne=true;
       break;
     }
   } //close iy loop
   if (!foundOne){
     xpoint[ix-1] = x;
     ypoint[ix-1] = 0;
   }
   //limit->Fill(x,ymin);
 }   //close ix loop
 
 //---------------------------------------------------------------
 // We book a TGraph to store the results
 //---------------------------------------------------------------
 TGraphErrors *gr = new TGraphErrors(npoints,xpoint,ypoint,ex,ey);
 gr->SetName(name);
 gr->SetTitle(name);

 return gr;
 
}
