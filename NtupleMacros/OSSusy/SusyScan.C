#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <iostream>
#include <iomanip>
#include <math.h>
#include <string>
#include <fstream>
#include "TLegendEntry.h"
#include "TMath.h"
#include "TStyle.h"
#include "zn.cc"
//#include "histtools.h"

bool debug = false;
const double pi = acos(-1.);
bool writefile = true; //write output root file

//---mSUGRA scan parameters---

// const int   nm0points    = 81;
// const float m0min        = 0.;
// const float m0max        = 4050.;
// const int   nm12points   = 26;
// const float m12min       = 100.;
// const float m12max       = 620.;

const int   nm0points    = 21;
const float m0min        = 0.;
const float m0max        = 1050.;
const int   nm12points   = 16;
const float m12min       = 100.;
const float m12max       = 420.;

//---Values for 100/pb---
// float bkg_syst_error    = 0.25;
// float upper_limit_100   = 62.3; //<<---- WRONG BUT NOT USED
// float upper_limit_175   = 11.627; 
// float zbi_ul_175        = 2.47;
// char* suffix            = "100pb";
// float lumi_scale_factor = 1.; 
// float ABCD_zbi_UL       = 2.33;

//---Values for 1/fb---
float bkg_syst_error    = 0.25;
//float upper_limit_100   = 62.3;  //<<---- WRONG BUT NOT USED
//float upper_limit_175   = 73.121; 
float zbi_ul_175        = 2.17;
char* suffix            = "1fb";
float lumi_scale_factor = 10.; 
float ABCD_zbi_UL       = 1.87;

//---Zoom in m12:m0 plot---
float xaxismin = 0.;
float xaxismax = 1000.;
float yaxismin = 100.;
float yaxismax = 400.;

//---Full range m12:m0 plot---
//float xaxismin = 0.;
//float xaxismax = 4000.;
//float yaxismin = 100.;
//float yaxismax = 600.;

TH2F* hpred100;
TH2F* hobs100; 
TH2F* hzbi100; 
TH2F* hzn100; 
TH2F* hexc100; 
TH2F* hpred175;
TH2F* hobs175; 
TH2F* hzbi175; 
TH2F* hzn175; 
TH2F* hexc175; 
TH2F* h_ABCD_pred;
TH2F* h_ABCD_obs;
TH2F* h_ABCD_zbi;
TH2F* h_ABCD_zn;
TH2F* h_ABCD_exc;

float getM0FromIndex(int index){

  float binsize = (m0max - m0min) / (float) nm0points;
  float m0      = index * binsize + m0min;
  return m0;

}

float getM12FromIndex(int index){

  float binsize = (m12max - m12min) / (float) nm12points;
  float m12      = index * binsize + m12min;
  return m12;

}

namespace hist {
  void add(const char* outHistName, const char* patORpfx0);
  TLegend* legend(TCanvas* canvas, Option_t* option = "lpf", Bool_t addColor = kFALSE, Int_t token = -1, 
		  Float_t xmin = 0.75, Float_t ymin = 0.75, Float_t xmax = 0.99, Float_t ymax = 0.99);
}

void loadHist(const char* filename, const char* directory = 0, const char* pfx = 0, 
	      const char* pat = "*", Bool_t doAdd = kFALSE);

void BookHists();

//double getZBi(double n_on , double n_off , double tau);
double getZBi(double n_on , double mu_b_hat , double sigma_b);

void drawPlot(TH2F* h , TCanvas *c , TH2F* hscan, bool printgif = false);

int getXBin(TH2F* h, float xval);
int getYBin(TH2F* h, float yval);

void SusyScan(bool printgif = false) {

  char* smfilename   = "root/victory_baseline_metgt50_sumjetptgt200_cand01_calo_tcmet_all_3x.root";
  //char* smfilename   = "/home/jribnik/devel/tas/CMS2/NtupleMacros/OSSusy/victory_baseline_metgt50_sumjetptgt200_cand01_calo_tcmet_3x.root";
  char* susyfilename = "root/victory_baseline_metgt50_sumjetptgt200_cand01_calo_tcmet_LMscan_3x.root";
  //char* susyfilename = "/tas03/home/jribnik/susyscan/OSSusy/root/victory_baseline_metgt50_sumjetptgt200_cand01_calo_tcmet_3x.root";

  cout << endl;
  cout << "Using SM file     : " << smfilename   << endl;
  cout << "Using SUSY file   : " << susyfilename << endl;
  cout << endl;

  //TFile *fsm = TFile::Open(smfilename);
  TFile *fsusy = TFile::Open(susyfilename);
  //TFile *fout;

  gROOT->SetStyle("Plain");
  gROOT->ProcessLine("gStyle->SetOptStat(0)");
  ofstream ofile(Form("SusyScan_%s.txt",suffix));
  ofstream ofile_ABCD(Form("SusyScan_ABCD_%s.txt",suffix));

  //Load and add histos--------------------------------------------
  loadHist(smfilename,   0, 0, "*dilPt_allj*",              kTRUE);
  loadHist(smfilename,   0, 0, "*dilPtSmeared_allj*",       kTRUE);
  loadHist(smfilename,   0, 0, "*tcmet_allj*",              kTRUE);
  loadHist(smfilename,   0, 0, "*sumJetPt_tcmetsqrtsumet*", kTRUE);
  loadHist(susyfilename, 0, 0, "*susy_hdilPt*",             kTRUE);
  loadHist(susyfilename, 0, 0, "*susy_htcmet*",             kTRUE);
  loadHist(susyfilename, 0, 0, "*susy_hmet_sumjetpt*",      kTRUE);

  std::cout << "adding histograms..." << std::endl;
  hist::add("sm_dilPt",        "^[^L].*_hdilPt_allj_all$");
  hist::add("sm_dilPtSmeared", "^[^L].*_hdilPtSmeared_allj_all$");  
  hist::add("sm_tcmet",        "^[^L].*_htcmet_allj_all$");  
  hist::add("sm_ABCD",         "^[^L].*_sumJetPt_tcmetsqrtsumet_allj_all$");  
  
  TH1F* sm_dilPt =        (TH1F*)gDirectory->Get("sm_dilPt");
  TH1F* sm_dilPtSmeared = (TH1F*)gDirectory->Get("sm_dilPtSmeared");
  TH1F* sm_tcmet =        (TH1F*)gDirectory->Get("sm_tcmet");
  TH1F* sm_ABCD =         (TH1F*)gDirectory->Get("sm_ABCD");
  
  if(sm_dilPt == 0){
    cout<<"Error unable to find histos sm_dilPt_allj_all"<<endl;
    return;
  }
  if(sm_dilPtSmeared == 0){
    cout<<"Error unable to find histos sm_dilPtSmeared_allj_all"<<endl;
    return;
  }
  if(sm_tcmet == 0){
    cout<<"Error unable to find histos sm_tcmet_allj_all"<<endl;
    return;
  }
  if(sm_ABCD == 0){
    cout<<"Error unable to find SM ABCD histos"<<endl;
    return;
  }

  if(debug) cout << "Loaded all histos" << endl;
  //   if(writefile){
  //     fout = new TFile(Form("SusyScan_%s.root",suffix),"RECREATE");
  //     fout->cd();
  //   }

  BookHists();
  int width = 7;
  
  ofile  << "|" << setw(width) << "m0"        << setw(width) << setw(4)  
         << "|" << setw(width) << "m12"       << setw(width) << setw(4)   
         << "|" << setw(width) << "pred100"   << setw(width)  
         << "|" << setw(width) << "obs100"    << setw(width)  
         << "|" << setw(width) << "ZBi100"    << setw(width)  
         << "|" << setw(width) << "ZN100"     << setw(width)  
         << "|" << setw(width) << "excl?"     << setw(width)  
         << "|" << setw(width) << "pred175"   << setw(width)  
         << "|" << setw(width) << "obs175"    << setw(width)  
         << "|" << setw(width) << "ZBi175"    << setw(width)  
         << "|" << setw(width) << "ZN175"     << setw(width)  
         << "|" << setw(width) << "excl?"     << setw(width) << "|" << endl;
  
  ofile_ABCD  << "|" << setw(width) << "m0"        << setw(width) << setw(4)  
              << "|" << setw(width) << "m12"       << setw(width) << setw(4)   
              << "|" << setw(width) << "pred"      << setw(width)  
              << "|" << setw(width) << "obs"       << setw(width)  
              << "|" << setw(width) << "ZBi"       << setw(width)  
              << "|" << setw(width) << "ZN"        << setw(width)  
              << "|" << setw(width) << "excl?"     << setw(width) << "|" << endl;  
           
  //loop over susy scan points
  if(debug) cout << "Begin loop over mSUGRA points" << endl;

  for(int im0 = 0 ; im0 < nm0points ; im0++){
    
    for(int im12 = 0 ; im12 < nm12points ; im12++){
  
      //---------------------------------------------------------------------------------
      //VICTORY METHOD
      //---------------------------------------------------------------------------------
      
      TH1F* hdilPt = (TH1F*)sm_dilPt->Clone(Form("hdilPt_m0_%i_m12_%i",im0,im12));
      TH1F* htcmet = (TH1F*)sm_tcmet->Clone(Form("htcmet_m0_%i_m12_%i",im0,im12));

      TH1F* hsusydilPt   = (TH1F*) fsusy->Get(Form("susy_hdilPt_m0_%i_m12_%i",im0,im12));
      TH1F* hsusytcmet   = (TH1F*) fsusy->Get(Form("susy_htcmet_m0_%i_m12_%i",im0,im12));
      
      hdilPt->Add(hsusydilPt);
      htcmet->Add(hsusytcmet);
      
      if(lumi_scale_factor > 1){
        hdilPt->    Scale( lumi_scale_factor );
        htcmet->    Scale( lumi_scale_factor );
        hsusydilPt->Scale( lumi_scale_factor );
        hsusytcmet->Scale( lumi_scale_factor );
      }

      //get predicted yield
      int bin50  = hdilPt->FindBin(50);
      int bin100 = hdilPt->FindBin(100);
      int bin175 = hdilPt->FindBin(175);
      double scale = hdilPt->Integral(bin50, 101) / hdilPt->Integral(0, 101);
      hdilPt->Scale( 1. / scale );
      double dilpt_prediction100 = hdilPt->Integral(bin100, 101);
      double dilpt_prediction175 = hdilPt->Integral(bin175, 101);

      //get observed yield
      bin100 = htcmet->FindBin(100);
      bin175 = htcmet->FindBin(175);
      
      double met_observed100 = htcmet->Integral(bin100, 101);
      double met_observed175 = htcmet->Integral(bin175, 101);

      //get observed mSUGRA yield
      //double susy_met_observed100 = hsusytcmet->Integral(bin100, 101);
      double susy_met_observed175 = hsusytcmet->Integral(bin175, 101);

      //get sparms
      float m0  = getM0FromIndex ( im0  );
      float m12 = getM12FromIndex( im12 );

      //get Z-values
      float ZBi100 = getZBi( met_observed100 , dilpt_prediction100 , bkg_syst_error * dilpt_prediction100  );
      float ZBi175 = getZBi( met_observed175 , dilpt_prediction175 , bkg_syst_error * dilpt_prediction175  );

      //float ZBi100 = getZBi( met_observed100 , dilpt_prediction100 * scale , scale );
      //float ZBi175 = getZBi( met_observed175 , dilpt_prediction175 * scale , scale );

      float ZN100 = getzn( met_observed100 - dilpt_prediction100 , dilpt_prediction100 , 0 , bkg_syst_error * 100);
      float ZN175 = getzn( met_observed175 - dilpt_prediction175 , dilpt_prediction175 , 0 , bkg_syst_error * 100);
      
      //determine if point is excluded
      //float excl100 = met_observed100 > upper_limit_100 ? 1 : 0;
      //float excl175 = met_observed175 > upper_limit_175 ? 1 : 0;

      //currently we don't perform exclusion for met > 100 GeV cut
      float excl100 = 0;

      //exclude point if ZBi significance (for met > 175 GeV) exceeds zbi_ul_175
      if( ZBi175 != ZBi175 )
        cout << "Found a nan!!! m0 " << m0 << " m12 " << m12 << endl;

      float excl175 = ZBi175 > zbi_ul_175 ? 1 : 0;
      
      //print text to file
      ofile << "|" 
            << setw(width) << setprecision(3) << m0 << setw(width) << "|" 
            << setw(width) << m12 << setw(width) << "|" 
            << setw(width) << dilpt_prediction100 << setw(width) << "|" 
            << setw(width) << met_observed100 << setw(width) << "|"
            << setw(width) << ZBi100 << setw(width) << "|" 
            << setw(width) << ZN100 << setw(width) << "|" 
            << setw(width) << excl100 << setw(width) << "|"
            << setw(width) << dilpt_prediction175 << setw(width) << "|" 
            << setw(width) << met_observed175 << setw(width) << "|"
            << setw(width) << ZBi175 << setw(width) << "|" 
            << setw(width) << ZN175 << setw(width) << "|" 
            << setw(width) << excl175 << setw(width) << "|" << endl;
    
      //fill histos
      //hpred100 -> Fill(m0,m12, dilpt_prediction100);
      //hobs100  -> Fill(m0,m12, met_observed100);
      //hzbi100  -> Fill(m0,m12, ZBi100);
      //hzn100   -> Fill(m0,m12, ZN100);
      //hexc100  -> Fill(m0,m12, excl100);
      hpred175 -> Fill(m0,m12, dilpt_prediction175);
      //hobs175  -> Fill(m0,m12, met_observed175);
      hobs175  -> Fill(m0,m12, susy_met_observed175);  //susy observed = total observed - SM observed
      hzbi175  -> Fill(m0,m12, ZBi175);
      hzn175   -> Fill(m0,m12, ZN175);
      hexc175  -> Fill(m0,m12, excl175);

      //---------------------------------------------------------------------------------
      //ABCD METHOD
      //---------------------------------------------------------------------------------
      if(debug) cout << "Do ABCD method m0 " << im0 << " m12 " << im12 << endl;
      
      TH2F* h_ABCD        = (TH2F*) sm_ABCD->Clone(Form("sm_ABCD_m0_%i_m12_%i",im0,im12));
      TH2F* h_ABCD_SUSY   = (TH2F*) fsusy->Get(Form("susy_hmet_sumjetpt_m0_%i_m12_%i",im0,im12));
    
      if( h_ABCD_SUSY == 0 ){
        cout << "ERROR CAN'T GET SUSY ABCD HISTOS" << endl;
        exit(0);
      }

      h_ABCD->Add(h_ABCD_SUSY);
      
      if(lumi_scale_factor > 1){
        h_ABCD->        Scale( lumi_scale_factor );
        h_ABCD_SUSY->   Scale( lumi_scale_factor );
      }
      
      //set ABCD regions
      float x1=100;
      float x2=250;
      float x3=250;
      float x4=1000;
      
      float y1 = 4;
      float y2 = 7;
      float y3 = 7;
      float y4 = 25;

      //get bins for 4 ABCD regions
      int ix1 = getXBin( h_ABCD , x1 );
      int ix2 = getXBin( h_ABCD , x2 );
      int ix3 = getXBin( h_ABCD , x3 );
      int ix4 = getXBin( h_ABCD , x4 );
      
      int iy1 = getYBin( h_ABCD , y1 );
      int iy2 = getYBin( h_ABCD , y2 );
      int iy3 = getYBin( h_ABCD , y3 );
      int iy4 = getYBin( h_ABCD , y4 );
      
      //get total yields in each region
      float A = h_ABCD->Integral( ix1+1, ix2, iy3+1, iy4 );
      float B = h_ABCD->Integral( ix1+1, ix2, iy1+1, iy2 );
      float C = h_ABCD->Integral( ix3+1, ix4, iy1+1, iy2 );
      float D = h_ABCD->Integral( ix3+1, ix4, iy3+1, iy4 );

      //get SUSY yields in each region
      float A_SUSY = h_ABCD->Integral( ix1+1, ix2, iy3+1, iy4 );
      float B_SUSY = h_ABCD->Integral( ix1+1, ix2, iy1+1, iy2 );
      float C_SUSY = h_ABCD->Integral( ix3+1, ix4, iy1+1, iy2 );
      float D_SUSY = h_ABCD->Integral( ix3+1, ix4, iy3+1, iy4 );

      if(debug) cout << "A B C D " << A << " " << B << " " << C << " " << D << endl;

      //get predicted/observed yields and significances
      float ABCD_pred = 1.5*A*C/B;
      float ABCD_obs  = D;
      //float ABCD_zbi  = getZBi(obs, pred, sqrt(pred + pow( syst_bkg_error * pred ,2) ) );
      float ABCD_zbi  = getZBi(ABCD_obs, ABCD_pred, bkg_syst_error * ABCD_pred );
      float ABCD_zn   = getzn( ABCD_obs - ABCD_pred, ABCD_pred , 0 , bkg_syst_error * 100);

      //exclude point if ZBi significance (for met > 175 GeV) exceeds zbi_ul_175
      if( ABCD_zbi != ABCD_zbi)
        cout << "ABCD method: found a nan!!! m0 " << m0 << " m12 " << m12 << endl;

      float ABCD_exc = ABCD_zbi > ABCD_zbi_UL ? 1 : 0;
      
      //print text to file
      ofile_ABCD << "|" 
                 << setw(width) << setprecision(3) << m0 << setw(width) << "|" 
                 << setw(width) << m12 << setw(width) << "|" 
                 << setw(width) << ABCD_pred << setw(width) << "|" 
                 << setw(width) << ABCD_obs << setw(width) << "|"
                 << setw(width) << ABCD_zbi << setw(width) << "|" 
                 << setw(width) << ABCD_zn << setw(width) << "|" 
                 << setw(width) << ABCD_exc << setw(width) << "|" << endl;
            
    
      //fill histos
      h_ABCD_pred -> Fill(m0,m12, ABCD_pred);
      //hobs175  -> Fill(m0,m12,  D);       //total observed
      h_ABCD_obs  -> Fill(m0,m12, D_SUSY);  //susy observed = total observed - SM observed
      h_ABCD_zbi  -> Fill(m0,m12, ABCD_zbi);
      h_ABCD_zn   -> Fill(m0,m12, ABCD_zn);
      h_ABCD_exc  -> Fill(m0,m12, ABCD_exc);








      if(debug) cout << "End of mSUGRA loop" << endl;
    }
      
  }
  
  if(debug) cout << "Drawing histos" << endl;

  //TH2 hscan is used to mask missing mSUGRA points 
  //TFile *fscan = TFile::Open("susyscan.root");  
  //TH2F*  hscan = (TH2F*) fscan->Get("h");

  TH2F* hscan;

  //TCanvas *cscan = new TCanvas("cscan","",800,600);
  //cscan->cd();
  //hscan->Draw("contz"); 

  const int ncan = 5;
  TCanvas *can[ncan];

  //draw victory method plots
  //drawPlot(hpred100 , can[0] , printgif);
  //drawPlot(hobs100  , can[1] , printgif);
  //drawPlot(hzbi100  , can[2] , printgif);
  //drawPlot(hzn100   , can[3] , printgif);
  //drawPlot(hexc100  , can[4] , printgif);

  //drawPlot(hpred175 , can[0] , hscan , printgif);
  //drawPlot(hobs175  , can[1] , hscan , printgif);
  //drawPlot(hzbi175  , can[2] , hscan , printgif);
  //drawPlot(hzn175   , can[3] , hscan , printgif);
  drawPlot(hexc175  , can[4] , hscan , printgif);

  //draw ABCD plots
  //drawPlot(h_ABCD_pred  , can[5] , hscan , printgif);
  //drawPlot(h_ABCD_obs   , can[6] , hscan , printgif);
  //drawPlot(h_ABCD_zbi   , can[7] , hscan , printgif);
  //drawPlot(h_ABCD_zn    , can[8] , hscan , printgif);
  drawPlot(h_ABCD_exc   , can[9] , hscan , printgif);

  ofile.close();

  //if(writefile){
  //fout->cd();
  //fout->Write();
  //fout->Close();
  //}
  
}

void drawPlot(TH2F* h , TCanvas *c , TH2F* hscan, bool printgif){
  c = new TCanvas(Form("%s_can",h->GetName()),Form("%s_can",h->GetName()),800,600);
  c -> Divide(2,1);

  char* drawmode = "colz";  

  /*
  if(strcmp(h->GetName(),"hexc175") == 0){
    h->SetContour(1);
    h->SetContourLevel(0, 0.5);
    drawmode= "cont4";
    h->SetMinimum(0.);
    h->SetMaximum(1.);

    //mask outliers
    if(strcmp(suffix,"1fb")==0){
      h->SetBinContent(13,1,0);
      h->SetBinContent(11,1,0);
      h->SetBinContent(11,2,0);
    }
  }
  if(strcmp(h->GetName(),"hzbi175") == 0){
    h->SetContour(2);
    h->SetContourLevel(0, 3.);
    h->SetContourLevel(1, 5.);
    drawmode= "cont4";
    h->SetMinimum(0.);
    h->SetMaximum(5.);

    //mask outliers
    if(strcmp(suffix,"1fb")==0){
      h->SetBinContent(11,2,0);
      h->SetBinContent(8,1,0);
      h->SetBinContent(8,6,0);

    }
  }
  if(strcmp(h->GetName(),"hzn175") == 0){
    h->SetContour(2);
    h->SetContourLevel(0, 3.);
    h->SetContourLevel(1, 5.);
    drawmode= "cont4";
    h->SetMinimum(0.);
    h->SetMaximum(5.);

    //mask outliers
    if(strcmp(suffix,"100pb")==0){
      h->SetBinContent(7,2,0);
    }

    //mask outliers
    if(strcmp(suffix,"1fb")==0){
      h->SetBinContent(13,1,0);
      h->SetBinContent(11,1,0);
      h->SetBinContent(11,2,0);
    }
  }
  if(strcmp(h->GetName(),"hobs175") == 0){
    h->SetMinimum(0.);
    //h->SetMaximum(5.);
  }
  if(strcmp(h->GetName(),"hpred175") == 0){
    h->SetMinimum(0.);
    //h->SetMaximum(5.);
  }
  */

  h->Draw(drawmode);
  //h->Draw("col");
  //h->Draw("cont3 same");
  h->GetXaxis()->SetRangeUser( xaxismin , xaxismax );
  h->GetYaxis()->SetRangeUser( yaxismin , yaxismax );


//---Draw boxes over missing mSUGRA points---
//   TBox box;
//   box.SetFillColor(11);
//   float xmin;
//   float xmax;
//   float ymin;
//   float ymax;

//   for(int ibinx = 0 ; ibinx < 20 ; ibinx++){
//     for(int ibiny = 0 ; ibiny < 15 ; ibiny++){
      
//       if(hscan->GetBinContent( ibinx+1 , ibiny+1 )==1){
       
//         xmin = ibinx*50;
//         xmax = (ibinx+1)*50;
//         ymin = ibiny*20 + 100;
//         ymax = (ibiny+1)*20 + 100;
//         box.DrawBox(xmin,ymin,xmax,ymax);
//         //cout<<"ibinx "<<ibinx<<" ibiny "<<ibiny<<endl;
//         //cout<<xmin<<" "<<xmax<<" "<<ymin<<" "<<ymax<<endl;
//       }
//     }
//   }

  if(printgif){
    c -> Modified();
    c -> Update();
    c -> Print(Form("plots/%s_%s.gif" , h->GetName() , suffix));
    c -> Print(Form("plots/%s_%s.pdf" , h->GetName() , suffix));
    c -> Print(Form("plots/%s_%s.eps" , h->GetName() , suffix));
    c -> Print(Form("plots/%s_%s.C"   , h->GetName() , suffix));
  }
}


double getZBi(double n_on , double mu_b_hat , double sigma_b){

  // total events in signal region (S+B)
  //double n_on     = 140.;
  // mean num of BG events expected in sig. region
  //double mu_b_hat = 83.33;
  // uncertainty of mu_b_hat
  //double sigma_b  = 8.333;       
  // scale factor to corresp. Noff/Non

  double tau      = mu_b_hat / (sigma_b*sigma_b);
  double n_off    = tau*mu_b_hat;
  double P_Bi     = TMath::BetaIncomplete(1./(1.+tau), n_on, n_off+1);
  //double Z_Bi     = sqrt(2.)*TMath::ErfInverse(1 - 2.*P_Bi);  
  
  if(P_Bi == 0){
    return 100; // set ZBi = 100 if P_Bi = 0
  }
  else if(P_Bi < 1.e-10){
    double mu = -2*log(P_Bi * sqrt( 2*pi ));
    return sqrt( mu - log(mu) );
  }
  else{
    return sqrt(2.)*TMath::ErfInverse(1 - 2.*P_Bi);
  }
 
}



void BookHists(){
  
  //book histos for predicted/observed yields and significance
  int   nm0points_    = 81;
  float m0min_        = -25.;
  float m0max_        = 4025.;
  int   nm12points_   = 26;
  float m12min_       = 90.;
  float m12max_       = 610.;

  //victory method histos
  hpred100 = new TH2F("hpred100" , "Predicted Yields (met>100)" , nm0points_,m0min_,m0max_,nm12points_,m12min_,m12max_);
  hobs100  = new TH2F("hobs100"  , "Observed Yields (met>100)"  , nm0points_,m0min_,m0max_,nm12points_,m12min_,m12max_);
  hzbi100  = new TH2F("hzbi100"  , "ZBI (met>100)"              , nm0points_,m0min_,m0max_,nm12points_,m12min_,m12max_);
  hzn100   = new TH2F("hzn100"   , "ZN (met>100)"               , nm0points_,m0min_,m0max_,nm12points_,m12min_,m12max_);
  hexc100  = new TH2F("hexc100"  , "Exclusion (met>100)"        , nm0points_,m0min_,m0max_,nm12points_,m12min_,m12max_);
  hpred175 = new TH2F("hpred175" , "Predicted Yields (met>175)" , nm0points_,m0min_,m0max_,nm12points_,m12min_,m12max_);
  hobs175  = new TH2F("hobs175"  , "Observed Yields (met>175)"  , nm0points_,m0min_,m0max_,nm12points_,m12min_,m12max_);
  hzbi175  = new TH2F("hzbi175"  , "ZBI (met>175)"              , nm0points_,m0min_,m0max_,nm12points_,m12min_,m12max_);
  hzn175   = new TH2F("hzn175"   , "ZN (met>175)"               , nm0points_,m0min_,m0max_,nm12points_,m12min_,m12max_);
  hexc175  = new TH2F("hexc175"  , "Exclusion (met>175)"        , nm0points_,m0min_,m0max_,nm12points_,m12min_,m12max_);  

  //ABCD method histos
  h_ABCD_pred = new TH2F("h_ABCD_pred" , "Predicted Yields (ABCD)" , nm0points_,m0min_,m0max_,nm12points_,m12min_,m12max_);
  h_ABCD_obs  = new TH2F("h_ABCD_obs"  , "Observed Yields (ABCD)"  , nm0points_,m0min_,m0max_,nm12points_,m12min_,m12max_);
  h_ABCD_zbi  = new TH2F("h_ABCD_zbi"  , "Z_{Bi} (ABCD)"           , nm0points_,m0min_,m0max_,nm12points_,m12min_,m12max_);
  h_ABCD_zn   = new TH2F("h_ABCD_zn"   , "Z_{N} (ABCD)"            , nm0points_,m0min_,m0max_,nm12points_,m12min_,m12max_);
  h_ABCD_exc  = new TH2F("h_ABCD_exc"  , "Exclusion (ABCD)"        , nm0points_,m0min_,m0max_,nm12points_,m12min_,m12max_);  

  hpred100   -> GetXaxis()->SetTitle("m_{0} (GeV)");
  hobs100    -> GetXaxis()->SetTitle("m_{0} (GeV)");
  hzbi100    -> GetXaxis()->SetTitle("m_{0} (GeV)");
  hzn100     -> GetXaxis()->SetTitle("m_{0} (GeV)");
  hexc100    -> GetXaxis()->SetTitle("m_{0} (GeV)");
  hpred175   -> GetXaxis()->SetTitle("m_{0} (GeV)");
  hobs175    -> GetXaxis()->SetTitle("m_{0} (GeV)");
  hzbi175    -> GetXaxis()->SetTitle("m_{0} (GeV)");
  hzn175     -> GetXaxis()->SetTitle("m_{0} (GeV)");
  hexc175    -> GetXaxis()->SetTitle("m_{0} (GeV)");
  h_ABCD_pred-> GetXaxis()->SetTitle("m_{0} (GeV)");
  h_ABCD_obs -> GetXaxis()->SetTitle("m_{0} (GeV)");
  h_ABCD_zbi -> GetXaxis()->SetTitle("m_{0} (GeV)");
  h_ABCD_zn  -> GetXaxis()->SetTitle("m_{0} (GeV)");
  h_ABCD_exc -> GetXaxis()->SetTitle("m_{0} (GeV)");

  hpred100   -> GetYaxis()->SetTitle("m_{1/2} (GeV)");
  hobs100    -> GetYaxis()->SetTitle("m_{1/2} (GeV)");
  hzbi100    -> GetYaxis()->SetTitle("m_{1/2} (GeV)");
  hzn100     -> GetYaxis()->SetTitle("m_{1/2} (GeV)");
  hexc100    -> GetYaxis()->SetTitle("m_{1/2} (GeV)");
  hpred175   -> GetYaxis()->SetTitle("m_{1/2} (GeV)");
  hobs175    -> GetYaxis()->SetTitle("m_{1/2} (GeV)");
  hzbi175    -> GetYaxis()->SetTitle("m_{1/2} (GeV)");
  hzn175     -> GetYaxis()->SetTitle("m_{1/2} (GeV)");
  hexc175    -> GetYaxis()->SetTitle("m_{1/2} (GeV)");
  h_ABCD_pred-> GetXaxis()->SetTitle("m_{1/2} (GeV)");
  h_ABCD_obs -> GetXaxis()->SetTitle("m_{1/2} (GeV)");
  h_ABCD_zbi -> GetXaxis()->SetTitle("m_{1/2} (GeV)");
  h_ABCD_zn  -> GetXaxis()->SetTitle("m_{1/2} (GeV)");
  h_ABCD_exc -> GetXaxis()->SetTitle("m_{1/2} (GeV)");

  hpred100->SetTitle("Predicted Yield (met > 100 GeV)");
  hobs100 ->SetTitle("Observed Yield (met > 100 GeV)");
  hzbi100 ->SetTitle("Z_{Bi} (met > 100 GeV)");
  hzn100  ->SetTitle("Z_{N} (met > 100 GeV)");
  hexc100 ->SetTitle("Exclusion (met > 100 GeV)");

  //hpred175->SetTitle("Predicted Yield (met > 175 GeV)");
  //hobs175 ->SetTitle("Observed Yield (met > 175 GeV)");

  if(strcmp(suffix,"100pb")==0){
    hpred175 ->SetTitle("Predicted Background Yield (Victory) (100 pb^{-1})");
    hobs175 ->SetTitle("Observed SUSY Yield (Victory) (100 pb^{-1})");
    hzbi175 ->SetTitle("Z_{Bi} (Victory) (100 pb^{-1})");
    hzn175 ->SetTitle("Z_{N} (Victory) (100 pb^{-1})");
    hexc175 ->SetTitle("95% CL Excluded Region (Victory) (100 pb^{-1})");

    h_ABCD_pred ->SetTitle("Predicted Background Yield (ABCD) (100 pb^{-1})");
    h_ABCD_obs  ->SetTitle("Observed SUSY Yield (ABCD) (100 pb^{-1})");
    h_ABCD_zbi  ->SetTitle("Z_{Bi} (ABCD) (100 pb^{-1})");
    h_ABCD_zn   ->SetTitle("Z_{N} (ABCD) (100 pb^{-1})");
    h_ABCD_exc  ->SetTitle("95% CL Excluded Region (ABCD) (100 pb^{-1})");
  }

  if(strcmp(suffix,"1fb")==0){
    hpred175 ->SetTitle("Predicted Background Yield (Victory) (1 fb^{-1})");
    hobs175 ->SetTitle("Observed SUSY Yield (Victory) (1 fb^{-1})");
    hzbi175 ->SetTitle("Z_{Bi} (Victory) (1 fb^{-1})");
    hzn175 ->SetTitle("Z_{N} (Victory) (1 fb^{-1})");
    hexc175 ->SetTitle("95% CL Excluded Region (Victory) (1 fb^{-1})");

    h_ABCD_pred ->SetTitle("Predicted Background Yield (ABCD) (1 fb^{-1})");
    h_ABCD_obs  ->SetTitle("Observed SUSY Yield (ABCD) (1 fb^{-1})");
    h_ABCD_zbi  ->SetTitle("Z_{Bi} (ABCD) (1 fb^{-1})");
    h_ABCD_zn   ->SetTitle("Z_{N} (ABCD) (1 fb^{-1})");
    h_ABCD_exc  ->SetTitle("95% CL Excluded Region (ABCD) (1 fb^{-1})");
  }


 
}



int getXBin(TH2F* h, float xval){
  
  float binsize = (h->GetXaxis()->GetXmax() - h->GetXaxis()->GetXmin())/h->GetNbinsX();
  int bin = (int) (xval - h->GetXaxis()->GetXmin()) / binsize;  
  
  return bin;

}

int getYBin(TH2F* h, float yval){
  
  float binsize = (h->GetYaxis()->GetXmax() - h->GetYaxis()->GetXmin())/h->GetNbinsY();
  int bin = (int) (yval - h->GetYaxis()->GetXmin()) / binsize; 
  
  return bin;
}
















/*
double getZBi(double n_on , double n_off , double tau){

  double P_Bi     = TMath::BetaIncomplete(1./(1.+tau), n_on, n_off+1);

  //use approximation for Z if p is very small
  if(P_Bi < 1.e-10){
    double mu = -2*log(P_Bi * sqrt( 2*pi ));
    return sqrt( mu - log(mu) );
  }
 
  else{
    return sqrt(2.)*TMath::ErfInverse(1 - 2.*P_Bi);
  }

}
*/
