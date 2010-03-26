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
//#include "histtools.h"

int nm0points  = 2;
int nm12points = 2;

TH2F* hpred100;
TH2F* hobs100; 
TH2F* hzbi100; 
TH2F* hpred175;
TH2F* hobs175; 
TH2F* hzbi175; 

namespace hist {
  void add(const char* outHistName, const char* patORpfx0);
  TLegend* legend(TCanvas* canvas, Option_t* option = "lpf", Bool_t addColor = kFALSE, Int_t token = -1, 
		  Float_t xmin = 0.75, Float_t ymin = 0.75, Float_t xmax = 0.99, Float_t ymax = 0.99);
}

void loadHist(const char* filename, const char* directory = 0, const char* pfx = 0, 
	      const char* pat = "*", Bool_t doAdd = kFALSE);

float getm0( int m0index ){
 
  if(m0index == 0) return 650.;
  if(m0index == 1) return 700.;

  return -9999.;
}

float getm12( int m12index ){
 
  if(m12index == 0) return 100.;
  if(m12index == 1) return 140.;

  return -9999.;
}

void BookHists();

double getZBi(double n_on , double mu_b_hat , double sigma_b);

void SusyScan(char* filename) {

  TFile *fin = TFile::Open(filename);
  gROOT->SetStyle("Plain");
  gROOT->ProcessLine("gStyle->SetOptStat(0)");
  ofstream ofile("SusyScan.txt");

  //Load and add histos--------------------------------------------
  loadHist(filename, 0, 0, "*dilPt_allj*",        kTRUE);
  loadHist(filename, 0, 0, "*susy_hdilPt*",       kTRUE);
  loadHist(filename, 0, 0, "*susy_htcmet*",       kTRUE);
  loadHist(filename, 0, 0, "*dilPtSmeared_allj*", kTRUE);
  loadHist(filename, 0, 0, "*tcmet_allj*",        kTRUE);
  
  std::cout << "adding histograms..." << std::endl;
  hist::add("sm_dilPt",        "^[^L].*_hdilPt_allj_all$");
  hist::add("sm_dilPtSmeared", "^[^L].*_hdilPtSmeared_allj_all$");  
  hist::add("sm_tcmet",        "^[^L].*_htcmet_allj_all$");  
  
  TH1F* sm_dilPt =        (TH1F*)gDirectory->Get("sm_dilPt");
  TH1F* sm_dilPtSmeared = (TH1F*)gDirectory->Get("sm_dilPtSmeared");
  TH1F* sm_tcmet =        (TH1F*)gDirectory->Get("sm_tcmet");
  
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
  
  //TFile *fout = new TFile("SusyScan.root","RECREATE");
  //fout->cd();
               
  BookHists();

  cout   << "|    m0     |   m12     |  pred100  |  obs100   |  pred175  |  obs175   |" << endl; 
  ofile  << "|    m0     |   m12     |  pred100  |  obs100   |  pred175  |  obs175   |" << endl; 

  //loop over susy scan points
  for(int im0 = 0 ; im0 < nm0points ; im0++){
    
    for(int im12 = 0 ; im12 < nm12points ; im12++){
  
      TH1F* hdilPt = (TH1F*)sm_dilPt->Clone(Form("hdilPt_m0_%i_m12_%i",im0,im12));
      TH1F* htcmet = (TH1F*)sm_tcmet->Clone(Form("htcmet_m0_%i_m12_%i",im0,im12));

      TH1F* hsusydilPt   = (TH1F*) fin->Get(Form("susy_hdilPt_m0%i_m12%i",im0,im12));
      TH1F* hsusytcmet   = (TH1F*) fin->Get(Form("susy_htcmet_m0%i_m12%i",im0,im12));
      
      hdilPt->Add(hsusydilPt);
      htcmet->Add(hsusytcmet);
      
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
      
      float m0  = getm0( im0 );
      float m12 = getm0( im12 );

      //print text to file
      ofile << "|" << setw(6) << setprecision(3) << m0 << setw(6) << "|" << setw(6) << m12 << setw(6) << "|" 
           << setw(6) << dilpt_prediction100 << setw(6) << "|" << setw(6) << met_observed100 << setw(6) << "|"
           << setw(6) << dilpt_prediction175 << setw(6) << "|" << setw(6) << met_observed175 << setw(6) << "|" << endl;
      
      cout << "|" << setw(6) << setprecision(3) << m0 << setw(6) << "|" << setw(6) << m12 << setw(6) << "|" 
           << setw(6) << dilpt_prediction100 << setw(6) << "|" << setw(6) << met_observed100 << setw(6) << "|"
           << setw(6) << dilpt_prediction175 << setw(6) << "|" << setw(6) << met_observed175 << setw(6) << "|" << endl;
    
      //fill histos
      hpred100 -> Fill(m0,m12,dilpt_prediction100);
      hobs100  -> Fill(m0,m12,met_observed100);
      hzbi100  -> Fill(m0,m12, getZBi( met_observed100 , dilpt_prediction100 , 0.2 * dilpt_prediction100  ) );
      hpred175 -> Fill(m0,m12,dilpt_prediction175);
      hobs175  -> Fill(m0,m12,met_observed175);
      hzbi175  -> Fill(m0,m12, getZBi( met_observed175 , dilpt_prediction175 , 0.2 * dilpt_prediction175 ) );

      //cout<<" ZBi100 "<<getZBi( met_observed100 , dilpt_prediction100 , 0.2 * dilpt_prediction100  )
      //    <<" ZBi175 "<<getZBi( met_observed175 , dilpt_prediction175 , 0.2 * dilpt_prediction175  )<<endl;
    }
      
  }
  
  TCanvas *c100 = new TCanvas("c100","met > 100 GeV canvas",1200,900);
  c100->Divide(2,2);
  c100->cd(1);
  hpred100->Draw("colz");
  c100->cd(2);
  hobs100->Draw("colz");
  c100->cd(3);
  hzbi100->Draw("colz");

  TCanvas *c175 = new TCanvas("c175","met > 175 GeV canvas",1200,900);
  c175->Divide(2,2);
  c175->cd(1);
  hpred175->Draw("colz");
  c175->cd(2);
  hobs175->Draw("colz");
  c175->cd(3);
  hzbi175->Draw("colz");

  
  ofile.close();
  //fout->cd();
  //fout->Write();
  //fout->Close();

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
  double Z_Bi     = sqrt(2.)*TMath::ErfInverse(1 - 2.*P_Bi);

  return Z_Bi;
}

void BookHists(){
  
  //book histos for predicted/observed yields and significance
  int   nm0    = 80;
  float m0min  = 0.;
  float m0max  = 4000.;
  int   nm12   = 25;
  float m12min = 100.;
  float m12max = 800.;

  hpred100 = new TH2F("hpred100" , "Predicted Yields (met>100)" , nm0,m0min,m0max,nm12,m12min,m12max);
  hobs100  = new TH2F("hobs100"  , "Observed Yields (met>100)"  , nm0,m0min,m0max,nm12,m12min,m12max);
  hzbi100  = new TH2F("hzbi100"  , "ZBI (met>100)"              , nm0,m0min,m0max,nm12,m12min,m12max);
  hpred175 = new TH2F("hpred175" , "Predicted Yields (met>175)" , nm0,m0min,m0max,nm12,m12min,m12max);
  hobs175  = new TH2F("hobs175"  , "Observed Yields (met>175)"  , nm0,m0min,m0max,nm12,m12min,m12max);
  hzbi175  = new TH2F("hzbi175"  , "ZBI (met>175)"              , nm0,m0min,m0max,nm12,m12min,m12max);
  
  hpred100->GetXaxis()->SetTitle("m_{0} (GeV)");
  hobs100 ->GetXaxis()->SetTitle("m_{0} (GeV)");
  hzbi100 ->GetXaxis()->SetTitle("m_{0} (GeV)");
  hpred175->GetXaxis()->SetTitle("m_{0} (GeV)");
  hobs175 ->GetXaxis()->SetTitle("m_{0} (GeV)");
  hzbi175 ->GetXaxis()->SetTitle("m_{0} (GeV)");

  hpred100->GetYaxis()->SetTitle("m_{1/2} (GeV)");
  hobs100 ->GetYaxis()->SetTitle("m_{1/2} (GeV)");
  hzbi100 ->GetYaxis()->SetTitle("m_{1/2} (GeV)");
  hpred175->GetYaxis()->SetTitle("m_{1/2} (GeV)");
  hobs175 ->GetYaxis()->SetTitle("m_{1/2} (GeV)");
  hzbi175 ->GetYaxis()->SetTitle("m_{1/2} (GeV)");

  hpred100->SetTitle("Predicted Yield (met > 100 GeV)");
  hobs100 ->SetTitle("Observed Yield (met > 100 GeV)");
  hzbi100 ->SetTitle("Z_{Bi} (met > 100 GeV)");
  hpred175->SetTitle("Predicted Yield (met > 175 GeV)");
  hobs175 ->SetTitle("Observed Yield (met > 175 GeV)");
  hzbi175 ->SetTitle("Z_{Bi} (met > 175 GeV)");

}

















