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
#include "TLegendEntry.h"
#include "THStack.h"
#include "TStyle.h"
#include <sstream>
#include "TLatex.h"
#include "TLine.h"
#include "TMath.h"
//#include "zn.cc"

using namespace std;

//-------------------------------------------------------

bool plotData         = true;  //overlay data
bool drawLegend       = true;  //draw the legend
bool drawYieldsLegend = false; //draw the legend with ABCD yields
bool drawall          = true;  //combine MC into single TH2
bool staterror        = true;
const char* data      = "dataskim";

const int width1      = 20;
const int width2      = 4;
const int prec        = 2;

//-------------------------------------------------------

const bool  makeLatexPlot = false;        //plot in latex style
char* pm         = " +/- ";
char* delim      = "|";
char* delimstart = "|";
char* delimend   = "|";
char* ee         = "ee";
char* mm         = "mm";
char* em         = "em";

//-------------------------------------------------------




inline double fround(double n, unsigned d)
{
  return floor(n * pow(10., d) + .5) / pow(10., d);
}

double getZBi(double n_on , double mu_b_hat , double sigma_b);
//int getXBin(TH2F* h, float xval);
//int getYBin(TH2F* h, float yval);
void drawSquare(float x1, float y1, float x2, float y2, int color = 1);
void getPrediction(char* filename, bool display = false, 
		   float x1 = 30, float x2 = 50,  float x3 = 60,  float x4 = 150,
		   float y1 = 76, float y2 = 106, float y3 = 116, float y4 = 200);

void printLine(){

  if( makeLatexPlot ){
    cout << "\\hline" << endl;
  }
  else{
    cout << "-------------------------------------------------------------------------"
         << "------------------------------------------------------------------------" << endl;
  }
}
 


float histError2D( TH2F* hist , int xbin1 , int xbin2 , int ybin1 , int ybin2 ){

  float err2 = 0;

  for( int i = xbin1 ; i <= xbin2 ; ++i ){
    for( int j = ybin1 ; j <= ybin2 ; ++j ){
      //cout << i << " " << j << endl;
      err2 += pow( hist->GetBinError( i , j ) , 2 );
    }
  }

  return sqrt( err2 );

}

void printRow( char* sample, float A, float B, float C, float D, float dA = 0, float dB = 0, float dC = 0, float dD = 0 );
// float A;
// float B;
// float C;
// float D;
// float pred;
// float obs;
// float prederr;
// float obserr;

void plotABCD( char* filename ){

  if( makeLatexPlot ){
    pm         = " $\\pm$ ";
    delim      = "&";
    delimstart = "";
    delimend   = "\\\\";
    ee         = "$ee$";
    mm         = "$\\mu\\mu$";
    em         = "$e\\mu$";
  }
  
  float x1=125;
  float x2=300;
  float x3=300;
  float x4=1500;

  float y1 = 4.5;
  float y2 = 8.5;
  float y3 = 8.5;
  float y4 = 30;

//   float x1=132;
//   float x2=315;
//   float x3=315;
//   float x4=1500;
  
//   float y1 = 4.4;
//   float y2 = 8.3;
//   float y3 = 8.3;
//   float y4 = 30;
  
  bool display = true;
  getPrediction(filename, display,x1,x2,x3,x4,y1,y2,y3,y4);

 
  
}

void getPrediction(char* filename, bool display, 
                   float x1, float x2, float x3, float x4,
		   float y1, float y2, float y3, float y4){
  
  //TFile *file      = TFile::Open("output/oct15th/ossusy_JPT_tcmet_bitmask.root");
  TFile *file = TFile::Open(filename);
  

  const float lumi_scale     = 1.;//30./11.06;
  const float syst_bkg_error = 0.25;

  int colors[]={2,4,6,8};
  //gROOT->SetStyle("Plain");
  //gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
    
  vector<char*> samples;
  vector<char*> samples_tex;
  samples.push_back("ttdil");         samples_tex.push_back("$t\\bar{t}\\rightarrow \\ell^{+}\\ell^{-}$");
  samples.push_back("ttotr");         samples_tex.push_back("$t\\bar{t}\\rightarrow \\mathrm{other}$");
  samples.push_back("DYall");         samples_tex.push_back("$Z^0 \\rightarrow \\ell^{+}\\ell^{-}$");
  //samples.push_back("DYee");          samples_tex.push_back("$Z^0 \\rightarrow e^{+}e^{-}$");
  //samples.push_back("DYmm");          samples_tex.push_back("$Z^0 \\rightarrow \\mu^{+}\\mu^{-}$");
  //samples.push_back("DYtautau");      samples_tex.push_back("$Z^0 \\rightarrow \\tau^{+}\\tau^{-}$");
  samples.push_back("wjets");         samples_tex.push_back("$W^{\\pm}$ + jets");
  samples.push_back("ww");            samples_tex.push_back("$W^+W^-$");
  samples.push_back("wz");            samples_tex.push_back("$W^{\\pm}Z^0$");
  samples.push_back("zz");            samples_tex.push_back("$Z^0Z^0$");            
  samples.push_back("tW");            samples_tex.push_back("single top");
  //samples.push_back("LM0");           samples_tex.push_back("LM0");
  //samples.push_back("LM1");           samples_tex.push_back("LM1");
  //samples.push_back("susy"); 
  //samples.push_back("Zjets");         samples_tex.push_back("$Z^0$ + jets");
  //samples.push_back("DYall");         samples_tex.push_back("$Z^0 \\rightarrow \\ell^{+}\\ell^{-}$");
  //samples.push_back("SMother");       samples_tex.push_back("SM other");
  //samples.push_back("SMtot");         samples_tex.push_back("SM total");

//   samples.push_back("ttdil");         samples_tex.push_back("$t\\bar{t}\\rightarrow \\ell^{+}\\ell^{-}$");
//   samples.push_back("DYall");         samples_tex.push_back("$Z^0 \\rightarrow \\ell^{+}\\ell^{-}$");
//   samples.push_back("SMother");       samples_tex.push_back("SM other");

  const unsigned int nsamples = samples.size();

  TH2F* h[nsamples];
  TH2F* hall = new TH2F();
  TH2F* hdata = new TH2F();

  float A[nsamples];
  float B[nsamples];
  float C[nsamples];
  float D[nsamples];
  
  float Atot = 0;
  float Btot = 0;
  float Ctot = 0;
  float Dtot = 0;

  float Aerr[nsamples];
  float Berr[nsamples];
  float Cerr[nsamples];
  float Derr[nsamples];
  
  float Atoterr2 = 0;
  float Btoterr2 = 0;
  float Ctoterr2 = 0;
  float Dtoterr2 = 0;

  int Adata = 0;
  int Bdata = 0;
  int Cdata = 0;
  int Ddata = 0;

  TCanvas *c1=new TCanvas("c1","c1",800,600);
  c1->cd();
  
  TLegend *leg = new TLegend(0.7,0.6,0.85,0.85);
  //char* var = "sumJetPt_tcmetsqrtsumet";
  //char* var = "habcd_nopresel";
  char* var = "habcd";
  //char* var = "sumJetPt_tcmetsumet";
  
  //loop over samples
  for( unsigned int isample = 0 ; isample < nsamples ; isample++ ){
    
    //add all SM backgrounds except ttdil and Zjets
    if(strcmp(samples.at(isample),"SMother")==0){
      
      h[isample]  = (TH2F*) file->Get(Form("ttotr_%s_allj_all",var));
      TH2* htemp1 = (TH2F*) file->Get(Form("ww_%s_allj_all",var));
      TH2* htemp2 = (TH2F*) file->Get(Form("wz_%s_allj_all",var));
      TH2* htemp3 = (TH2F*) file->Get(Form("zz_%s_allj_all",var));
      TH2* htemp4 = (TH2F*) file->Get(Form("wjets_%s_allj_all",var));
      TH2* htemp5 = (TH2F*) file->Get(Form("tW_%s_allj_all",var));
      //TH2* htemp6 = (TH2F*) file->Get(Form("Zjets_%s_allj_all",var));
      h[isample]->Add(htemp1);
      h[isample]->Add(htemp2);
      h[isample]->Add(htemp3);
      h[isample]->Add(htemp4);
      h[isample]->Add(htemp5);
      //h[isample]->Add(htemp6);
    }
    else if(strcmp(samples.at(isample),"DYall")==0){
      h[isample]  = (TH2F*) file->Get(Form("DYee_%s_allj_all",var));
      TH2* htemp1 = (TH2F*) file->Get(Form("DYmm_%s_allj_all",var));
      TH2* htemp2 = (TH2F*) file->Get(Form("DYtautau_%s_allj_all",var));
      h[isample]->Add(htemp1);
      h[isample]->Add(htemp2);
    }
    else if( strcmp(samples.at(isample) , "susy" ) == 0){
      cout << "susy scan not currently implemented, quitting" << endl;
      exit(0);
      //h[isample]     = (TH2F*) susyfile->Get("susy_hmet_sumjetpt_m0_1_m12_1");
    }else{
      h[isample]   = (TH2F*) file->Get(Form("%s_%s_allj_all",samples.at(isample),var));
    } 
    
             
    h[isample] -> Scale( lumi_scale );
  

    //if(strcmp(samples.at(isample),"Zjets")==0) h[isample]->Scale(3);
      

    //h[isample]->Scale( h[isample]->GetEntries() / h[isample]->Integral() );
    
    if(isample == 0){
      hall = (TH2F*) h[isample]->Clone();
    }else{
      TH2F* htemp = (TH2F*) h[isample]->Clone();
      hall->Add(htemp);
    }
    
    //get yields in 4 ABCD regions
    //int ix1 = getXBin( h[isample] , x1 );
    //int ix2 = getXBin( h[isample] , x2 );
    //int ix3 = getXBin( h[isample] , x3 );
    //int ix4 = getXBin( h[isample] , x4 );
    
    //int iy1 = getYBin( h[isample] , y1 );
    //int iy2 = getYBin( h[isample] , y2 );
    //int iy3 = getYBin( h[isample] , y3 );
    //int iy4 = getYBin( h[isample] , y4 );

    //int ix1 = 151;
    //int ix2 = 325;
    //int ix3 = 326;
    //int ix4 = 100000;
    
    //int iy1 = 46;
    //int iy2 = 70;
    //int iy3 = 71;
    //int iy4 = 100000;

    int ix1 = (int) ( x1 + 1 );
    int ix2 = (int) ( x2     );
    int ix3 = (int) ( x3 + 1 );
    int ix4 = (int) ( 2000 );

    int iy1 = (int) ( 10 * y1 + 1 );
    int iy2 = (int) ( 10 * y2     );
    int iy3 = (int) ( 10 * y2 + 1 );
    int iy4 = (int) ( 2000      );

    //cout << ix1 << " " << ix2 << " " << ix3 << " " << ix4 << endl;
    //cout << iy1 << " " << iy2 << " " << iy3 << " " << iy4 << endl;

    //yield for each sample
    A[isample] = h[isample]->Integral( ix1, ix2, iy3, iy4 );
    B[isample] = h[isample]->Integral( ix1, ix2, iy1, iy2 );
    C[isample] = h[isample]->Integral( ix3, ix4, iy1, iy2 );
    D[isample] = h[isample]->Integral( ix3, ix4, iy3, iy4 );

    Aerr[isample] = histError2D( h[isample] , ix1, ix2, iy3, iy4 );
    Berr[isample] = histError2D( h[isample] , ix1, ix2, iy1, iy2 );
    Cerr[isample] = histError2D( h[isample] , ix3, ix4, iy1, iy2 );
    Derr[isample] = histError2D( h[isample] , ix3, ix4, iy3, iy4 );

    //total yield
    Atot += h[isample]->Integral( ix1, ix2, iy3, iy4 );
    Btot += h[isample]->Integral( ix1, ix2, iy1, iy2 );
    Ctot += h[isample]->Integral( ix3, ix4, iy1, iy2 );
    Dtot += h[isample]->Integral( ix3, ix4, iy3, iy4 );
    
    Atoterr2 += pow( histError2D( h[isample] , ix1, ix2, iy3, iy4 ) , 2);
    Btoterr2 += pow( histError2D( h[isample] , ix1, ix2, iy1, iy2 ) , 2);
    Ctoterr2 += pow( histError2D( h[isample] , ix3, ix4, iy1, iy2 ) , 2);
    Dtoterr2 += pow( histError2D( h[isample] , ix3, ix4, iy3, iy4 ) , 2);
    
    if( plotData ){

      hdata  = (TH2F*) file->Get(Form("%s_%s_allj_all",data,var));
      
//       int ix1 = 151;
//       int ix2 = 325;
//       int ix3 = 326;
//       int ix4 = 100000;
      
//       int iy1 = 46;
//       int iy2 = 70;
//       int iy3 = 71;
//       int iy4 = 100000;
      
      //yield for each sample
      Adata = (int) hdata->Integral( ix1, ix2, iy3, iy4 );
      Bdata = (int) hdata->Integral( ix1, ix2, iy1, iy2 );
      Cdata = (int) hdata->Integral( ix3, ix4, iy1, iy2 );
      Ddata = (int) hdata->Integral( ix3, ix4, iy3, iy4 );
      
    }
    
    //draw histos
    if(display){
      
      leg->AddEntry(h[isample] , samples.at(isample) );
      h[isample]->SetLineColor( colors[isample] );
      h[isample]->SetMarkerColor( colors[isample] );
      h[isample]->SetMarkerSize(0.1);
      h[isample]->SetFillColor( 0 );

      h[isample]->RebinX(10);
      h[isample]->RebinY(2);

      if(isample==0){
        h[isample]->Draw("box");
        h[isample]->SetTitle("");
        h[isample]->GetXaxis()->SetTitle("tcmet");
        h[isample]->GetXaxis()->SetRangeUser(0,1500);
        h[isample]->GetYaxis()->SetRangeUser(0,30);
        //h[isample]->GetYaxis()->SetTitle("#slash{E}_{T}  /  #sqrt{H_{T}} (GeV^{1/2})");
        h[isample]->GetYaxis()->SetTitle("y (GeV^{1/2})");
        h[isample]->GetYaxis()->SetTitleOffset(1);
        h[isample]->GetXaxis()->SetTitle("H_{T} (GeV)");
      }else{
        h[isample]->Draw("boxsame");
      }
    }
  }

  TH2F *hcut = new TH2F();

  if(drawall){

    hcut = (TH2F*) hall->Clone();
    //hall->Draw("colz");
    hall->Draw("box");
    hall->SetLineColor(4);
    hall->SetMarkerColor(0);
    hall->SetFillColor(0);
    hall->RebinX(20);
    hall->RebinY(4);
    hall->SetTitle("");
    hall->GetXaxis()->SetTitle("tcmet");
    hall->GetXaxis()->SetRangeUser(0,1500);
    hall->GetYaxis()->SetRangeUser(0,30);
    //hall->GetYaxis()->SetTitle("#slash{E}_{T}  /  #sqrt{H_{T}} (GeV^{1/2})");
    hall->GetYaxis()->SetTitle("y (GeV^{1/2})");
    hall->GetYaxis()->SetTitleOffset(1);
    hall->GetXaxis()->SetTitle("H_{T} (GeV)");
    hdata->SetMarkerColor(2);
    hdata->SetMarkerSize(1);
    //hdata->SetMarkerStyle(20);
    hdata->Draw("same");

    if( drawLegend ){
      TLegend *legall = new TLegend(0.7,0.7,0.9,0.9);
      legall->AddEntry(hall,"SM MC","f");
      legall->AddEntry(hdata,"Data","p");
      legall->SetFillColor(0);
      legall->SetBorderSize(1);
      legall->Draw();
    }
  }
  
  //drawSquare(50,0,250,76);
  leg->SetBorderSize(1);
  leg->SetFillColor(0);
  //leg->Draw();
  
  drawSquare(x1,y1,x2,y2);
  drawSquare(x3,y1,x4,y2);
  drawSquare(x1,y3,x2,y4);
  drawSquare(x3,y3,x4,y4);
  
  TLatex *t=new TLatex();
  t->SetTextSize(0.05);
  t->SetTextColor(1);
  t->DrawLatex(x1+50,y1+1.5,"B");
  t->DrawLatex(x1+50,y3+10,"A");
  t->DrawLatex(x3+50,y1+1.5,"C");
  t->DrawLatex(x3+50,y3+10,"D");

  t->SetTextSize(0.04);
  t->DrawLatex(350,27,"CMS");
  t->DrawLatex(350,25,"34.0 pb^{-1} at #sqrt{s} = 7 TeV");
  t->DrawLatex(350,23,"Events with ee/#mu#mu/e#mu");
  
  //TH2F *hcut=(TH2F*) h[0]->Clone();

  for( int ibinx = 1 ; ibinx <= hcut->GetXaxis()->GetNbins() ; ibinx++ ){
    for( int ibiny = 1 ; ibiny <= hcut->GetYaxis()->GetNbins() ; ibiny++ ){

      if( hcut->GetXaxis()->GetBinLowEdge( ibinx ) < x1 ) hcut->SetBinContent( ibinx , ibiny , 0 );
      if( hcut->GetYaxis()->GetBinLowEdge( ibiny ) < y1 ) hcut->SetBinContent( ibinx , ibiny , 0 );

    }
  }



//   hcut->GetXaxis()->SetRangeUser(x1,x4);
//   hcut->GetYaxis()->SetRangeUser(y1,y4);
  double cor = hcut->GetCorrelationFactor(1,2);  
  cout << "cor " << cor << endl;

  float pred = Atot*Ctot/Btot;
  float obs  = Dtot;
  float zbi  = getZBi(obs, pred, sqrt(pred + pow( syst_bkg_error * pred ,2) ) );
  //float zn   = getzn( obs - pred, pred , 0 , syst_bkg_error * 100);
  //t->DrawLatex(700,14,  Form("Observed  = %f",   obs));
  //t->DrawLatex(700,12,  Form("Predicted = %f",   pred));
  //t->DrawLatex(700,10,  Form("Z_{Bi} = %f",      zbi));
  //t->DrawLatex(700,8,   Form("Z_{N} = %f",       zn));
  //t->DrawLatex(400,23,  Form("cor. factor = %f", cor));

  if( drawYieldsLegend ){
    TBox *box = new TBox();
    box->SetFillColor(0);
    box->SetLineColor(1);
    box->SetLineWidth(2);
    box->DrawBox(690,15,1000,25);
    drawSquare(690,15,1000,25,1);
    
    float xtext = 700;
    float ytext = 24;
    //t->SetTextColor(2);
    t->SetTextSize(0.03);
    
    t->DrawLatex(xtext,ytext,   Form("N_{A} = %.2f",Atot));
    t->DrawLatex(xtext,ytext-2, Form("N_{B} = %.2f",Btot));
    t->DrawLatex(xtext,ytext-4, Form("N_{C} = %.2f",Ctot));
    t->DrawLatex(xtext,ytext-6, Form("N_{D} = %.2f",Dtot));
    t->DrawLatex(xtext,ytext-8, Form("N_{A}N_{C}/N_{B} = %.2f",Atot*Ctot/Btot));
    //   t->DrawLatex(xtext,ytext-10,Form("corr = %f",cor));
  }

  //print output to screen
  cout << setprecision(4) << endl;


  cout << "Observed  " << obs  << endl;
  cout << "Predicted " << pred << endl;
  cout << "ZBi       " << zbi  << endl;
  //cout << "ZN        " << zn   << endl << endl;

  printLine();

  cout  << delimstart << setw(width1) << "sample" << setw(width2) 
        << delim      << setw(width1) << "A"      << setw(width2) 
        << delim      << setw(width1) << "B"      << setw(width2) 
        << delim      << setw(width1) << "C"      << setw(width2) 
        << delim      << setw(width1) << "D"      << setw(width2) 
        << delim      << setw(width1) << "PRED"   << setw(width2) 
        << delimend << endl;


  printLine();
  
  for( unsigned int isample = 0 ; isample < nsamples ; isample++ ){

    char* sample = samples.at(isample);
    if( makeLatexPlot ) sample = samples_tex.at(isample);
    
    printRow( sample , A[isample], B[isample], C[isample], D[isample], 
              Aerr[isample], Berr[isample], Cerr[isample], Derr[isample] );
  }

  printLine();
  
  printRow( "total SM MC", Atot, Btot, Ctot, Dtot, sqrt(Atoterr2), sqrt(Btoterr2), sqrt(Ctoterr2), sqrt(Dtoterr2) );
  
  if( plotData ){
    
    printLine();
 
    printRow( "data", Adata, Bdata, Cdata, Ddata, sqrt(Adata), sqrt(Bdata), sqrt(Cdata), sqrt(Ddata) );

  }

  printLine();

  float ratio = Dtot / ( Atot * Ctot / Btot );

  float ratioerr = ratio * sqrt( Atoterr2/(Atot*Atot) + Btoterr2/(Btot*Btot) + Ctoterr2/(Ctot*Ctot) + Dtoterr2/(Dtot*Dtot) );
  cout << "obs/pred " << Form("%.2f",ratio) << " \\pm " << Form("%.2f",ratioerr) << endl;
 
}


void printRow( char* sample, float A, float B, float C, float D, float dA, float dB, float dC, float dD ){
  
  float pred    = B > 0 ? A * C / B : 0.;
  float prederr = (A > 0 && B > 0 && D > 0 ) ? pred * sqrt( pow( dA/A , 2 ) + pow( dB/B , 2) + pow( dC/C , 2) ) : 0.;
  
  stringstream sA;
  stringstream sB;
  stringstream sC;
  stringstream sD;
  stringstream sPred;

//   sA    << Form("%.2f", A     );
//   sB    << Form("%.2f", B     );
//   sC    << Form("%.2f", C     );
//   sD    << Form("%.2f", D     );
//   sPred << Form("%.2f %s %.2f", pred , pm , prederr );
  
  
  if( sample == "data"){
    sA    << A;
    sB    << B;
    sC    << C;
    sD    << D;
    sPred << Form("%.2f %s %.2f", pred , pm , prederr );
  }
  else{
    sA    << Form("%.2f %s %.2f", A    , pm , dA      );
    sB    << Form("%.2f %s %.2f", B    , pm , dB      );
    sC    << Form("%.2f %s %.2f", C    , pm , dC      );
    sD    << Form("%.2f %s %.2f", D    , pm , dD      );
    sPred << Form("%.2f %s %.2f", pred , pm , prederr );
  }
  







//   if( sample == "data"){
//     sA    << A   ;
//     sB    << B   ;
//     sC    << C   ;
//     sD    << D   ;
//     //sA    << fround( A    , prec) << pm << fround( dA      , prec);
//     //sB    << fround( B    , prec) << pm << fround( dB      , prec);
//     //sC    << fround( C    , prec) << pm << fround( dC      , prec);
//     //sD    << fround( D    , prec) << pm << fround( dD      , prec);
//     sPred << fround( pred , prec) << pm << fround( prederr , prec);
//   }
//   else{
//     sA    << Form("%.4f", A    );
//     sB    << Form("%.4f", B    );
//     sC    << Form("%.4f", C    );
//     sD    << Form("%.4f", D    );
//     sPred << Form("%.4f", pred );
//   }
   

  
  cout  << delimstart << setw(width1) << sample << setw(width2)
        << delim << setw(width1) << sA.str()    << setw(width2)
        << delim << setw(width1) << sB.str()    << setw(width2)
        << delim << setw(width1) << sC.str()     << setw(width2)
        << delim << setw(width1) << sD.str()     << setw(width2)
        << delim << setw(width1) << sPred.str()  << setw(width2)  << delimend << endl;
}



/*
int getXBin(TH2F* h, float xval){
  
  float binsize = (h->GetXaxis()->GetXmax() - h->GetXaxis()->GetXmin())/h->GetNbinsX();
  int bin = (int) ( (xval - h->GetXaxis()->GetXmin()) / binsize );  
  
  return bin;

}

int getYBin(TH2F* h, float yval){
  
  float binsize = (h->GetYaxis()->GetXmax() - h->GetYaxis()->GetXmin())/h->GetNbinsY();
  int bin = (int) ( (yval - h->GetYaxis()->GetXmin()) / binsize ); 
  
  return bin;
}
*/
void drawSquare(float x1, float y1, float x2, float y2, int color){

  TLine *line = new TLine();
  line->SetLineColor(color);
  line->SetLineWidth(2);

  line->DrawLine(x1,y1,x2,y1);
  line->DrawLine(x2,y1,x2,y2);
  line->DrawLine(x2,y2,x1,y2);
  line->DrawLine(x1,y2,x1,y1);

  delete line;
  
}

double getZBi(double n_on , double mu_b_hat , double sigma_b){

  double tau      = mu_b_hat / (sigma_b*sigma_b);
  double n_off    = tau*mu_b_hat;
  double P_Bi     = TMath::BetaIncomplete(1./(1.+tau), n_on, n_off+1);
  
  if(P_Bi == 0){
    return 100; // set ZBi = 100 if P_Bi = 0
  }
  else if(P_Bi < 1.e-10){
    double mu = -2*log(P_Bi * sqrt( 2*TMath::Pi() ));
    return sqrt( mu - log(mu) );
  }
  else{
    return sqrt(2.)*TMath::ErfInverse(1 - 2.*P_Bi);
  }
 
}






 /*
  static const int ncuts = 10;
  float cut[ncuts];
  float predArray[ncuts];
  float obsArray[ncuts];
  float baratio[ncuts];

   
  for(int i = 0; i<10; i++){
    getPrediction(false,x1,x2,x3,x4,y1,y2,y3,y4);
    cout<<"-----------------"<<endl;
    cout<<"met > "<<x3<<" GeV"<<endl;
    cout<<"predicted "<<pred<<endl; //" +/- "<<prederr<<endl;
    cout<<"observed  "<<obs <<endl; //" +/- "<<obserr<< endl;
    
    cut[i]=x3;
    predArray[i]=pred;
    obsArray[i]=obs;
    baratio[i] = B/A;
    x2 += 10;
    x3 += 10;
  }

  cout<<"-----------------"<<endl;
  cout<<"predicted "<<pred<<endl;
  cout<<"observed  "<<obs<<endl;

  TCanvas *c2=new TCanvas("c2","",800,600);
  c2->cd();

  TGraph *grpred  = new TGraph(ncuts,cut,predArray);
  TGraph *grobs   = new TGraph(ncuts,cut,obsArray);
  TGraph *grba    = new TGraph(ncuts,cut,baratio);
  grpred->SetMarkerColor(2);
  grobs ->SetMarkerColor(4);
  grba  ->SetMarkerColor(4);
  grpred->SetMarkerStyle(20);
  grobs-> SetMarkerStyle(20);
  grba -> SetMarkerStyle(20);
  grpred->SetMarkerSize(1);
  grobs-> SetMarkerSize(1);
  grba -> SetMarkerSize(1);
  grpred->Draw("AP");
  grpred->GetXaxis()->SetTitle("MET cut (GeV)");
  grpred->GetYaxis()->SetTitle("N(predicted)");
  grpred->SetTitle("Z+jets");
  //grobs->Draw("sameP");
  grba  ->Draw("sameP");
  grpred->SetMaximum(60);
  */
