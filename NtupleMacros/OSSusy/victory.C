#include "TROOT.h"
#include "TH1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TLegend.h"
#include <iostream>
#include <iomanip>
#include <math.h>
#include <string>
#include <sstream>
#include "TLegendEntry.h"
#include "TTree.h"
#include "TChain.h"
#include "TLatex.h"
#include "TLine.h"
//#include "histtools.h"

using namespace std;

//-------------------------------------------------------

//const float scale_factor = 2;
const int   width1      = 20;
const int   width2      = 4;

//const string plottitle  = "signal region (H_{T} > 300 GeV)";
//const string plottitle  = "SM MC";

//const char* iter        = "output/nov5th_v7";
//const char* iter        = "output_38X/nov5th_v2";
const char* iter        = "output_38X/nov5th_v5_skim";
//const char* iter        = "oct15th_v2";

const bool plotData     = true;
const int  nprec        = 2;
const char* data        = "dataskim";

const bool issignal       = true; //signal/control?

string plottitle  = "control region (125 < H_{T} < 300 GeV)";


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

// Scaling factor 1.46848
// Warning in <TCanvas::Constructor>: Deleting canvas with same name: c1

// |           | Predicted |  Observed |  Obs/Pred |
// |        MC |    1.3631 |    1.8767 |    1.3767 |

//-------------------------------------------------------

void plotHist( TH1F* hpred , TH1F* hobs , TH1F* hpred_data , TH1F* hobs_data,
               string title, string xtitle , int rebin, float metcut);

inline double fround(double n, double d){
  return floor(n * pow(10., d) + .5) / pow(10., d);
}
void printRow( string sample , TH1F* hobs , TH1F* hpred , float cut );
void printLine(){

  if( makeLatexPlot ){
    cout << "\\hline" << endl;
  }
  else{
    cout << "-----------------------------------------------------------------------------" << endl;
        
  }
}
void printHeader(){
  
  cout << delimstart << setw(width1) << "Sample"        << setw(width2)
       << delim      << setw(width1) << "Predicted"     << setw(width2)
       << delim      << setw(width1) << "Observed"      << setw(width2)
       << delim      << setw(width1) << "Obs/pred"      << setw(width2)
       << delimend   << endl;
}

void victory( bool printgif =false ){

  if( issignal ) plottitle  = "signal region (H_{T} > 300 GeV)";

  if( makeLatexPlot ){
    pm         = " $\\pm$ ";
    delim      = "&";
    delimstart = "";
    delimend   = "\\\\";
    ee         = "$ee$";
    mm         = "$\\mu\\mu$";
    em         = "$e\\mu$";
  }

  int nbins  = 40;
  float xmin = 0;
  float xmax = 21.25;
  
//   int nbins  = 40;
//   float xmin = 0;
//   float xmax = 20.;

  //declare MC samples to include
  vector <string> samples;
  samples.push_back("ttdil");
  samples.push_back("ttotr");
  //samples.push_back("DYee");
  //samples.push_back("DYmm");
  samples.push_back("DYtautau");
  samples.push_back("wjets");
  samples.push_back("ww");
  samples.push_back("wz");
  samples.push_back("zz");
  samples.push_back("tW");
  //samples.push_back("LM1");
  const unsigned int nsamples = samples.size();
  
  //create MC chains and histos
  TChain *ch[nsamples];
  TChain *ch_mctot = new TChain("t");
  
  TH1F* hmety[nsamples];
  TH1F* hdilpty[nsamples];
  TH1F* hdilpt[nsamples];

  TH1F* hmety_mctot     = new TH1F("hmety_mctot"    , "hmety_mctot"   , nbins,xmin,xmax);     
  TH1F* hdilpty_mctot   = new TH1F("hdilpty_mctot"  , "hdilpty_mctot" , nbins,xmin,xmax);     
  TH1F* hdilpt_mctot    = new TH1F("hdilpt_mctot"   , "hdilpt_mctot"  , 300,0,300);           

  hmety_mctot->Sumw2();    
  hdilpty_mctot->Sumw2();    
  hdilpt_mctot->Sumw2();       

  for( unsigned int i = 0 ; i < nsamples ; ++i ){
    ch[i] = new TChain("t");
    ch[i]->Add(Form( "%s/%s_smallTree.root" , iter , samples.at(i).c_str() ));
    ch_mctot->Add(Form( "%s/%s_smallTree.root" , iter , samples.at(i).c_str() ));
 
    hmety[i]    =  new TH1F(Form("%s_hmety"    , samples.at(i).c_str()) , 
                            Form("%s_hmety"    , samples.at(i).c_str()) , nbins,xmin,xmax);

    hdilpty[i]  =  new TH1F(Form("%s_hdilpty"  , samples.at(i).c_str()) , 
                            Form("%s_hmety"    , samples.at(i).c_str()) , nbins,xmin,xmax); 
    
    hdilpt[i]   =  new TH1F(Form("%s_hdilpt"   , samples.at(i).c_str()) , 
                            Form("%s_hmety"    , samples.at(i).c_str()) , 300,0,300);           
    
    hdilpty[i]->Sumw2();
    hmety[i]->Sumw2();
    hdilpt[i]->Sumw2();
  }
  
  //create data chain and histos
  TChain *chdata = new TChain("t");
  if( plotData )
    chdata->Add(Form("%s/%s_smallTree.root",iter,data));
 
  TH1F* hmety_data          =  new TH1F("hmety_data",   "", nbins,xmin,xmax);
  TH1F* hdilpty_data        =  new TH1F("hdilpty_data", "", nbins,xmin,xmax); 
  TH1F* hdilpt_data         =  new TH1F("hdilpt_data",  "", 300,0,300); 

  hdilpty_data->Sumw2();
  hmety_data->Sumw2();

  //declare cuts for event selection
  char* jetcutstring = "sumjetpt>125&&sumjetpt<300&&njets>1";
  if( issignal ) jetcutstring = "sumjetpt>300&&njets>1";
  TCut jetcut(jetcutstring);
  //TCut jetcut("sumjetpt>300&&njets>1");
  //TCut jetcut("sumjetpt>125&&sumjetpt<300&&njets>1");
  TCut zcut("passz==0");  
  TCut metcut("tcmet>50");
  TCut weight("weight");

  TCut fcut  = ( jetcut + zcut +metcut ) * weight;

  //TCut nullcut("ptl1>-9999");  
  //TCut weight("1");
  //TCut etacut("abs(etall)<2.5&&abs(etalt)<2.5");
  //TCut pt2010cut("ptl1>20&&ptl2>10");
  //TCut pt2020cut("ptl1>20&&ptl2>20");
  //TCut dphicut("dildphi<1");
  //TCut mucut("mull==1&&mult==1");  
  //TCut nlepcut("nlep==2");  
  //TCut weight("costhetaweight");
  //TCut weightcut("weight>0");

  TCanvas *ctemp =new TCanvas();
  ctemp->cd();

  //total MC
  ch_mctot->Draw("TMath::Min(tcmet/sqrt(sumjetpt),19.9)>>hmety_mctot",     fcut);
  ch_mctot->Draw("TMath::Min(dilpt/sqrt(sumjetpt),19.9)>>hdilpty_mctot",   fcut);
  ch_mctot->Draw("TMath::Min(dilpt,299.9)>>hdilpt_mctot",                  fcut);

 

  //apply dilepton pt scale factor
  float k = hdilpt_mctot->Integral() /  hdilpt_mctot->Integral( hdilpt_mctot->FindBin(50) , 100000 );
  cout << "Scaling factor " << k << endl;
  hdilpty_mctot->Scale( k );

  
  //data
  if( plotData ){
    chdata->Draw("TMath::Min(tcmet/sqrt(sumjetpt),19.9)>>hmety_data",     fcut);
    chdata->Draw("TMath::Min(dilpt/sqrt(sumjetpt),19.9)>>hdilpty_data",   fcut);
    chdata->Draw("TMath::Min(dilpt,299.9)>>hdilpt_data",                  fcut);
  }
  delete ctemp;

  float kdata = 1;
  if( plotData ){
    kdata = k;
    cout << "Using data scale factor from MC " << k << endl;
    hdilpty_data->Scale( kdata );
  }

  
  //MC
  TH1F* hmety_temp     = new TH1F("hmety_temp"    , "" , nbins,xmin,xmax);     
  TH1F* hdilpty_temp   = new TH1F("hdilpty_temp"  , "" , nbins,xmin,xmax);     
  TH1F* hdilpt_temp    = new TH1F("hdilpt_temp"   , "" , 300 , 0 , 300);     
  
  hmety_temp->Sumw2();
  hdilpty_temp->Sumw2();

  for( unsigned int i = 0 ; i < nsamples ; ++i ){

    ch[i]->Draw("TMath::Min(tcmet/sqrt(sumjetpt),19.9)>>hmety_temp",     fcut);
    ch[i]->Draw("TMath::Min(dilpt/sqrt(sumjetpt),19.9)>>hdilpty_temp",   fcut);
    ch[i]->Draw("TMath::Min(dilpt,299.9)>>hdilpt_temp",                  fcut);

    hmety[i]   = (TH1F*) hmety_temp->Clone(   Form("%s_hmety"    , samples.at(i).c_str()) );
    hdilpty[i] = (TH1F*) hdilpty_temp->Clone( Form("%s_hdilpty"  , samples.at(i).c_str()) );
    hdilpt[i]  = (TH1F*) hdilpt_temp->Clone(  Form("%s_hdilpt"   , samples.at(i).c_str()) );

    hdilpty[i]->Scale( k );
    hmety_temp->Reset();
    hdilpty_temp->Reset();

  }
  
  /*
  TH1F *hmety_mcsum   =  (TH1F*) hmety[0]->Clone();
  TH1F *hdilpty_mcsum =  (TH1F*) hdilpty[0]->Clone();
  TH1F *hdilpt_mcsum  =  (TH1F*) hdilpt[0]->Clone();
  
  for( unsigned int i = 1 ; i < nsamples ; ++i ){

    hmety[i]->Scale( scale_factor );
    hdilpty[i]->Scale( scale_factor );
    hdilpt[i]->Scale( scale_factor );
    
    hmety_mcsum->Add( hmety[i] );
    hdilpty_mcsum->Add( hdilpty[i] );
    hdilpt_mcsum->Add( hdilpt[i] );

  }

  float knew = hdilpt_mcsum->Integral() /  hdilpt_mcsum->Integral( hdilpt_mcsum->FindBin(50) , 100000 );
  cout << "knew " << knew << endl;
  hdilpty_mcsum->Scale(knew);
  printRow( "mcsum" , hmety_mcsum , hdilpty_mcsum , 8.5 );
  */

  //dineutrino pt
  TCanvas *c1=new TCanvas("c1","",800,600);
  c1->cd();  
  
  plotHist( hdilpty_mctot , hmety_mctot , hdilpty_data , hmety_data , 
            plottitle , "y   (GeV^{1/2})" , 4 , 8.5 );

  if( printgif ) c1->Print("plots/victory.gif");




  cout << endl << endl;

  printLine();
  printHeader();
  printLine();

  for( unsigned int i = 0 ; i < nsamples ; ++i ){
    printRow( samples.at(i) , hmety[i] , hdilpty[i] , 8.5 );
  }

  printLine();
  printRow( "total MC" , hmety_mctot , hdilpty_mctot , 8.5 );
  printLine();
  if( plotData ){
    printRow( "data" , hmety_data , hdilpty_data , 8.5 );
    printLine();
  }









  /*

  //TCanvas *c1=new TCanvas("c1","",800,1200);
  //c1->Divide(1,2);
  //c1->cd(1);

  cout << hmety->GetBinContent(17) << " " << hmety->GetBinError(17) << endl;
  cout << hdilpty->GetBinContent(17) << " " << hdilpty->GetBinError(17) << endl;
  cout << hmety->GetBinContent(17) / hdilpty->GetBinContent(17) << endl;

  //TLatex t;
  //t.SetNDC();
  //t.DrawLatex(0.2,0.2,Form("k = %.3f",k));

  //TCanvas *c2 = new TCanvas("c2","",800,600);
  //c2->cd();
  c1->cd(2);

  TH1F* hratio = (TH1F*) hmety->Clone("ratio");
  hratio->Divide( hdilpty );
  //hratio->Rebin(2);
  //hratio->Scale(0.5);

  hratio->Draw();
  hratio->SetMinimum(-1);
  hratio->SetMaximum(5);
  hratio->GetXaxis()->SetTitle("tcmet / #sqrt{sumJetPt} (GeV^{1/2})");
  hratio->GetYaxis()->SetTitle("Observed / Predicted");
  cout << hratio->GetBinContent(17) << " " << hratio->GetBinError(17) << endl;

  TLine line;
  line.SetLineStyle(2);
  line.DrawLine(0,1.38,21.25,1.38);
  line.DrawLine(8.5,-1,8.5,5);
  */




  //t->Draw("TMath::Min(nupt ,299)>>hnupt_jets","sumjetpt>200 && njets>1");
  //t->Draw("TMath::Min(dilpt,299)>>hdilpt_jets","sumjetpt>200 && njets>1");
  //t->Draw("TMath::Min(dilpt,299)>>hdilpt_nupt50","nupt>50");

  /*
  //dineutrino pt (nupt > 50 GeV)
  float integral   =  hdilpt_nupt50->Integral();
  float integral50 =  hdilpt_nupt50->Integral( hdilpt_nupt50->FindBin(50) , 100000 );
  float k          =  integral / integral50;
  cout << "Scaling factor " << k << endl;
  hdilpt_nupt50->Scale( k );

  TCanvas *c2=new TCanvas("c2","",800,600);
  c2->cd();
  plotHist( hdilpt_nupt50 , hnupt , "" , "neutrino p_{T} (GeV)" , 1);

  //dineutrino pt (with jet cuts)
  TCanvas *c3=new TCanvas("c3","",800,600);
  c3->cd();
  plotHist(  hdilpt_jets , hnupt_jets , "sumJetPt > 200 GeV, nJets #geq 2" , "neutrino p_{T} (GeV)" , 1);
  */


  
  //ee 
  //t->Draw("TMath::Min(tcmet,299)>>hmet","(sumjetpt>200&&passz==0&&njets>1&&leptype==0)");
  //t->Draw("TMath::Min(dilpt,299)>>hdilpt","(tcmet>50&&sumjetpt>200&&passz==0&&njets>1&&leptype==0)");

  //mm
  //t->Draw("TMath::Min(tcmet,299)>>hmet","(sumjetpt>200&&passz==0&&njets>1&&leptype==1)");
  //t->Draw("TMath::Min(dilpt,299)>>hdilpt","(tcmet>50&&sumjetpt>200&&passz==0&&njets>1&&leptype==1)");

  //em
  //t->Draw("TMath::Min(tcmet,299)>>hmet","(sumjetpt>200&&passz==0&&njets>1&&leptype==2)");
  //t->Draw("TMath::Min(dilpt,299)>>hdilpt","(tcmet>50&&sumjetpt>200&&passz==0&&njets>1&&leptype==2)");



  //float integral   =  hmet->Integral();
  //float integral50 =  hmet->Integral( hdilpt->FindBin(50) , 100000 );
 
}


float calculateHistError( TH1F* h , int minbin , int maxbin ){

  float err2 = 0;
  for( int i = minbin ; i <= maxbin ; ++i ){
    err2 += pow( h->GetBinError(i) , 2 );
  }

  return sqrt(err2);
}



void printRow( string sample , TH1F* hobs , TH1F* hpred , float cut ){

  stringstream sobs;
  stringstream spred;
  stringstream sratio;

  int   minbin     = hobs->FindBin( cut );
  int   maxbin     = 10000;

  float pred    = hpred->Integral( minbin , maxbin);
  float prederr = calculateHistError( hpred , minbin , maxbin );
  
  float obs     = hobs->Integral(  minbin , maxbin );
  float obserr  = calculateHistError( hobs , minbin , maxbin );

  float ratio    = pred > 0 ? obs / pred : 0;
  float ratioerr = obs > 0 && pred > 0 ? ratio * sqrt( pow( prederr / pred , 2 ) + pow( obserr / obs , 2 ) ) : 0;
  
  sobs   << Form( "%.2f" , obs   ) << pm << Form( "%.2f" , obserr   ); 
  spred  << Form( "%.2f" , pred  ) << pm << Form( "%.2f" , prederr  ); 
  sratio << Form( "%.2f" , ratio ) << pm << Form( "%.2f" , ratioerr ); 

  cout << delimstart << setw(width1) << sample        << setw(width2)
       << delim      << setw(width1) << spred.str()   << setw(width2)
       << delim      << setw(width1) << sobs.str()    << setw(width2)
       << delim      << setw(width1) << sratio.str()  << setw(width2)
       << delimend   << endl;

}




void plotHist( TH1F* hpred , TH1F* hobs , TH1F* hpred_data, TH1F* hobs_data,
               string title, string xtitle , int rebin , float metcut){

  hpred->SetLineColor(4);
  hobs ->SetLineColor(2);
  hpred->SetMarkerColor(4);
  hobs ->SetMarkerColor(2);
  hpred->SetMarkerStyle(23);
  hobs->SetMarkerStyle(4);
  hpred_data->SetLineColor(4);
  hobs_data ->SetLineColor(2);
  hpred_data->SetMarkerColor(4);
  hobs_data ->SetMarkerColor(2);
  hpred_data->SetMarkerStyle(23);
  hobs_data->SetMarkerStyle(4);
  hpred_data->SetMarkerSize(2);
  hobs_data->SetMarkerSize(2);
  hpred->SetFillColor(0);
  hobs ->SetFillColor(0);
  hpred->SetTitle( title.c_str() );
  hpred->GetXaxis()->SetTitle( xtitle.c_str() );
  hpred->GetXaxis()->SetTitleSize(0.055);
  hpred->GetYaxis()->SetTitle( "Events   " );
  hpred->GetYaxis()->SetTitleSize(0.055);
  hpred->GetYaxis()->SetTitleOffset(1);
  //hpred->GetYaxis()->SetTitle( "Entries / 2.125 GeV" );

  int metbin = hobs->FindBin( metcut );
  float pred = hpred->Integral( metbin , 100000000);
  float obs  = hobs->Integral(  metbin , 100000000);

  float pred_data = hpred_data->Integral( metbin , 100000000);
  float obs_data  = hobs_data->Integral(  metbin , 100000000);

  //cout << "Predicted (met>" << metcut << ") " << pred << endl;
  //cout << "Observed  (met>" << metcut << ") " << obs  << endl;
  //cout << "obs/pred            " << obs/pred << endl;

  //int width1 = 10;
  //int width2 = 2;

  cout << endl;
  cout << "|" << setw(width1) << ""          << setw(width2)
       << "|" << setw(width1) << "Predicted" << setw(width2)
       << "|" << setw(width1) << "Observed"  << setw(width2)
       << "|" << setw(width1) << "Obs/Pred"  << setw(width2) << "|" << endl;

  cout << "|" << setw(width1) << "MC"      << setw(width2)
       << "|" << setw(width1) << fround(pred,nprec)      << setw(width2)
       << "|" << setw(width1) << fround(obs,nprec)       << setw(width2)
       << "|" << setw(width1) << fround(obs/pred,nprec)  << setw(width2) << "|" << endl;

  if( plotData ){
    cout << "|" << setw(width1) << "data"      << setw(width2)
         << "|" << setw(width1) << fround(pred_data,nprec)      << setw(width2)
         << "|" << setw(width1) << fround(obs_data,nprec)       << setw(width2)
         << "|" << setw(width1) << fround(obs_data/pred_data,nprec)  << setw(width2) << "|" << endl;
  }
  cout << endl;


  stringstream spred;
  spred << "Pred: " << pred ;
  stringstream sobs;
  sobs << "Obs:  " << obs ;

  if( rebin > 1 ){
    hobs->Rebin ( rebin );
    hpred->Rebin( rebin );
    hobs_data->Rebin ( rebin );
    hpred_data->Rebin( rebin );
  }
 
  if( plotData ){
    hpred->Draw("hist");
    hobs->Draw("samehist");
    hpred_data->Draw("sameE1");
    hobs_data->Draw("sameE1");
    hpred->SetLineWidth(2);
    hobs->SetLineWidth(2);
  }else{
    hpred->Draw("E1");
    hobs->Draw("sameE1");
  } 
  gPad->SetLogy(1);

  //TLatex t;
  //t.SetNDC();
  //t.DrawLatex(0.55,0.65,Form("obs/pred = %.2f",obs/pred));

  TLegend *leg = new TLegend(0.66,0.67,0.98,0.95);
  //leg->AddEntry(hpred, Form("pred>%.1f GeV = %.3f",metcut,pred) ,"l");
  //leg->AddEntry(hobs,  Form("obs>%.1f GeV = %.3f",metcut,obs)  ,"l");
  leg->AddEntry(hpred, "MC predicted" ,"l");
  leg->AddEntry(hobs,  "MC observed"  ,"l");
  if( plotData ){
    leg->AddEntry(hpred_data, "data predicted");
    leg->AddEntry(hobs_data,  "data observed" );
  }
  leg->SetFillColor(0);
  leg->SetBorderSize(1);
  leg->Draw();

  TLine line;
  line.SetLineWidth(2);
  line.SetLineStyle(2);
  //line.DrawLine( metcut , hpred->GetMinimum() , metcut , 2 * hpred->GetMaximum() );
  line.DrawLine( metcut , 0.5 * TMath::Min( hobs->GetMinimum() , hpred->GetMinimum() ),
                 metcut , 2 *   TMath::Max( hobs->GetMaximum() , hpred->GetMaximum() ) );

  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextSize(0.04);
  t->DrawLatex(0.18,0.32,"CMS");
  t->DrawLatex(0.18,0.27,"34.0 pb^{-1} at #sqrt{s} = 7 TeV");
  t->DrawLatex(0.18,0.22,"Events with ee/#mu#mu/e#mu");

  if( issignal ) t->DrawLatex(0.18,0.17,"H_{T} > 300 GeV (signal)");
  else           t->DrawLatex(0.18,0.17,"125 < H_{T} < 300 GeV (control)");
 
}
