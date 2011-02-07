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
#include "TRandom3.h"
#include <iomanip>

using namespace std;

inline double fround(double n, double d){
  return floor(n * pow(10., d) + .5) / pow(10., d);
}

//-------------------------------------------------------

float scale             = 1.;//35./34.85;
const int   width1      = 20;
const int   width2      = 4;
const bool  printdata   = true;     //print data yields?
const bool  showError   = true;    //show errors for MC?
const bool  noLumiScale = false;     //for MC, show event counts
const char* yield       = "yield";
const char* data        = "dataskim";

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

void print( TH1F* h , string label , bool error );
void printLine();

void printYields(string filename){

  if( makeLatexPlot ){
    pm         = " $\\pm$ ";
    delim      = "&";
    delimstart = "";
    delimend   = "\\\\";
    ee         = "$ee$";
    mm         = "$\\mu\\mu$";
    em         = "$e\\mu$";
  }

  TFile *f=TFile::Open( filename.c_str() );

  vector<char*> mcsamples;
  vector<char*> mcsamples_tex;
  //mcsamples.push_back("ttall"); 
  mcsamples.push_back("ttdil");      mcsamples_tex.push_back("$t\\bar{t}\\rightarrow \\ell^{+}\\ell^{-}$");
  mcsamples.push_back("ttotr");      mcsamples_tex.push_back("$t\\bar{t}\\rightarrow \\mathrm{other}$");
  //mcsamples.push_back("DYee");       mcsamples_tex.push_back("Z^0 \\rightarrow e^{+}e^{-}");
  //mcsamples.push_back("DYmm");       mcsamples_tex.push_back("Z^0 \\rightarrow \\mu^{+}\\mu^{-}");
  //mcsamples.push_back("DYtautau");   mcsamples_tex.push_back("Z^0 \\rightarrow \\tau^{+}\\tau^{-}");
  //mcsamples.push_back("DYtot");      mcsamples_tex.push_back("$Z^0 \\rightarrow \\ell^{+}\\ell^{-}$");  
  mcsamples.push_back("DYall");      mcsamples_tex.push_back("$Z^0 \\rightarrow \\ell^{+}\\ell^{-}$");  
  //mcsamples.push_back("Zjets");      mcsamples_tex.push_back("$Z^0$ + jets");
  mcsamples.push_back("wjets");      mcsamples_tex.push_back("$W^{\\pm}$ + jets");
  
  //   //mcsamples.push_back("vv");         mcsamples_tex.push_back("$VV$");
  mcsamples.push_back("ww");         mcsamples_tex.push_back("$W^+W^-$");
  mcsamples.push_back("wz");         mcsamples_tex.push_back("$W^{\\pm}Z^0$");
  mcsamples.push_back("zz");         mcsamples_tex.push_back("$Z^0Z^0$");
  mcsamples.push_back("tW");         mcsamples_tex.push_back("single top");
  //   //mcsamples.push_back("ttpythia");         mcsamples_tex.push_back("tt");
  //   //mcsamples.push_back("ttpythiapileup");  mcsamples_tex.push_back("tt (PU)");
  
  const unsigned int nmcsamples = mcsamples.size();
  
  vector<char*> susysamples;
  susysamples.push_back("LM0");
  susysamples.push_back("LM1");
  const unsigned int nsusysamples = susysamples.size();
 
  TH1F* h      = new TH1F();  
  TH1F* hmctot = new TH1F();

  printLine();

  //print header
  cout << delimstart << setw(width1) << "Sample"    << setw(width2)
       << delim      << setw(width1) << ee        << setw(width2)
       << delim      << setw(width1) << mm        << setw(width2)
       << delim      << setw(width1) << em        << setw(width2)
       << delim      << setw(width1) << "tot"       << setw(width2) 
       << delimend   << endl;
  
  printLine();

  //print SM MC samples
  for(unsigned int imcsample = 0 ; imcsample < nmcsamples ; imcsample++){

    if( strcmp( mcsamples.at(imcsample) , "DYall" ) == 0 ){
      TH1F* hee = (TH1F*) f->Get(Form("DYee_%s",yield));
      TH1F* hmm = (TH1F*) f->Get(Form("DYmm_%s",yield));
      TH1F* htt = (TH1F*) f->Get(Form("DYtautau_%s",yield));
      h = (TH1F*) hee->Clone();
      h->Add(hmm);
      h->Add(htt);
    }
    else{
      h = (TH1F*) f->Get(Form("%s_%s",mcsamples.at(imcsample),yield));
    }

    if(h==0) continue;

    if( imcsample == 0 ) hmctot = (TH1F*) h->Clone();
    else                 hmctot->Add(h);
    
    if( makeLatexPlot )  print( h , mcsamples_tex[imcsample] , showError );
    else                 print( h , mcsamples[imcsample]     , showError );

  }

  printLine();
  
  //print sum of SM MC samples
  print( hmctot , "total SM MC" , showError );

  printLine();

  if( printdata ){
    h = (TH1F*) f->Get(Form("%s_%s",data,yield));
    if( h != 0 ){
      print( h , "data" , false );
      printLine();
    }
  }

  //print SUSY MC samples  
  if( nsusysamples > 0 ){
    
    for(unsigned int isusysample = 0 ; isusysample < nsusysamples ; isusysample++){
      
      h = (TH1F*) f->Get(Form("%s_%s",susysamples.at(isusysample),yield));
      if(h==0) continue;
      
      print( h , susysamples[isusysample] , showError );
      
    }
    printLine();
  }

}

void printLine(){

  if( makeLatexPlot ){
    cout << "\\hline" << endl;
  }
  else{
    cout << "------------------------------------------------"
         << "------------------------------------------------" << endl;
  }
}


void print( TH1F* h , string label , bool error ){

  stringstream see;
  stringstream smm;
  stringstream sem;
  stringstream stot;

  float lumiscale = 1;

  if( noLumiScale ) lumiscale = h->GetEntries() / h->Integral();

  if( label == "data" ){
    if( error ){
      see  << Form("%.0f",lumiscale*scale*h->GetBinContent(2)) << pm << Form("%.0f",lumiscale*scale*h->GetBinError(2));
      smm  << Form("%.0f",lumiscale*scale*h->GetBinContent(3)) << pm << Form("%.0f",lumiscale*scale*h->GetBinError(3));
      sem  << Form("%.0f",lumiscale*scale*h->GetBinContent(4)) << pm << Form("%.0f",lumiscale*scale*h->GetBinError(4));
      stot << Form("%.0f",lumiscale*scale*h->GetBinContent(1)) << pm << Form("%.0f",lumiscale*scale*h->GetBinError(1));
    }else{
      see  << Form("%.0f",lumiscale*scale*h->GetBinContent(2));
      smm  << Form("%.0f",lumiscale*scale*h->GetBinContent(3));
      sem  << Form("%.0f",lumiscale*scale*h->GetBinContent(4));
      stot << Form("%.0f",lumiscale*scale*h->GetBinContent(1));
    }
  }else{
    if( error ){
      see  << Form("%.2f",lumiscale*scale*h->GetBinContent(2)) << pm << Form("%.2f",lumiscale*scale*h->GetBinError(2));
      smm  << Form("%.2f",lumiscale*scale*h->GetBinContent(3)) << pm << Form("%.2f",lumiscale*scale*h->GetBinError(3));
      sem  << Form("%.2f",lumiscale*scale*h->GetBinContent(4)) << pm << Form("%.2f",lumiscale*scale*h->GetBinError(4));
      stot << Form("%.2f",lumiscale*scale*h->GetBinContent(1)) << pm << Form("%.2f",lumiscale*scale*h->GetBinError(1));
    }else{
      see  << Form("%.0f",lumiscale*scale*h->GetBinContent(2));
      smm  << Form("%.0f",lumiscale*scale*h->GetBinContent(3));
      sem  << Form("%.0f",lumiscale*scale*h->GetBinContent(4));
      stot << Form("%.0f",lumiscale*scale*h->GetBinContent(1));
    }
  }

  cout << delimstart << setw(width1) << label      << setw(width2)
       << delim      << setw(width1) << see.str()  << setw(width2)
       << delim      << setw(width1) << smm.str()  << setw(width2)
       << delim      << setw(width1) << sem.str()  << setw(width2)
       << delim      << setw(width1) << stot.str() << setw(width2)
       << delimend   << endl;
  
  
}

 
