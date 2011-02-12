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
const float scale     = 1;//100./2.77;//30/2.77; //apply scaling
const int   nprec     = 2;       //ndigits precision
const int   width1    = 15;
const int   width2    = 4;
const bool  printdata = true;     //print data yields?
const bool  showError = false;    //show errors for MC?
//-------------------------------------------------------

void print( TH1F* h , string label , bool error );

void printYields(string filename){

  TFile *f=TFile::Open( filename.c_str() );

  vector<char*> mcsamples;
  mcsamples.push_back("ttdil"); 
  mcsamples.push_back("ttotr"); 
  mcsamples.push_back("Zjets"); 
  mcsamples.push_back("wjets"); 
  mcsamples.push_back("qcdpt15"); 
  mcsamples.push_back("qcdpt30"); 
  mcsamples.push_back("ww"); 
  mcsamples.push_back("wz"); 
  mcsamples.push_back("zz"); 
  mcsamples.push_back("tW"); 
  const unsigned int nmcsamples = mcsamples.size();

  vector<char*> susysamples;
  susysamples.push_back("LM0");
  susysamples.push_back("LM1");
  const unsigned int nsusysamples = susysamples.size();
 
  TH1F* h      = new TH1F();  
  TH1F* hmctot = new TH1F();

  cout << endl << endl
       << "------------------------------------------------" 
       << "------------------------------------------------" << endl;

  //print header


  cout << "|" << setw(width1) << "Sample"    << setw(width2)
       << "|" << setw(width1) << "ee"        << setw(width2)
       << "|" << setw(width1) << "mm"        << setw(width2)
       << "|" << setw(width1) << "em"        << setw(width2)
       << "|" << setw(width1) << "tot"       << setw(width2) 
       << "|" << endl;
  

  //print SM MC samples
  for(unsigned int imcsample = 0 ; imcsample < nmcsamples ; imcsample++){

    h = (TH1F*) f->Get(Form("%s_yield",mcsamples.at(imcsample)));
    if(h==0) continue;

    if( imcsample == 0 ) hmctot = (TH1F*) h->Clone();
    else                 hmctot->Add(h);
    
    print( h , mcsamples[imcsample] , showError );

  }

  cout << "------------------------------------------------" 
       << "------------------------------------------------" << endl;

  //print sum of SM MC samples
  print( hmctot , "tot SM MC" , showError );

  cout << "------------------------------------------------" 
       << "------------------------------------------------" << endl;

  if( printdata ){
    h = (TH1F*) f->Get("data_yield");
    if( h != 0 ){
      print( h , "data" , false );
      cout << "------------------------------------------------" 
           << "------------------------------------------------" << endl;
    }
  }

  //print SUSY MC samples  
  if( nsusysamples > 0 ){
    
    for(unsigned int isusysample = 0 ; isusysample < nsusysamples ; isusysample++){
      
      h = (TH1F*) f->Get(Form("%s_yield",susysamples.at(isusysample)));
      if(h==0) continue;
      
      print( h , susysamples[isusysample] , showError );
      
    }
    
    cout << "------------------------------------------------" 
         << "------------------------------------------------" << endl;
  }
  

}


void print( TH1F* h , string label , bool error ){

  stringstream see;
  stringstream smm;
  stringstream sem;
  stringstream stot;

  if( error ){
    see  << fround(scale*h->GetBinContent(2),nprec) << " +/- " << fround(scale*h->GetBinError(2),nprec);
    smm  << fround(scale*h->GetBinContent(3),nprec) << " +/- " << fround(scale*h->GetBinError(3),nprec);
    sem  << fround(scale*h->GetBinContent(4),nprec) << " +/- " << fround(scale*h->GetBinError(4),nprec);
    stot << fround(scale*h->GetBinContent(1),nprec) << " +/- " << fround(scale*h->GetBinError(1),nprec);
  }else{
    see  << fround(scale*h->GetBinContent(2),nprec);
    smm  << fround(scale*h->GetBinContent(3),nprec);
    sem  << fround(scale*h->GetBinContent(4),nprec);
    stot << fround(scale*h->GetBinContent(1),nprec);
  }

  cout << "|" << setw(width1) << label      << setw(width2)
       << "|" << setw(width1) << see.str()  << setw(width2)
       << "|" << setw(width1) << smm.str()  << setw(width2)
       << "|" << setw(width1) << sem.str()  << setw(width2)
       << "|" << setw(width1) << stot.str() << setw(width2)
       << "|" << endl;
  
  
}
