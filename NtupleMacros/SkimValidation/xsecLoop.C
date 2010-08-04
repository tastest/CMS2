#define xsecLoop_cxx
#include "xsecLoop.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <fstream>
#include <strstream>
#ifndef __CINT__

using namespace std;

void xsecLoop::Loop(std::vector<int>& numLepVsRun, Bool_t doMu)
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if(doMu) {
        // do the muon loop
        for( int mu = 0; mu < nmu; ++mu ) {
          if( 
             musid->at( mu )  &&
             musiso->at( mu ) < 0.15
             ) {
            numLepVsRun.at( run ) += 1;
            //             std::cout<<"muId  "<< musid->at( mu )  <<std::endl;
            //             std::cout<<"muIso "<< musiso->at( mu )  <<std::endl;
          }
        }
      }
      else {
        //        std::cout<<"now for the electrons "<<nel<<std::endl;
        // do the electron loop
        for( int el = 0; el < nel; ++el ) {
          //             std::cout<<"elId  "<< elsid->at( el )  <<std::endl;
          //             std::cout<<"elIso "<< elsiso->at( el )  <<std::endl;
          if( 
             elsid->at( el )  &&
             elsiso->at( el ) < 0.15
             ) {
            numLepVsRun.at( run ) += 1;
          }
        }
      }
   }
}

void xsecLoop::readFile(TString infile, std::vector<Double_t>& varVsRun) {

  TString current_file = "";
  TString line = "";
  float temp = 0;

  TString histoname = "";

  current_file.Append(infile);
  ifstream stream(current_file);

  Int_t linecounter = 0;
  Double_t var1 = 0.;
  Double_t var2 = 0.;

  Double_t prev_var1 = 0.;
  Double_t prev_var2 = 0.;

  Int_t final_cluster_count = 0;

  while ( line.ReadLine(stream) ) {
    istrstream stream5(line.Data());
    stream5 >>  var1 ;       
    stream5 >>  var2 ;       

    if ( var1 == prev_var1 ) continue;

    // save variables as previous and go to next line:
    //    prev_iRun     = iRun  ;
    prev_var1     = var1  ;
    prev_var2     = var2  ;
    
    linecounter++;

    if( 42 != 42 ) {
      cout
        << linecounter << "  "
        << var1   <<"  "   
        << var2   <<"  "   <<endl;   
    }
     
     varVsRun[var1] = var2;

  }
  //  std::cout<<"lines read: "<<linecounter<<std::endl;
}



#endif
