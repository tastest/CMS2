//-------------------------------------------------------------------------------
// Use this utility to reweight MC to match the nvtx distribution in data.
// This requires a root file with a histogram named "hratio" whose bin contents 
// specify the weight for each nvtx bin. This weight is determined by plotting
// the nvtx distribution for a data sample and a MC sample and taking the ratio.
//
// A sample file can be found at:
// /tas/benhoob/vtxreweight/vtxreweight_Spring11MC_153pb_Zselection.root
//
// This root file was made with this macro: 
// /tas/benhoob/vtxreweight/make_vtxreweight_Spring11MC_153pb_Zselection.cc
//
// and was produced after applying a Z selection.
//
// You can make your own root file by modifying this macro
//
// 
// usage:
//   
//    //include header
//    #include "../Tools/vtxreweight.cc"
//
//    //initialize
//    bool verbose = true;
//    char* vtxfile = "/tas/benhoob/vtxreweight/vtxreweight_Spring11MC_153pb_Zselection.root";
//    set_vtxreweight_rootfile( vtxfile , verbose );
//
//    //in the event loop
//    float myvtxweight = vtxweight();
//
// 
// PLEASE NOTE: ALWAYS CHECK THAT THE WEIGHTING HAS BEEN DONE PROPERLY
// BY COMPARING THE DATA NVTX DIST WITH THE WEIGHTED MC NVTX DISTRIBUTION,
// WHERE NVTX IS THE NUMBER OF DA VERTICES PASSING isGoodDAVertex()
//------------------------------------------------------------------------------

// $Id: msugraCrossSection.cc,v 1.1 2011/07/07 12:11:51 benhoob Exp $

// CINT is allowed to see this, but nothing else:
#include "msugraCrossSection.h"

#ifndef __CINT__

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <set>
#include <string>

//#include "../CORE/CMS2.h"
//#include "../CORE/eventSelections.h"

bool loaded_file = false;

void set_msugra_file ( const char* filename , bool verbose ){

  ifile.open(filename);

  if( !ifile.is_open() ){
    cout << "msugraCrossSection.cc: error, couldn't open file : " << filename << endl;
    cout << "Quitting!" << endl;
    exit(0);
  }

  if( verbose ){
    char line[200];

    while( !ifile.eof()){
      ifile.getline(line,200);
      cout << line << endl;
    }
  }

  loaded_file = true;

}

double getMsugraCrossSection( double my_m0 , double my_m12, double my_tanb , bool verbose ){

  if( !loaded_file ){
    cout << "musgraCrossSection.cc: You need to do"              << endl;
    cout << "set_msugra_file( filename )"                        << endl;
    cout << "before calling getCrossSection()"                   << endl;
    cout << "a sample file can be found at"                      << endl;
    cout << "/tas/benhoob/msugra/goodModelNames_tanbeta10.txt"   << endl;
    cout << "now, quitting"                                      << endl;
    exit(2);
  }
  
  ifile.clear();
  ifile.seekg(0);

  double m0, m12, tanb, A0, mu=1.0;
  double signMu, xsec;

  string line;
  string xsecstring;

  bool foundPoint = false;
  double my_xsec = -1;

  while( !ifile.eof() && !foundPoint ){

    ifile >> line;
    ifile >> xsecstring;

    TString model_params(line);
    if (!model_params.Contains("msugra"))
      continue;

    TObjArray* tokens = model_params.Tokenize("_");
    m0      = ((TObjString*)tokens->At(1))->GetString().Atof();
    m12     = ((TObjString*)tokens->At(2))->GetString().Atof();
    tanb    = ((TObjString*)tokens->At(3))->GetString().Atof();
    A0      = ((TObjString*)tokens->At(4))->GetString().Atof();
    mu      = ((TObjString*)tokens->At(5))->GetString().Atof();

    xsec = TString(xsecstring).Atof();
   
    if( fabs(m0-my_m0) < 0.1 && fabs(m12-my_m12) < 0.1 && fabs(tanb-my_tanb) < 0.1 ){

      if( verbose ){
	cout << "Found musgra point:" << endl;
	cout << Form("m0 = %.0f m1/2 = %.0f tanb = %.0f A = %.0f mu = %.0f xsec = %.10f",m0,m12,tanb,A0,mu,xsec) << endl;
      }

      foundPoint = true;
      my_xsec = xsec;
    }
  }

  if( verbose && my_xsec < 0 ){
    cout << "Didn't find point:" << endl;
    cout << Form("m0 = %.0f m1/2 = %.0f tanb = %.0f",my_m0,my_m12,my_tanb) << endl;
  }

  return 1E9 * my_xsec;

}


#endif // __CUNT__

