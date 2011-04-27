//-------------------------------------------------------------------------------
// Use this utility to reweight MC to match the nvtx distribution in data.
// This requires a root file with a histogram named "hratio" whose bin contents 
// specify the weight for each nvtx bin. A sample file can be found at:
// /tas/benhoob/vtxreweight/vtxreweight_Spring11MC_23pbPR.root
//
// You can produce your own root file using the macro:
//
// usage:
//   
//    //include header
//    #include "../Tools/vtxreweight.cc"
//
//    //initialize
//    bool verbose = true;
//    set_vtxreweight_rootfile("vtxreweight.root",verbose);
//
//    //in the event loop
//    float myvtxweight = vtxweight();
//
// 


// $Id: vtxreweight.cc,v 1.1 2011/04/27 18:28:20 benhoob Exp $

// CINT is allowed to see this, but nothing else:
#include "vtxreweight.h"

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

#include "../CORE/CMS2.h"
#include "../CORE/eventSelections.h"

bool loaded_vtxreweight_hist = false;

float vtxweight( bool isData ){

  if( isData ) return 1;

  if( !loaded_vtxreweight_hist ){
    cout << "vtxreweight.cc: You need to do"                              << endl;
    cout << "set_vtxreweight_rootfile( filename )"                        << endl;
    cout << "before calling vtxweight()"                                  << endl;
    cout << "a sample vtxreweight file can be found at"                   << endl;
    cout << "/tas/benhoob/vtxreweight/vtxreweight_Spring11MC_23pbPR.root" << endl;
    cout << "now, quitting"                                               << endl;
    exit(2);
  }

  int ndavtx = 0;
    
  for (size_t v = 0; v < cms2.davtxs_position().size(); ++v){
    if(isGoodDAVertex(v)) ++ndavtx;
  }

  if( ndavtx == 0 ){
    cout << "vtxreweight.cc: warning 0 good vertices found, returning 0 weight" << endl;
    return 0;
  }

  if( ndavtx > vtxreweight_hist->GetNbinsX() ) ndavtx = vtxreweight_hist->GetNbinsX();

  return vtxreweight_hist->GetBinContent(ndavtx);

}

void set_vtxreweight_rootfile ( const char* filename , bool verbose ){
  TFile* file = TFile::Open(filename);

  if( file == 0 ){
    cout << "vtxreweight.cc: error, couldn't open file : " << filename << endl;
    cout << "Quitting!" << endl;
    exit(0);
  }

  vtxreweight_hist = (TH1F*) file->Get("hratio");

  if( vtxreweight_hist == 0 ){
    cout << "vtxreweight.cc: error, couldn't retrieve hist hratio from file " << filename << endl;
    cout << "Quitting" << endl;
    exit(1);
  }

  if( verbose ){
    cout << "Opened vtx reweighting file " << filename << endl;
    
    cout << "|" << setw(10) << "nvtx"   << setw(4) 
	 << "|" << setw(10) << "weight" << setw(4) << "|" << endl;

    for(unsigned int ibin = 1 ; ibin <= vtxreweight_hist->GetNbinsX() ; ++ibin ){

      cout << "|" << setw(10) << ibin                                   << setw(4) 
	   << "|" << setw(10) << vtxreweight_hist->GetBinContent(ibin) << setw(4) << "|" << endl;

    }

  loaded_vtxreweight_hist = true;
  }
}

#endif // __CUNT__

