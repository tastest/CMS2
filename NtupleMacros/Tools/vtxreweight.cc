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

// $Id: vtxreweight.cc,v 1.4 2011/05/20 14:24:22 warren Exp $

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

float vtxweight( bool isData , bool useDAVertices ){

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

  int nvtx = 0;
  
  //---------------------
  // count DA vertices
  //---------------------

  if( useDAVertices ){
    //cout << "Counting DA vertices" << endl;
    for (size_t v = 0; v < cms2.davtxs_position().size(); ++v){
      if(isGoodDAVertex(v)) ++nvtx;
    }
  }

  //---------------------
  // count DA vertices
  //---------------------

  else{
    //cout << "Counting standard vertices" << endl;
    for (size_t v = 0; v < cms2.vtxs_position().size(); ++v){
      if(isGoodVertex(v)) ++nvtx;
    }
  }

  if( nvtx == 0 ){
    cout << "vtxreweight.cc: warning 0 good vertices found, returning 0 weight" << endl;
    return 0;
  }

  if( nvtx > vtxreweight_hist->GetNbinsX() ) nvtx = vtxreweight_hist->GetNbinsX();

  //cout << "nvtx " << nvtx << " weight " << vtxreweight_hist->GetBinContent(nvtx) << endl;
  return vtxreweight_hist->GetBinContent(nvtx);

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

    for(unsigned int ibin = 1 ; ibin <= (unsigned int)vtxreweight_hist->GetNbinsX() ; ++ibin ){

      cout << "|" << setw(10) << ibin                                   << setw(4) 
	   << "|" << setw(10) << vtxreweight_hist->GetBinContent(ibin) << setw(4) << "|" << endl;

    }
  }
  loaded_vtxreweight_hist = true;

}

#endif // __CUNT__

