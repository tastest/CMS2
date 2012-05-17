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
// WHERE NVTX IS THE NUMBER OF VERTICES PASSING isGoodVertex()
//------------------------------------------------------------------------------

// $Id: vtxreweight.cc,v 1.8 2012/05/17 10:33:26 benhoob Exp $

// CINT is allowed to see this, but nothing else:
#include "vtxreweight.h"

#ifndef __CINT__
#include <iostream>
#include <iomanip>
#include "../CORE/CMS2.h"
#include "../CORE/eventSelections.h"
#include "../CORE/trackSelections.h"
#include "TH1.h"
#include "TFile.h"


TH1F* vtxreweight_hist;
TH1F* vtxreweight_hist_alt;
//note: this hist should have the first bin be 1 vertex so that getnbinsx corresponds exactly to nvertices

bool loaded_vtxreweight_hist = false;
bool loaded_vtxreweight_hist_alt = false;

float vtxweight_n( const int nvertices , const bool isData , const bool usealt ) {

  if( isData ) return 1;

  if( (!usealt && !loaded_vtxreweight_hist) || 
	  (usealt && !loaded_vtxreweight_hist_alt) ){
    std::cout << "vtxreweight.cc: You need to do"                              << std::endl;
    std::cout << "set_vtxreweight_rootfile( filename )"                        << std::endl;
    std::cout << "before calling vtxweight()"                                  << std::endl;
    std::cout << "a sample vtxreweight file can be found at"                   << std::endl;
    std::cout << "/tas/benhoob/vtxreweight/vtxreweight_Spring11MC_23pbPR.root" << std::endl;
    std::cout << "now, quitting"                                               << std::endl;
    exit(2);
  }

  int nvtx = nvertices;
  TH1F* hist = (usealt ? vtxreweight_hist_alt : vtxreweight_hist );
  if( nvtx > hist->GetNbinsX() )
	nvtx = hist->GetNbinsX();

  float weight = 0;
  weight = hist->GetBinContent( hist->FindBin(nvtx) );
  if( weight <= 0 ) //we don't want to kill events bc they have no weight
	weight = 1.;
  //cout << "nvtx " << nvtx << " weight " << weight << endl;
  return weight;

}


//this is the original version--it used to count vertices and then get the weight
//move the weight to the above (new) version
float vtxweight( const bool isData , const bool usealt ){

  if( isData ) return 1;

  if( (!usealt && !loaded_vtxreweight_hist) || 
	  (usealt && !loaded_vtxreweight_hist_alt) ){
    std::cout << "vtxreweight.cc: You need to do"                              << std::endl;
    std::cout << "set_vtxreweight_rootfile( filename )"                        << std::endl;
    std::cout << "before calling vtxweight()"                                  << std::endl;
    std::cout << "a sample vtxreweight file can be found at"                   << std::endl;
    std::cout << "/tas/benhoob/vtxreweight/vtxreweight_Spring11MC_23pbPR.root" << std::endl;
    std::cout << "now, quitting"                                               << std::endl;
    exit(2);
  }

  int nvtx = 0;
  
  //cout << "Counting standard vertices" << endl;
  for (size_t v = 0; v < cms2.vtxs_position().size(); ++v){
    if(isGoodVertex(v)) ++nvtx;
  }

  if( nvtx == 0 ){
    cout << "vtxreweight.cc: warning 0 good vertices found, returning 0 weight" << endl;
    return 0;
  }

  return vtxweight_n( nvtx, isData, usealt );
}

void set_vtxreweight_rootfile ( const char* filename , const bool verbose , const bool usealt ){
  TFile* file = TFile::Open(filename);

  if( file == 0 ){
    cout << "vtxreweight.cc: error, couldn't open file : " << filename << endl;
    cout << "Quitting!" << endl;
    exit(0);
  }

  if( !usealt )
	vtxreweight_hist = (TH1F*) file->Get("hratio");
  else 
	vtxreweight_hist_alt = (TH1F*) file->Get("hratio_alt");
  TH1F* hist = (usealt ? vtxreweight_hist_alt : vtxreweight_hist );

  if( hist == 0 ){
    cout << "vtxreweight.cc: error, couldn't retrieve hist hratio"
		 << (usealt ? "_alt" : "") << " from file " << filename << endl;
    cout << "Quitting" << endl;
    exit(1);
  }

  if( verbose ){
    cout << "Opened vtx reweighting file " << filename << endl;
    
    cout << "|" << std::setw(10) << "nvtx"   << std::setw(4) 
		 << "|" << std::setw(10) << "weight" << std::setw(4) << "|" << endl;

    for(unsigned int ibin = 1 ; ibin <= (unsigned int)hist->GetNbinsX() ; ++ibin ){

      cout << "|" << std::setw(10) << ibin << std::setw(4) 
		   << "|" << std::setw(10) << hist->GetBinContent(ibin) << std::setw(4) << "|" << endl;

    }
  }

  if( !usealt )
	loaded_vtxreweight_hist = true;
  else
	loaded_vtxreweight_hist_alt = true;

}

#endif // __CUNT__

