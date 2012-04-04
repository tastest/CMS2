
#ifndef VTXREWEIGHT_H
#define VTXREWEIGHT_H

#include "TH1.h"

float vtxweight( const bool isData = false , const bool useDAVertices = true , const bool usealt = false );
float vtxweight_n( const int nvtx , const bool isData = false , const bool usealt = false );
void set_vtxreweight_rootfile ( const char* filename = "vtxreweight.root" , bool verbose = false , const bool usealt = false );
//note: this hist should have the first bin be 1 vertex so that getnbinsx corresponds exactly to nvertices
TH1F* vtxreweight_hist;
TH1F* vtxreweight_hist_alt;

#endif

