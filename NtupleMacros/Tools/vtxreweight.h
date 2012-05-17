#ifndef VTXREWEIGHT_H
#define VTXREWEIGHT_H

float vtxweight( const bool isData = false , const bool usealt = false );
float vtxweight_n( const int nvtx , const bool isData = false , const bool usealt = false );
void set_vtxreweight_rootfile ( const char* filename = "vtxreweight.root" , bool verbose = false , const bool usealt = false );

#endif

