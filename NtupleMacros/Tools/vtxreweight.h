
#ifndef VTXREWEIGHT_H
#define VTXREWEIGHT_H

float vtxweight( bool isData = false );
void set_vtxreweight_rootfile ( const char* filename = "vtxreweight.root", bool verbose = false );
TH1F* vtxreweight_hist;

#endif

