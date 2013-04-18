#ifndef PLOTUTILITIES_H
#define PLOTUTILITIES_H

#include "TH1F.h"

typedef TH1F H;

// for interactive stuff
// comes from old "histtools.C"

void saveHist(const char* filename, const char* pat="*");
void deleteHistos();

#endif

