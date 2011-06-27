
#ifndef HISTTOOLS_H
#define HISTTOOLS_H

#include "TH1.h"
#include "TGraph.h"

typedef TH1F H;

H cumulate (const H &in, bool increasing); 
TGraph eff_rej (const H &signal, H &background, bool normalize, bool increasing);

#endif

