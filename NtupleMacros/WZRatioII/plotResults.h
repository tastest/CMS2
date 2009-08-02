
#ifndef PLOTRESULTS_H
#define PLOTRESULTS_H

#include "TString.h"
class HistogramUtilities;

void test();
void plotEffVar(HistogramUtilities &h1, TString name, Int_t rebin);


void plotEff(HistogramUtilities &h1, TString name, TString det, bool ascending);
void plotResults(TString det);



#endif

