
#ifndef PLOTRESULTS_H
#define PLOTRESULTS_H

#include "TString.h"
class HistogramUtilities;

void plotEff(HistogramUtilities &h1, TString name, TString det, bool ascending);
void plotResults(TString det);
void test();

#endif

