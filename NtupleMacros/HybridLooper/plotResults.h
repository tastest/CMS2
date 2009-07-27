
#ifndef PLOTRESULTS_H
#define PLOTRESULTS_H

#include "TString.h"
class HistogramUtilities;

void plotEff(HistogramUtilities &h1, TString name, TString det);
void plotResults(TString det);

#endif

