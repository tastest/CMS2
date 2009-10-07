
#ifndef PLOTRESULTS_H
#define PLOTRESULTS_H

#include "TString.h"
#include "TH2F.h"

class HistogramUtilities;


void plotResults();

void plotResultsDilep(TString hyp);
void plotResultsLep(TString hyp);

double abcd(TH2F* h, double x1=0., double x2=0., double x3=0., double x4=0., double y1=0., double y2=0., double y3=0., double y4=0.);
double integrateTH2F(TH2F* h, double xlow=0., double xhgh=0., double ylow=0., double yhgh=0.);
void printTH2F(TH2F* h);
void compareTH2F(TH2F* h1, TH2F* h2);


#endif

