
#ifndef PLOTRESULTS_H
#define PLOTRESULTS_H

#include "TString.h"

class TCanvas;
class HistogramUtilities;
class THStack;
class TArrow;

void plotEff(HistogramUtilities &h1, TString name, TString saveName, TString det, bool ascending, int rebin = 1, bool legendOnRight = true, float cutValEB = -1.0, float cutValEE = -1.0);
void plotStack(HistogramUtilities &h1, TString name, TString saveName, TString det, int rebin = 1, float cutValEB = -1.0, float cutValEE = -1.0);


TArrow *getArrow(TString det, THStack *st,float cutValEB, float cutValEE);


void plotAllResults();
void plotResults(TString det);

void plotAllResultsW();
void plotResultsW(TString det, TString fileStamp);


void test();

#endif

