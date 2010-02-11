
#ifndef PLOTRESULTS_H
#define PLOTRESULTS_H

#include "TString.h"

class TCanvas;
class HistogramUtilities;
class THStack;
class TArrow;


void plotNormOverlay(HistogramUtilities &h1, TString name, TString saveName, TString det, int rebin, bool legendOnRight, float cutValEB, float cutValEE);
void plotEff(HistogramUtilities &h1, TString name, TString saveName, TString det, bool ascending, int rebin = 1, bool legendOnRight = true, float cutValEB = -1.0, float cutValEE = -1.0);
void plotStack(HistogramUtilities &h1, TString name, TString titleX, TString saveName, TString det, int rebin = 1, float cutValEB = -1.0, float cutValEE = -1.0);
TArrow *getArrow(THStack *st, TString det, float cutValEB, float cutValEE);
void plot2DSB(HistogramUtilities &h1, TString name, TString xTitle, TString yTitle, TString saveName, TString det);

void plotResultsID(TString det, TString hyp, TString fileStamp, TString saveName);

#endif

