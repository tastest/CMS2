
#ifndef PLOTRESULTS_H
#define PLOTRESULTS_H

#include "TString.h"

class TCanvas;
class HistogramUtilities;
class THStack;
class TArrow;

void plotStack(HistogramUtilities &h1, TString name, TString titleX, TString saveName, TString det, TString norm, int rebin = 1, float cutValEB = -1.0, float cutValEE = -1.0);
void plotDataRefOverlayStack(HistogramUtilities &hData, HistogramUtilities &hRef, TString name, TString titleX, TString saveName, TString det, TString norm, int rebin = 1,  float cutValEB = -1.0, float cutValEE = -1.0);
TArrow *getArrow(THStack *st, TString det, float cutValEB, float cutValEE, float max = -1.0);

void plotResults(TString hyp, TString version, TString fileStamp, TString refFileStamp, TString norm, const float &luminorm);

void makeTables(TString fileStamp, TString refFileStamp, TString norm, const float &luminorm);

#endif

