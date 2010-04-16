
#ifndef PLOTRESULTS_H
#define PLOTRESULTS_H

#include "TString.h"

class TCanvas;
class HistogramUtilities;
class THStack;
class TArrow;
class TH1F;

void plotEff(HistogramUtilities &h1, TString name, TString saveName, TString det, bool ascending, int rebin = 1, bool legendOnRight = true, float cutValEB = -1.0, float cutValEE = -1.0);
void plotStack(HistogramUtilities &h1, TString name, TString titleX, TString saveName, TString det, int rebin = 1, float cutValEB = -1.0, float cutValEE = -1.0);


TArrow *getArrow(THStack *st, TString det, float cutValEB, float cutValEE);


void plot2DSB(HistogramUtilities &h1, TString name, TString xTitle, TString yTitle, TString saveName, TString det);

void plotResultsW(TString det, TString fileStamp, TString version);

void setError(TH1F *h1_numer, TH1F *h1_denom, TH1F *h1_eff);

void plotValidationOverlay(HistogramUtilities &h1, TString name_before, TString name_after, TString sn1, TString sn2, TString saveName, TString det, int rebin = 1,  bool plotDist = true, bool plotEff = true);

#endif

