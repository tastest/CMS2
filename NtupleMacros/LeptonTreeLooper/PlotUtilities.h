#ifndef PLOTUTILITIES_H
#define PLOTUTILITIES_H

#include "TH1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TFile.h"

typedef TH1F H;

H cumulate (const H &in, bool increasing);
TGraph eff_rej (const H &signal, H &background, bool normalize, bool increasing);
void saveHist(const char* filename, const char* pat="*");
void deleteHistos();
TCanvas *ComparePlots(TFile *f, const char *hist1, const char *hist2, const char *label1, const char *label2, unsigned int rebin, bool norm);

#endif

