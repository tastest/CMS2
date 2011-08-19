#ifndef GATHER_H
#define GATHER_H

#include "TH1F.h"
#include "TCut.h"
#include "TChain.h"
#include "TCanvas.h"
#include "BabySample.h"

#include <vector>

class TGraphAsymmErrors;

//
// utility functions
//

bool sortHistsByIntegral(TH1* h1, TH1* h2);
void addToLegend(TLegend *leg, TH1F *hist, TString opt);
void makeStack(std::vector<TH1F*> &v_hists);
TH1F* slideIntegrated(TH1F* h1);
void PreselectBabies(std::vector<BabySample*> bss, TCut cut);
void PrintBins(const TH1F *h1);
void PrintBins(const TGraphAsymmErrors *h1);
void ZeroBinError(TH1F *h1);
void ComputeEfficiency(TH1F *gr_eff_data, const TH1F *h1_data_denom, const TH1F *h1_data_numer);

//
// luminosity functions
//
float GetAllLumi(TChain *c, float zPerPb);
float GetNewLumi(TChain *c, float zPerPb);
float GetZPerPb(TChain *c, float lumi);

float GetIntLumi(TChain *c, float lumi, int brun, int bls, int erun, int els);
float GetIntLumi(BabySample *bs, float lumi);

//
// drawing functions
//

TCanvas* TagAndProbe(const char *savename, TCut evsel, TCut var1, TCut var2,
        TCut tag1, TCut tag2, TCut probe1, TCut probe2, TCut sel1, TCut sel2,
        float intlumipb, unsigned int nbins, float xlo, float xhi, bool integrated, std::vector<BabySample*> bss);


TCanvas* TriggerMonitor(const char *savename, TCut var, TCut sel, TCut trig, float intlumipb, unsigned int nbins, float xlo, float xhi, bool integrated, 
            BabySample *bs);

TCanvas* DrawAll(TCut field, const char *savename, TCut sel, float intlumipb, unsigned int nbins, float xlo, float xhi, bool integrated,
           std::vector<BabySample*> bss);

TH1F* Plot(BabySample *bs, TCut var, TCut selection, float intlumipb,
        unsigned int nbins, float xlo, float xhi, bool integrated, unsigned int isfx);

static unsigned int gDrawAllCount = 0;

#endif
