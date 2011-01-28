#ifndef GATHER_H
#define GATHER_H

#include "TH1F.h"
#include "TCut.h"
#include "TChain.h"
#include "TCanvas.h"
#include "BabySample.h"

#include <vector>

//
// utility functions
//

bool sortHistsByIntegral(TH1* h1, TH1* h2);
void addToLegend(TLegend *leg, TH1F *hist, TString opt);
void makeStack(std::vector<TH1F*> &v_hists);
TH1F* slideIntegrated(TH1F* h1);

//
// luminosity functions
//

float GetIntLumi(TChain *c, float lumi, int brun, int bls, int erun, int els);
float GetIntLumi(BabySample *bs, float lumi);

//
// drawing functions
//

TCanvas* DrawAll(TCut field, const char *savename, TCut sel, TCut presel, float intlumipb, unsigned int nbins, float xlo, float xhi, bool integrated,
           std::vector<BabySample*> bss);

TH1F* Plot(const BabySample *bs, TCut var, TCut selection, float intlumipb,
        unsigned int nbins, float xlo, float xhi, bool integrated, unsigned int isfx);

static unsigned int gDrawAllCount = 0;

#endif
