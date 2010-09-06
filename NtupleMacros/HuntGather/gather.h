#ifndef GATHER_H
#define GATHER_H

class TChain;
class TH1F;

float GetIntLumi(float lumi, int brun, int bls, int erun, int els);
float GetIntLumi(float lumi);
TH1F* Plot(const char *pfx, TChain *chain, const char *field, TCut sel, TCut presel, float intlumifb, float kfactor,
           unsigned int nbins, float xlo, float xhi);
// These should only be used with data, where intlumifb need not be specified and 
// most likely you aren't scaling
TH1F* Plot(const char *pfx, TChain *chain, const char *field, TCut sel, TCut presel, float kfactor,
           unsigned int nbins, float xlo, float xhi);
TH1F* Plot(const char *pfx, TChain *chain, const char *field, TCut sel, TCut presel,
           unsigned int nbins, float xlo, float xhi);

#endif
