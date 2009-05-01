
#include "EffH2F.h"

EffH2F::EffH2F( char* name, char* title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup)
{
  numer = new TH2F( Form("%s_%s", name, "numer"), Form("%s_%s", name, "numer"), nbinsx, xlow, xup, nbinsy, ylow, yup );
  denom = new TH2F( Form("%s_%s", name, "denom"), Form("%s_%s", name, "denom"), nbinsx, xlow, xup, nbinsy, ylow, yup );
  eff = new TH2F( Form("%s_%s", name, "eff"), Form("%s_%s", name, "eff"), nbinsx, xlow, xup, nbinsy, ylow, yup );

  numer->Sumw2();
  denom->Sumw2();
  eff->Sumw2();
}

EffH2F::~EffH2F() {
  delete numer;
  delete denom;
  delete eff;
}

void EffH2F::MakeEff() {
  eff->Divide( numer, denom );
}

