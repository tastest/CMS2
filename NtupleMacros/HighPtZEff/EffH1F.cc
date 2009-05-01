
#include "EffH1F.h"



EffH1F::EffH1F( char* name, char* title, Int_t nbinsx, Double_t xlow, Double_t xup )
{
  numer = new TH1F( Form("%s_%s", name, "numer"), Form("%s_%s", name, "numer"), nbinsx, xlow, xup );
  denom = new TH1F( Form("%s_%s", name, "denom"), Form("%s_%s", name, "denom"), nbinsx, xlow, xup );
  eff = new TH1F( Form("%s_%s", name, "eff"), Form("%s_%s", name, "eff"), nbinsx, xlow, xup );

  numer->Sumw2();
  denom->Sumw2();
  eff->Sumw2();
}

EffH1F::~EffH1F() {
  delete numer;
  delete denom;
  delete eff;
}

void EffH1F::MakeEff() {
  eff->Divide( numer, denom );
}

