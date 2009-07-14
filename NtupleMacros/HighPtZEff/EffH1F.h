// -*- C++ -*-    
//my (warren) attempt to make an efficiency histo class

#include "TH1F.h"
#include "TGraphAsymmErrors.h"

//using namespace std;

#ifndef EFFH1F_H
#define EFFH1F_H

class EffH1F {
  
public:
  EffH1F() { }
  EffH1F( char* name, char* title, Int_t nbinsx, Double_t xlow, Double_t xup);
  ~EffH1F ( );

  void MakeEff(const double ymin = 0.7, const bool rebin = false, const Float_t n = 500);

  //data members  
  TH1F* numer;
  TH1F* denom;
  TH1F* eff;
  TGraphAsymmErrors *gr_eff ;

};

#endif
