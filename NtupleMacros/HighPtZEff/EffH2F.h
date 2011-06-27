// -*- C++ -*-    
//my (warren) attempt to make an efficiency histo class

#include "TH2F.h"
//#include <stdio.h>

#ifndef EFFH2F_H
#define EFFH2F_H

class EffH2F {
  
public:
  EffH2F() { }
  EffH2F( char* name, char* title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup);
  ~EffH2F ( );

  void MakeEff();

  //data members  
  TH2F* numer;
  TH2F* denom;
  TH2F* eff;

};

#endif
