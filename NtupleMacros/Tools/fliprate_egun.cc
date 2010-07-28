//-------------------------------------------------------------------
// Returns the sign flip rate from the electron gun
// Note that the electron gun only goes up to 100 GeV
// For electron pt>100 returns the flip rate in the 90-100 Gev bin
// For electron pt<10  returns zero
// For electron abs(eta)>2.4 returns zero & prints out an error
//
// Usage:
// double flipRate   = getSingleEleFlipRate(el_pt, el_eta)
// double fliRateErr = getSingleEleFlipRateError(el_pt, el_eta)
//
// Claudio & Derek 23 July 2009
//-------------------------------------------------------------
#include "fliprate_egun.h"
#include <iostream>
#include <stdio.h>
#include <math.h>

using namespace std;

double getSingleEleNum(double el_pt, double el_eta) {

  el_eta = fabs(el_eta);

  if( el_eta < 1.28 ){ 
    if( el_pt > 70 ) return 611;
    if( el_pt > 50 ) return 390;
    if( el_pt > 40 ) return 163;
    if( el_pt > 30 ) return 140;
    if( el_pt > 10 ) return 123;
    return 0.0;
  }
  if( el_eta < 1.56 ){ 
    if( el_pt > 70 ) return 397;
    if( el_pt > 50 ) return 226;
    if( el_pt > 40 ) return 85;
    if( el_pt > 30 ) return 61;
    if( el_pt > 10 ) return 55;
    return 0.0;
  }
  if( el_eta < 1.84 ){ 
    if( el_pt > 70 ) return 709;
    if( el_pt > 50 ) return 380;
    if( el_pt > 40 ) return 153;
    if( el_pt > 30 ) return 113;
    if( el_pt > 10 ) return 69;
    return 0.0;
  }
  if( el_eta < 2.12 ){ 
    if( el_pt > 70 ) return 841;
    if( el_pt > 50 ) return 447;
    if( el_pt > 40 ) return 196;
    if( el_pt > 30 ) return 143;
    if( el_pt > 10 ) return 106;
    return 0.0;
  }
  if( el_eta < = 2.5 ){ 
    if( el_pt > 70 ) return 1237;
    if( el_pt > 50 ) return 636;
    if( el_pt > 40 ) return 301;
    if( el_pt > 30 ) return 235;
    if( el_pt > 10 ) return 294;
    return 0.0;
  }
  std::cout << "Error: eta > 2.5 value found" << endl;
  return 0.0;
}


double getSingleEleDenom(double el_pt, double el_eta) {

  el_eta = fabs(el_eta);

  if( el_eta < 1.28 ){ 
    if( el_pt > 70 ) return 118827;
    if( el_pt > 50 ) return 79313;
    if( el_pt > 40 ) return 39510;
    if( el_pt > 30 ) return 39205;
    if( el_pt > 10 ) return 69763;
    return 0.0;
  }
  if( el_eta < 1.56 ){ 
    if( el_pt > 70 ) return 19732;
    if( el_pt > 50 ) return 13022;
    if( el_pt > 40 ) return 6391;
    if( el_pt > 30 ) return 6060;
    if( el_pt > 10 ) return 9686;
    return 0.0;
  }
  if( el_eta < 1.84 ){ 
    if( el_pt > 70 ) return 22468;
    if( el_pt > 50 ) return 15126;
    if( el_pt > 40 ) return 7377;
    if( el_pt > 30 ) return 7148;
    if( el_pt > 10 ) return 10356;
    return 0.0;
  }
  if( el_eta < 2.12 ){ 
    if( el_pt > 70 ) return 23300;
    if( el_pt > 50 ) return 15541;
    if( el_pt > 40 ) return 7841;
    if( el_pt > 30 ) return 7594;
    if( el_pt > 10 ) return 12488;
    return 0.0;
  }
  if( el_eta < = 2.5 ){ 
    if( el_pt > 70 ) return 28736;
    if( el_pt > 50 ) return 18739;
    if( el_pt > 40 ) return 9310;
    if( el_pt > 30 ) return 8982;
    if( el_pt > 10 ) return 16488;
    return 0.0;
  }
  std::cout << "Error: eta > 2.5 value found" << endl;
  return 0.0;
}


double getSingleEleFlipRate(double el_pt, double el_eta) {
  if( el_pt < 10.0 || fabs(el_eta) > 2.5 ){
    std::cout << "Error in 'getSingleEleFlipRate': pt or eta value found out of range" << endl;
    return 0.0;
  }
  return getSingleEleNum(el_pt, fabs(el_eta))/getSingleEleDenom(el_pt, fabs(el_eta));
}

double getSingleEleFlipRateError(double el_pt, double el_eta) {
  //the binomial error
  if( el_pt < 10.0 || fabs(el_eta) > 2.5 ){
    std::cout << "Error in 'getSingleEleFlipRate': pt or eta value found out of range" << endl;
    return 0.0;
  }
  double num   = getSingleEleNum(el_pt,   fabs(el_eta));
  double denom = getSingleEleDenom(el_pt, fabs(el_eta));
  double p = num/denom;
  return sqrt(p*(1-p)/denom);
}
