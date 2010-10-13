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
    if( el_pt > 70 ) return 31;
    if( el_pt > 50 ) return 20;
    if( el_pt > 40 ) return 5;
    if( el_pt > 30 ) return 9;
    if( el_pt > 10 ) return 4;
    return 0.0;
  }
  if( el_eta < 1.56 ){ 
    if( el_pt > 70 ) return 48;
    if( el_pt > 50 ) return 30;
    if( el_pt > 40 ) return 12;
    if( el_pt > 30 ) return 8;
    if( el_pt > 10 ) return 5;
    return 0.0;
  }
  if( el_eta < 1.84 ){ 
    if( el_pt > 70 ) return 92;
    if( el_pt > 50 ) return 67;
    if( el_pt > 40 ) return 35;
    if( el_pt > 30 ) return 20;
    if( el_pt > 10 ) return 13;
    return 0.0;
  }
  if( el_eta < 2.12 ){ 
    if( el_pt > 70 ) return 72;
    if( el_pt > 50 ) return 29;
    if( el_pt > 40 ) return 16;
    if( el_pt > 30 ) return 13;
    if( el_pt > 10 ) return 6;
    return 0.0;
  }
  if( el_eta < = 2.5 ){ 
    if( el_pt > 70 ) return 93;
    if( el_pt > 50 ) return 42;
    if( el_pt > 40 ) return 17;
    if( el_pt > 30 ) return 9;
    if( el_pt > 10 ) return 20;
    return 0.0;
  }
  std::cout << "Error: eta > 2.5 value found" << endl;
  return 0.0;
}


double getSingleEleDenom(double el_pt, double el_eta) {

  el_eta = fabs(el_eta);

  if( el_eta < 1.28 ){ 
    if( el_pt > 70 ) return 117286;
    if( el_pt > 50 ) return 78397;
    if( el_pt > 40 ) return 39090;
    if( el_pt > 30 ) return 38858;
    if( el_pt > 10 ) return 69323;
    return 0.0;
  }
  if( el_eta < 1.56 ){ 
    if( el_pt > 70 ) return 18452;
    if( el_pt > 50 ) return 12276;
    if( el_pt > 40 ) return 6087;
    if( el_pt > 30 ) return 5793;
    if( el_pt > 10 ) return 9374;
    return 0.0;
  }
  if( el_eta < 1.84 ){ 
    if( el_pt > 70 ) return 20720;
    if( el_pt > 50 ) return 14174;
    if( el_pt > 40 ) return 7021;
    if( el_pt > 30 ) return 6852;
    if( el_pt > 10 ) return 10145;
    return 0.0;
  }
  if( el_eta < 2.12 ){ 
    if( el_pt > 70 ) return 21436;
    if( el_pt > 50 ) return 14537;
    if( el_pt > 40 ) return 7401;
    if( el_pt > 30 ) return 7241;
    if( el_pt > 10 ) return 12178;
    return 0.0;
  }
  if( el_eta < = 2.5 ){ 
    if( el_pt > 70 ) return 25029;
    if( el_pt > 50 ) return 16802;
    if( el_pt > 40 ) return 8462;
    if( el_pt > 30 ) return 8277;
    if( el_pt > 10 ) return 15646;
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
