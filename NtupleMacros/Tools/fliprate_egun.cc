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

  if( el_eta < 0.5 ){ 
    if( el_pt > 90) return 0;
    if( el_pt > 80) return 2;
    if( el_pt > 70) return 2;
    if( el_pt > 60) return 3;
    if( el_pt > 50) return 2;
    if( el_pt > 40) return 1;
    if( el_pt > 30) return 3;
    if( el_pt > 20) return 0;
    if( el_pt > 10) return 0;
    return 0.0;
  }
  if( el_eta < 1 ){ 
    if( el_pt > 90) return 8;
    if( el_pt > 80) return 3;
    if( el_pt > 70) return 3;
    if( el_pt > 60) return 4;
    if( el_pt > 50) return 5;
    if( el_pt > 40) return 2;
    if( el_pt > 30) return 6;
    if( el_pt > 20) return 2;
    if( el_pt > 10) return 1;
    return 0.0;
  }
  if( el_eta < 1.28 ){ 
    if( el_pt > 90) return 12;
    if( el_pt > 80) return 12;
    if( el_pt > 70) return 8;
    if( el_pt > 60) return 15;
    if( el_pt > 50) return 10;
    if( el_pt > 40) return 6;
    if( el_pt > 30) return 10;
    if( el_pt > 20) return 2;
    if( el_pt > 10) return 1;
    return 0.0;
  }
  if( el_eta < 1.56 ){ 
    if( el_pt > 90) return 32;
    if( el_pt > 80) return 32;
    if( el_pt > 70) return 22;
    if( el_pt > 60) return 22;
    if( el_pt > 50) return 19;
    if( el_pt > 40) return 17;
    if( el_pt > 30) return 9;
    if( el_pt > 20) return 6;
    if( el_pt > 10) return 1;
    return 0.0;
  }
  if( el_eta < 1.84 ){ 
    if( el_pt > 90) return 56;
    if( el_pt > 80) return 56;
    if( el_pt > 70) return 50;
    if( el_pt > 60) return 62;
    if( el_pt > 50) return 41;
    if( el_pt > 40) return 38;
    if( el_pt > 30) return 18;
    if( el_pt > 20) return 12;
    if( el_pt > 10) return 0;
    return 0.0;
  }
  if( el_eta < 2.12 ){ 
    if( el_pt > 90) return 50;
    if( el_pt > 80) return 35;
    if( el_pt > 70) return 36;
    if( el_pt > 60) return 36;
    if( el_pt > 50) return 22;
    if( el_pt > 40) return 26;
    if( el_pt > 30) return 13;
    if( el_pt > 20) return 6;
    if( el_pt > 10) return 4;
    return 0.0;
  }
  if( el_eta < 2.5 ){ 
    if( el_pt > 90) return 47;
    if( el_pt > 80) return 56;
    if( el_pt > 70) return 48;
    if( el_pt > 60) return 37;
    if( el_pt > 50) return 34;
    if( el_pt > 40) return 29;
    if( el_pt > 30) return 9;
    if( el_pt > 20) return 21;
    if( el_pt > 10) return 4;
    return 0.0;
  }
  std::cout << "Error: eta > 2.5 value found" << endl; 
  return 0.0;
}


double getSingleEleDenom(double el_pt, double el_eta) {

  el_eta = fabs(el_eta);

  if( el_eta < 0.5 ){ 
    if( el_pt > 90) return 29598;
    if( el_pt > 80) return 29458;
    if( el_pt > 70) return 29894;
    if( el_pt > 60) return 29800;
    if( el_pt > 50) return 29817;
    if( el_pt > 40) return 29646;
    if( el_pt > 30) return 29468;
    if( el_pt > 20) return 28424;
    if( el_pt > 10) return 23816;
    return 0.0;
  }
  if( el_eta < 1 ){ 
    if( el_pt > 90) return 29226;
    if( el_pt > 80) return 29247;
    if( el_pt > 70) return 29497;
    if( el_pt > 60) return 29223;
    if( el_pt > 50) return 29944;
    if( el_pt > 40) return 29770;
    if( el_pt > 30) return 29486;
    if( el_pt > 20) return 29112;
    if( el_pt > 10) return 24737;
    return 0.0;
  }
  if( el_eta < 1.28 ){ 
    if( el_pt > 90) return 16042;
    if( el_pt > 80) return 16359;
    if( el_pt > 70) return 16036;
    if( el_pt > 60) return 16218;
    if( el_pt > 50) return 16093;
    if( el_pt > 40) return 16086;
    if( el_pt > 30) return 15891;
    if( el_pt > 20) return 14635;
    if( el_pt > 10) return 11896;
    return 0.0;
  }
  if( el_eta < 1.56 ){ 
    if( el_pt > 90) return 12192;
    if( el_pt > 80) return 12205;
    if( el_pt > 70) return 11799;
    if( el_pt > 60) return 12129;
    if( el_pt > 50) return 11795;
    if( el_pt > 40) return 11776;
    if( el_pt > 30) return 11147;
    if( el_pt > 20) return 10052;
    if( el_pt > 10) return 7545;
    return 0.0;
  }
  if( el_eta < 1.84 ){ 
    if( el_pt > 90) return 12990;
    if( el_pt > 80) return 13375;
    if( el_pt > 70) return 13431;
    if( el_pt > 60) return 13415;
    if( el_pt > 50) return 13127;
    if( el_pt > 40) return 12808;
    if( el_pt > 30) return 12308;
    if( el_pt > 20) return 10435;
    if( el_pt > 10) return 6272;
    return 0.0;
  }
  if( el_eta < 2.12 ){ 
    if( el_pt > 90) return 13452;
    if( el_pt > 80) return 13478;
    if( el_pt > 70) return 13617;
    if( el_pt > 60) return 13479;
    if( el_pt > 50) return 13585;
    if( el_pt > 40) return 13490;
    if( el_pt > 30) return 13059;
    if( el_pt > 20) return 12067;
    if( el_pt > 10) return 8339;
    return 0.0;
  }
  if( el_eta < 2.5 ){ 
    if( el_pt > 90) return 14884;
    if( el_pt > 80) return 15043;
    if( el_pt > 70) return 15072;
    if( el_pt > 60) return 14891;
    if( el_pt > 50) return 14931;
    if( el_pt > 40) return 15031;
    if( el_pt > 30) return 14592;
    if( el_pt > 20) return 14307;
    if( el_pt > 10) return 12050;
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
