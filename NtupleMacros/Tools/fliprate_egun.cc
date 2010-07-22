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
    if( el_pt > 90.0 )	return 26.0;
    if( el_pt > 80.0 )	return 12.0;
    if( el_pt > 70.0 )	return 9.0;
    if( el_pt > 60.0 )	return 12.0;
    if( el_pt > 50.0 )	return  6.0;
    if( el_pt > 40.0 )	return  5.0;
    if( el_pt > 30.0 )	return  4.0;
    if( el_pt > 20.0 )	return  2.0;
    if( el_pt > 10.0 )	return  2.0;
    return 0.0;
  }
  if( el_eta < 1.0 ){
    if( el_pt > 90.0 )  return 32.0;
    if( el_pt > 80.0 )  return 31.0;
    if( el_pt > 70.0 )  return 22.0;
    if( el_pt > 60.0 )  return 19.0;
    if( el_pt > 50.0 )  return 17.0;
    if( el_pt > 40.0 )  return 10.0;
    if( el_pt > 30.0 )  return 12.0;
    if( el_pt > 20.0 )  return  3.0;
    if( el_pt > 10.0 )  return  1.0;
    return 0.0;
  }
  if( el_eta < 1.28 ){
    if( el_pt > 90.0 )  return 41.0;
    if( el_pt > 80.0 )  return 39.0;
    if( el_pt > 70.0 )  return 34.0;
    if( el_pt > 60.0 )  return 30.0;
    if( el_pt > 50.0 )  return 26.0;
    if( el_pt > 40.0 )  return 16.0;
    if( el_pt > 30.0 )  return 13.0;
    if( el_pt > 20.0 )  return  6.0;
    if( el_pt > 10.0 )  return  1.0;
    return 0.0;
  }
  if( el_eta < 1.56 ){
    if( el_pt > 90.0 )  return 74.0;
    if( el_pt > 80.0 )  return 88.0;
    if( el_pt > 70.0 )  return 90.0;
    if( el_pt > 60.0 )  return 64.0;
    if( el_pt > 50.0 )  return 57.0;
    if( el_pt > 40.0 )  return 42.0;
    if( el_pt > 30.0 )  return 29.0;
    if( el_pt > 20.0 )  return 11.0;
    if( el_pt > 10.0 )  return 2.0;
    return 0.0;
  }
  if( el_eta < 1.84 ){
    if( el_pt > 90.0 )  return 228.0;
    if( el_pt > 80.0 )  return 240.0;
    if( el_pt > 70.0 )  return 215.0;
    if( el_pt > 60.0 )  return 195.0;
    if( el_pt > 50.0 )  return 124.0;
    if( el_pt > 40.0 )  return 100.0;
    if( el_pt > 30.0 )  return 55.0;
    if( el_pt > 20.0 )  return 16.0;
    if( el_pt > 10.0 )  return 9.0;
    return 0.0;
  }
  if( el_eta < 2.12 ){
    if( el_pt > 90.0 )  return 165.0;
    if( el_pt > 80.0 )  return 199.0;
    if( el_pt > 70.0 )  return 166.0;
    if( el_pt > 60.0 )  return 160.0;
    if( el_pt > 50.0 )  return 112.0;
    if( el_pt > 40.0 )  return 89.0;
    if( el_pt > 30.0 )  return 60.0;
    if( el_pt > 20.0 )  return 28.0;
    if( el_pt > 10.0 )  return 6.0;
    return 0.0;
  }
  if( el_eta <= 2.5 ){
    if( el_pt > 90.0 )  return 250.0;
    if( el_pt > 80.0 )  return 250.0;
    if( el_pt > 70.0 )  return 227.0;
    if( el_pt > 60.0 )  return 211.0;
    if( el_pt > 50.0 )  return 187.0;
    if( el_pt > 40.0 )  return 168.0;
    if( el_pt > 30.0 )  return 124.0;
    if( el_pt > 20.0 )  return 64.0;
    if( el_pt > 10.0 )  return 32.0;
    return 0.0;
  }
  std::cout << "Error: eta > 2.5 value found" << endl; 
  return 0.0;
}


double getSingleEleDenom(double el_pt, double el_eta) {

  el_eta = fabs(el_eta);

  if( el_eta < 0.5 ){
    if( el_pt > 90.0 )	return 30479.0;
    if( el_pt > 80.0 )	return 31026.0;
    if( el_pt > 70.0 )	return 30869.0;
    if( el_pt > 60.0 )	return 30163.0;
    if( el_pt > 50.0 )	return 30171.0;
    if( el_pt > 40.0 )	return 30690.0;
    if( el_pt > 30.0 )	return 30450.0;
    if( el_pt > 20.0 )	return 29890.0;
    if( el_pt > 10.0 )	return 26416.0;
    return 0.0;
  }
  if( el_eta < 1.0 ){
    if( el_pt > 90.0 )  return 31752.0;
    if( el_pt > 80.0 )  return 32158.0;
    if( el_pt > 70.0 )  return 32041.0;
    if( el_pt > 60.0 )  return 32160.0;
    if( el_pt > 50.0 )  return 31503.0;
    if( el_pt > 40.0 )  return 31799.0;
    if( el_pt > 30.0 )  return 31506.0;
    if( el_pt > 20.0 )  return 29683.0;
    if( el_pt > 10.0 )  return 26905.0;
    return 0.0;
  }
  if( el_eta < 1.28 ){
    if( el_pt > 90.0 )  return 17909.0;
    if( el_pt > 80.0 )  return 18537.0;
    if( el_pt > 70.0 )  return 18780.0;
    if( el_pt > 60.0 )  return 18594.0;
    if( el_pt > 50.0 )  return 18140.0;
    if( el_pt > 40.0 )  return 18025.0;
    if( el_pt > 30.0 )  return 17392.0;
    if( el_pt > 20.0 )  return 16069.0;
    if( el_pt > 10.0 )  return 12778.0;
    return 0.0;
  }
  if( el_eta < 1.56 ){
    if( el_pt > 90.0 )  return 13680.0;
    if( el_pt > 80.0 )  return 14568.0;
    if( el_pt > 70.0 )  return 15104.0;
    if( el_pt > 60.0 )  return 14733.0;
    if( el_pt > 50.0 )  return 14777.0;
    if( el_pt > 40.0 )  return 14499.0;
    if( el_pt > 30.0 )  return 13572.0;
    if( el_pt > 20.0 )  return 12063.0;
    if( el_pt > 10.0 )  return 9166.0;
    return 0.0;
  }
  if( el_eta < 1.84 ){
    if( el_pt > 90.0 )  return 16990.0;
    if( el_pt > 80.0 )  return 17608.0;
    if( el_pt > 70.0 )  return 17624.0;
    if( el_pt > 60.0 )  return 17817.0;
    if( el_pt > 50.0 )  return 17688.0;
    if( el_pt > 40.0 )  return 17822.0;
    if( el_pt > 30.0 )  return 16628.0;
    if( el_pt > 20.0 )  return 14626.0;
    if( el_pt > 10.0 )  return 8920.0;
    return 0.0;
  }
  if( el_eta < 2.12 ){
    if( el_pt > 90.0 )  return 16554.0;
    if( el_pt > 80.0 )  return 17241.0;
    if( el_pt > 70.0 )  return 17441.0;
    if( el_pt > 60.0 )  return 17629.0;
    if( el_pt > 50.0 )  return 17659.0;
    if( el_pt > 40.0 )  return 17627.0;
    if( el_pt > 30.0 )  return 17311.0;
    if( el_pt > 20.0 )  return 16194.0;
    if( el_pt > 10.0 )  return 12052.0;
    return 0.0;
  }
  if( el_eta <= 2.5 ){
    if( el_pt > 90.0 )  return 18590.0;
    if( el_pt > 80.0 )  return 19566.0;
    if( el_pt > 70.0 )  return 19833.0;
    if( el_pt > 60.0 )  return 19576.0;
    if( el_pt > 50.0 )  return 20099.0;
    if( el_pt > 40.0 )  return 20491.0;
    if( el_pt > 30.0 )  return 20188.0;
    if( el_pt > 20.0 )  return 20162.0;
    if( el_pt > 10.0 )  return 17844.0;
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

