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
    if( el_pt > 80.0 )	return 13.0;
    if( el_pt > 70.0 )	return 10.0;
    if( el_pt > 60.0 )	return 12.0;
    if( el_pt > 50.0 )	return  8.0;
    if( el_pt > 40.0 )	return  6.0;
    if( el_pt > 30.0 )	return  5.0;
    if( el_pt > 20.0 )	return  3.0;
    if( el_pt > 10.0 )	return  2.0;
    return 0.0;
  }
  if( el_eta < 1.0 ){
    if( el_pt > 90.0 )  return 32.0;
    if( el_pt > 80.0 )  return 34.0;
    if( el_pt > 70.0 )  return 23.0;
    if( el_pt > 60.0 )  return 20.0;
    if( el_pt > 50.0 )  return 17.0;
    if( el_pt > 40.0 )  return 13.0;
    if( el_pt > 30.0 )  return 14.0;
    if( el_pt > 20.0 )  return  4.0;
    if( el_pt > 10.0 )  return  6.0;
    return 0.0;
  }
  if( el_eta < 1.28 ){
    if( el_pt > 90.0 )  return 42.0;
    if( el_pt > 80.0 )  return 41.0;
    if( el_pt > 70.0 )  return 35.0;
    if( el_pt > 60.0 )  return 33.0;
    if( el_pt > 50.0 )  return 27.0;
    if( el_pt > 40.0 )  return 17.0;
    if( el_pt > 30.0 )  return 14.0;
    if( el_pt > 20.0 )  return  9.0;
    if( el_pt > 10.0 )  return  2.0;
    return 0.0;
  }
  if( el_eta < 1.56 ){
    if( el_pt > 90.0 )  return 82.0;
    if( el_pt > 80.0 )  return 93.0;
    if( el_pt > 70.0 )  return 96.0;
    if( el_pt > 60.0 )  return 69.0;
    if( el_pt > 50.0 )  return 66.0;
    if( el_pt > 40.0 )  return 50.0;
    if( el_pt > 30.0 )  return 34.0;
    if( el_pt > 20.0 )  return 16.0;
    if( el_pt > 10.0 )  return 8.0;
    return 0.0;
  }
  if( el_eta < 1.84 ){
    if( el_pt > 90.0 )  return 249.0;
    if( el_pt > 80.0 )  return 268.0;
    if( el_pt > 70.0 )  return 241.0;
    if( el_pt > 60.0 )  return 228.0;
    if( el_pt > 50.0 )  return 151.0;
    if( el_pt > 40.0 )  return 131.0;
    if( el_pt > 30.0 )  return 87.0;
    if( el_pt > 20.0 )  return 42.0;
    if( el_pt > 10.0 )  return 27.0;
    return 0.0;
  }
  if( el_eta < 2.12 ){
    if( el_pt > 90.0 )  return 183.0;
    if( el_pt > 80.0 )  return 229.0;
    if( el_pt > 70.0 )  return 188.0;
    if( el_pt > 60.0 )  return 177.0;
    if( el_pt > 50.0 )  return 144.0;
    if( el_pt > 40.0 )  return 118.0;
    if( el_pt > 30.0 )  return 92.0;
    if( el_pt > 20.0 )  return 48.0;
    if( el_pt > 10.0 )  return 17.0;
    return 0.0;
  }
  if( el_eta <= 2.5 ){
    if( el_pt > 90.0 )  return 288.0;
    if( el_pt > 80.0 )  return 298.0;
    if( el_pt > 70.0 )  return 290.0;
    if( el_pt > 60.0 )  return 262.0;
    if( el_pt > 50.0 )  return 251.0;
    if( el_pt > 40.0 )  return 227.0;
    if( el_pt > 30.0 )  return 193.0;
    if( el_pt > 20.0 )  return 108.0;
    if( el_pt > 10.0 )  return 63.0;
    return 0.0;
  }
  std::cout << "Error: eta > 2.5 value found" << endl; 
  return 0.0;
}


double getSingleEleDenom(double el_pt, double el_eta) {

  el_eta = fabs(el_eta);

  if( el_eta < 0.5 ){
    if( el_pt > 90.0 )	return 30492.0;
    if( el_pt > 80.0 )	return 31046.0;
    if( el_pt > 70.0 )	return 30888.0;
    if( el_pt > 60.0 )	return 30185.0;
    if( el_pt > 50.0 )	return 30199.0;
    if( el_pt > 40.0 )	return 30718.0;
    if( el_pt > 30.0 )	return 30480.0;
    if( el_pt > 20.0 )	return 29931.0;
    if( el_pt > 10.0 )	return 26461.0;
    return 0.0;
  }
  if( el_eta < 1.0 ){
    if( el_pt > 90.0 )  return 31770.0;
    if( el_pt > 80.0 )  return 32174.0;
    if( el_pt > 70.0 )  return 32055.0;
    if( el_pt > 60.0 )  return 32183.0;
    if( el_pt > 50.0 )  return 31528.0;
    if( el_pt > 40.0 )  return 31828.0;
    if( el_pt > 30.0 )  return 31533.0;
    if( el_pt > 20.0 )  return 29714.0;
    if( el_pt > 10.0 )  return 26948.0;
    return 0.0;
  }
  if( el_eta < 1.28 ){
    if( el_pt > 90.0 )  return 17922.0;
    if( el_pt > 80.0 )  return 18561.0;
    if( el_pt > 70.0 )  return 18797.0;
    if( el_pt > 60.0 )  return 18620.0;
    if( el_pt > 50.0 )  return 18176.0;
    if( el_pt > 40.0 )  return 18054.0;
    if( el_pt > 30.0 )  return 17420.0;
    if( el_pt > 20.0 )  return 16103.0;
    if( el_pt > 10.0 )  return 12802.0;
    return 0.0;
  }
  if( el_eta < 1.56 ){
    if( el_pt > 90.0 )  return 13708.0;
    if( el_pt > 80.0 )  return 14604.0;
    if( el_pt > 70.0 )  return 15142.0;
    if( el_pt > 60.0 )  return 14776.0;
    if( el_pt > 50.0 )  return 14817.0;
    if( el_pt > 40.0 )  return 14547.0;
    if( el_pt > 30.0 )  return 13621.0;
    if( el_pt > 20.0 )  return 12108.0;
    if( el_pt > 10.0 )  return 9209.0;
    return 0.0;
  }
  if( el_eta < 1.84 ){
    if( el_pt > 90.0 )  return 17060.0;
    if( el_pt > 80.0 )  return 17710.0;
    if( el_pt > 70.0 )  return 17726.0;
    if( el_pt > 60.0 )  return 17935.0;
    if( el_pt > 50.0 )  return 17803.0;
    if( el_pt > 40.0 )  return 17937.0;
    if( el_pt > 30.0 )  return 16760.0;
    if( el_pt > 20.0 )  return 14729.0;
    if( el_pt > 10.0 )  return 8994.0;
    return 0.0;
  }
  if( el_eta < 2.12 ){
    if( el_pt > 90.0 )  return 16625.0;
    if( el_pt > 80.0 )  return 17345.0;
    if( el_pt > 70.0 )  return 17555.0;
    if( el_pt > 60.0 )  return 17726.0;
    if( el_pt > 50.0 )  return 17809.0;
    if( el_pt > 40.0 )  return 17785.0;
    if( el_pt > 30.0 )  return 17459.0;
    if( el_pt > 20.0 )  return 16361.0;
    if( el_pt > 10.0 )  return 12200.0;
    return 0.0;
  }
  if( el_eta <= 2.5 ){
    if( el_pt > 90.0 )  return 18820.0;
    if( el_pt > 80.0 )  return 19841.0;
    if( el_pt > 70.0 )  return 20168.0;
    if( el_pt > 60.0 )  return 19906.0;
    if( el_pt > 50.0 )  return 20512.0;
    if( el_pt > 40.0 )  return 20889.0;
    if( el_pt > 30.0 )  return 20635.0;
    if( el_pt > 20.0 )  return 20646.0;
    if( el_pt > 10.0 )  return 18432.0;
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

