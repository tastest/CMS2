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
    if( el_pt > 90.0 )	return 28.0;
    if( el_pt > 80.0 )	return 12.0;
    if( el_pt > 70.0 )	return 11.0;
    if( el_pt > 60.0 )	return 17.0;
    if( el_pt > 50.0 )	return 12.0;
    if( el_pt > 40.0 )	return 12.0;
    if( el_pt > 30.0 )	return 8.0;
    if( el_pt > 20.0 )	return  7.0;
    if( el_pt > 10.0 )	return  5.0;
    return 0.0;
  }
  if( el_eta < 1.0 ){
    if( el_pt > 90.0 )  return 34.0;
    if( el_pt > 80.0 )  return 34.0;
    if( el_pt > 70.0 )  return 34.0;
    if( el_pt > 60.0 )  return 25.0;
    if( el_pt > 50.0 )  return 27.0;
    if( el_pt > 40.0 )  return 19.0;
    if( el_pt > 30.0 )  return 22.0;
    if( el_pt > 20.0 )  return 10.0;
    if( el_pt > 10.0 )  return  9.0;
    return 0.0;
  }
  if( el_eta < 1.28 ){
    if( el_pt > 90.0 )  return 61.0;
    if( el_pt > 80.0 )  return 37.0;
    if( el_pt > 70.0 )  return 55.0;
    if( el_pt > 60.0 )  return 39.0;
    if( el_pt > 50.0 )  return 34.0;
    if( el_pt > 40.0 )  return 29.0;
    if( el_pt > 30.0 )  return 16.0;
    if( el_pt > 20.0 )  return 18.0;
    if( el_pt > 10.0 )  return  5.0;
    return 0.0;
  }
  if( el_eta < 1.56 ){
    if( el_pt > 90.0 )  return 121.0;
    if( el_pt > 80.0 )  return 111.0;
    if( el_pt > 70.0 )  return 103.0;
    if( el_pt > 60.0 )  return 91.0;
    if( el_pt > 50.0 )  return 76.0;
    if( el_pt > 40.0 )  return 67.0;
    if( el_pt > 30.0 )  return 50.0;
    if( el_pt > 20.0 )  return 29.0;
    if( el_pt > 10.0 )  return 13.0;
    return 0.0;
  }
  if( el_eta < 1.84 ){
    if( el_pt > 90.0 )  return 309.0;
    if( el_pt > 80.0 )  return 298.0;
    if( el_pt > 70.0 )  return 278.0;
    if( el_pt > 60.0 )  return 250.0;
    if( el_pt > 50.0 )  return 189.0;
    if( el_pt > 40.0 )  return 174.0;
    if( el_pt > 30.0 )  return 113.0;
    if( el_pt > 20.0 )  return 76.0;
    if( el_pt > 10.0 )  return 32.0;
    return 0.0;
  }
  if( el_eta < 2.12 ){
    if( el_pt > 90.0 )  return 236.0;
    if( el_pt > 80.0 )  return 235.0;
    if( el_pt > 70.0 )  return 209.0;
    if( el_pt > 60.0 )  return 196.0;
    if( el_pt > 50.0 )  return 158.0;
    if( el_pt > 40.0 )  return 164.0;
    if( el_pt > 30.0 )  return 122.0;
    if( el_pt > 20.0 )  return 79.0;
    if( el_pt > 10.0 )  return 45.0;
    return 0.0;
  }
  if( el_eta <= 2.5 ){
    if( el_pt > 90.0 )  return 335.0;
    if( el_pt > 80.0 )  return 317.0;
    if( el_pt > 70.0 )  return 309.0;
    if( el_pt > 60.0 )  return 260.0;
    if( el_pt > 50.0 )  return 280.0;
    if( el_pt > 40.0 )  return 258.0;
    if( el_pt > 30.0 )  return 219.0;
    if( el_pt > 20.0 )  return 157.0;
    if( el_pt > 10.0 )  return 113.0;
    return 0.0;
  }
  std::cout << "Error: eta > 2.5 value found" << endl; 
  return 0.0;
}


double getSingleEleDenom(double el_pt, double el_eta) {

  el_eta = fabs(el_eta);

  if( el_eta < 0.5 ){
    if( el_pt > 90.0 )	return 31489.0;
    if( el_pt > 80.0 )	return 31338.0;
    if( el_pt > 70.0 )	return 30862.0;
    if( el_pt > 60.0 )	return 30196.0;
    if( el_pt > 50.0 )	return 30340.0;
    if( el_pt > 40.0 )	return 30928.0;
    if( el_pt > 30.0 )	return 30644.0;
    if( el_pt > 20.0 )	return 30311.0;
    if( el_pt > 10.0 )	return 27229.0;
    return 0.0;
  }
  if( el_eta < 1.0 ){
    if( el_pt > 90.0 )  return 32646.0;
    if( el_pt > 80.0 )  return 32375.0;
    if( el_pt > 70.0 )  return 32108.0;
    if( el_pt > 60.0 )  return 32410.0;
    if( el_pt > 50.0 )  return 31849.0;
    if( el_pt > 40.0 )  return 32202.0;
    if( el_pt > 30.0 )  return 32030.0;
    if( el_pt > 20.0 )  return 30543.0;
    if( el_pt > 10.0 )  return 28259.0;
    return 0.0;
  }
  if( el_eta < 1.28 ){
    if( el_pt > 90.0 )  return 19137.0;
    if( el_pt > 80.0 )  return 18866.0;
    if( el_pt > 70.0 )  return 19075.0;
    if( el_pt > 60.0 )  return 18773.0;
    if( el_pt > 50.0 )  return 18566.0;
    if( el_pt > 40.0 )  return 18434.0;
    if( el_pt > 30.0 )  return 17999.0;
    if( el_pt > 20.0 )  return 17082.0;
    if( el_pt > 10.0 )  return 14523.0;
    return 0.0;
  }
  if( el_eta < 1.56 ){
    if( el_pt > 90.0 )  return 15522.0;
    if( el_pt > 80.0 )  return 14949.0;
    if( el_pt > 70.0 )  return 15503.0;
    if( el_pt > 60.0 )  return 15170.0;
    if( el_pt > 50.0 )  return 15152.0;
    if( el_pt > 40.0 )  return 14957.0;
    if( el_pt > 30.0 )  return 14488.0;
    if( el_pt > 20.0 )  return 13212.0;
    if( el_pt > 10.0 )  return 11309.0;
    return 0.0;
  }
  if( el_eta < 1.84 ){
    if( el_pt > 90.0 )  return 17898.0;
    if( el_pt > 80.0 )  return 17889.0;
    if( el_pt > 70.0 )  return 17892.0;
    if( el_pt > 60.0 )  return 18417.0;
    if( el_pt > 50.0 )  return 18063.0;
    if( el_pt > 40.0 )  return 18635.0;
    if( el_pt > 30.0 )  return 17798.0;
    if( el_pt > 20.0 )  return 16213.0;
    if( el_pt > 10.0 )  return 10044.0;
    return 0.0;
  }
  if( el_eta < 2.12 ){
    if( el_pt > 90.0 )  return 17604.0;
    if( el_pt > 80.0 )  return 17472.0;
    if( el_pt > 70.0 )  return 17612.0;
    if( el_pt > 60.0 )  return 17919.0;
    if( el_pt > 50.0 )  return 18161.0;
    if( el_pt > 40.0 )  return 18108.0;
    if( el_pt > 30.0 )  return 18085.0;
    if( el_pt > 20.0 )  return 17223.0;
    if( el_pt > 10.0 )  return 13206.0;
    return 0.0;
  }
  if( el_eta <= 2.5 ){
    if( el_pt > 90.0 )  return 20035.0;
    if( el_pt > 80.0 )  return 19722.0;
    if( el_pt > 70.0 )  return 20278.0;
    if( el_pt > 60.0 )  return 19938.0;
    if( el_pt > 50.0 )  return 20545.0;
    if( el_pt > 40.0 )  return 21077.0;
    if( el_pt > 30.0 )  return 20963.0;
    if( el_pt > 20.0 )  return 21234.0;
    if( el_pt > 10.0 )  return 19485.0;
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

