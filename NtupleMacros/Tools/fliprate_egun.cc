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

double getSingleEleNum(double el_pt, double el_eta) {

  el_eta = fabs(el_eta);

  if( el_eta < 0.5 ){
    if( el_pt > 90.0 )	return 26.0;
    if( el_pt > 80.0 )	return 24.0;
    if( el_pt > 70.0 )	return 29.0;
    if( el_pt > 60.0 )	return 22.0;
    if( el_pt > 50.0 )	return  9.0;
    if( el_pt > 40.0 )	return 25.0;
    if( el_pt > 30.0 )	return 12.0;
    if( el_pt > 20.0 )	return  7.0;
    if( el_pt > 10.0 )	return  5.0;
    return 0.0;
  }
  if( el_eta < 1.0 ){
    if( el_pt > 90.0 )  return 31.0;
    if( el_pt > 80.0 )  return 51.0;
    if( el_pt > 70.0 )  return 43.0;
    if( el_pt > 60.0 )  return 39.0;
    if( el_pt > 50.0 )  return 37.0;
    if( el_pt > 40.0 )  return 24.0;
    if( el_pt > 30.0 )  return 20.0;
    if( el_pt > 20.0 )  return  7.0;
    if( el_pt > 10.0 )  return  3.0;
    return 0.0;
  }
  if( el_eta < 1.28 ){
    if( el_pt > 90.0 )  return 47.0;
    if( el_pt > 80.0 )  return 56.0;
    if( el_pt > 70.0 )  return 40.0;
    if( el_pt > 60.0 )  return 50.0;
    if( el_pt > 50.0 )  return 39.0;
    if( el_pt > 40.0 )  return 37.0;
    if( el_pt > 30.0 )  return 22.0;
    if( el_pt > 20.0 )  return 10.0;
    if( el_pt > 10.0 )  return  8.0;
    return 0.0;
  }
  if( el_eta < 1.56 ){
    if( el_pt > 90.0 )  return 96.0;
    if( el_pt > 80.0 )  return 103.0;
    if( el_pt > 70.0 )  return 98.0;
    if( el_pt > 60.0 )  return 97.0;
    if( el_pt > 50.0 )  return 88.0;
    if( el_pt > 40.0 )  return 50.0;
    if( el_pt > 30.0 )  return 59.0;
    if( el_pt > 20.0 )  return 22.0;
    if( el_pt > 10.0 )  return  8.0;
    return 0.0;
  }
  if( el_eta < 1.84 ){
    if( el_pt > 90.0 )  return 240.0;
    if( el_pt > 80.0 )  return 256.0;
    if( el_pt > 70.0 )  return 241.0;
    if( el_pt > 60.0 )  return 210.0;
    if( el_pt > 50.0 )  return 180.0;
    if( el_pt > 40.0 )  return 146.0;
    if( el_pt > 30.0 )  return 94.0;
    if( el_pt > 20.0 )  return 45.0;
    if( el_pt > 10.0 )  return 16.0;
    return 0.0;
  }
  if( el_eta < 2.12 ){
    if( el_pt > 90.0 )  return 198.0;
    if( el_pt > 80.0 )  return 232.0;
    if( el_pt > 70.0 )  return 201.0;
    if( el_pt > 60.0 )  return 206.0;
    if( el_pt > 50.0 )  return 162.0;
    if( el_pt > 40.0 )  return 113.0;
    if( el_pt > 30.0 )  return 77.0;
    if( el_pt > 20.0 )  return 43.0;
    if( el_pt > 10.0 )  return 10.0;
    return 0.0;
  }
  if( el_eta <= 2.4 ){
    if( el_pt > 90.0 )  return 180.0;
    if( el_pt > 80.0 )  return 217.0;
    if( el_pt > 70.0 )  return 220.0;
    if( el_pt > 60.0 )  return 188.0;
    if( el_pt > 50.0 )  return 172.0;
    if( el_pt > 40.0 )  return 174.0;
    if( el_pt > 30.0 )  return 137.0;
    if( el_pt > 20.0 )  return 100.0;
    if( el_pt > 10.0 )  return  39.0;
    return 0.0;
  }
  std::cout << "Error: eta > 2.4 value found" << endl; 
  return 0.0;
}


double getSingleEleDenom(double el_pt, double el_eta) {

  el_eta = fabs(el_eta);

  if( el_eta < 0.5 ){
    if( el_pt > 90.0 )	return 38019.0;
    if( el_pt > 80.0 )	return 39008.0;
    if( el_pt > 70.0 )	return 39668.0;
    if( el_pt > 60.0 )	return 39292.0;
    if( el_pt > 50.0 )	return 39538.0;
    if( el_pt > 40.0 )	return 39462.0;
    if( el_pt > 30.0 )	return 39498.0;
    if( el_pt > 20.0 )	return 38681.0;
    if( el_pt > 10.0 )	return 36077.0;
    return 0.0;
  }
  if( el_eta < 1.0 ){
    if( el_pt > 90.0 )  return 38665.0;
    if( el_pt > 80.0 )  return 39954.0;
    if( el_pt > 70.0 )  return 39946.0;
    if( el_pt > 60.0 )  return 39874.0;
    if( el_pt > 50.0 )  return 40589.0;
    if( el_pt > 40.0 )  return 40002.0;
    if( el_pt > 30.0 )  return 39947.0;
    if( el_pt > 20.0 )  return 39383.0;
    if( el_pt > 10.0 )  return 36947.0;
    return 0.0;
  }
  if( el_eta < 1.28 ){
    if( el_pt > 90.0 )  return 20181.0;
    if( el_pt > 80.0 )  return 21839.0;
    if( el_pt > 70.0 )  return 21684.0;
    if( el_pt > 60.0 )  return 21629.0;
    if( el_pt > 50.0 )  return 21948.0;
    if( el_pt > 40.0 )  return 21431.0;
    if( el_pt > 30.0 )  return 21777.0;
    if( el_pt > 20.0 )  return 21046.0;
    if( el_pt > 10.0 )  return 18560.0;
    return 0.0;
  }
  if( el_eta < 1.56 ){
    if( el_pt > 90.0 )  return 15161.0;
    if( el_pt > 80.0 )  return 17173.0;
    if( el_pt > 70.0 )  return 17489.0;
    if( el_pt > 60.0 )  return 17561.0;
    if( el_pt > 50.0 )  return 17794.0;
    if( el_pt > 40.0 )  return 17651.0;
    if( el_pt > 30.0 )  return 17281.0;
    if( el_pt > 20.0 )  return 16145.0;
    if( el_pt > 10.0 )  return 13391.0;
    return 0.0;
  }
  if( el_eta < 1.84 ){
    if( el_pt > 90.0 )  return 18861.0;
    if( el_pt > 80.0 )  return 20511.0;
    if( el_pt > 70.0 )  return 20640.0;
    if( el_pt > 60.0 )  return 20627.0;
    if( el_pt > 50.0 )  return 20522.0;
    if( el_pt > 40.0 )  return 20085.0;
    if( el_pt > 30.0 )  return 19473.0;
    if( el_pt > 20.0 )  return 18245.0;
    if( el_pt > 10.0 )  return 13404.0;
    return 0.0;
  }
  if( el_eta < 2.12 ){
    if( el_pt > 90.0 )  return 18737.0;
    if( el_pt > 80.0 )  return 19583.0;
    if( el_pt > 70.0 )  return 19915.0;
    if( el_pt > 60.0 )  return 19982.0;
    if( el_pt > 50.0 )  return 19891.0;
    if( el_pt > 40.0 )  return 19915.0;
    if( el_pt > 30.0 )  return 19552.0;
    if( el_pt > 20.0 )  return 19276.0;
    if( el_pt > 10.0 )  return 17462.0;
    return 0.0;
  }
  if( el_eta <= 2.4 ){
    if( el_pt > 90.0 )  return 15319.0;
    if( el_pt > 80.0 )  return 16760.0;
    if( el_pt > 70.0 )  return 17310.0;
    if( el_pt > 60.0 )  return 16995.0;
    if( el_pt > 50.0 )  return 17632.0;
    if( el_pt > 40.0 )  return 17522.0;
    if( el_pt > 30.0 )  return 17377.0;
    if( el_pt > 20.0 )  return 17227.0;
    if( el_pt > 10.0 )  return 16620.0;
    return 0.0;
  }
  std::cout << "Error: eta > 2.4 value found" << endl; 
  return 0.0;
}


double getSingleEleFlipRate(double el_pt, double el_eta) { 
  if( el_pt < 10.0 || fabs(el_eta) > 2.4 ){
    std::cout << "Error in 'getSingleEleFlipRate': pt or eta value found out of range" << endl;  
    return 0.0;
  } 
  return getSingleEleNum(el_pt, fabs(el_eta))/getSingleEleDenom(el_pt, fabs(el_eta));
}

double getSingleEleFlipRateError(double el_pt, double el_eta) {
  //the binomial error
  if( el_pt < 10.0 || fabs(el_eta) > 2.4 ){
    std::cout << "Error in 'getSingleEleFlipRate': pt or eta value found out of range" << endl;
    return 0.0;
  }
  double num   = getSingleEleNum(el_pt,   fabs(el_eta));
  double denom = getSingleEleDenom(el_pt, fabs(el_eta));
  double p = num/denom;
  return sqrt(p*(1-p)/denom);
}

