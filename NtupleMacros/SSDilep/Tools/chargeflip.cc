#include "chargeflip.h"
#include <iostream>
#include <stdio.h>
#include <math.h>



using namespace std;

/////////////////////////////////////////////////////////////////////
double getZtoEEMCNum(double el_pt, double el_eta) {

  el_eta = fabs(el_eta);

  if(el_eta < 1.0 && el_pt < 30.0) 
    return 3.0;
   

  if(el_eta < 0.5) {
    if(el_pt >=30.0 && el_pt < 40.0 ) 
      return 4.0;
    
    if(el_pt >= 40.0 && el_pt < 50.0 )
      return 8.0;
  }


  if(el_eta >= 0.5 && el_eta < 1.0) {
    if(el_pt >= 30.0 && el_pt < 40.0)
      return 4.0;

    if(el_pt >= 40.0 && el_pt < 50.0 )
      return 18.0;
  }

  if(el_eta < 1.0 && el_pt >= 50.0)
    return 10.0;
  
  
  if(el_eta >=1.0 && el_eta < 1.56) {
    
    if(el_pt < 30.0)
      return 7.0;
  }
   
  if(el_eta >=1.0 && el_eta < 1.28) {
    if(el_pt >=30.0 &&  el_pt < 40.0)
      return 10.0;

    if(el_pt >= 40.0 && el_pt < 50.0)
      return 11.0;
    
    if(el_pt >=50.0)
      return 6.0;

  }

  if(el_eta >= 1.28 && el_eta < 1.56) {
    
    
    if(el_pt >= 30.0 && el_pt < 40.0)
      return 18.;

    if(el_pt >= 40 && el_pt < 50.)
      return 24.;

    if(el_pt >= 50)
      return 6.;

  }
      
  if(el_eta >= 1.56 && el_eta < 1.84) {

    if(el_pt < 30.0)
      return 5.0;
    
    if(el_pt >= 30.0 && el_pt < 40.0)
      return 34.;

    if(el_pt >= 40 && el_pt < 50.)
      return 55.;

    if(el_pt >= 50)
      return 19.;

  }

  if(el_eta >= 1.84 && el_eta < 2.12) {

    if(el_pt < 30.0)
      return 11.;

    if(el_pt >= 30.0 && el_pt < 40.0)
      return 24.;

    if(el_pt >= 40.0 && el_pt < 50)
      return 29.;
    
    if(el_pt >= 50)
      return 4.;

  }


  if(el_eta >= 2.12) {
    
    if(el_pt < 30.0)
      return 11.;

    if(el_pt >= 30.0 && el_pt < 40.0)
      return 26.;

    if(el_pt >= 40.0 && el_pt < 50.0)
      return 44.;

    if(el_pt >= 50.0)
      return 5.;
  
  }

  cout << "Should never get here!" << __FILE__ << __LINE__ << endl; 
  return 0.0;

  
}


/////////////////////////////////////////////////////////////////////
double getZtoEEMCDenom(double el_pt, double el_eta) {

  el_eta = fabs(el_eta);

  if(el_eta < 1.0 && el_pt < 30.0) 
    return 16216. + getZtoEEMCNum(el_pt, el_eta);
   

  if(el_eta < 0.5) {
    if(el_pt >=30.0 && el_pt < 40.0 ) 
      return 19863.0 + getZtoEEMCNum(el_pt, el_eta);
    
    if(el_pt >= 40.0 && el_pt < 50.0 )
      return 23842.0 + getZtoEEMCNum(el_pt, el_eta);
  }

  if(el_eta >= 0.5 && el_eta < 1.0) {
    if(el_pt >= 30.0 && el_pt < 40.0)
      return 18119.0 + getZtoEEMCNum(el_pt, el_eta);
    
    if(el_pt >= 40.0 && el_pt < 50.0 )
      return 23277.0 + getZtoEEMCNum(el_pt, el_eta);
  }
  
  
  if(el_eta < 1.0 && el_pt >= 50.0)
    return 13580.0 + getZtoEEMCNum(el_pt, el_eta);

  
  if(el_eta >=1.0 && el_eta < 1.56) {
    if(el_pt < 30.0)
      return 7596.0 + getZtoEEMCNum(el_pt, el_eta);
  }
  

  if(el_eta >=1.0 && el_eta < 1.28) {
    
    if(el_pt >=30.0 &&  el_pt < 40.0)
      return 8769.0 + getZtoEEMCNum(el_pt, el_eta);

    if(el_pt >= 40.0 && el_pt < 50.0)
      return 11459.0 + getZtoEEMCNum(el_pt, el_eta);
    
    if(el_pt >=50)
      return 3243.0 + getZtoEEMCNum(el_pt, el_eta);

  }

  if(el_eta >= 1.28 && el_eta < 1.56) {
            
    if(el_pt >= 30.0 && el_pt < 40.0)
      return 6262. + getZtoEEMCNum(el_pt, el_eta);

    if(el_pt >= 40 && el_pt < 50.)
      return 8182. + getZtoEEMCNum(el_pt, el_eta);

    if(el_pt >= 50)
      return 2235. + getZtoEEMCNum(el_pt, el_eta);

  }
      
  if(el_eta > 1.56 && el_eta < 1.84) {

    if(el_pt < 30.0)
      return 3703. + getZtoEEMCNum(el_pt, el_eta);
    
    if(el_pt >= 30.0 && el_pt < 40.0)
      return 6193. + getZtoEEMCNum(el_pt, el_eta);

    if(el_pt >= 40 && el_pt < 50.)
      return 8476. + getZtoEEMCNum(el_pt, el_eta);

    if(el_pt >= 50)
      return 2456. + getZtoEEMCNum(el_pt, el_eta);

  }

  if(el_eta >= 1.84 && el_eta < 2.12) {

    if(el_pt < 30.0)
      return 4031. + getZtoEEMCNum(el_pt, el_eta);

    if(el_pt >= 30.0 && el_pt < 40.0)
      return 5816. + getZtoEEMCNum(el_pt, el_eta);

    if(el_pt >= 40.0 && el_pt < 50)
      return 7198. + getZtoEEMCNum(el_pt, el_eta);
    
    if(el_pt >= 50)
      return 2049. + getZtoEEMCNum(el_pt, el_eta);

  }


  if(el_eta >= 2.12) {
    
    if(el_pt < 30.0)
      return 3644. + getZtoEEMCNum(el_pt, el_eta);

    if(el_pt >= 30.0 && el_pt < 40.0)
      return 4667. + getZtoEEMCNum(el_pt, el_eta);

    if(el_pt >= 40. && el_pt < 50.0)
      return 5010. + getZtoEEMCNum(el_pt, el_eta);

    if(el_pt >= 50.0)
      return 1304. + getZtoEEMCNum(el_pt, el_eta);
  
  }

  cout << "Should never get here! " << __FILE__ << __LINE__ << endl;
  return 0.0;

}

/////////////////////////////////////////////////////////////////////

double getZtoEEMCFlipRate(double el_pt, double el_eta) { 
  
  return getZtoEEMCNum(el_pt, fabs(el_eta))/getZtoEEMCDenom(el_pt, fabs(el_eta));

}

/////////////////////////////////////////////////////////////////////
double getZtoEEMCFlipRateError(double el_pt, double el_eta) {


  //the binomial error
  double num   = getZtoEEMCNum(el_pt,   fabs(el_eta));
  double denom = getZtoEEMCDenom(el_pt, fabs(el_eta));
  double p = num/denom;
  
  return sqrt(p*(1-p)/denom);
}

