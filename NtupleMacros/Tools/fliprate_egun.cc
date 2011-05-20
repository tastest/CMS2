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
    if( el_pt > 95 ) return 9;
    if( el_pt > 90 ) return 5;
    if( el_pt > 85 ) return 4;
    if( el_pt > 80 ) return 6;
    if( el_pt > 75 ) return 7;
    if( el_pt > 70 ) return 5;
    if( el_pt > 65 ) return 8;
    if( el_pt > 60 ) return 4;
    if( el_pt > 55 ) return 4;
    if( el_pt > 50 ) return 3;
    if( el_pt > 45 ) return 6;
    if( el_pt > 40 ) return 6;
    if( el_pt > 35 ) return 6;
    if( el_pt > 30 ) return 3;
    if( el_pt > 25 ) return 5;
    if( el_pt > 20 ) return 0;
    if( el_pt > 10 ) return 1;
    return 0.0;
  }
  if( el_eta < 1 ){ 
    if( el_pt > 95 ) return 17;
    if( el_pt > 90 ) return 25;
    if( el_pt > 85 ) return 16;
    if( el_pt > 80 ) return 9;
    if( el_pt > 75 ) return 18;
    if( el_pt > 70 ) return 10;
    if( el_pt > 65 ) return 17;
    if( el_pt > 60 ) return 15;
    if( el_pt > 55 ) return 13;
    if( el_pt > 50 ) return 21;
    if( el_pt > 45 ) return 14;
    if( el_pt > 40 ) return 7;
    if( el_pt > 35 ) return 8;
    if( el_pt > 30 ) return 9;
    if( el_pt > 25 ) return 12;
    if( el_pt > 20 ) return 4;
    if( el_pt > 10 ) return 2;
    return 0.0;
  }
  if( el_eta < 1.479 ){ 
    if( el_pt > 95 ) return 88;
    if( el_pt > 90 ) return 76;
    if( el_pt > 85 ) return 95;
    if( el_pt > 80 ) return 73;
    if( el_pt > 75 ) return 84;
    if( el_pt > 70 ) return 72;
    if( el_pt > 65 ) return 75;
    if( el_pt > 60 ) return 75;
    if( el_pt > 55 ) return 77;
    if( el_pt > 50 ) return 65;
    if( el_pt > 45 ) return 45;
    if( el_pt > 40 ) return 39;
    if( el_pt > 35 ) return 48;
    if( el_pt > 30 ) return 33;
    if( el_pt > 25 ) return 25;
    if( el_pt > 20 ) return 19;
    if( el_pt > 10 ) return 8;
    return 0.0;
  }
  if( el_eta < 1.8 ){ 
    if( el_pt > 95 ) return 144;
    if( el_pt > 90 ) return 138;
    if( el_pt > 85 ) return 129;
    if( el_pt > 80 ) return 116;
    if( el_pt > 75 ) return 122;
    if( el_pt > 70 ) return 103;
    if( el_pt > 65 ) return 111;
    if( el_pt > 60 ) return 108;
    if( el_pt > 55 ) return 90;
    if( el_pt > 50 ) return 93;
    if( el_pt > 45 ) return 91;
    if( el_pt > 40 ) return 69;
    if( el_pt > 35 ) return 46;
    if( el_pt > 30 ) return 44;
    if( el_pt > 25 ) return 30;
    if( el_pt > 20 ) return 16;
    if( el_pt > 10 ) return 9;
    return 0.0;
  }
  if( el_eta < 2 ){ 
    if( el_pt > 95 ) return 70;
    if( el_pt > 90 ) return 60;
    if( el_pt > 85 ) return 70;
    if( el_pt > 80 ) return 60;
    if( el_pt > 75 ) return 56;
    if( el_pt > 70 ) return 74;
    if( el_pt > 65 ) return 58;
    if( el_pt > 60 ) return 54;
    if( el_pt > 55 ) return 59;
    if( el_pt > 50 ) return 60;
    if( el_pt > 45 ) return 51;
    if( el_pt > 40 ) return 43;
    if( el_pt > 35 ) return 25;
    if( el_pt > 30 ) return 28;
    if( el_pt > 25 ) return 17;
    if( el_pt > 20 ) return 10;
    if( el_pt > 10 ) return 6;
    return 0.0;
  }
  if( el_eta < 2.1 ){ 
    if( el_pt > 95 ) return 40;
    if( el_pt > 90 ) return 33;
    if( el_pt > 85 ) return 33;
    if( el_pt > 80 ) return 29;
    if( el_pt > 75 ) return 28;
    if( el_pt > 70 ) return 23;
    if( el_pt > 65 ) return 27;
    if( el_pt > 60 ) return 21;
    if( el_pt > 55 ) return 27;
    if( el_pt > 50 ) return 15;
    if( el_pt > 45 ) return 18;
    if( el_pt > 40 ) return 17;
    if( el_pt > 35 ) return 11;
    if( el_pt > 30 ) return 12;
    if( el_pt > 25 ) return 8;
    if( el_pt > 20 ) return 8;
    if( el_pt > 10 ) return 6;
    return 0.0;
  }
  if( el_eta < 2.2 ){ 
    if( el_pt > 95 ) return 35;
    if( el_pt > 90 ) return 35;
    if( el_pt > 85 ) return 27;
    if( el_pt > 80 ) return 29;
    if( el_pt > 75 ) return 36;
    if( el_pt > 70 ) return 20;
    if( el_pt > 65 ) return 24;
    if( el_pt > 60 ) return 21;
    if( el_pt > 55 ) return 19;
    if( el_pt > 50 ) return 32;
    if( el_pt > 45 ) return 21;
    if( el_pt > 40 ) return 19;
    if( el_pt > 35 ) return 12;
    if( el_pt > 30 ) return 10;
    if( el_pt > 25 ) return 13;
    if( el_pt > 20 ) return 5;
    if( el_pt > 10 ) return 9;
    return 0.0;
  }
  if( el_eta <= 2.4 ){ 
    if( el_pt > 95 ) return 64;
    if( el_pt > 90 ) return 84;
    if( el_pt > 85 ) return 68;
    if( el_pt > 80 ) return 64;
    if( el_pt > 75 ) return 56;
    if( el_pt > 70 ) return 57;
    if( el_pt > 65 ) return 54;
    if( el_pt > 60 ) return 47;
    if( el_pt > 55 ) return 38;
    if( el_pt > 50 ) return 46;
    if( el_pt > 45 ) return 39;
    if( el_pt > 40 ) return 29;
    if( el_pt > 35 ) return 23;
    if( el_pt > 30 ) return 25;
    if( el_pt > 25 ) return 14;
    if( el_pt > 20 ) return 20;
    if( el_pt > 10 ) return 14;
    return 0.0;
  }
  std::cout << "Error: eta > 2.4 value found" << endl;
  return 0.0;
}


double getSingleEleDenom(double el_pt, double el_eta) {

  el_eta = std::fabs(el_eta);

  if( el_eta < 0.5 ){ 
    if( el_pt > 95 ) return 80174;
    if( el_pt > 90 ) return 79910;
    if( el_pt > 85 ) return 79820;
    if( el_pt > 80 ) return 79722;
    if( el_pt > 75 ) return 79890;
    if( el_pt > 70 ) return 79653;
    if( el_pt > 65 ) return 79892;
    if( el_pt > 60 ) return 80077;
    if( el_pt > 55 ) return 80011;
    if( el_pt > 50 ) return 79597;
    if( el_pt > 45 ) return 80053;
    if( el_pt > 40 ) return 79263;
    if( el_pt > 35 ) return 79239;
    if( el_pt > 30 ) return 78273;
    if( el_pt > 25 ) return 77894;
    if( el_pt > 20 ) return 75100;
    if( el_pt > 10 ) return 114108;
    return 0.0;
  }
  if( el_eta < 1 ){ 
    if( el_pt > 95 ) return 79315;
    if( el_pt > 90 ) return 79561;
    if( el_pt > 85 ) return 79406;
    if( el_pt > 80 ) return 79842;
    if( el_pt > 75 ) return 79654;
    if( el_pt > 70 ) return 79671;
    if( el_pt > 65 ) return 79640;
    if( el_pt > 60 ) return 79599;
    if( el_pt > 55 ) return 79540;
    if( el_pt > 50 ) return 79343;
    if( el_pt > 45 ) return 79944;
    if( el_pt > 40 ) return 79623;
    if( el_pt > 35 ) return 78557;
    if( el_pt > 30 ) return 78303;
    if( el_pt > 25 ) return 77378;
    if( el_pt > 20 ) return 75091;
    if( el_pt > 10 ) return 122322;
    return 0.0;
  }
  if( el_eta < 1.479 ){ 
    if( el_pt > 95 ) return 71523;
    if( el_pt > 90 ) return 70668;
    if( el_pt > 85 ) return 71662;
    if( el_pt > 80 ) return 71434;
    if( el_pt > 75 ) return 71467;
    if( el_pt > 70 ) return 71204;
    if( el_pt > 65 ) return 70951;
    if( el_pt > 60 ) return 71044;
    if( el_pt > 55 ) return 70664;
    if( el_pt > 50 ) return 70559;
    if( el_pt > 45 ) return 70218;
    if( el_pt > 40 ) return 69528;
    if( el_pt > 35 ) return 68507;
    if( el_pt > 30 ) return 67672;
    if( el_pt > 25 ) return 65776;
    if( el_pt > 20 ) return 61672;
    if( el_pt > 10 ) return 93412;
    return 0.0;
  }
  if( el_eta < 1.8 ){ 
    if( el_pt > 95 ) return 37802;
    if( el_pt > 90 ) return 38002;
    if( el_pt > 85 ) return 37636;
    if( el_pt > 80 ) return 38035;
    if( el_pt > 75 ) return 38432;
    if( el_pt > 70 ) return 37862;
    if( el_pt > 65 ) return 37989;
    if( el_pt > 60 ) return 37686;
    if( el_pt > 55 ) return 37738;
    if( el_pt > 50 ) return 37321;
    if( el_pt > 45 ) return 37092;
    if( el_pt > 40 ) return 36538;
    if( el_pt > 35 ) return 35582;
    if( el_pt > 30 ) return 34886;
    if( el_pt > 25 ) return 32819;
    if( el_pt > 20 ) return 30064;
    if( el_pt > 10 ) return 36945;
    return 0.0;
  }
  if( el_eta < 2 ){ 
    if( el_pt > 95 ) return 27683;
    if( el_pt > 90 ) return 27689;
    if( el_pt > 85 ) return 27736;
    if( el_pt > 80 ) return 27917;
    if( el_pt > 75 ) return 27924;
    if( el_pt > 70 ) return 27763;
    if( el_pt > 65 ) return 28018;
    if( el_pt > 60 ) return 27977;
    if( el_pt > 55 ) return 28120;
    if( el_pt > 50 ) return 27549;
    if( el_pt > 45 ) return 27622;
    if( el_pt > 40 ) return 27251;
    if( el_pt > 35 ) return 26803;
    if( el_pt > 30 ) return 26001;
    if( el_pt > 25 ) return 25559;
    if( el_pt > 20 ) return 23249;
    if( el_pt > 10 ) return 29668;
    return 0.0;
  }
  if( el_eta < 2.1 ){ 
    if( el_pt > 95 ) return 13372;
    if( el_pt > 90 ) return 13612;
    if( el_pt > 85 ) return 13518;
    if( el_pt > 80 ) return 13317;
    if( el_pt > 75 ) return 13707;
    if( el_pt > 70 ) return 13562;
    if( el_pt > 65 ) return 13637;
    if( el_pt > 60 ) return 13714;
    if( el_pt > 55 ) return 13635;
    if( el_pt > 50 ) return 13628;
    if( el_pt > 45 ) return 13481;
    if( el_pt > 40 ) return 13384;
    if( el_pt > 35 ) return 13437;
    if( el_pt > 30 ) return 13341;
    if( el_pt > 25 ) return 13083;
    if( el_pt > 20 ) return 12750;
    if( el_pt > 10 ) return 17362;
    return 0.0;
  }
  if( el_eta < 2.2 ){ 
    if( el_pt > 95 ) return 12545;
    if( el_pt > 90 ) return 12730;
    if( el_pt > 85 ) return 12914;
    if( el_pt > 80 ) return 12765;
    if( el_pt > 75 ) return 12899;
    if( el_pt > 70 ) return 12897;
    if( el_pt > 65 ) return 12758;
    if( el_pt > 60 ) return 12824;
    if( el_pt > 55 ) return 13093;
    if( el_pt > 50 ) return 12773;
    if( el_pt > 45 ) return 12724;
    if( el_pt > 40 ) return 12958;
    if( el_pt > 35 ) return 12976;
    if( el_pt > 30 ) return 12877;
    if( el_pt > 25 ) return 12626;
    if( el_pt > 20 ) return 12046;
    if( el_pt > 10 ) return 17226;
    return 0.0;
  }
  if( el_eta <= 2.4 ){ 
    if( el_pt > 95 ) return 23896;
    if( el_pt > 90 ) return 24146;
    if( el_pt > 85 ) return 24045;
    if( el_pt > 80 ) return 24087;
    if( el_pt > 75 ) return 24197;
    if( el_pt > 70 ) return 24369;
    if( el_pt > 65 ) return 24486;
    if( el_pt > 60 ) return 24379;
    if( el_pt > 55 ) return 24598;
    if( el_pt > 50 ) return 24470;
    if( el_pt > 45 ) return 24531;
    if( el_pt > 40 ) return 24214;
    if( el_pt > 35 ) return 24295;
    if( el_pt > 30 ) return 24178;
    if( el_pt > 25 ) return 23611;
    if( el_pt > 20 ) return 23022;
    if( el_pt > 10 ) return 31178;
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
