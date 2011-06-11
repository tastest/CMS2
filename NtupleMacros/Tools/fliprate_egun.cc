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
    if( el_pt > 95 ) return 13;
    if( el_pt > 90 ) return 6;
    if( el_pt > 85 ) return 7;
    if( el_pt > 80 ) return 13;
    if( el_pt > 75 ) return 9;
    if( el_pt > 70 ) return 6;
    if( el_pt > 65 ) return 11;
    if( el_pt > 60 ) return 11;
    if( el_pt > 55 ) return 3;
    if( el_pt > 50 ) return 8;
    if( el_pt > 45 ) return 10;
    if( el_pt > 40 ) return 9;
    if( el_pt > 35 ) return 11;
    if( el_pt > 30 ) return 12;
    if( el_pt > 25 ) return 3;
    if( el_pt > 20 ) return 3;
    if( el_pt > 10 ) return 2;
    return 0.0;
  }
  if( el_eta < 1 ){ 
    if( el_pt > 95 ) return 20;
    if( el_pt > 90 ) return 20;
    if( el_pt > 85 ) return 26;
    if( el_pt > 80 ) return 26;
    if( el_pt > 75 ) return 20;
    if( el_pt > 70 ) return 17;
    if( el_pt > 65 ) return 20;
    if( el_pt > 60 ) return 25;
    if( el_pt > 55 ) return 16;
    if( el_pt > 50 ) return 24;
    if( el_pt > 45 ) return 19;
    if( el_pt > 40 ) return 23;
    if( el_pt > 35 ) return 13;
    if( el_pt > 30 ) return 10;
    if( el_pt > 25 ) return 6;
    if( el_pt > 20 ) return 11;
    if( el_pt > 10 ) return 7;
    return 0.0;
  }
  if( el_eta < 1.479 ){ 
    if( el_pt > 95 ) return 177;
    if( el_pt > 90 ) return 155;
    if( el_pt > 85 ) return 147;
    if( el_pt > 80 ) return 163;
    if( el_pt > 75 ) return 145;
    if( el_pt > 70 ) return 159;
    if( el_pt > 65 ) return 152;
    if( el_pt > 60 ) return 142;
    if( el_pt > 55 ) return 115;
    if( el_pt > 50 ) return 121;
    if( el_pt > 45 ) return 112;
    if( el_pt > 40 ) return 99;
    if( el_pt > 35 ) return 90;
    if( el_pt > 30 ) return 64;
    if( el_pt > 25 ) return 54;
    if( el_pt > 20 ) return 41;
    if( el_pt > 10 ) return 13;
    return 0.0;
  }
  if( el_eta < 1.8 ){ 
    if( el_pt > 95 ) return 290;
    if( el_pt > 90 ) return 302;
    if( el_pt > 85 ) return 260;
    if( el_pt > 80 ) return 274;
    if( el_pt > 75 ) return 280;
    if( el_pt > 70 ) return 251;
    if( el_pt > 65 ) return 193;
    if( el_pt > 60 ) return 223;
    if( el_pt > 55 ) return 209;
    if( el_pt > 50 ) return 158;
    if( el_pt > 45 ) return 163;
    if( el_pt > 40 ) return 144;
    if( el_pt > 35 ) return 138;
    if( el_pt > 30 ) return 98;
    if( el_pt > 25 ) return 69;
    if( el_pt > 20 ) return 46;
    if( el_pt > 10 ) return 11;
    return 0.0;
  }
  if( el_eta < 2 ){ 
    if( el_pt > 95 ) return 166;
    if( el_pt > 90 ) return 153;
    if( el_pt > 85 ) return 141;
    if( el_pt > 80 ) return 146;
    if( el_pt > 75 ) return 148;
    if( el_pt > 70 ) return 157;
    if( el_pt > 65 ) return 117;
    if( el_pt > 60 ) return 133;
    if( el_pt > 55 ) return 129;
    if( el_pt > 50 ) return 107;
    if( el_pt > 45 ) return 117;
    if( el_pt > 40 ) return 119;
    if( el_pt > 35 ) return 96;
    if( el_pt > 30 ) return 74;
    if( el_pt > 25 ) return 49;
    if( el_pt > 20 ) return 31;
    if( el_pt > 10 ) return 20;
    return 0.0;
  }
  if( el_eta < 2.1 ){ 
    if( el_pt > 95 ) return 70;
    if( el_pt > 90 ) return 72;
    if( el_pt > 85 ) return 69;
    if( el_pt > 80 ) return 66;
    if( el_pt > 75 ) return 73;
    if( el_pt > 70 ) return 66;
    if( el_pt > 65 ) return 60;
    if( el_pt > 60 ) return 50;
    if( el_pt > 55 ) return 58;
    if( el_pt > 50 ) return 58;
    if( el_pt > 45 ) return 64;
    if( el_pt > 40 ) return 46;
    if( el_pt > 35 ) return 33;
    if( el_pt > 30 ) return 30;
    if( el_pt > 25 ) return 25;
    if( el_pt > 20 ) return 26;
    if( el_pt > 10 ) return 4;
    return 0.0;
  }
  if( el_eta < 2.2 ){ 
    if( el_pt > 95 ) return 77;
    if( el_pt > 90 ) return 65;
    if( el_pt > 85 ) return 83;
    if( el_pt > 80 ) return 73;
    if( el_pt > 75 ) return 81;
    if( el_pt > 70 ) return 58;
    if( el_pt > 65 ) return 56;
    if( el_pt > 60 ) return 54;
    if( el_pt > 55 ) return 49;
    if( el_pt > 50 ) return 46;
    if( el_pt > 45 ) return 43;
    if( el_pt > 40 ) return 37;
    if( el_pt > 35 ) return 37;
    if( el_pt > 30 ) return 30;
    if( el_pt > 25 ) return 22;
    if( el_pt > 20 ) return 24;
    if( el_pt > 10 ) return 22;
    return 0.0;
  }
  if( el_eta <= 2.4 ){ 
    if( el_pt > 95 ) return 151;
    if( el_pt > 90 ) return 135;
    if( el_pt > 85 ) return 138;
    if( el_pt > 80 ) return 116;
    if( el_pt > 75 ) return 119;
    if( el_pt > 70 ) return 123;
    if( el_pt > 65 ) return 99;
    if( el_pt > 60 ) return 107;
    if( el_pt > 55 ) return 89;
    if( el_pt > 50 ) return 73;
    if( el_pt > 45 ) return 85;
    if( el_pt > 40 ) return 66;
    if( el_pt > 35 ) return 63;
    if( el_pt > 30 ) return 55;
    if( el_pt > 25 ) return 46;
    if( el_pt > 20 ) return 41;
    if( el_pt > 10 ) return 25;
    return 0.0;
  }
  std::cout << "Error: eta > 2.4 value found" << endl;
  return 0.0;
}


double getSingleEleDenom(double el_pt, double el_eta) {

  el_eta = std::fabs(el_eta);

  if( el_eta < 0.5 ){ 
    if( el_pt > 95 ) return 146039;
    if( el_pt > 90 ) return 145996;
    if( el_pt > 85 ) return 145655;
    if( el_pt > 80 ) return 145579;
    if( el_pt > 75 ) return 145586;
    if( el_pt > 70 ) return 146280;
    if( el_pt > 65 ) return 145877;
    if( el_pt > 60 ) return 145949;
    if( el_pt > 55 ) return 145744;
    if( el_pt > 50 ) return 145034;
    if( el_pt > 45 ) return 145704;
    if( el_pt > 40 ) return 145591;
    if( el_pt > 35 ) return 144881;
    if( el_pt > 30 ) return 143408;
    if( el_pt > 25 ) return 142485;
    if( el_pt > 20 ) return 138025;
    if( el_pt > 10 ) return 211781;
    return 0.0;
  }
  if( el_eta < 1 ){ 
    if( el_pt > 95 ) return 144904;
    if( el_pt > 90 ) return 144773;
    if( el_pt > 85 ) return 144301;
    if( el_pt > 80 ) return 144859;
    if( el_pt > 75 ) return 144819;
    if( el_pt > 70 ) return 144986;
    if( el_pt > 65 ) return 145269;
    if( el_pt > 60 ) return 144663;
    if( el_pt > 55 ) return 144649;
    if( el_pt > 50 ) return 144620;
    if( el_pt > 45 ) return 145117;
    if( el_pt > 40 ) return 144704;
    if( el_pt > 35 ) return 143959;
    if( el_pt > 30 ) return 143506;
    if( el_pt > 25 ) return 141976;
    if( el_pt > 20 ) return 137829;
    if( el_pt > 10 ) return 226752;
    return 0.0;
  }
  if( el_eta < 1.479 ){ 
    if( el_pt > 95 ) return 130698;
    if( el_pt > 90 ) return 130670;
    if( el_pt > 85 ) return 131265;
    if( el_pt > 80 ) return 130788;
    if( el_pt > 75 ) return 130542;
    if( el_pt > 70 ) return 130913;
    if( el_pt > 65 ) return 130703;
    if( el_pt > 60 ) return 131326;
    if( el_pt > 55 ) return 130205;
    if( el_pt > 50 ) return 130500;
    if( el_pt > 45 ) return 128679;
    if( el_pt > 40 ) return 127828;
    if( el_pt > 35 ) return 126136;
    if( el_pt > 30 ) return 124081;
    if( el_pt > 25 ) return 120124;
    if( el_pt > 20 ) return 113099;
    if( el_pt > 10 ) return 171132;
    return 0.0;
  }
  if( el_eta < 1.8 ){ 
    if( el_pt > 95 ) return 67948;
    if( el_pt > 90 ) return 67662;
    if( el_pt > 85 ) return 68025;
    if( el_pt > 80 ) return 67695;
    if( el_pt > 75 ) return 67964;
    if( el_pt > 70 ) return 67763;
    if( el_pt > 65 ) return 68007;
    if( el_pt > 60 ) return 67484;
    if( el_pt > 55 ) return 67503;
    if( el_pt > 50 ) return 67047;
    if( el_pt > 45 ) return 66825;
    if( el_pt > 40 ) return 66613;
    if( el_pt > 35 ) return 64770;
    if( el_pt > 30 ) return 63052;
    if( el_pt > 25 ) return 60261;
    if( el_pt > 20 ) return 55350;
    if( el_pt > 10 ) return 68121;
    return 0.0;
  }
  if( el_eta < 2 ){ 
    if( el_pt > 95 ) return 49578;
    if( el_pt > 90 ) return 50082;
    if( el_pt > 85 ) return 50358;
    if( el_pt > 80 ) return 50193;
    if( el_pt > 75 ) return 50708;
    if( el_pt > 70 ) return 50537;
    if( el_pt > 65 ) return 50322;
    if( el_pt > 60 ) return 50720;
    if( el_pt > 55 ) return 50076;
    if( el_pt > 50 ) return 50555;
    if( el_pt > 45 ) return 50261;
    if( el_pt > 40 ) return 49947;
    if( el_pt > 35 ) return 49041;
    if( el_pt > 30 ) return 48351;
    if( el_pt > 25 ) return 46704;
    if( el_pt > 20 ) return 43626;
    if( el_pt > 10 ) return 55808;
    return 0.0;
  }
  if( el_eta < 2.1 ){ 
    if( el_pt > 95 ) return 24458;
    if( el_pt > 90 ) return 24811;
    if( el_pt > 85 ) return 24682;
    if( el_pt > 80 ) return 24646;
    if( el_pt > 75 ) return 25254;
    if( el_pt > 70 ) return 25106;
    if( el_pt > 65 ) return 25073;
    if( el_pt > 60 ) return 25313;
    if( el_pt > 55 ) return 25093;
    if( el_pt > 50 ) return 25157;
    if( el_pt > 45 ) return 25310;
    if( el_pt > 40 ) return 25348;
    if( el_pt > 35 ) return 24817;
    if( el_pt > 30 ) return 24641;
    if( el_pt > 25 ) return 24511;
    if( el_pt > 20 ) return 23796;
    if( el_pt > 10 ) return 33407;
    return 0.0;
  }
  if( el_eta < 2.2 ){ 
    if( el_pt > 95 ) return 23971;
    if( el_pt > 90 ) return 24140;
    if( el_pt > 85 ) return 24021;
    if( el_pt > 80 ) return 24273;
    if( el_pt > 75 ) return 23960;
    if( el_pt > 70 ) return 24085;
    if( el_pt > 65 ) return 24132;
    if( el_pt > 60 ) return 24303;
    if( el_pt > 55 ) return 24450;
    if( el_pt > 50 ) return 24511;
    if( el_pt > 45 ) return 24325;
    if( el_pt > 40 ) return 24191;
    if( el_pt > 35 ) return 24301;
    if( el_pt > 30 ) return 24104;
    if( el_pt > 25 ) return 23691;
    if( el_pt > 20 ) return 22936;
    if( el_pt > 10 ) return 33778;
    return 0.0;
  }
  if( el_eta <= 2.4 ){ 
    if( el_pt > 95 ) return 45660;
    if( el_pt > 90 ) return 46147;
    if( el_pt > 85 ) return 46462;
    if( el_pt > 80 ) return 46245;
    if( el_pt > 75 ) return 45964;
    if( el_pt > 70 ) return 45963;
    if( el_pt > 65 ) return 46431;
    if( el_pt > 60 ) return 46485;
    if( el_pt > 55 ) return 46067;
    if( el_pt > 50 ) return 46443;
    if( el_pt > 45 ) return 46323;
    if( el_pt > 40 ) return 46183;
    if( el_pt > 35 ) return 46043;
    if( el_pt > 30 ) return 45490;
    if( el_pt > 25 ) return 44732;
    if( el_pt > 20 ) return 43611;
    if( el_pt > 10 ) return 60727;
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
