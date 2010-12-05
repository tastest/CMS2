#ifndef MCSUSYKFACTOR_H
#define MCSUSYKFACTOR_H
//----------------------------------------------------
// k-factors for LMx samples 
//
// Usage:
//       kfactorSUSY("lm0");
//
//--------------------------------------------------------

#include "CMS2.h"

float kfactorSUSY(string sample);
double lmdata(int ipart1, int ipart2, string prefix);
int sfinalState(int ipart1, int ipart2);

#endif

