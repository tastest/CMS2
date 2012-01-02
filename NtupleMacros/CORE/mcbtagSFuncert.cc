#include <math.h>
#include <stdlib.h>
#include <set>
#include "TDatabasePDG.h"
#include "Math/VectorUtil.h"
#include "CMS2.h"
#include "Math/LorentzVector.h"
#include "TSystem.h"
#include "TFile.h"
#include "TAxis.h"
#include <iostream>

//------------------------------------------------------------------------
//------------------------------------------------------------------------
// Code by Claudio - Jan 1, 2012
// The function btagEventWeight returns the "event weight" for an event with 
// two or more btags due to the "btagging scale factors". 
// This event_weight is meant to multiply the evt_scale1fb variable and
// whatever other scale factors (eg: lepton efficiency scale factors)
// Note: this only applies if the btag scale factors are < 1
//
// The function btagEventUncertainty returns the "event uncertainty" due to 
// btag scale factors.  It is a relative uncertainty.
// It should be used as follows:
// - Call this function for each event which passes the full selection, including
//   the >=2 btag requirement
// - Multiply the "event uncertainty" by the "event weight" defined above as
//   well evt_scale1fb, etc.  Let's call this quantity "scaled event uncertainty"
//   which is an absolute uncertainty.
// - The sum of "scaled event uncertainty" over all events passing the cuts
//   is the systematic uncertainty on the yield in 1 inv fb.
//
//-----------------------------------------------------------------------
//----------------------------------------------------------------------
// The btagging scale factor and its uncertainty as a function of pt is 
// hardwired in the functions btagScaleFactor and btagScaleFactorError  
// THE VALUES ARE NOT QUITE RIGHT.  PLEASE PUT THE RIGHT ONES IN
float btagScaleFactor(float jetpt) {
  if (jetpt < 240.) return 0.96;
  return 0.96;
}
float btagScaleFactorError(float jetpt) {
  if (jetpt < 240.) return 0.04;
  return 0.144;
}
//
//--------------------------------------------------------------------------
// In order to calculate the "event uncertainty" we need the actual values
// of the btagging efficiencies (for data).  These come from some database 
// or some plots or something.  
// THE FUNCTION BELOW IS JUST A PLACE HOLDER.  PLEASE FIX IT.
// Note: these are meant to be the efficiencies for jets in the fiducial region,
// i.e., something like abs(eta)<2.5.
// The btag efficiency does not need to be perfect, since it is only used for
// calculating uncertainties
float btagEff(float jetpt) {
  return 0.643;
}
//------------------------------------------------------------------------
// Here comes btagEventWeight
// Inputs:
// nbjets = number of reconstructed tagged jets (must be between 2 and 4; if
//          there are 5 or more btag jets, set nbjets=4 and only pass the
//          forst four to this function)
// pt1    = pt of the first  btagged jet
// pt2    = pt of the second btagged jet
// pt3    = pt of the third  btagged jet
// pt4    = pt of the fourth btegged jet
// (note: these do not need to be truth matched)
// Returns a negative number if something goes wrong
float btagEventWeight(int nbjets, float pt1, float pt2, float pt3=0., float pt4=0.){

  // protect against bad input values
  if (nbjets < 2 || nbjets > 4) {
    std::cout << "Illegal nbjets = " << nbjets << " in btagEventWeight" <<endl;
    return -1.;
  }

  // get the scale factors for these jets;
  float sf1 = btagScaleFactor(pt1);
  float sf2 = btagScaleFactor(pt2);
  float sf3 = btagScaleFactor(pt3);
  float sf4 = btagScaleFactor(pt4);

  // calculate the event weight for the two jet case
  if (nbjets == 2) return sf1*sf2;

  // this is the 3 jet case
  if (nbjets == 3) {
    float temp = sf1*sf2+sf1*sf3+sf2*sf3-2*sf1*sf2*sf3;
    if (temp < 0.) {
      std::cout << "Something wrong in btagEventWeight (3 jet case)" << std::endl;
      return -1.;
    }
    return temp;
  }

  // we get here if there are 4 jets
  float temp = sf1*sf2 + sf1*sf3 + sf1*sf4 + sf2*sf3 + sf2*sf4 + sf3*sf4;
  temp = temp + 3*sf1*sf2*sf3*sf4;
  temp = temp - 2*(sf1*sf2*sf3 + sf1*sf2*sf4 + sf1*sf3*sf4 + sf2*sf3*sf4);
  if (temp < 0.) {
    std::cout << "Something wrong in btagEventWeight (4 jet case)" << std::endl;
    return -1.;
  }
  return temp;
}
//------------------------------------------------------------------------
//------------------------------------------------------------------------
// Here comes btagEventUncertainty.
// Note this is quite approximate, but should be good enough 
// as an uncertainty.
// Inputs:
// nbjets = number of b quarks at status = 3 (must be between 2 and 4)
// pt1    = pt of the first  b quark
// pt2    = pt of the second b quark
// pt3    = pt of the third  b quark
// pt4    = pt of the fourth b quark
// eta1   = eta of the first  b quark
// eta2   = eta of the second b quark
// eta3   = eta of the third  b quark
// eta4   = eta of the fourth b quark
float btagEventUncertainty(int nbjets, float pt1, float eta1, float pt2, float eta2, float pt3=0., float eta3=0., float pt4=0., float eta4=0.) {
  // protect against bad input values
  if (nbjets < 2 || nbjets > 4) {
    std::cout << "Illegal nbjets = " << nbjets << " in btagEventUncertainty" <<endl;
    return -1.;
  }

  // Count fiducial jets and load arrays of efficiencies and errors

  float minpt  = 40.;
  float etacut = 2.5;
  int mynbjet  = 0;
  float eff[4];
  float effErr[4];
  if (pt1 > minpt && fabs(eta1) < etacut) {
    eff[mynbjet]    =  btagEff(pt1);
    effErr[mynbjet] =  eff[mynbjet]*btagScaleFactorError(pt1)/btagScaleFactor(pt1);
    mynbjet++;
  }
  if (pt2 > minpt && fabs(eta2) < etacut) {
    eff[mynbjet]    =  btagEff(pt2);
    effErr[mynbjet] =  eff[mynbjet]*btagScaleFactorError(pt2)/btagScaleFactor(pt2);
    mynbjet++;
  }
  if (nbjets >= 3 && pt3 > minpt && fabs(eta3) < etacut) {
    eff[mynbjet]    =  btagEff(pt3);
    effErr[mynbjet] =  eff[mynbjet]*btagScaleFactorError(pt3)/btagScaleFactor(pt3);
    mynbjet++;
  }
  if (nbjets == 4 && pt4 > minpt && fabs(eta4) < etacut) {
    eff[mynbjet]    =  btagEff(pt4);
    effErr[mynbjet] =  eff[mynbjet]*btagScaleFactorError(pt4)/btagScaleFactor(pt4);
    mynbjet++;
  }

  // If we find < 2 bquarks, we must be dealing with edge effects 
  // We will pretend that we have 2 bjets at threshold
  if (mynbjet < 2) {
    eff[1]    = btagEff(minpt);
    effErr[1] = eff[1]*btagScaleFactorError(minpt)/btagScaleFactor(minpt);
    if (mynbjet == 0) {
      eff[0]    = eff[1];
      effErr[1] = effErr[0];
    }
    mynbjet = 2;
  }

  // Here for two jets
  if (mynbjet == 2) return (eff[0]*effErr[1]+eff[1]*effErr[0])/(eff[0]*eff[1]);

  // Here for three jets
  if (mynbjet == 3) {
    float eps = eff[0]*eff[1] + eff[0]*eff[2] + eff[1]*eff[2] - 2*eff[0]*eff[1]*eff[2];
    if (eps < 0) {
      std::cout << "Negative eps in btagEventUncertainty (3 jet case)" << std::endl;
      return 0.;
    }
    float temp = (eff[1]+eff[2]-2*eff[1]*eff[2])*effErr[0];
    temp       = (eff[0]+eff[2]-2*eff[0]*eff[2])*effErr[1] + temp; 
    temp       = (eff[0]+eff[1]-2*eff[0]*eff[1])*effErr[2] + temp; 
    return fabs(temp)/eps;
  }
  
  // Here for four jets
  if (mynbjet == 4) {
    float eps = 0.;
    for (int i=0; i<4; i++) {
      for (int j=i+1; j<4; j++) {
	eps = eps + eff[i]*eff[j];
      }
    }
    eps = eps + 3*eff[0]*eff[1]*eff[2]*eff[3];
    eps = eps - 2*(eff[0]*eff[1]*eff[2] + eff[0]*eff[1]*eff[3] + 
                   eff[0]*eff[2]*eff[3] + eff[1]*eff[2]*eff[3]);
    if (eff < 0) {
      std::cout << "Negative eps in btagEventUncertainty (4 jet case)" << std::endl;
      return 0.;
    }
    float temp =  effErr[0]*(eff[1] + eff[2] + eff[3] + 3*eff[1]*eff[2]*eff[3]);
    temp = temp + effErr[1]*(eff[0] + eff[2] + eff[3] + 3*eff[0]*eff[2]*eff[3]);
    temp = temp + effErr[2]*(eff[0] + eff[1] + eff[3] + 3*eff[0]*eff[1]*eff[3]);
    temp = temp + effErr[3]*(eff[0] + eff[1] + eff[2] + 3*eff[0]*eff[1]*eff[2]);
    temp = temp - 2 * effErr[0] * (eff[1]*eff[2] + eff[1]*eff[3] + eff[2]*eff[3]);
    temp = temp - 2 * effErr[1] * (eff[0]*eff[2] + eff[0]*eff[3] + eff[2]*eff[3]);
    temp = temp - 2 * effErr[2] * (eff[0]*eff[1] + eff[0]*eff[3] + eff[1]*eff[3]);
    temp = temp - 2 * effErr[3] * (eff[0]*eff[1] + eff[0]*eff[2] + eff[1]*eff[2]);
    return fabs(temp)/eps;
  }

  // We should never get here!!!
  std::cout << "Something wrong in btagEventUncertainty: mynbjet = " << mynbjet << std::endl; 
  return 0;
}

