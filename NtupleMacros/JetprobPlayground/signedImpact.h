// -*- C++ -*-

#ifndef SIGNED_IMPACT_H
#define SIGNED_IMPACT_H

#include <vector>
#include "Math/LorentzVector.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

std::pair<double, double> signedImpact (const LorentzVector &jet, unsigned int trkidx);
void getSignedImpact(unsigned int         trkJetIndex,
                     std::vector<double>& signedImpactVec,
                     std::vector<double>& signedImpactErrVec,
		     std::vector<unsigned int>&    trackIndexVec);

double getJetProb(std::vector<double>& signedImpactVec,
                  std::vector<double>& signedImpactErrVec,
                  double               maxSignedImpact,
                  double               minSDSignificance);

#endif
