#ifndef babymakercommon_h
#define babymakercommon_h

#include <vector>

#include "Math/LorentzVector.h"
#include "Math/PxPyPzE4D.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;
typedef std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > VofP4;

float deltaPhi (float, float);
double dRbetweenVectors (const LorentzVector&, const LorentzVector&);
bool sortByPt (const LorentzVector&, const LorentzVector&);
// Only use this for pfjets as it is currently written!
bool sortByPFJetPt (const unsigned int&, const unsigned int&);
bool isGoodPFJet(unsigned int pfjeti);

#endif
