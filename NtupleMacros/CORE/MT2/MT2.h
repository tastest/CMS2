#include "MT2Utility.cc"
#include "Math/LorentzVector.h"

enum enum_mt2_method { BISECT, GRID };

// MT2 declaration
double MT2(
  const float,
  const float,
  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >,
  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >,
  float = 0.0,
  bool = false
);

// MT2J declaration
double MT2J(
  const float,
  const float,
  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >,
  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >,
  const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > >,
  float = 0.0,
  enum_mt2_method = BISECT,
  bool = false
);

vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > CleanJets(
  const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > >,
  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >,
  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >,
  const float,
  const float,
  const float
);


