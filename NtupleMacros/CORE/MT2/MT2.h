#include "MT2Utility.cc"
#include "Math/LorentzVector.h"

// leptonic (WW) mt2 declaration
double MT2(	
	const float,
	const float,	
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >,
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >
);

// N jets (tt) mt2 declaration
double MT2J(
  const float,
  const float,
  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >,
  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >,
	const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > >
);


