#ifndef UTILITIES_H
#define UTILITIES_H

#include "TH1F.h"
#include "TH2D.h"
#include "TVector3.h"
#include <algorithm>
#include <set>
#include "Math/VectorUtil.h"

TVector3 correctMETforTracks() ;
double trkIsolation(int trk_index);
double dRbetweenVectors(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > vec1, 
			ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > vec2 );

extern TH2D *rfhist;

#endif
