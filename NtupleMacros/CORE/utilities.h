#ifndef UTILITIES_H
#define UTILITIES_H

#include "TH1F.h"
#include "TH2D.h"
#include "TVector3.h"
#include <algorithm>
#include <set>
#include "Math/VectorUtil.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

//TVector3 correctMETforTracks() ;
double trkIsolation(int trk_index);
double dRbetweenVectors(const LorentzVector & vec1, 
			const LorentzVector & vec2 );
bool   sortByPt(const LorentzVector & vec1, const LorentzVector & vec2 );

class DorkyEvent;
class DorkyEventIdentifier;

extern TH2D *rfhist;

#endif
