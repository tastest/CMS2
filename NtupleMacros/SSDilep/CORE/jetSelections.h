#include "TLorentzVector.h"
#include "Math/LorentzVector.h"
#include <vector>
#include <iostream>
#include <utility>

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

bool comparePt(const LorentzVector &lv1,
               const LorentzVector &lv2);

std::vector<LorentzVector> getCaloJets(int i_hyp); 

