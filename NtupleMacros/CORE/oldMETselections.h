#ifndef SELECTION_H
#define SELECTION_H
#include "TLorentzVector.h"
#include "Math/LorentzVector.h"
#include <vector>
#include <iostream>
#include <utility>

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

bool pass2Met(int index, const class TVector3& corr);
double MetSpecial(double MET, double MetPhi, int index);
bool pass4Met(int index, const class TVector3& corr);
bool met10(int index, const class TVector3& corr);
bool met1(int index, const class TVector3& corr);
bool sumEt10(double sumEt);
bool sumEt1(double sumEt);
bool passMet_OF20_SF30(int hypIdx, bool useTcMet);
bool metSimple (float threshold, const TVector3& corr);
bool metBalance (int i_hyp, const TVector3& corr);
bool metProjected (int i_hyp, const TVector3& corr);
bool pass5Met (int i_hyp, const TVector3& corr);
bool passMetVJets09(float value, bool useTcMet);

#endif
