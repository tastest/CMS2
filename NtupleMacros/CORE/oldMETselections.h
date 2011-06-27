#ifndef OLD_MET_SELECTION_H
#define OLD_MET_SELECTION_H
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


//---------------------------------------------
// function to calculate latest tcMET
//---------------------------------------------
metStruct correctedTCMET(bool usePV = false, bool useElectronVetoCone = false, bool useHFcleaning = true, 
                         bool useHCALcleaning = true, bool useECALcleaning = false, bool printout = false, ostream& ostr = std::cout);

//---------------------------------------------
// calorimeter cleaning for tcMET
//---------------------------------------------
metStruct cleanTCMET(float met_x, float met_y, float sumet, bool useHFcleaning = true, bool useHCALcleaning = true, bool useECALcleaning = false);

#endif
