#ifndef GETTCMETFROMCALOMET_H
#define GETTCMETFROMCALOMET_H
#include "TLorentzVector.h"
#include "Math/LorentzVector.h"
#include <vector>
#include "TH2.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

struct metStruct{
  metStruct() : met(-999.), metphi(-999.), metx(-999.), mety(-999.), sumet(-999.)  {}
  float met;
  float metphi;
  float metx;
  float mety;
  float sumet;
};

metStruct getTcmetFromCaloMet(TH2F* rf);
bool  isMuon               ( int index );
bool  isElectron           ( int index );
bool  isGoodTrack          ( int index );
void  correctMETforMuon    ( LorentzVector p4, int index );
void  correctMETforPion    ( int index , TH2F* rf);
void  correctSumEtForMuon  ( LorentzVector p4, int index );
void  correctSumEtForPion  ( int index , TH2F* rf);

//float tcmet_x;
//float tcmet_y;
//float tcsumet;
//TH2F* rf;

#endif
