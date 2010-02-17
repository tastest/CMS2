#ifndef METSELECTIONS_H
#define METSELECTIONS_H

#include "../Tools/tcmet_looper/getTcmetFromCaloMet.h"

//
// met selections
//

//---------------------------------------------
// function to calculate projected MET.
// takes three parameters as input:
//
// 1. met
// 2. met phi
// 3. hypothesis index
//---------------------------------------------
float projectedMET( float met, float metPhi, int hyp_index );

//---------------------------------------------
// function to correct tcMET for electron bug
//---------------------------------------------
metStruct correctedTCMET();

//---------------------------------------------
// utility function find deltaPhi between met
// and nearest hypothesis lepton
//---------------------------------------------
float nearestHypLeptonPhi( flaot metPhi, int hyp_index);

#endif

