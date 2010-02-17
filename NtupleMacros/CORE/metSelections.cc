//
// met selections
//

#include "metSelections.h"
#include "CMS2.h"
#include "TMath.h"
#include "../Tools/tcmet_looper/getResponseFunction_fit.C"
#include "selections.h"

//---------------------------------------------
// function to calculate projected MET.
// takes three parameters as input:
//
// 1. met
// 2. met phi
// 3. hypothesis index
//---------------------------------------------
float projectedMET( float met, float metPhi, int hyp_index ) {

     float deltaPhi = nearestHypLeptonPhi(metPhi, hyp_index);

     return (deltaPhi < TMath::Pi() / 2.)) ? met * sin(deltaPhi) : met;
}

//---------------------------------------------
// function to correct tcMET for electron bug
//---------------------------------------------
metStruct correctedTCMET() {

     TH2F* rf = getResponseFunction_fit();

     metStruct met = getTcmetFromCaloMet(rf);

     return met;
}

//---------------------------------------------
// utility function find deltaPhi between met
// and nearest hypothesis lepton
//---------------------------------------------
float nearestHypLeptonPhi( float metPhi, int hyp_index ) {

     //WARNING!  This was designed to work in a dilepton environment - NOT a trilepton 
     float tightDPhi = min(fabs(cms2.hyp_lt_p4()[hyp_index].phi() - metPhi), 2 * TMath::Pi() - fabs(cms2.hyp_lt_p4()[i_hyp].phi() - metPhi));
     float looseDPhi = min(fabs(cms2.hyp_ll_p4()[hyp_index].phi() - metPhi), 2 * TMath::Pi() - fabs(cms2.hyp_ll_p4()[i_hyp].phi() - metPhi));
     
     return min(tightDPhi, looseDPhi);

}
