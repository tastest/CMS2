//
// met selections
//

#include <math.h>
#include <algorithm>
#include "TMath.h"
#include "CMS2.h"
#include "trackSelections.h"
#include "metSelections.h"
#include "Math/LorentzVector.h"


typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;


//---------------------------------------------------
// Function that checks whether met (or tcmet) was 
// corrected for a given muon.  This uses the value maps
//---------------------------------------------------
bool wasMetCorrectedForThisMuon(int imu, whichMetType type) {
  bool answer=true;
  switch(type) {
    case usingTcMet:
      if (cms2.mus_tcmet_flag().at(imu) == 0 || 
          cms2.mus_tcmet_flag().at(imu) == 4) answer = false;
      break;
    case usingCaloMet:
      if (cms2.mus_met_flag().at(imu) == 0) answer = false;
      break;
    default:
      std::cout << "Illegal call to wasMetCorrectedForThisMuon" <<std::endl;
  }
}

//-----------------------------------------------------------
// Function that corrects the met (or tcmet) for a given
// muon in case this was not done in reco.  Uses value maps
//-------------------------------------------------------------
void fixMetForThisMuon(int imu, float& metX, float& metY, whichMetTypr type) {
  bool wasItCorrected = wasMetCorrectedForThisMuon(imu, type);
  if (!wasItCorrected) {
    switch(type) {

      case usingTcMet:
	if (cms2.mus_tcmet_flag()[imu] == 0) {
	  metX += cms2.mus_met_deltax()[imu] - cms2.mus_p4()[imu].x();
	  metY += cms2.mus_met_deltay()[imu] - cms2.mus_p4()[imu].y();
	} else if (cms2.mus_tcmet_flag()[imu] == 4) {
	  metX += - mus_tcmet_deltax()[imu] + mus_met_deltax()[imu] - mus_p4()[iMu].x(); 
	  metY += - mus_tcmet_deltay()[imu] + mus_met_deltay()[imu] - mus_p4()[iMu].y(); 
	}
	break;

      case usingCaloMet:
	metX += cms2.mus_met_deltax()[imu] - cms2.mus_p4()[imu].x();
	metY += cms2.mus_met_deltay()[imu] - cms2.mus_p4()[imu].y();
	break;
    }
  }
}

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
     return ((deltaPhi < TMath::Pi() / 2.) ? met * sin(deltaPhi) : met);
}

//---------------------------------------------
// as above but simpler for single lepton events
//---------------------------------------------
float projectedMETW( float met, float metPhi, float leptonPhi) {
     float deltaPhi = acos(cos(metPhi - leptonPhi));
     return ((deltaPhi < TMath::Pi() / 2.) ? met * sin(deltaPhi) : met);
}


//---------------------------------------------
// function to calculate latest tcMET
//---------------------------------------------
#include "tcmet/getTcmetFromCaloMet.icc"
#include "tcmet/getResponseFunction_fit.icc"
metStruct correctedTCMET(bool usePV, bool useElectronVetoCone, bool useHFcleaning, bool useHCALcleaning, bool useECALcleaning, bool printout, ostream& ostr) 
{
     // static because we only want to get the response function once
     static TH2F* rf = getResponseFunction_fit();
     return getTcmetFromCaloMet(rf, useElectronVetoCone, usePV, useHFcleaning, useHCALcleaning, useECALcleaning, printout,ostr);
}

//---------------------------------------------
// calorimeter cleaning for tcMET
//---------------------------------------------
#include "tcmet/cleanTcmet.icc"
metStruct cleanTCMET (float met_x, float met_y, float sumet, bool useHFcleaning, bool useHCALcleaning, bool useECALcleaning)
{
     metStruct tcmetStruct;
     tcmetStruct.met    = sqrt(met_x * met_x + met_y * met_y);
     tcmetStruct.metphi = atan2(met_y, met_x);
     tcmetStruct.metx   = met_x;
     tcmetStruct.mety   = met_y;
     tcmetStruct.sumet  = sumet;
     
     if (useHFcleaning)
	  tcmetStruct = cleanTCMETforHFspikes(tcmetStruct.metx, tcmetStruct.mety, tcmetStruct.sumet);

     if (useHCALcleaning)
	  tcmetStruct = cleanTCMETforHCALnoise(tcmetStruct.metx, tcmetStruct.mety, tcmetStruct.sumet);

     if (useECALcleaning)
	  tcmetStruct = cleanTCMETforECALspikes(tcmetStruct.metx, tcmetStruct.mety, tcmetStruct.sumet);

     return tcmetStruct;
}


//---------------------------------------------
// utility function find deltaPhi between met
// and nearest hypothesis lepton
//---------------------------------------------
float nearestHypLeptonPhi( float metPhi, int hyp_index ) {

     //WARNING!  This was designed to work in a dilepton environment - NOT a trilepton 
     float tightDPhi = min(fabs(cms2.hyp_lt_p4()[hyp_index].phi() - metPhi), (float)(2 * TMath::Pi()) - fabs(cms2.hyp_lt_p4()[hyp_index].phi() - metPhi));
     float looseDPhi = min(fabs(cms2.hyp_ll_p4()[hyp_index].phi() - metPhi), (float)(2 * TMath::Pi()) - fabs(cms2.hyp_ll_p4()[hyp_index].phi() - metPhi));
     
     return min(tightDPhi, looseDPhi);

}
