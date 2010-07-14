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

  return answer;
}

//-----------------------------------------------------------
// Function that corrects the met (or tcmet) for a given
// muon in case this was not done in reco.  Uses value maps
//-------------------------------------------------------------
void fixMetForThisMuon(int imu, float& metX, float& metY, whichMetType type) {
  bool wasItCorrected = wasMetCorrectedForThisMuon(imu, type);
  if (!wasItCorrected) {
    switch(type) {

      case usingTcMet:
	if (cms2.mus_tcmet_flag()[imu] == 0) {
	  metX += cms2.mus_met_deltax()[imu] - cms2.mus_p4()[imu].x();
	  metY += cms2.mus_met_deltay()[imu] - cms2.mus_p4()[imu].y();
	} else if (cms2.mus_tcmet_flag()[imu] == 4) {
	     metX += - cms2.mus_tcmet_deltax()[imu] + cms2.trks_trk_p4()[cms2.mus_trkidx()[imu]].px() // undo the pion correction
		  + cms2.mus_met_deltax()[imu] - cms2.mus_p4()[imu].x(); // perform the muon correction
	     metY += - cms2.mus_tcmet_deltay()[imu] + cms2.trks_trk_p4()[cms2.mus_trkidx()[imu]].py() // undo the pion correction
		  + cms2.mus_met_deltay()[imu] - cms2.mus_p4()[imu].y(); // perform the muon correction
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

//---------------------------------------------
// correct tcMET for any hypothesis muons
// that have not been corrected for
//---------------------------------------------
metStruct correctTCMETforHypMuons (int hyp_index, float met_x, float met_y, float sumet)
{
     metStruct tcmetStruct;
     tcmetStruct.met     = sqrt(met_x * met_x + met_y * met_y);
     tcmetStruct.metphi  = atan2(met_y, met_x);
     tcmetStruct.metx    = met_x;
     tcmetStruct.mety    = met_y;
     tcmetStruct.sumet   = sumet; 

     if (cms2.hyp_type()[hyp_index] ==3)
	  return tcmetStruct;

     unsigned int i_lt = cms2.hyp_lt_index()[hyp_index];
     unsigned int i_ll = cms2.hyp_ll_index()[hyp_index];

     if (abs(cms2.hyp_lt_id()[hyp_index]) == 13)
     {
	  if(cms2.mus_tcmet_flag()[i_lt] == 0)
	  {
	       met_x += cms2.mus_met_deltax()[i_lt] - cms2.mus_p4()[i_lt].x();
	       met_y += cms2.mus_met_deltay()[i_lt] - cms2.mus_p4()[i_lt].y();
	  }
	  else if (cms2.mus_tcmet_flag()[i_lt] == 4)
	  {
	       met_x += -cms2.mus_tcmet_deltax()[i_lt] + cms2.mus_met_deltax()[i_lt] - cms2.mus_p4()[i_lt].x() + cms2.trks_trk_p4()[cms2.mus_trkidx()[i_lt]].x(); 
	       met_y += -cms2.mus_tcmet_deltay()[i_lt] + cms2.mus_met_deltay()[i_lt] - cms2.mus_p4()[i_lt].y() + cms2.trks_trk_p4()[cms2.mus_trkidx()[i_lt]].y(); 
	  }
     }
     if (abs(cms2.hyp_ll_id()[hyp_index]) == 13)
     {
	  if(cms2.mus_tcmet_flag()[i_ll] == 0)
	  { 
	       met_x += cms2.mus_met_deltax()[i_ll] - cms2.mus_p4()[i_ll].x(); 
	       met_y += cms2.mus_met_deltay()[i_ll] - cms2.mus_p4()[i_ll].y(); 
	  }
	  else if (cms2.mus_tcmet_flag()[i_ll] == 4)
	  { 
	       met_x += -cms2.mus_tcmet_deltax()[i_ll] + cms2.mus_met_deltax()[i_ll] - cms2.mus_p4()[i_ll].x() + cms2.trks_trk_p4()[cms2.mus_trkidx()[i_ll]].x();  
	       met_y += -cms2.mus_tcmet_deltay()[i_ll] + cms2.mus_met_deltay()[i_ll] - cms2.mus_p4()[i_ll].y() + cms2.trks_trk_p4()[cms2.mus_trkidx()[i_ll]].y();  
	  } 
     }

     tcmetStruct.met     = sqrt(met_x * met_x + met_y * met_y);
     tcmetStruct.metphi  = atan2(met_y, met_x);
     tcmetStruct.metx    = met_x;
     tcmetStruct.mety    = met_y;
     tcmetStruct.sumet   = sumet; 

     return tcmetStruct;
}
