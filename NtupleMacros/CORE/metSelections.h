#ifndef METSELECTIONS_H
#define METSELECTIONS_H

//
// met selections
//

struct metStruct{
  metStruct() : met(-999.), metphi(-999.), metx(-999.), mety(-999.), sumet(-999.)  {}
  float met;
  float metphi;
  float metx;
  float mety;
  float sumet;
};


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
// as above but simpler for single lepton events
//---------------------------------------------
float projectedMETW( float met, float metPhi, float leptonPhi);

//---------------------------------------------
// function to calculate latest tcMET
//---------------------------------------------
metStruct correctedTCMET(bool usePV = false, bool useElectronVetoCone = false, bool useHFcleaning = true, 
                         bool useHCALcleaning = true, bool useECALcleaning = false, bool printout = false, ostream& ostr = cout);

//---------------------------------------------
// calorimeter cleaning for tcMET
//---------------------------------------------
metStruct cleanTCMET(float met_x, float met_y, float sumet, bool useHFcleaning = true, bool useHCALcleaning = true, bool useECALcleaning = false);

//---------------------------------------------
// utility function find deltaPhi between met
// and nearest hypothesis lepton
//---------------------------------------------
float nearestHypLeptonPhi( float metPhi, int hyp_index);

#endif

