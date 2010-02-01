#ifndef MUON_SELECTIONS_H
#define MUON_SELECTIONS_H
#include "CORE/CMS2.h"
//
// muon id
// see http://www.t2.ucsd.edu/tastwiki/bin/view/CMS/MuonID
//

enum ReferencePointType { 
  ZeroZeroZero = 0, 
  BeamSpot = 1, 
  PrimaryVertex = 2 
};

enum SelectionType { 
  FakableObject = 0, 
  Nominal = 1, 
  Tight = 2 
}; 

//
// desined for high Pt muons (above 10 GeV). Nominal selection
// is tight enough for most if not all analysis. Tight 
// version is a backup.
//
bool muonId(unsigned int index, SelectionType type = Nominal);

//
// In order to caclulate the impact parameter we need to define
// the reference point. 
//
bool muond0(unsigned int index, ReferencePointType type = BeamSpot, double cut = 0.2);

//
// relative isolation
//
double muonIsoValue(unsigned int index);
bool muonIso(unsigned int index, double cut = 0.1){
  return muonIsoValue(index) < cut;
}

//
// Trigger
// for now we don't match muons to triggered objects
//
bool passedMuonTriggerRequirements(){
  return cms2.passHLTTrigger("HLT_Mu9");
}

//
// cosmics rejection
//
bool isCosmics(unsigned int index);

#endif

