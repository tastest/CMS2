#include "CMS2.h"
#include "muonSelections.h"
#include <iostream>
//------------------------------------------------------------------
// Apply various muon identification requirements
//------------------------------------------------------------------
bool muonId(unsigned int index, SelectionType type){

  float isovalue;

  switch(type) {
  case Nominal:
    isovalue = 0.1;
    break;
  case NominalTTbar:
    isovalue = 0.15;
    break;
  default:
    std::cout << "muonID ERROR: requested muon type is not defined. Abort." << std::endl;
    return false;
  } 
  return 
    muonIdNotIsolated( index, type ) && 
    muonIsoValue(index) < isovalue;          // Isolation cut
}
bool muonIdNotIsolated(unsigned int index, SelectionType type){
  switch (type) {

  case Nominal:
    if ( cms2.mus_p4()[index].pt() < 5.) {
      std::cout << "muonID ERROR: requested muon is too low pt,  Abort." << std::endl;
      return false;
    }
    if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.5)  return false; // eta cut
    if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) return false; //glb fit chisq
    if (((cms2.mus_type().at(index)) & (1<<1)) == 0)    return false; // global muon
    if (((cms2.mus_type().at(index)) & (1<<2)) == 0)    return false; // tracker muon
    if (cms2.mus_validHits().at(index) < 11)            return false; // # of tracker hits
    if (cms2.mus_iso_ecalvetoDep().at(index) > 4)       return false; // ECalE < 4 
    if (cms2.mus_iso_hcalvetoDep().at(index) > 6)       return false; // HCalE < 6 
    if (cms2.mus_gfit_validSTAHits().at(index) == 0)    return false; // Glb fit must have hits in mu chambers
    if (TMath::Abs(cms2.mus_d0corr().at(index)) > 0.02) return false; // d0 from beamspot
    return true;
    break;

  case NominalTTbar:
    if ( cms2.mus_p4()[index].pt() < 5.) {
      std::cout << "muonID ERROR: requested muon is too low pt,  Abort." << std::endl;
      return false;
    }
    if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.5)  return false; // eta cut
    if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) return false; //glb fit chisq
    if (((cms2.mus_type().at(index)) & (1<<1)) == 0)    return false; // global muon
    if (((cms2.mus_type().at(index)) & (1<<2)) == 0)    return false; // tracker muon
    if (cms2.mus_validHits().at(index) < 11)            return false; // # of tracker hits
    if (TMath::Abs(cms2.mus_d0corr().at(index)) > 0.02) return false; // d0 from beamspot
    return true;

  default:
    std::cout << "muonID ERROR: requested muon type is not defined. Abort." << std::endl;
    return false;
  }
}
//------------------------------------------------------------------------------
// Calculate relative muon isolation variable.  Below 20 GeV, special treatment
//-------------------------------------------------------------------------------
double muonIsoValue(unsigned int index){
  double sum =  cms2.mus_iso03_sumPt().at(index) +
    cms2.mus_iso03_emEt().at(index)  +
    cms2.mus_iso03_hadEt().at(index);
  double pt  = cms2.mus_p4().at(index).pt();
  return sum/max(pt,20.);
}
//-----------------------------------------------
// Nothing here at the moment?
//-----------------------------------------------
bool isCosmics(unsigned int index){
  return false;
}

bool passedMuonTriggerRequirements()
{
  return cms2.passHLTTrigger("HLT_Mu9");
}
