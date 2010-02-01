// CMS2 includes
#include "muonSelections.h"
#include "CMS2.h"
#include <iostream>

bool muonId(unsigned int index, SelectionType type){
  switch (type) {
  case Nominal:
    if (((cms2.mus_type().at(index)) & (1<<1)) == 0) return false; // global muon
    if (((cms2.mus_type().at(index)) & (1<<2)) == 0) return false; // tracker muon
    if (cms2.mus_validHits().at(index) < 11)    return false; 
    if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) return false; 
    if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.4) return false;
    return true;
    break;
  default:
    std::cout << "ERROR: requested muon type is not defined. Abort." << std::endl;
    return false;
  }
}

bool muond0(unsigned int index, ReferencePointType type, double cut){
  switch (type){
  case BeamSpot:
    return TMath::Abs(cms2.mus_d0corr().at(index)) > cut;
    break;
  case ZeroZeroZero:
    return TMath::Abs(cms2.mus_d0().at(index)) > cut;
    break;
  default:
    std::cout << "ERROR: unsupported d0 selection. Abort." << std::endl;
    return false;
  }
}

double muonIsoValue(unsigned int index, double cut){
  // ECAL veto cone 
  if (cms2.mus_iso_ecalvetoDep().at(index) > 4) return false; // ECalE < 4 
  // HCAL veto cone
  if (cms2.mus_iso_hcalvetoDep().at(index) > 6) return false; // HCalE < 6 
  double sum =  cms2.mus_iso03_sumPt().at(index) +
    cms2.mus_iso03_emEt().at(index)  +
    cms2.mus_iso03_hadEt().at(index);
  double pt  = cms2.mus_p4().at(index).pt();
  return sum/max(pt,20.);
}

bool isCosmics(unsigned int index){
  return false;
}


