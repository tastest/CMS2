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
        case muonSelectionFO_mu_v1:
            isovalue = 0.40;
            break;
        case muonSelectionFO_mu_ttbar:
            isovalue = 0.40;
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

    if ( cms2.mus_p4()[index].pt() < 5.0) {
        std::cout << "muonID ERROR: requested muon is too low pt,  Abort." << std::endl;
        return false;
    }


    switch (type) {

        case Nominal:
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
            if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.5)  return false; // eta cut
            if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) return false; //glb fit chisq
            if (((cms2.mus_type().at(index)) & (1<<1)) == 0)    return false; // global muon
            if (((cms2.mus_type().at(index)) & (1<<2)) == 0)    return false; // tracker muon
            if (cms2.mus_validHits().at(index) < 11)            return false; // # of tracker hits
            if (TMath::Abs(cms2.mus_d0corr().at(index)) > 0.02) return false; // d0 from beamspot
            return true;

        case muonSelectionFO_mu_v1:
            if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.5)  return false; // eta cut
            if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 50) return false;
            if (((cms2.mus_type().at(index)) & (1<<1)) == 0)    return false; // global muon
            if (((cms2.mus_type().at(index)) & (1<<2)) == 0)    return false; // tracker muon
            if (cms2.mus_validHits().at(index) < 11)            return false; // # of tracker hits
            if (cms2.mus_iso_hcalvetoDep().at(index) > 12)       return false;
            if (cms2.mus_gfit_validSTAHits().at(index) == 0)    return false; // Glb fit must have hits in mu chambers
            if (TMath::Abs(cms2.mus_d0corr().at(index)) > 0.2) return false; // d0 from beamspot
            return true;
            break;

        case muonSelectionFO_mu_ttbar:
            if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.5)  return false; // eta cut
            if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 50) return false; //glb fit chisq
            if (((cms2.mus_type().at(index)) & (1<<1)) == 0)    return false; // global muon
            if (((cms2.mus_type().at(index)) & (1<<2)) == 0)    return false; // tracker muon
            if (cms2.mus_validHits().at(index) < 11)            return false; // # of tracker hits
            if (TMath::Abs(cms2.mus_d0corr().at(index)) > 0.2) return false; // d0 from beamspot
            return true;
            break;


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

//--------------------------------------------------------------
//  Remove cosmics by looking for back-to-back muon-track pairs
//--------------------------------------------------------------
// http://indico.cern.ch/contributionDisplay.py?contribId=2&confId=86834
bool isCosmics(unsigned int index){
    for (int itrk=0; itrk < int(cms2.trks_trk_p4().size()); ++itrk) {
        const LorentzVector& mu_p4  = cms2.mus_trk_p4().at(index);
        const LorentzVector& trk_p4 = cms2.trks_trk_p4().at(itrk);
        double sprod = mu_p4.px()*trk_p4.px()+mu_p4.py()*trk_p4.py()+mu_p4.pz()*trk_p4.pz();
        if ( acos( -(sprod/trk_p4.P()/mu_p4.P()) ) < 0.01 &&
                fabs(trk_p4.pt()-mu_p4.pt())/mu_p4.pt() < 0.05 )
            return true;
    }
    return false;
}

bool passedMuonTriggerRequirements()
{
    return cms2.passHLTTrigger("HLT_Mu9");
}
