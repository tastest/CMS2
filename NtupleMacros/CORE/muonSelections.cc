#include "CMS2.h"
#include "muonSelections.h"
#include "trackSelections.h"
#include "utilities.h"
#include <iostream>
#include "eventSelections.h"
#include "trackSelections.h"

//------------------------------------------------------------------
// Apply various muon identification requirements
//------------------------------------------------------------------
bool muonId(unsigned int index, SelectionType type, int vertex_index){

    float isovalue;
    bool truncated = true;

    switch(type) {
    case Nominal:
        isovalue = 0.1;
        break;
    case NominalTTbar:
        isovalue = 0.15;
        break;
    case NominalTTbarV2:
        isovalue = 0.15;
        break;
    case muonSelectionFO_mu_v1:
        isovalue = 0.40;
        break;
    case muonSelectionFO_mu_ttbar:
        isovalue = 0.40;
        break;
    case muonSelectionFO_mu_ttbar_iso10:
        isovalue = 1.0;
        break;
    case NominalWWV0:
    case NominalWWV1:
        isovalue = 0.15;
        break;
    case muonSelectionFO_mu_wwV1:
    case muonSelectionFO_mu_ww:
        isovalue = 0.40;
        break;
    case muonSelectionFO_mu_smurf_04:
        if (!muonIdNotIsolated( index, type )) return false;
        return muonIsoValuePF(index,0) < 0.40;
        break;
    case muonSelectionFO_mu_wwV1_iso10_d0:
    case muonSelectionFO_mu_wwV1_iso10:
    case muonSelectionFO_mu_ww_iso10:
        isovalue = 1.0;
        break;
    case muonSelectionFO_mu_smurf_10:
        if (!muonIdNotIsolated( index, type )) return false;
        return muonIsoValuePF(index,0) < 1.0;
        break;
    case NominalSS:
    case NominalSSv2:
        isovalue = 0.1;
        break;
    case NominalSSd0PV:
        isovalue = 0.1;
        break;
    case NominalSSnod0:
        isovalue = 0.1;
        break;
    case muonSelectionFO_mu_ss:
    case muonSelectionFO_mu_ssV2:
        isovalue = 0.40;
        break;
    case muonSelectionFO_mu_ss_iso10:
    case muonSelectionFO_mu_ssV2_iso10:
        isovalue = 1.0;
        break;
    case OSGeneric_v1:
        isovalue = 0.15;
        break;
    case OSGeneric_v2:
        isovalue  = 0.15;
        break;
    case OSGeneric_v2_FO:
        isovalue  = 0.4;
        break;
    case OSGeneric_v3:
	    truncated = false;
	    isovalue  = 0.15;
        break;
    case OSGeneric_v3_FO:
	    truncated = false;
        isovalue  = 0.4;
        break;
    case OSZ_v1:
        isovalue = 0.15;
        break;
    case OSZ_v2:
        isovalue = 0.15;
        break;
    case NominalTTbar_pass6:
        isovalue = 0.15;
        break;
    case muonSelectionFO_mu_ttbar_pass6:
        isovalue = 0.4;
        break; 
    case NominalSmurfV3:
        if (!muonIdNotIsolated( index, type )) return false;
        if (cms2.mus_p4().at(index).pt()<20) 
            return muonIsoValue(index,false) < 0.1;
        else
            return muonIsoValue(index,false) < 0.15;
        break;
    case NominalSSv3:
        if (!muonIdNotIsolated(index, type, vertex_index))
            return false;
        return (muonIsoValue(index, false) < 0.15);
        break;
    case muonSelectionFO_ssV3:
        if (!muonIdNotIsolated(index, type, vertex_index))
            return false;
        return (muonIsoValue(index, false) < 0.40);
        break;
    case NominalSSv4:
        isovalue = 0.4;
        break;
    case muonSelectionFO_ssV4:
        if (!muonIdNotIsolated(index, type))
            return false;
        return (muonIsoValue(index, false) < 0.15);
        break;
    case NominalSmurfV4:
        if (!muonIdNotIsolated( index, type )) return false;
        if (cms2.mus_p4().at(index).pt()>20) {
            if (TMath::Abs(cms2.mus_p4()[index].eta())<1.479) return muonIsoValuePF(index,0) < 0.22;
            else return muonIsoValuePF(index,0) < 0.20;
        } else {
            return muonIsoValuePF(index,0) < 0.11;
        }
        break;
    case NominalSmurfV5:
        if (!muonIdNotIsolated( index, type )) return false;
        if (cms2.mus_p4().at(index).pt()>20) {
            if (TMath::Abs(cms2.mus_p4()[index].eta())<1.479) return muonIsoValuePF(index,0,0.3) < 0.13;
            else return muonIsoValuePF(index,0,0.3) < 0.09;
        } else {
            if (TMath::Abs(cms2.mus_p4()[index].eta())<1.479) return muonIsoValuePF(index,0,0.3) < 0.06;
            else return muonIsoValuePF(index,0,0.3) < 0.05;
        }
        break;
    default:
        std::cout << "muonID ERROR: requested muon type is not defined. Abort." << std::endl;
        exit(1);
        return false;
    } 
    return 
        muonIdNotIsolated( index, type, vertex_index ) && 
        muonIsoValue(index,truncated) < isovalue;          // Isolation cut
}
bool muonIdNotIsolated(unsigned int index, SelectionType type, int vertex_index){

    if ( cms2.mus_p4()[index].pt() < 5.0) {
        std::cout << "muonID ERROR: requested muon is too low pt,  Abort." << std::endl;
        return false;
    }

    int vtxidx = firstGoodDAvertex();    

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

    case NominalTTbarV2:
        if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.5)  return false; // eta cut
        if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) return false; //glb fit chisq
        if (((cms2.mus_type().at(index)) & (1<<1)) == 0)    return false; // global muon
        if (((cms2.mus_type().at(index)) & (1<<2)) == 0)    return false; // tracker muon
        if (cms2.mus_validHits().at(index) < 11)            return false; // # of tracker hits
        if (cms2.mus_gfit_validSTAHits().at(index) == 0)    return false; // Glb fit must have hits in mu chambers
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

    case muonSelectionFO_mu_ttbar_iso10:
        if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.5)  return false; // eta cut
        if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 50) return false; //glb fit chisq
        if (((cms2.mus_type().at(index)) & (1<<1)) == 0)    return false; // global muon
        if (((cms2.mus_type().at(index)) & (1<<2)) == 0)    return false; // tracker muon
        if (cms2.mus_validHits().at(index) < 11)            return false; // # of tracker hits
        if (TMath::Abs(cms2.mus_d0corr().at(index)) > 0.2) return false; // d0 from beamspot
        return true;
        break;

    case NominalWWV0:
        if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.4)  return false; // eta cut
        if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) return false; //glb fit chisq
        if (((cms2.mus_type().at(index)) & (1<<1)) == 0)    return false; // global muon
        if (((cms2.mus_type().at(index)) & (1<<2)) == 0)    return false; // tracker muon
        if (cms2.mus_validHits().at(index) < 11)            return false; // # of tracker hits  
        if (cms2.mus_gfit_validSTAHits().at(index)==0 )     return false; // Glb fit must have hits in mu chambers
        if (TMath::Abs(mud0PV(index)) >= 0.02)              return false; // d0 from pvtx
        return true;
        break;

    case muonSelectionFO_mu_wwV1:
    case muonSelectionFO_mu_wwV1_iso10:
    case NominalWWV1:
        if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.4)  return false; // eta cut
        if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) return false; //glb fit chisq
        if (((cms2.mus_type().at(index)) & (1<<1)) == 0)    return false; // global muon
        if (((cms2.mus_type().at(index)) & (1<<2)) == 0)    return false; // tracker muon
        if (cms2.mus_validHits().at(index) < 11)            return false; // # of tracker hits  
        if (cms2.mus_gfit_validSTAHits().at(index)==0 )     return false; // Glb fit must have hits in mu chambers
        if (TMath::Abs(mud0PV_wwV1(index)) >= 0.02)         return false; // d0 from pvtx
        if (TMath::Abs(mudzPV_wwV1(index)) >= 1.0)          return false; // dz from pvtx
        if (cms2.mus_ptErr().at(index)/cms2.mus_p4().at(index).pt()>0.1) return false;
        if (cms2.trks_valid_pixelhits().at(cms2.mus_trkidx().at(index))==0) return false;
        if (cms2.mus_nmatches().at(index)<2) return false;
        return true;
        break;

    case muonSelectionFO_mu_wwV1_iso10_d0: // same as muonSelectionFO_mu_wwV1_iso10 but with looser d0 cut
        if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.4)  return false; // eta cut
        if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) return false; //glb fit chisq
        if (((cms2.mus_type().at(index)) & (1<<1)) == 0)    return false; // global muon
        if (((cms2.mus_type().at(index)) & (1<<2)) == 0)    return false; // tracker muon
        if (cms2.mus_validHits().at(index) < 11)            return false; // # of tracker hits  
        if (cms2.mus_gfit_validSTAHits().at(index)==0 )     return false; // Glb fit must have hits in mu chambers
        if (TMath::Abs(mud0PV_wwV1(index)) >= 0.2)         return false; // d0 from pvtx
        if (TMath::Abs(mudzPV_wwV1(index)) >= 1.0)          return false; // dz from pvtx
        if (cms2.mus_ptErr().at(index)/cms2.mus_p4().at(index).pt()>0.1) return false;
        if (cms2.trks_valid_pixelhits().at(cms2.mus_trkidx().at(index))==0) return false;
        if (cms2.mus_nmatches().at(index)<2) return false;
        return true;
        break;

    case muonSelectionFO_mu_ww:
        if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.4)  return false; // eta cut
        if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) return false; //glb fit chisq
        if (((cms2.mus_type().at(index)) & (1<<1)) == 0)    return false; // global muon
        if (((cms2.mus_type().at(index)) & (1<<2)) == 0)    return false; // tracker muon
        if (cms2.mus_validHits().at(index) < 11)            return false; // # of tracker hits  
        if (cms2.mus_gfit_validSTAHits().at(index)==0 )     return false; // Glb fit must have hits in mu chambers
        if (TMath::Abs(mud0PV(index)) >= 0.02)              return false; // d0 from pvtx
        return true;

    case muonSelectionFO_mu_ww_iso10:
        if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.4)  return false; // eta cut
        if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) return false; //glb fit chisq
        if (((cms2.mus_type().at(index)) & (1<<1)) == 0)    return false; // global muon
        if (((cms2.mus_type().at(index)) & (1<<2)) == 0)    return false; // tracker muon
        if (cms2.mus_validHits().at(index) < 11)            return false; // # of tracker hits  
        if (cms2.mus_gfit_validSTAHits().at(index)==0 )     return false; // Glb fit must have hits in mu chambers
        if (TMath::Abs(mud0PV(index)) >= 0.02)              return false; // d0 from pvtx
        return true;


    case NominalSS:
        if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.4)  return false; // eta cut
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

    case NominalSSd0PV:
        if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.4)  return false; // eta cut
        if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) return false; //glb fit chisq
        if (((cms2.mus_type().at(index)) & (1<<1)) == 0)    return false; // global muon
        if (((cms2.mus_type().at(index)) & (1<<2)) == 0)    return false; // tracker muon
        if (cms2.mus_validHits().at(index) < 11)            return false; // # of tracker hits
        if (cms2.mus_iso_ecalvetoDep().at(index) > 4)       return false; // ECalE < 4 
        if (cms2.mus_iso_hcalvetoDep().at(index) > 6)       return false; // HCalE < 6 
        if (cms2.mus_gfit_validSTAHits().at(index) == 0)    return false; // Glb fit must have hits in mu chambers
        if (TMath::Abs(mud0PV(index)) > 0.02) return false; // d0 from beamspot
        return true;
        break;

    case NominalSSnod0:
        if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.4)  return false; // eta cut
        if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) return false; //glb fit chisq
        if (((cms2.mus_type().at(index)) & (1<<1)) == 0)    return false; // global muon
        if (((cms2.mus_type().at(index)) & (1<<2)) == 0)    return false; // tracker muon
        if (cms2.mus_validHits().at(index) < 11)            return false; // # of tracker hits
        if (cms2.mus_iso_ecalvetoDep().at(index) > 4)       return false; // ECalE < 4 
        if (cms2.mus_iso_hcalvetoDep().at(index) > 6)       return false; // HCalE < 6 
        if (cms2.mus_gfit_validSTAHits().at(index) == 0)    return false; // Glb fit must have hits in mu chambers
        return true;
        break;

    case muonSelectionFO_mu_ss:
        if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.4)  return false; // eta cut
        if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 50) return false; //glb fit chisq
        if (((cms2.mus_type().at(index)) & (1<<1)) == 0)    return false; // global muon
        if (((cms2.mus_type().at(index)) & (1<<2)) == 0)    return false; // tracker muon
        if (cms2.mus_validHits().at(index) < 11)            return false; // # of tracker hits
        if (TMath::Abs(cms2.mus_d0corr().at(index)) > 0.2) return false; // d0 from beamspot
        return true;
        break;

    case muonSelectionFO_mu_ss_iso10:
        if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.4)  return false; // eta cut
        if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 50) return false; //glb fit chisq
        if (((cms2.mus_type().at(index)) & (1<<1)) == 0)    return false; // global muon
        if (((cms2.mus_type().at(index)) & (1<<2)) == 0)    return false; // tracker muon
        if (cms2.mus_validHits().at(index) < 11)            return false; // # of tracker hits
        if (TMath::Abs(cms2.mus_d0corr().at(index)) > 0.2) return false; // d0 from beamspot
        return true;
        break;

    case OSGeneric_v1:
        //this selector is NominalTTbarV2 + eta < 2.4 + dpt/pt < 0.1
        if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.4)                       return false; // eta cut
        if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) return false; //glb fit chisq
        if (((cms2.mus_type().at(index)) & (1<<1)) == 0)                         return false; // global muon
        if (((cms2.mus_type().at(index)) & (1<<2)) == 0)                         return false; // tracker muon
        if (cms2.mus_validHits().at(index) < 11)                                 return false; // # of tracker hits
        if (cms2.mus_gfit_validSTAHits().at(index) == 0)                         return false; // Glb fit must have hits in mu chambers
        if (TMath::Abs(cms2.mus_d0corr().at(index)) > 0.02)                      return false; // d0 from beamspot
        if (cms2.mus_ptErr().at(index)/cms2.mus_p4().at(index).pt()>0.1)         return false; // dpt/pt 
        return true;

    case OSGeneric_v2:
        //baseline selector for 2011 OS analysis
        if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.4)                       return false; // eta cut
        if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) return false; // glb fit chisq
        if (((cms2.mus_type().at(index)) & (1<<1)) == 0)                         return false; // global muon
        if (((cms2.mus_type().at(index)) & (1<<2)) == 0)                         return false; // tracker muon
        if (cms2.mus_validHits().at(index) < 11)                                 return false; // # of tracker hits
        if (cms2.mus_gfit_validSTAHits().at(index) == 0)                         return false; // Glb fit must have hits in mu chambers
        if (TMath::Abs(mud0PV_smurfV3(index)) > 0.02)                            return false; // d0(PV) < 0.02 cm
        if (TMath::Abs(mudzPV_smurfV3(index)) > 1  )                             return false; // dz(PV) < 1 cm
        if (cms2.mus_ptErr().at(index)/cms2.mus_p4().at(index).pt()>0.1)         return false; // dpt/pt < 0.1
        return true;

    case OSGeneric_v3:
        // copy of OSGeneric_v2 ID, with different iso requirement
	    return muonIdNotIsolated( index, OSGeneric_v2, -1 );

    case OSGeneric_v2_FO:
	    //baseline FO for 2011 OS analysis
	    //relax cuts: chi2/ndf < 50, d0 < 0.2, reliso < 0.4
        if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.4)                       return false; // eta cut
        if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 50) return false; // glb fit chisq
        if (((cms2.mus_type().at(index)) & (1<<1)) == 0)                         return false; // global muon
        if (((cms2.mus_type().at(index)) & (1<<2)) == 0)                         return false; // tracker muon
        if (cms2.mus_validHits().at(index) < 11)                                 return false; // # of tracker hits
        if (cms2.mus_gfit_validSTAHits().at(index) == 0)                         return false; // Glb fit must have hits in mu chambers
        if (TMath::Abs(mud0PV_smurfV3(index)) > 0.2)                             return false; // d0(PV) < 0.2 cm
        if (TMath::Abs(mudzPV_smurfV3(index)) > 1  )                             return false; // dz(PV) < 1 cm
        if (cms2.mus_ptErr().at(index)/cms2.mus_p4().at(index).pt()>0.1)         return false; // dpt/pt < 0.1
        return true;

    case OSGeneric_v3_FO:
        // copy of OSGeneric_v2_FO ID, with different iso requirement
	    return muonIdNotIsolated( index, OSGeneric_v2_FO, -1 );

    case OSZ_v1:
        //this selector is Nominal + eta < 2.4 + dpt/pt < 0.1
        if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.4)                       return false; // eta cut
        if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) return false; //glb fit chisq
        if (((cms2.mus_type().at(index)) & (1<<1)) == 0)                         return false; // global muon
        if (((cms2.mus_type().at(index)) & (1<<2)) == 0)                         return false; // tracker muon
        if (cms2.mus_validHits().at(index) < 11)                                 return false; // # of tracker hits
        if (cms2.mus_iso_ecalvetoDep().at(index) > 4)                            return false; // ECalE < 4 
        if (cms2.mus_iso_hcalvetoDep().at(index) > 6)                            return false; // HCalE < 6 
        if (cms2.mus_gfit_validSTAHits().at(index) == 0)                         return false; // Glb fit must have hits in mu chambers
        if (TMath::Abs(cms2.mus_d0corr().at(index)) > 0.02)                      return false; // d0 from beamspot
        if (cms2.mus_ptErr().at(index)/cms2.mus_p4().at(index).pt()>0.1)         return false; // dpt/pt 
        return true;
        break;

    case OSZ_v2:
        // baseline selector for 2011 Z+MET analysis
        if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.4)                       return false; // eta cut
        if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) return false; // glb fit chisq
        if (((cms2.mus_type().at(index)) & (1<<1)) == 0)                         return false; // global muon
        if (((cms2.mus_type().at(index)) & (1<<2)) == 0)                         return false; // tracker muon
        if (cms2.mus_validHits().at(index) < 11)                                 return false; // # of tracker hits
        if (cms2.mus_iso_ecalvetoDep().at(index) > 4)                            return false; // ECalE < 4 
        if (cms2.mus_iso_hcalvetoDep().at(index) > 6)                            return false; // HCalE < 6 
        if (cms2.mus_gfit_validSTAHits().at(index) == 0)                         return false; // Glb fit must have hits in mu chambers
        if (TMath::Abs(cms2.mus_d0corr().at(index)) > 0.02)                      return false; // d0 from beamspot
        if (cms2.mus_ptErr().at(index)/cms2.mus_p4().at(index).pt()>0.1)         return false; // dpt/pt 
        if (cms2.mus_iso_ecalvetoDep().at(index) > 4)                            return false; // ECalE < 4 
        if (cms2.mus_iso_hcalvetoDep().at(index) > 6)                            return false; // HCalE < 6 
        if (!isPFMuon(index,true,1.0))                                           return false; // require muon is pfmuon with same pt
        return true;
        break;

    case NominalSSv2:
        if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.4)  return false; // eta cut
        if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) return false; //glb fit chisq
        if (((cms2.mus_type().at(index)) & (1<<1)) == 0)    return false; // global muon
        if (((cms2.mus_type().at(index)) & (1<<2)) == 0)    return false; // tracker muon
        if (cms2.mus_validHits().at(index) < 11)            return false; // # of tracker hits
        if (cms2.mus_iso_ecalvetoDep().at(index) > 4)       return false; // ECalE < 4 
        if (cms2.mus_iso_hcalvetoDep().at(index) > 6)       return false; // HCalE < 6 
        if (cms2.mus_gfit_validSTAHits().at(index) == 0)    return false; // Glb fit must have hits in mu chambers
        if (TMath::Abs(cms2.mus_d0corr().at(index)) > 0.02) return false; // d0 from beamspot
        if (cms2.mus_ptErr().at(index)/cms2.mus_p4().at(index).pt()>0.1) return false;
        return true;
        break;

    case muonSelectionFO_mu_ssV2:
    case muonSelectionFO_mu_ssV2_iso10:
        if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.4)  return false; // eta cut
        if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 50) return false; //glb fit chisq
        if (((cms2.mus_type().at(index)) & (1<<1)) == 0)    return false; // global muon
        if (((cms2.mus_type().at(index)) & (1<<2)) == 0)    return false; // tracker muon
        if (cms2.mus_validHits().at(index) < 11)            return false; // # of tracker hits
        if (TMath::Abs(cms2.mus_d0corr().at(index)) > 0.2) return false; // d0 from beamspot
        if (cms2.mus_ptErr().at(index)/cms2.mus_p4().at(index).pt()>0.1) return false;
        return true;
        break;
    case NominalTTbar_pass6:
        if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.4)  return false; // eta cut
        if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) return false; //glb fit chisq
        if (((cms2.mus_type().at(index)) & (1<<1)) == 0)    return false; // global muon
        if (((cms2.mus_type().at(index)) & (1<<2)) == 0)    return false; // tracker muon
        if (cms2.mus_validHits().at(index) < 11)            return false; // # of tracker hits
        if (cms2.mus_gfit_validSTAHits().at(index) == 0)    return false; // Glb fit must have hits in mu chambers
        if (TMath::Abs(cms2.mus_d0corr().at(index)) > 0.02) return false; // d0 from beamspot
        return true;
        break;
    case muonSelectionFO_mu_ttbar_pass6:
        if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.4)  return false; // eta cut
        if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 50) return false; //glb fit chisq
        if (((cms2.mus_type().at(index)) & (1<<1)) == 0)    return false; // global muon
        if (((cms2.mus_type().at(index)) & (1<<2)) == 0)    return false; // tracker muon
        if (cms2.mus_validHits().at(index) < 11)            return false; // # of tracker hits
        if (TMath::Abs(cms2.mus_d0corr().at(index)) > 0.2) return false; // d0 from beamspot
        return true;
        break;
    case NominalSmurfV3:
    case NominalSmurfV4:
    case NominalSmurfV5:
    case muonSelectionFO_mu_smurf_04:
    case muonSelectionFO_mu_smurf_10:
        if (type == NominalSmurfV3 || type == NominalSmurfV4 || type == NominalSmurfV5){
            if (cms2.mus_p4().at(index).pt()<20){
                if (TMath::Abs(mud0PV_smurfV3(index)) >= 0.01)    return false; // d0 from pvtx
            } else {
                if (TMath::Abs(mud0PV_smurfV3(index)) >= 0.02)    return false; // d0 from pvtx
            }
        } else {
            if (TMath::Abs(mud0PV_smurfV3(index)) >= 0.2)    return false; // d0 from pvtx
        }
        if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.4)  return false; // eta cut
        if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) return false; //glb fit chisq
        if (((cms2.mus_type().at(index)) & (1<<1)) == 0)    return false; // global muon
        if (((cms2.mus_type().at(index)) & (1<<2)) == 0)    return false; // tracker muon
        if (cms2.mus_validHits().at(index) < 11)            return false; // # of tracker hits  
        if (cms2.mus_gfit_validSTAHits().at(index)==0 )     return false; // Glb fit must have hits in mu chambers
        if (TMath::Abs(mudzPV_smurfV3(index)) >= 0.1)       return false; // dz from pvtx
        if (cms2.mus_ptErr().at(index)/cms2.mus_p4().at(index).pt()>0.1) return false;
        if (cms2.trks_valid_pixelhits().at(cms2.mus_trkidx().at(index))==0) return false;
        if (cms2.mus_nmatches().at(index)<2) return false;
        return true;
        break;
        //baseline selector for 2011 SS analysis
    case NominalSSv3:
        if (fabs(cms2.mus_p4().at(index).eta()) > 2.4)                           return false; // eta cut
        if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) return false; // glb fit chisq
        if (((cms2.mus_type().at(index)) & (1<<1)) == 0)                         return false; // global muon
        if (((cms2.mus_type().at(index)) & (1<<2)) == 0)                         return false; // tracker muon
        if (cms2.mus_validHits().at(index) < 11)                                 return false; // # of tracker hits
        if (cms2.mus_gfit_validSTAHits().at(index) == 0)                         return false; // Glb fit must have hits in mu chambers
        if (cms2.mus_ptErr().at(index)/cms2.mus_p4().at(index).pt() > 0.1)       return false; // dpt/pt < 0.1
        if (cms2.mus_iso_ecalvetoDep().at(index) > 4)                            return false; // ECalE < 4 
        if (cms2.mus_iso_hcalvetoDep().at(index) > 6)                            return false; // HCalE < 6 
        // now, cut on dz, d0 using supplied vertex index
        //if (vertex_index < 0) return false; // are our hypothesis leptons consistent with a vertex
        // note, the d0 computed here is w.r.t. a DA primary vertex
        // NOT a standard primary vertex
        // if not consistent vertex is found, computer d0 w.r.t to the beamSpot
        // if the muon doesn't have an associated track. get d0 w.r.t. beamSpot
        if (vertex_index < 0 || cms2.mus_trkidx().at(index) < 0) {
            if (fabs(cms2.mus_d0corr().at(index)) > 0.02)
                return false;
        }
        else {
            if (fabs(trks_d0_pv(cms2.mus_trkidx().at(index), vertex_index, true).first) > 0.02)
                return false;
        }
        return true;
        break;
        //baseline selector for 2011 SS analysis
    case muonSelectionFO_ssV3:
        if (fabs(cms2.mus_p4().at(index).eta()) > 2.4)                           return false; // eta cut
        if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 50) return false; // glb fit chisq
        if (((cms2.mus_type().at(index)) & (1<<1)) == 0)                         return false; // global muon
        if (((cms2.mus_type().at(index)) & (1<<2)) == 0)                         return false; // tracker muon
        if (cms2.mus_validHits().at(index) < 11)                                 return false; // # of tracker hits
        if (cms2.mus_gfit_validSTAHits().at(index) == 0)                         return false; // Glb fit must have hits in mu chambers
        if (cms2.mus_ptErr().at(index)/cms2.mus_p4().at(index).pt() > 0.1)       return false; // dpt/pt < 0.1
        // now, cut on dz, d0 using supplied vertex index
        //if (vertex_index < 0) return false; // are our hypothesis leptons consistent with a vertex
        // note, the d0 computed here is w.r.t. a DA primary vertex
        // NOT a standard primary vertex
        // if not consistent vertex is found, computer d0 w.r.t to the beamSpot
        // if the muon doesn't have an associated track. get d0 w.r.t. beamSpot
        if (vertex_index < 0 || cms2.mus_trkidx().at(index) < 0) {
            if (fabs(cms2.mus_d0corr().at(index)) > 0.2)
                return false;
        }
        else {
            if (fabs(trks_d0_pv(cms2.mus_trkidx().at(index), vertex_index, true).first) > 0.2)
                return false;
        }
        return true;
        break;
    case NominalSSv4:
        if (fabs(cms2.mus_p4().at(index).eta()) > 2.4)                           return false; // eta cut
        if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) return false; // glb fit chisq
        if (((cms2.mus_type().at(index)) & (1<<1)) == 0)                         return false; // global muon
        if (((cms2.mus_type().at(index)) & (1<<2)) == 0)                         return false; // tracker muon
        if (cms2.mus_validHits().at(index) < 11)                                 return false; // # of tracker hits
        if (cms2.mus_gfit_validSTAHits().at(index) == 0)                         return false; // Glb fit must have hits in mu chambers
        if (cms2.mus_ptErr().at(index)/cms2.mus_p4().at(index).pt() > 0.1)       return false; // dpt/pt < 0.1
        if (cms2.mus_iso_ecalvetoDep().at(index) > 4)                            return false; // ECalE < 4 
        if (cms2.mus_iso_hcalvetoDep().at(index) > 6)                            return false; // HCalE < 6
        // cut on d0, dz using first good DA vertex
        // if there isn't a good vertex, use the beamSpot
        if ((vertex_index < 0 && vtxidx < 0) || cms2.mus_trkidx().at(index) < 0) {
            if (fabs(cms2.mus_d0corr().at(index)) > 0.02)
                return false;
        }
        else if (vertex_index >= 0) {
            if (fabs(trks_d0_pv(cms2.mus_trkidx().at(index), vertex_index, true).first) > 0.02)
                return false;
        }
        else if (vtxidx >= 0) {
            if (fabs(trks_d0_pv(cms2.mus_trkidx().at(index), vtxidx, true).first) > 0.02)
                return false;            
        }
        return true;
        break;
        //baseline selector for 2011 SS analysis
    case muonSelectionFO_ssV4:
        if (fabs(cms2.mus_p4().at(index).eta()) > 2.4)                           return false; // eta cut
        if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 50) return false; // glb fit chisq
        if (((cms2.mus_type().at(index)) & (1<<1)) == 0)                         return false; // global muon
        if (((cms2.mus_type().at(index)) & (1<<2)) == 0)                         return false; // tracker muon
        if (cms2.mus_validHits().at(index) < 11)                                 return false; // # of tracker hits
        if (cms2.mus_gfit_validSTAHits().at(index) == 0)                         return false; // Glb fit must have hits in mu chambers
        if (cms2.mus_ptErr().at(index)/cms2.mus_p4().at(index).pt() > 0.1)       return false; // dpt/pt < 0.1
        // cut on d0, dz using first good DA vertex
        // if there isn't a good vertex, use the beamSpot
        if ((vertex_index < 0 && vtxidx < 0) || cms2.mus_trkidx().at(index) < 0) {
            if (fabs(cms2.mus_d0corr().at(index)) > 0.02)
                return false;
        }
        else if (vertex_index >= 0) {
            if (fabs(trks_d0_pv(cms2.mus_trkidx().at(index), vertex_index, true).first) > 0.02)
                return false;
        }
        else if (vtxidx >= 0) {
            if (fabs(trks_d0_pv(cms2.mus_trkidx().at(index), vtxidx, true).first) > 0.02)
                return false;            
        }
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
double muonIsoValueOriginal(unsigned int index){
    double sum =  cms2.mus_iso03_sumPt().at(index) + cms2.mus_iso03_emEt().at(index) + cms2.mus_iso03_hadEt().at(index);
    double pt  = cms2.mus_p4().at(index).pt();
    return sum/max(pt,20.);
}
double muonIsoValue(unsigned int index, bool truncated ){
    return ( muonIsoValue_TRK( index, truncated ) + muonIsoValue_ECAL( index, truncated ) + muonIsoValue_HCAL( index, truncated ) );
}
double muonIsoValue_FastJet(unsigned int index, bool truncated ){
    return ( muonIsoValue_TRK( index, truncated ) + TMath::Max( muonIsoValue_ECAL( index, truncated ) + muonIsoValue_HCAL( index, truncated ) - mu_fastjet_rel_offset(index,truncated) , 0.0 ) );
}
double mu_fastjet_rel_offset(unsigned int index, bool truncated ){
    double pt        = cms2.mus_p4().at(index).pt();
    if(truncated) pt = max( pt, 20.0 );
    double offset = TMath::Pi() * pow( 0.3 , 2 ) * cms2.evt_rho();
    return offset / pt;
}
double muonIsoValue_TRK(unsigned int index, bool truncated ){
    double pt        = cms2.mus_p4().at(index).pt();
    if(truncated) pt = max( pt, 20.0 );
    return cms2.mus_iso03_sumPt().at(index) / pt;
}
double muonIsoValue_ECAL(unsigned int index, bool truncated ){
    double pt  = cms2.mus_p4().at(index).pt();
    if(truncated) pt = max( pt, 20.0 );
    return cms2.mus_iso03_emEt().at(index) / pt;
}
double muonIsoValue_HCAL(unsigned int index, bool truncated){
    double pt  = cms2.mus_p4().at(index).pt();
    if(truncated) pt = max( pt, 20.0 );
    return cms2.mus_iso03_hadEt().at(index) / pt;
}

double muonIsoValuePF( unsigned int imu, unsigned int idavtx, float coner, float minptn, float dzcut){
    float pfciso = 0;
    float pfniso = 0;
    int mutkid = cms2.mus_trkidx().at(imu);
    float mudz = mutkid>=0 ? trks_dz_dapv(mutkid,idavtx).first : cms2.mus_sta_z0corr().at(imu);
    for (unsigned int ipf=0; ipf<cms2.pfcands_p4().size(); ++ipf){
        float dR = dRbetweenVectors(cms2.pfcands_p4().at(ipf),cms2.mus_p4().at(imu));
        if (dR>coner) continue;
        float pfpt = cms2.pfcands_p4().at(ipf).pt();
        if (cms2.pfcands_charge().at(ipf)==0) {
            //neutrals
            if (pfpt>minptn) pfniso+=pfpt;
        } else {
            //charged
            //avoid double counting of muon itself
            int pftkid = cms2.pfcands_trkidx().at(ipf);
            if (mutkid>=0 && pftkid>=0 && mutkid==pftkid) continue;
            //first check electrons with gsf track
            if (abs(cms2.pfcands_particleId().at(ipf))==11 && cms2.pfcands_pfelsidx().at(ipf)>=0 && cms2.pfels_elsidx().at(cms2.pfcands_pfelsidx().at(ipf))>=0) {
                int gsfid = cms2.els_gsftrkidx().at(cms2.pfels_elsidx().at(cms2.pfcands_pfelsidx().at(ipf))); 
                if (gsfid>=0) { 
                    if(fabs(gsftrks_dz_dapv( gsfid,idavtx ).first - mudz )<dzcut) {//dz cut
                        pfciso+=pfpt;
                    }   
                    continue;//and avoid double counting
                }
            }
            //then check anything that has a ctf track
            if (cms2.pfcands_trkidx().at(ipf)>=0) {//charged (with a ctf track)
                if(fabs( trks_dz_dapv(cms2.pfcands_trkidx().at(ipf),idavtx).first - mudz )<dzcut) {//dz cut
                    pfciso+=pfpt;
                }
            } 
        }
    } 
    return (pfciso+pfniso)/cms2.mus_p4().at(imu).pt();
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

//------------------------------------------------------------------------------
// Muon d0 corrected by PV
//-------------------------------------------------------------------------------

double mud0PV(unsigned int index){
    if ( cms2.vtxs_sumpt().empty() ) return 9999.;
    unsigned int iMax = 0;
    double sumPtMax = cms2.vtxs_sumpt().at(0);
    for ( unsigned int i = iMax+1; i < cms2.vtxs_sumpt().size(); ++i )
        if ( cms2.vtxs_sumpt().at(i) > sumPtMax ){
            iMax = i;
            sumPtMax = cms2.vtxs_sumpt().at(i);
        }
    double dxyPV = cms2.mus_d0()[index]-
        cms2.vtxs_position()[iMax].x()*sin(cms2.mus_trk_p4()[index].phi())+
        cms2.vtxs_position()[iMax].y()*cos(cms2.mus_trk_p4()[index].phi());
    return dxyPV;
}


double mud0PV_wwV1(unsigned int index){
    if ( cms2.vtxs_sumpt().empty() ) return 9999.;
    double sumPtMax = -1;
    int iMax = -1;
    for ( unsigned int i = 0; i < cms2.vtxs_sumpt().size(); ++i ){
        // if (!isGoodVertex(i)) continue;
        // Copied from eventSelections.cc 
        if (cms2.vtxs_isFake()[i]) continue;
        if (cms2.vtxs_ndof()[i] < 4.) continue;
        if (cms2.vtxs_position()[i].Rho() > 2.0) continue;
        if (fabs(cms2.vtxs_position()[i].Z()) > 24.0) continue;
        if ( cms2.vtxs_sumpt().at(i) > sumPtMax ){
            iMax = i;
            sumPtMax = cms2.vtxs_sumpt().at(i);
        }
    }
    if (iMax<0) return 9999.;
    double dxyPV = cms2.mus_d0()[index]-
        cms2.vtxs_position()[iMax].x()*sin(cms2.mus_trk_p4()[index].phi())+
        cms2.vtxs_position()[iMax].y()*cos(cms2.mus_trk_p4()[index].phi());
    return dxyPV;
}

double mud0PV_smurfV3(unsigned int index){
    int vtxIndex = 0;
    double dxyPV = cms2.mus_d0()[index]-
        cms2.davtxs_position()[vtxIndex].x()*sin(cms2.mus_trk_p4()[index].phi())+
        cms2.davtxs_position()[vtxIndex].y()*cos(cms2.mus_trk_p4()[index].phi());
    return dxyPV;
}

double dzPV_mu(const LorentzVector& vtx, const LorentzVector& p4, const LorentzVector& pv){
    return (vtx.z()-pv.z()) - ((vtx.x()-pv.x())*p4.x()+(vtx.y()-pv.y())*p4.y())/p4.pt() * p4.z()/p4.pt();
}

double mudzPV_smurfV3(unsigned int index){
    int vtxIndex = 0;
    double dzpv = dzPV_mu(cms2.mus_vertex_p4()[index], cms2.mus_trk_p4()[index], cms2.davtxs_position()[vtxIndex]);
    return dzpv;
}

double mudzPV_wwV1(unsigned int index){
    if ( cms2.vtxs_sumpt().empty() ) return 9999.;
    double sumPtMax = -1;
    int iMax = -1;
    for ( unsigned int i = 0; i < cms2.vtxs_sumpt().size(); ++i ){
        // if (!isGoodVertex(i)) continue;
        // Copied from eventSelections.cc 
        if (cms2.vtxs_isFake()[i]) continue;
        if (cms2.vtxs_ndof()[i] < 4.) continue;
        if (cms2.vtxs_position()[i].Rho() > 2.0) continue;
        if (fabs(cms2.vtxs_position()[i].Z()) > 24.0) continue;
        if ( cms2.vtxs_sumpt().at(i) > sumPtMax ){
            iMax = i;
            sumPtMax = cms2.vtxs_sumpt().at(i);
        }
    }
    if (iMax<0) return 9999.;
    // double dzpv = cms2.mus_z0corr()[index]-cms2.vtxs_position()[iMax].z();
    const LorentzVector& vtx = cms2.mus_vertex_p4()[index];
    const LorentzVector& p4 = cms2.mus_trk_p4()[index];
    const LorentzVector& pv = cms2.vtxs_position()[iMax];
    return (vtx.z()-pv.z()) - ((vtx.x()-pv.x())*p4.x()+(vtx.y()-pv.y())*p4.y())/p4.pt() * p4.z()/p4.pt(); 
    /* directly from NtupleMacros/WW/doAnalysis.cc
       double dzpv = dzPV(cms2.mus_vertex_p4()[index], cms2.mus_trk_p4()[index], cms2.vtxs_position()[iMax]);
       double dzPV(const LorentzVector& vtx, const LorentzVector& p4, const LorentzVector& pv){
       return (vtx.z()-pv.z()) - ((vtx.x()-pv.x())*p4.x()+(vtx.y()-pv.y())*p4.y())/p4.pt() * p4.z()/p4.pt();
       }*/
}


bool isPFMuon( int index , bool requireSamePt , float dpt_max ){

    int ipf = cms2.mus_pfmusidx().at( index );

    //--------------------------
    // require matched pfmuon
    //--------------------------

    if( ipf >= int(cms2.pfmus_p4().size()) || ipf < 0 ) return false;

    //----------------------------------------------------
    // require PFMuon pt = reco muon pt (within dpt_max)
    //----------------------------------------------------

    if( requireSamePt ){

        float pt_pf = cms2.pfmus_p4().at(ipf).pt();
        float pt    = cms2.mus_p4().at(index).pt();

        if( fabs( pt_pf - pt ) > dpt_max ) return false;

    }

    return true;

}

