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
	 case muonSelectionFO_mu_wwV1_iso10:
	 case muonSelectionFO_mu_ww_iso10:
		  isovalue = 1.0;
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
	 case OSZ_v1:
		  isovalue = 0.15;
		  break;
	 default:
		  std::cout << "muonID ERROR: requested muon type is not defined. Abort." << std::endl;
		  exit(1);
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
		  if ( cms2.mus_p4()[index].pt() < 20.0 )             return false; // pt cut
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
		  if ( cms2.mus_p4()[index].pt() < 20.0 )             return false; // pt cut
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

	 case muonSelectionFO_mu_ww:
		  if ( cms2.mus_p4()[index].pt() < 20.0 )             return false; // pt cut
		  if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.4)  return false; // eta cut
		  if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) return false; //glb fit chisq
		  if (((cms2.mus_type().at(index)) & (1<<1)) == 0)    return false; // global muon
		  if (((cms2.mus_type().at(index)) & (1<<2)) == 0)    return false; // tracker muon
		  if (cms2.mus_validHits().at(index) < 11)            return false; // # of tracker hits  
		  if (cms2.mus_gfit_validSTAHits().at(index)==0 )     return false; // Glb fit must have hits in mu chambers
		  if (TMath::Abs(mud0PV(index)) >= 0.02)              return false; // d0 from pvtx
		  return true;
	    
	 case muonSelectionFO_mu_ww_iso10:
		  if ( cms2.mus_p4()[index].pt() < 20.0 )             return false; // pt cut
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
