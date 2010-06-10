#include "CMS2.h"
#include "photonSelections.h"
#include <iostream>
#include "TSystem.h"

//-----------
// Photon ID
//-----------
bool photonId( const unsigned int iPhoton, PhotonSelectionType type /* default defined in header */ ){

//  float drMin    = 0.2;
//  bool  isPhoton = true;
//  LorentzVector photon_sc_p4;
//  LorentzVector els_sc_p4;
//  float els_scet;
//  float els_sce;
//  int   els_typ;

  float dR_outer = 0.4;
  float dR_inner = 0.05;
  float sumHollow = 0;

  switch (type) {

    case Yuri:

      if( fabs( cms2.photons_p4().at(iPhoton).pt() ) < 10 ) return false;                                               // Pt


      if( fabs( cms2.photons_p4().at(iPhoton).eta() ) > 1.479 ) return false;                                           // Eta
      if( cms2.photons_ecalIso().at(iPhoton) >= ( 4.2 + .004*cms2.photons_p4().at(iPhoton).pt() ) ) return false;       // ECAL Isolation
      if( cms2.photons_hcalIso().at(iPhoton) >= ( 2.2 + .001*cms2.photons_p4().at(iPhoton).pt() ) ) return false;       // HCAL Isolation
      if( cms2.photons_hOverE().at(iPhoton)  >= 0.05 ) return false;                                                    // H over E
      if( cms2.photons_sigmaIEtaIEta().at(iPhoton) >= 0.013 ) return false;                                             // Eta width
      //if( cms2.photons_tkIsoHollow().at(iPhoton)  >= ( 2 + 0.001*cms2.photons_p4().at(iPhoton).pt() ) ) return false;   // Hollow

      // Hollow Cone Isolation
      for( unsigned int iTrack=0; iTrack < cms2.trks_trk_p4().size(); iTrack++ ){
        float dR = ROOT::Math::VectorUtil::DeltaR( cms2.photons_p4().at(iPhoton), cms2.trks_trk_p4().at(iTrack) );
        if( dR < dR_outer && dR > dR_inner ){
          sumHollow += cms2.trks_trk_p4().at(iTrack).pt();         
        }
      }
      if( sumHollow  >= ( 2.0 + 0.001*cms2.photons_p4().at(iPhoton).pt() ) ) return false;   // Hollow


//      //
//      photon_sc_p4 = cms2.scs_pos_p4().at( cms2.photons_scindex().at(iPhoton) );
//      for( unsigned int iEl = 0; iEl < cms2.els_p4().size(); iEl++ ){
//        els_sce = cms2.els_eSC().at(iEl);
//        els_scet = (cms2.els_eSC().at(iEl) / cosh(cms2.els_etaSC().at(iEl)) );
//        els_typ = cms2.els_type().at(iEl) & (1<<2);
//	      if( cms2.els_scindex().at(iEl) < 0 ){
//	        cout << "els_scindex() is negative for (els_sce, els_scet, els_pt, els_Et, els_type & (1<<2) ): (" << els_sce << ", " 
//               << els_scet << ", " << cms2.els_p4().at(iEl).pt() << ", " << cms2.els_p4().at(iEl).Et() << ", " << els_typ << ")" << endl;
//	      }
//        //els_sc_p4 = cms2.scs_pos_p4().at( cms2.els_scindex().at(iEl) );
//        if( ROOT::Math::VectorUtil::DeltaR( photon_sc_p4, els_sc_p4 ) < drMin ) isPhoton = false;
//      }
//      //if( isPhoton == false ) return false;



      if( isSpikePhoton(iPhoton) ) return false;                                                                        // Spike Removal

      return true;

      break;

    default:
      std::cout << "photonID ERROR: requested photon type is not defined. Abort." << std::endl;
      gSystem->Exit(1);
      return false;

  }
}

bool isSpikePhoton( const unsigned int index ) {

    const int scidx = cms2.photons_scindex()[index];
    bool isSpike = false;
    if (scidx != -1) {
        //subtract twice max since max is in both 1x3 and 3x1, and we want neither
        const float r4 = (cms2.scs_e1x3()[scidx] + cms2.scs_e3x1()[scidx] - 2*cms2.scs_eMax()[scidx])/cms2.scs_eMax()[scidx];
        if (r4 < 0.05) isSpike = true;
    }

    return isSpike;

}

