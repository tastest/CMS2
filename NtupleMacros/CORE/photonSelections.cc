#include "photonSelections.h"
#include "jetSelections.cc"
#include "utilities.cc"
#include <iostream>
#include "TSystem.h"
#include "Math/VectorUtil.h"
#include "CMS2.h"

//-----------
// Photon ID
//-----------
bool photonId( const unsigned int iPhoton, PhotonSelectionType type /* default defined in header */ ){

  float dR_outer = 0.4;
  float dR_inner = 0.05;
  float sumHollow = 0;

  switch (type) {

    case Yuri:

      if( fabs( cms2.photons_p4().at(iPhoton).pt() ) < 10 ) return false;                                               // Pt
      if( fabs( cms2.photons_p4().at(iPhoton).eta() ) > 1.479 ) return false;                                           // Eta
      if( cms2.photons_ecalIso03().at(iPhoton) >= ( 4.2 + .004*cms2.photons_p4().at(iPhoton).pt() ) ) return false;       // ECAL Isolation
      if( cms2.photons_hcalIso03().at(iPhoton) >= ( 2.2 + .001*cms2.photons_p4().at(iPhoton).pt() ) ) return false;       // HCAL Isolation
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
      if( sumHollow  >= ( 2.0 + 0.001*cms2.photons_p4().at(iPhoton).pt() ) ) return false;                              // Hollow

      if( isSpikePhoton(iPhoton) ) return false;                                                                        // Spike Removal

      return true;

      break;

    default:
      std::cout << "photonID ERROR: requested photon type is not defined. Abort." << std::endl;
      gSystem->Exit(1);
      return false;

  }
}

//----------------------------
// Spike rejection for photons
// ( following electrons )
//----------------------------
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

//-----------------------------------------------------------------
//function to select a good EM object for met templates analysis
//function returns -1 if the photon fails the selection
//otherwise, function returns index of pfjet matched to EM object
//this pfjet must be excluded from the njets, sumJetPt summation
//USAGE:
// int pfjet_index = isGoodEMObject( photon_index );
// if( pfjet_index < 0 ) continue;
//-----------------------------------------------------------------
 

int isGoodEMObject( const unsigned int index ){

  //Minimum cut for making into babies. Tighter cut (0.95 was previous operating point) will be applied at template making.
  const float neutralemfcut = 0.7;

  //apply this cut at template creation time
  //if ( fabs( photons_p4().at(index).eta() ) > 1   )     return -1; //eta < 1
  if ( photons_p4().at(index).pt() < 22           )     return -1; //pt > 22 GeV
  if ( photons_hOverE().at(index) > 0.1           )     return -1; //h/e < 0.1

  //spike cleaning
  if( isSpikePhoton( index ) )                          return -1;

  //if photon survives to this point, find pfjet nearest photon
  //require pt > 10 GeV pfjet, eta < 2.5 within dr < 0.3 of photon

  float drmin       = 100;
  int   iMatchedJet = -1;

  for (unsigned int ijet = 0 ; ijet < pfjets_p4().size() ; ijet++) {
          
    LorentzVector vjet = pfjets_cor().at(ijet) * pfjets_p4().at(ijet);
    LorentzVector vg   = photons_p4().at(index);
          
    if( vjet.pt()  < 10  )       continue;
    if( fabs(vjet.eta()) > 2.5 )       continue;
    
    float dr = dRbetweenVectors(vjet, vg);
    
    if( dr < drmin ){
      drmin       = dr;
      iMatchedJet = ijet;
    }
  }

  if( iMatchedJet < 0 ) return -2; //change -1 to -2, etc, so i can keep track in my looper of the cut which fails
  if( drmin > 0.3 )     return -3;

  //require pfjet neutral EM fraction > 0.95
  float emfrac = pfjets_neutralEmE().at(iMatchedJet) / pfjets_p4().at(iMatchedJet).energy();
  if( emfrac < neutralemfcut )               return -4; 
  //too tight
  //if( !passesPFJetID(iMatchedJet) )                 return -5;

  return iMatchedJet;

}
