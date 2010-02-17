#include <algorithm>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>

#include "TChain.h"
#include "TChainElement.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TH2F.h"
#include "TMath.h"
#include "Math/LorentzVector.h"
#include "TRandom3.h"

#include "getTcmetFromCaloMet.h"
#include "CMS2_3x.h"
//#include "getResponseFunction_fit.C"

using namespace std;
using namespace tas;

enum TrackQuality { undefQuality=-1, loose=0, tight=1, highPurity=2, confirmed=3, goodIterative=4, qualitySize=5};

//--------------------------------------------------------------------

bool isTrackQuality( int index, int cuts ) {
  return ( ( cms2.trks_qualityMask().at(index) & cuts ) == cuts );
}

//--------------------------------------------------------------------

metStruct getTcmetFromCaloMet(TH2F* rf){

  //rf = getResponseFunction_fit();

  //initialize to calomet values
  tcmet_x = evt_met() * cos( evt_metPhi() );
  tcmet_y = evt_met() * sin( evt_metPhi() );
  tcsumet = evt_sumet();
  
  //correct calomet for muons
  for( unsigned int i = 0; i < mus_p4().size(); i++ ) {
    
    int flag = mus_tcmet_flag().at(i);
    
    if( flag == 0 ) continue;
    
    else if( flag == 1 ){ 
      correctMETforMuon  ( mus_gfit_p4().at(i), i );
      correctSumEtForMuon( mus_gfit_p4().at(i), i );
    }
    else if( flag == 2 ){
      correctMETforMuon  ( mus_trk_p4().at(i), i );
      correctSumEtForMuon( mus_trk_p4().at(i), i );
    }
    else if( flag == 3 || flag==5 ){
      correctMETforMuon  ( mus_p4().at(i), i );
      correctSumEtForMuon( mus_p4().at(i), i );
    }
    else if( flag == 4 ){
      correctMETforPion  ( mus_trkidx().at(i) , rf);
      correctSumEtForPion( mus_trkidx().at(i) , rf);
    }
    
    else
      std::cout << "Muon has an undefined flag!\n";
  }
  
  //correct muon-corrected calomet for pions
  for( unsigned int i = 0; i < trks_trk_p4().size(); i++ ) {
    
    if( isMuon( i ) ) 
      continue;
    
    if( isElectron( i ) ) 
      continue;
    
    if( !isGoodTrack( i ) ) 
      continue;
    
    correctMETforPion( i , rf);
    correctSumEtForPion( i , rf);
  }
  
  float tcmet    = TMath::Sqrt( tcmet_x * tcmet_x + tcmet_y * tcmet_y );
  float tcmetPhi = atan2(tcmet_y,tcmet_x); 
  
  metStruct tcmetStruct;
  tcmetStruct.met    = tcmet;
  tcmetStruct.metphi = tcmetPhi;
  tcmetStruct.metx   = tcmet_x;
  tcmetStruct.mety   = tcmet_y;
  tcmetStruct.sumet  = tcsumet;
  
  return tcmetStruct;
}

//--------------------------------------------------------------------

bool isGoodTrack( int index ) {
  
  if( trks_trk_p4().at(index).pt() > 100 )                           return false;
  if( fabs( trks_trk_p4().at(index).eta() ) > 2.65 )                 return false;
  if( fabs( trks_d0corr().at(index) ) > 0.1 )                        return false;
  if( trks_validHits().at(index) < 6 )                               return false;
  if( trks_chi2().at(index) / trks_ndof().at(index) > 5 )            return false;
  if( trks_ptErr().at(index) / trks_trk_p4().at(index).pt() > 0.20 ) return false;
  if( !isTrackQuality( index, (1 << highPurity) ) )                  return false;
  if( fabs( trks_outer_p4().at(index).eta() ) > 5 )                  return false;  
  if( fabs( trks_outer_p4().at(index).phi() ) > 2 * TMath::Pi() )    return false; 
  
  return true;
}

//--------------------------------------------------------------------

void correctMETforMuon( LorentzVector p4, int index ) {

  float deltax = mus_met_deltax().at(index);
  float deltay = mus_met_deltay().at(index);

  tcmet_x -= ( p4.px() - deltax );
  tcmet_y -= ( p4.py() - deltay );
}

//--------------------------------------------------------------------

void correctMETforPion( int index , TH2F* rf) {

  float deltax = 0;
  float deltay = 0;
  
  if( trks_trk_p4().at(index).pt() > 1 ) {
    
    int bin =        rf->FindBin(trks_trk_p4().at(index).eta(), trks_trk_p4().at(index).pt() );
    float response = rf->GetBinContent( bin );
    
    deltax = response * trks_trk_p4().at(index).P() * 
      sin( trks_outer_p4().at(index).Theta() ) * cos( trks_outer_p4().at(index).phi() );
    deltay = response * trks_trk_p4().at(index).P() * 
      sin( trks_outer_p4().at(index).Theta() ) * sin( trks_outer_p4().at(index).phi() );
  }
  
  tcmet_x -= ( trks_trk_p4().at(index).px() - deltax );
  tcmet_y -= ( trks_trk_p4().at(index).py() - deltay );
}

//--------------------------------------------------------------------

void correctSumEtForMuon( LorentzVector p4, int index ){
  
  float deltax = mus_met_deltax().at(index);
  float deltay = mus_met_deltay().at(index);

  tcsumet += ( p4.pt() - TMath::Sqrt( deltax * deltax + deltay * deltay ) );
}

//--------------------------------------------------------------------

void correctSumEtForPion( int index , TH2F* rf) {
  
  if( trks_trk_p4().at(index).pt() < 1) {
    tcsumet += trks_trk_p4().at(index).pt();
  }
  
  else {
    int bin = rf->FindBin( trks_trk_p4().at(index).eta() , trks_trk_p4().at(index).pt()); 
    double fracTrackEnergy = rf->GetBinContent( bin); 
    tcsumet += ( 1 - fracTrackEnergy ) * trks_trk_p4().at(index).pt();
  }
}

//--------------------------------------------------------------------

bool isMuon( int index ) {

  for( unsigned int i = 0; i < mus_p4().size(); i++ ) {

    if( mus_trkidx().at(i) == index ) return true;
  }

  return false;
}

//--------------------------------------------------------------------

bool isElectron( int index ) {

  for( unsigned int i = 0; i < els_p4().size(); i++ ) {

    if( els_trkidx().at(i) == index && els_hOverE().at(i) < 0.1 ) return true;
  }

  return false;
}

//--------------------------------------------------------------------
