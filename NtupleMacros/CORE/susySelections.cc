#include <assert.h>
#include <algorithm>
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TDatabasePDG.h"
#include "triggerUtils.h"
#include "CMS2.h"

#include "susySelections.h"
#include "triggerUtils.cc"
#include "electronSelections.h"
#include "muonSelections.h"

using namespace tas;

/*****************************************************************************************/
//print event info
/*****************************************************************************************/
void printEventInfo(){
  cout << cms2.evt_dataset() << endl;
  cout << cms2.evt_run() << " " << cms2.evt_lumiBlock() << " " << cms2.evt_event() << endl;
}


/*****************************************************************************************/
//veto Z->mumugamma events                                                                
/*****************************************************************************************/
bool vetoZmumuGamma( unsigned int hypIdx , float emax , float minmass , float maxmass ){

  //we only care about mumu hyp types
  if( cms2.hyp_type().at(hypIdx) != 0 ) return false; 

  //we only want to veto events that started *outside* of Z window 
  //and moved inside after including ecal deposit
  if( cms2.hyp_p4().at(hypIdx).mass() > minmass && cms2.hyp_p4().at(hypIdx).mass() < maxmass )
    return false;

  //get ecal deposits for each muon
  int ill = cms2.hyp_ll_index().at( hypIdx );
  int ilt = cms2.hyp_lt_index().at( hypIdx );

  //energy = ET * cosh(eta)
  float ell = cms2.mus_iso_ecalvetoDep().at(ill) * cosh( cms2.mus_p4().at(ill).eta() ) ;
  float elt = cms2.mus_iso_ecalvetoDep().at(ilt) * cosh( cms2.mus_p4().at(ilt).eta() ) ;

  //don't veto event if ecal deposits are both less than emax
  if( ell < emax && elt < emax ) return false;

  LorentzVector vll = cms2.hyp_ll_p4().at(hypIdx);
  LorentzVector vlt = cms2.hyp_lt_p4().at(hypIdx);

  //if ecal deposit is greater than emax, scale muon momentum up
  if( ell > emax ) vll = vll * ( 1 + ell / vll.energy() );
  if( elt > emax ) vlt = vlt * ( 1 + elt / vlt.energy() );

  float mass = ( vll + vlt ).mass();

  //check if dimuon mass, including extra ecal energy, is inside Z mass window
  if( mass > minmass && mass < maxmass ){
    //cout << "Veto Z->mumugamma event!" << endl;
    //cout << cms2.evt_run() << " " << cms2.evt_lumiBlock() << " " << cms2.evt_event() << endl;
    return true;
  }

  return false;
}

/*****************************************************************************************/
//generalized Z veto
/*****************************************************************************************/
bool ZVetoGeneral( float ptcut , float minmass ,  float maxmass , SelectionType mutype ){

  for(unsigned int i = 0; i < hyp_p4().size(); ++i) {

    //check that hyp leptons come from same vertex
    if( !hypsFromSameVtx( i ) )    continue;

    //opposite-sign, same-flavor
    if( cms2.hyp_lt_id().at(i)     *  cms2.hyp_ll_id().at(i) > 0 ) continue;
    if( cms2.hyp_type().at(i) == 1 || cms2.hyp_type().at(i) == 2 ) continue;

    //check that both lepton pt's > ptcut
    if( cms2.hyp_ll_p4().at(i).pt() < ptcut ) continue; 
    if( cms2.hyp_lt_p4().at(i).pt() < ptcut ) continue; 
          
    //muon ID
    if ( abs( cms2.hyp_ll_id().at(i) ) == 13  && !( muonId( cms2.hyp_ll_index().at(i) , mutype ) ) )   continue;
    if ( abs( cms2.hyp_lt_id().at(i) ) == 13  && !( muonId( cms2.hyp_lt_index().at(i) , mutype ) ) )   continue;
          
    //electron ID
    if ( abs( cms2.hyp_ll_id().at(i) ) == 11  && !( pass_electronSelection( cms2.hyp_ll_index().at(i) , electronSelection_el_OSV1 , false , false ))) continue;
    if ( abs( cms2.hyp_lt_id().at(i) ) == 11  && !( pass_electronSelection( cms2.hyp_lt_index().at(i) , electronSelection_el_OSV1 , false , false ))) continue;
          
    if( cms2.hyp_p4().at(i).mass() > minmass && cms2.hyp_p4().at(i).mass() < maxmass ){
      //cout << "General Z veto: mass " << cms2.hyp_p4().at(i).mass() << endl;
      //printEventInfo();
      return true; 
    }
  }

  return false;
}



/*****************************************************************************************/
//passes the OS SUSY trigger selection
/*****************************************************************************************/
bool passSUSYTrigger_v1( bool isData , int hypType ) {

  //int run = cms2.evt_run();

  //currently do not require any triggers for MC
  if( !isData ) return true;

  //mumu
  if( hypType == 0 ){

    // This is overkill, as Mu15 should be a susbset of Mu11 (for example)
    // But why not.
  
    if( passUnprescaledHLTTrigger("HLT_DoubleMu3_v2") ) return true;
    if( passUnprescaledHLTTrigger("HLT_DoubleMu5_v1") )   return true;
    if( passUnprescaledHLTTrigger("HLT_Mu11") )           return true;
    if( passUnprescaledHLTTrigger("HLT_Mu13_v1") )        return true;
    if( passUnprescaledHLTTrigger("HLT_Mu15_v1") )        return true;
    if( passUnprescaledHLTTrigger("HLT_Mu17_v1") )        return true; //<<<---Added 2e32
    if( passUnprescaledHLTTrigger("HLT_Mu19_v1") )        return true; //<<<---Added 2e32
    
    //if( run <= 147116 ){
    if( passUnprescaledHLTTrigger("HLT_DoubleMu3") )    return true; //136033-147116
    if( passUnprescaledHLTTrigger("HLT_Mu9") )          return true; //136033-147116
    //}
    //if( run <= 144114 ){
    if( passUnprescaledHLTTrigger("HLT_Mu7") )          return true; //140116-144114
    //}
    //if( run <= 141882 ){
    if( passUnprescaledHLTTrigger("HLT_Mu5") )          return true; //136033-141882
    //}
    
  }
 


  //ee
  else if( hypType == 3 ){

    //Added for 2e32 single ele trigs:
    if( passUnprescaledHLTTrigger("HLT_Ele17_SW_TighterEleIdIsol_L1R_v1")) return true;
    if( passUnprescaledHLTTrigger("HLT_Ele17_SW_TighterEleIdIsol_L1R_v2")) return true;
    if( passUnprescaledHLTTrigger("HLT_Ele17_SW_TighterEleIdIsol_L1R_v3")) return true;

    if( passUnprescaledHLTTrigger("HLT_Ele22_SW_TighterEleId_L1R_v2"))     return true;
    if( passUnprescaledHLTTrigger("HLT_Ele22_SW_TighterEleId_L1R_v3")) return true;
    if( passUnprescaledHLTTrigger("HLT_Ele22_SW_TighterCaloIdIsol_L1R_v2")) return true;

    if( passUnprescaledHLTTrigger("HLT_Ele27_SW_TightCaloEleIdTrack_L1R_v1")) return true;

    if( passUnprescaledHLTTrigger("HLT_Ele32_SW_TightCaloEleIdTrack_L1R_v1")) return true;
    if( passUnprescaledHLTTrigger("HLT_Ele32_SW_TighterEleId_L1R_v2")) return true;

   // This is a family of never prescaled Ele17 triggers
    // Some of them are probably prescaled at 2e32?
    //got prescaled at 148058...
    if( passUnprescaledHLTTrigger("HLT_Ele17_SW_TightEleId_L1R") )          return true; 
    if( passUnprescaledHLTTrigger("HLT_Ele17_SW_TighterEleId_L1R_v1") )     return true;
    
    // These are unprescaled double triggers
    if( passUnprescaledHLTTrigger("HLT_DoubleEle15_SW_L1R_v1") )                 return true;
    if( passUnprescaledHLTTrigger("HLT_DoubleEle17_SW_L1R_v1") )                 return true; //    
           
    if( passUnprescaledHLTTrigger("HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v1") ) return true; // 147390-->
    if( passUnprescaledHLTTrigger("HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v2") ) return true; // 147390-->
    if( passUnprescaledHLTTrigger("HLT_Ele17_SW_TightCaloEleId_SC8HE_L1R_v1") )  return true; // 147196-148058 

    // These are double triggers that became prescaled at some point
    //if( run <= 147116 ){
    if( passUnprescaledHLTTrigger("HLT_DoubleEle10_SW_L1R") ) return true; //141956-147116
    //}
    //if( run <= 141882 ){
    if( passUnprescaledHLTTrigger("HLT_DoubleEle5_SW_L1R") )  return true; //136033-141882
    //}
      
      
    //if( run <= 147116 ){
    // This is a family of Ele17 triggers that are not
    // there anymore or are not unprescaled anymore.
    // The loosest one is indicated as master
    // In principle, the other two could be skipped
    if( passUnprescaledHLTTrigger("HLT_Ele17_SW_CaloEleId_L1R") )                 return true; //146428-147116 <---- master
    if( passUnprescaledHLTTrigger("HLT_Ele17_SW_EleId_L1R") )                     return true; //146428-147116
    if( passUnprescaledHLTTrigger("HLT_Ele17_SW_LooseEleId_L1R") )                return true; //146428-147116
    //}

    //if( run <= 144114 ){
    // This is a family of Ele15 triggers that are not
    // there anymore or are not unprescaled anymore.
    // The loosest one is indicated as master
    // In principle, the other one could be skipped
    if( passUnprescaledHLTTrigger("HLT_Ele15_SW_CaloEleId_L1R") ) return true; //141956-144114 <---- master
    if( passUnprescaledHLTTrigger("HLT_Ele15_SW_EleId_L1R") )     return true; //141956-144114
    //}
    
    // These are Ele15 that got prescaled in a different
    // time frame as the previous family...use both
    //if( run <= 143962 ){
    if( passUnprescaledHLTTrigger("HLT_Ele15_SW_L1R") ) return true; //140058-143962
    //}
    //if( run <= 141882 ){
    if( passUnprescaledHLTTrigger("HLT_Ele15_LW_L1R") ) return true; //136033-141882
    //}
    
    // This is a bonus trigger, should not be needed, but why not?
    //if( run <= 144114 ){
    if( passUnprescaledHLTTrigger("HLT_Ele20_SW_L1R") ) return true; //140058-144114
    //}

    // This is the Ele10 family that came and went at ramdom times
    // Use all (although in principle the one with EleID could be
    // skipped in the subset of runs where triggers with no EleId
    // were present)
    //if( run <= 144114 ){
    if( passUnprescaledHLTTrigger("HLT_Ele10_SW_EleId_L1R") ) return true;      //141956-144114
    //}
    //if( run <= 141882 ){
    if( passUnprescaledHLTTrigger("HLT_Ele10_LW_EleId_L1R") ) return true;      //136033-141882
    //} 
    //if( run <= 139980 ){
    if( passUnprescaledHLTTrigger("HLT_Ele10_LW_L1R") )       return true;      //136033-139980
    if( passUnprescaledHLTTrigger("HLT_Ele10_SW_L1R") )       return true;      //139195-139980
    //}

  }

  //emu
  else if( hypType == 1 || hypType == 2 ){
    
    //---------------------------
    // single muon triggers
    //---------------------------
    if( passUnprescaledHLTTrigger("HLT_Mu11") )           return true;
    if( passUnprescaledHLTTrigger("HLT_Mu13_v1") )        return true;
    if( passUnprescaledHLTTrigger("HLT_Mu15_v1") )        return true;
    
    //if( run <= 147116 ){
    if( passUnprescaledHLTTrigger("HLT_Mu9") )          return true; //136033-147116
    //}
    //if( run <= 144114 ){
    if( passUnprescaledHLTTrigger("HLT_Mu7") )          return true; //140116-144114
    //}
    //if( run <= 141882 ){
    if( passUnprescaledHLTTrigger("HLT_Mu5") )          return true; //136033-141882
    //}

    if( passUnprescaledHLTTrigger("HLT_Mu17_v1") )        return true; //<<<---Added 2e32
    if( passUnprescaledHLTTrigger("HLT_Mu19_v1") )        return true; //<<<---Added 2e32
 
    
    //---------------------------
    // single electron triggers
    //---------------------------
    
    //Added 2e32
    if( passUnprescaledHLTTrigger("HLT_Ele17_SW_TighterEleIdIsol_L1R_v1")) return true;
    if( passUnprescaledHLTTrigger("HLT_Ele17_SW_TighterEleIdIsol_L1R_v2")) return true;
    if( passUnprescaledHLTTrigger("HLT_Ele17_SW_TighterEleIdIsol_L1R_v3")) return true;

    if( passUnprescaledHLTTrigger("HLT_Ele22_SW_TighterEleId_L1R_v2"))     return true;
    if( passUnprescaledHLTTrigger("HLT_Ele22_SW_TighterEleId_L1R_v3")) return true;
    if( passUnprescaledHLTTrigger("HLT_Ele22_SW_TighterCaloIdIsol_L1R_v2")) return true;

    if( passUnprescaledHLTTrigger("HLT_Ele27_SW_TightCaloEleIdTrack_L1R_v1")) return true;

    if( passUnprescaledHLTTrigger("HLT_Ele32_SW_TightCaloEleIdTrack_L1R_v1")) return true;
    if( passUnprescaledHLTTrigger("HLT_Ele32_SW_TighterEleId_L1R_v2")) return true;

    
    // This is a family of never prescaled Ele17 triggers
    // Some of them are probably prescaled at 2e32?
    // they did get prescaled 148058... never say never...
    if( passUnprescaledHLTTrigger("HLT_Ele17_SW_TightEleId_L1R") )          return true;
    if( passUnprescaledHLTTrigger("HLT_Ele17_SW_TighterEleId_L1R_v1") )     return true;
    
    //if( run <= 147116 ){
    // This is a family of Ele17 triggers that are not
    // there anymore or are not unprescaled anymore.
    // The loosest one is indicated as master
    // In principle, the other two could be skipped
    if( passUnprescaledHLTTrigger("HLT_Ele17_SW_CaloEleId_L1R") )                 return true; //146428-147116
    if( passUnprescaledHLTTrigger("HLT_Ele17_SW_EleId_L1R") )                     return true; //146428-147116
    if( passUnprescaledHLTTrigger("HLT_Ele17_SW_LooseEleId_L1R") )                return true; //146428-147116
    //}
    
    //if( run <= 144114 ){
    // This is a family of Ele15 triggers that are not
    // there anymore or are not unprescaled anymore.
    // The loosest one is indicated as master
    // In principle, the other one could be skipped
    if( passUnprescaledHLTTrigger("HLT_Ele15_SW_CaloEleId_L1R") ) return true; //141956-144114 
    if( passUnprescaledHLTTrigger("HLT_Ele15_SW_EleId_L1R") )     return true; //141956-144114
    //}
    
    // These are Ele15 that got prescaled in a different
    // time frame as the previous family...use both
    //if( run <= 143962 ){
    if( passUnprescaledHLTTrigger("HLT_Ele15_SW_L1R") ) return true; //140058-143962
    //}
    //if( run <= 141882 ){
    if( passUnprescaledHLTTrigger("HLT_Ele15_LW_L1R") ) return true; //136033-141882
    //}
    
    // This is a bonus trigger, should not be needed, but why not?
    //if( run <= 144114 ){
    if( passUnprescaledHLTTrigger("HLT_Ele20_SW_L1R") ) return true; //140058-144114
    //}
    
    // This is the Ele10 family that came and went at ramdom times
    // Use all (although in principle the one with EleID could be
    // skipped in the subset of runs where triggers with no EleId
    // were present)
    //if( run <= 144114 ){
    if( passUnprescaledHLTTrigger("HLT_Ele10_SW_EleId_L1R") ) return true;      //141956-144114
    //}
    //if( run <= 141882 ){
    if( passUnprescaledHLTTrigger("HLT_Ele10_LW_EleId_L1R") ) return true;      //136033-141882
    //} 
    //if( run <= 139980 ){
    if( passUnprescaledHLTTrigger("HLT_Ele10_LW_L1R") )       return true;      //136033-139980
    if( passUnprescaledHLTTrigger("HLT_Ele10_SW_L1R") )       return true;      //139195-139980
    //}
    
    //---------------------------
    // e-mu cross triggers
    //---------------------------
    
    if( passUnprescaledHLTTrigger("HLT_Mu5_Ele5_v1") )          return true;
    if( passUnprescaledHLTTrigger("HLT_Mu5_Ele9_v1") )          return true;
    
    //Added 2e32
    if( passUnprescaledHLTTrigger("HLT_Mu11_Ele8_v1") )     return true;    // Mu3 L1 Seed
    if( passUnprescaledHLTTrigger("HLT_Mu8_Ele8_v1") )      return true;    // Mu3 L1 Seed
    if( passUnprescaledHLTTrigger("HLT_Mu5_Ele13_v2") )     return true;    // EG8 L1 Seed
    if( passUnprescaledHLTTrigger("HLT_Mu5_Ele17_v1") )     return true;    // EG8 L1 Seed & open mu Seed
    if( passUnprescaledHLTTrigger("HLT_Mu5_Ele17_v2") )     return true;    // EG8 L1 Seed & open mu Seed

  }
  
  else{
    cout << "ERROR: unrecognized hyptype " << hypType << ", quitting" << endl;
    exit(0);
  }

  return false;
}


//----------------------------------------------------------
//this is a simplified version of the SUSY triggers
//----------------------------------------------------------


bool passSimpleSUSYTrigger_v1( bool isData ) {

  //currently do not require any triggers for MC
  if( !isData ) return true;

  if( passUnprescaledHLTTrigger("HLT_Mu9") )                                    return true;            
  if( passUnprescaledHLTTrigger("HLT_Mu15_v1 147196-148862") )                  return true;
  if( passUnprescaledHLTTrigger("HLT_Ele10_LW_L1R") )                           return true;             
  if( passUnprescaledHLTTrigger("HLT_Ele10_LW_EleId_L1R") )                     return true;      
  if( passUnprescaledHLTTrigger("HLT_Ele15_SW_L1R") )                           return true; 
  if( passUnprescaledHLTTrigger("HLT_Ele15_SW_CaloEleId_L1R") )                 return true;  
  if( passUnprescaledHLTTrigger("HLT_Ele17_SW_CaloEleId_L1R") )                 return true; 
  if( passUnprescaledHLTTrigger("HLT_Ele17_SW_TightEleId_L1R") )                return true; 
  if( passUnprescaledHLTTrigger("HLT_DoubleEle10_SW_L1R") )                     return true;     
  if( passUnprescaledHLTTrigger("HLT_DoubleEle15_SW_L1R_v1") )                  return true; 
  if( passUnprescaledHLTTrigger("HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v1") )  return true; 
  if( passUnprescaledHLTTrigger("HLT_Mu5_Ele5_v1") )                            return true; 
  if( passUnprescaledHLTTrigger("HLT_Mu8_Ele8_v1") )                            return true; 
  if( passUnprescaledHLTTrigger("HLT_Mu5_Ele17_v1") )                           return true;

  return false;
}
/*****************************************************************************************/
//hypothesis disambiguation. Returns the hypothesis that has mass closest to MZ
/*****************************************************************************************/
unsigned int selectBestZHyp(const vector<unsigned int> &v_goodHyps) {
  
  float mindeltam         = 100000;
  unsigned int bestHypIdx = 0;
  for(unsigned int i = 0; i < v_goodHyps.size(); i++) {
    unsigned int index = v_goodHyps.at(i);
    float deltam = fabs( hyp_p4()[index].mass() - 91. );
    if( deltam < mindeltam ) {
      mindeltam = deltam;
      bestHypIdx = index;
    }
  }

  return bestHypIdx;

}


