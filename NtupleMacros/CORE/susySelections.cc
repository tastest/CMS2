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

using namespace tas;

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

    // This is a family of never prescaled Ele17 triggers
    // Some of them are probably prescaled at 2e32?
    if( passUnprescaledHLTTrigger("HLT_Ele17_SW_TightEleId_L1R") )          return true;
    if( passUnprescaledHLTTrigger("HLT_Ele17_SW_TighterEleId_L1R_v1") )     return true;
    
    // These are unprescaled double triggers
    if( passUnprescaledHLTTrigger("HLT_DoubleEle15_SW_L1R_v1") )                 return true;
    if( passUnprescaledHLTTrigger("HLT_DoubleEle17_SW_L1R_v1") )                 return true; //<---Added 2e32
    if( passUnprescaledHLTTrigger("HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v1") ) return true;
    if( passUnprescaledHLTTrigger("HLT_Ele17_SW_TightCaloEleId_SC8HE_L1R_v1") )  return true;

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

    //Added 2e32
    if( passUnprescaledHLTTrigger("HLT_Ele17_SW_TighterEleIdIsol_L1R_v2")) return true;
    if( passUnprescaledHLTTrigger("HLT_Ele22_SW_TighterEleId_L1R_v2"))     return true;
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
    if( passUnprescaledHLTTrigger("HLT_Ele17_SW_TighterEleIdIsol_L1R_v2")) return true;
    if( passUnprescaledHLTTrigger("HLT_Ele22_SW_TighterEleId_L1R_v2"))     return true;
    
    // This is a family of never prescaled Ele17 triggers
    // Some of them are probably prescaled at 2e32?
    if( passUnprescaledHLTTrigger("HLT_Ele17_SW_TightEleId_L1R") )          return true;
    if( passUnprescaledHLTTrigger("HLT_Ele17_SW_TighterEleId_L1R_v1") )     return true;
    
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
    
    //---------------------------
    // e-mu cross triggers
    //---------------------------
    
    if( passUnprescaledHLTTrigger("HLT_Mu5_Ele5_v1") )          return true;
    if( passUnprescaledHLTTrigger("HLT_Mu5_Ele9_v1") )          return true;
    
    //Added 2e32
    if( passUnprescaledHLTTrigger("HLT_Mu11_Ele8_v1") )     return true;    // Mu3 L1 Seed
    if( passUnprescaledHLTTrigger("HLT_Mu8_Ele8_v1") )      return true;    // Mu3 L1 Seed
    if( passUnprescaledHLTTrigger("HLT_Mu5_Ele13_v2") )     return true;    // EG8 L1 Seed
    if( passUnprescaledHLTTrigger("HLT_Mu5_Ele17_v1") )     return true;    // EG8 L1 Seed

  }
  
  else{
    cout << "ERROR: unrecognized hyptype " << hypType << ", quitting" << endl;
    exit(0);
  }

  return false;
}

/*****************************************************************************************/
//hypothesis disambiguation. Returns the hypothesis that has mass closest to MZ
/*****************************************************************************************/
unsigned int selectBestZHyp(const vector<unsigned int> &v_goodHyps) {
  
  float mindeltam         = 100;
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


