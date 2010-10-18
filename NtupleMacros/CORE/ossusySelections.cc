#include <assert.h>
#include <algorithm>
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TDatabasePDG.h"
#include "triggerUtils.h"
#include "CMS2.h"

#include "ossusySelections.h"

using namespace tas;

/*****************************************************************************************/
//passes the OS SUSY trigger selection
/*****************************************************************************************/
bool passTrigger_OSV1( bool isData , int hypType ) {

  int run = cms2.evt_run();

  //currently do not require any triggers for MC
  if( !isData ) return true;

  //mumu
  if( hypType == 0 ){

    if( run < 147117 ){
      if(nHLTObjects("HLT_Mu9") != 0)            return true; 
      if(nHLTObjects("HLT_DoubleMu3") != 0)      return true;
    }

    //note here the overlap between HLT_Mu11 and HLT_Mu15_v1: is HLT_Mu15_v1 needed?
    //also note that for 147117 - 147196 the lowest muon pt cut is > 10 GeV
    if( run >= 147117 ) { 
      if(nHLTObjects("HLT_Mu11") > 0) return true;            
    }

    if( run >= 147196 ){
      if(nHLTObjects("HLT_Mu15_v1") != 0)        return true;
      if(nHLTObjects("HLT_DoubleMu5_v1") != 0)   return true;
    }

  }
 
  //ee
  else if( hypType == 3 ){
    
    if (run < 138000) {
      if(nHLTObjects("HLT_Ele10_LW_L1R") != 0 )      return true;
      if(nHLTObjects("HLT_Ele15_LW_L1R") != 0 )	     return true;
      if(nHLTObjects("HLT_DoubleEle5_SW_L1R") != 0 ) return true;
    }

    if (run >= 138000 && run < 141900) {
      if(nHLTObjects("HLT_Ele15_LW_L1R") != 0 )       return true;
      if(nHLTObjects("HLT_Ele10_LW_EleId_L1R") != 0 ) return true;
      if(nHLTObjects("HLT_DoubleEle5_SW_L1R") != 0 )  return true;
    }

    if (run >= 141900 && run < 141956) {
      if(nHLTObjects("HLT_Ele15_SW_L1R") != 0 )       return true;
      //need 10 GeV trigger here??
    }
    
    if( run >= 141956 && run < 144115){
      if(nHLTObjects("HLT_Ele10_SW_EleId_L1R") != 0 )                    return true;
      if(nHLTObjects("HLT_Ele15_SW_CaloEleId_L1R") != 0 )                return true;
      if(nHLTObjects("HLT_Ele15_SW_EleId_L1R") != 0 )                    return true;   
    }

    if( run >= 146428 && run < 147116 ){
      if(nHLTObjects("HLT_Ele17_SW_LooseEleId_L1R") != 0 )               return true;
      if(nHLTObjects("HLT_Ele17_SW_CaloEleId_L1R") != 0 )                return true;
      if(nHLTObjects("HLT_Ele17_SW_EleId_L1R") != 0 )                    return true;
      if(nHLTObjects("HLT_DoubleEle10_SW_L1R") != 0 )                    return true;   
    }

    if( run >= 147196 ){
      if(nHLTObjects("HLT_Ele17_SW_TightCaloEleId_SC8HE_L1R_v1") != 0 )  return true;
      if(nHLTObjects("HLT_Ele17_SW_TightEleIdIsol_L1R_v1") != 0 )        return true;
      if(nHLTObjects("HLT_DoubleEle15_SW_L1R_v1") != 0 )                 return true;
    }
      
  }

  //emu
  else if( hypType == 1 || hypType == 2 ){

    //single muon triggers
    if( run < 147117 ){
      if(nHLTObjects("HLT_Mu9") != 0)            return true; 
    }

    if( run >= 147117 ) { 
      if(nHLTObjects("HLT_Mu11") > 0) return true;            
    }

    if( run >= 147196 ){
      if(nHLTObjects("HLT_Mu15_v1") != 0)        return true;
    }

    //single electron triggers
    if (run < 138000) {
      if(nHLTObjects("HLT_Ele10_LW_L1R") != 0 )      return true;
      if(nHLTObjects("HLT_Ele15_LW_L1R") != 0 )	     return true;
    }

    if (run >= 138000 && run < 141900) {
      if(nHLTObjects("HLT_Ele15_LW_L1R") != 0 )       return true;
      if(nHLTObjects("HLT_Ele10_LW_EleId_L1R") != 0 ) return true;
    }

    if (run >= 141900 && run < 141956) {
      if(nHLTObjects("HLT_Ele15_SW_L1R") != 0 )       return true;
      //need 10 GeV trigger here??
    }
    
    if( run >= 141956 && run < 144115){
      if(nHLTObjects("HLT_Ele10_SW_EleId_L1R") != 0 )                    return true;
      if(nHLTObjects("HLT_Ele15_SW_CaloEleId_L1R") != 0 )                return true;
      if(nHLTObjects("HLT_Ele15_SW_EleId_L1R") != 0 )                    return true;   
    }

    if( run >= 146428 && run < 147116 ){
      if(nHLTObjects("HLT_Ele17_SW_LooseEleId_L1R") != 0 )               return true;
      if(nHLTObjects("HLT_Ele17_SW_CaloEleId_L1R") != 0 )                return true;
      if(nHLTObjects("HLT_Ele17_SW_EleId_L1R") != 0 )                    return true;
    }

    if( run >= 147196 ){
      if(nHLTObjects("HLT_Ele17_SW_TightCaloEleId_SC8HE_L1R_v1") != 0 )  return true;
      if(nHLTObjects("HLT_Ele17_SW_TightEleIdIsol_L1R_v1") != 0 )        return true;
    }
    
    //e-mu cross triggers
    if( run >= 147196 ){
      if(nHLTObjects("HLT_Mu5_Ele9_v1") != 0 )        return true;
      //if(nHLTObjects("HLT_Mu5_Ele15_v1") != 0 )       return true; //trigger not found???
      //if(nHLTObjects("HLT_Mu8_Ele8_v1") != 0 )        return true; //trigger not found???
    }
    
  }

  else{
    cout << "ERROR: unrecognized hyptype " << hypType << ", quitting" << endl;
    exit(0);
  }


  return false;
}

/*****************************************************************************************/
//hypothesis disambiguation. Returns the hypothesis that has mass closent to MZ
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
