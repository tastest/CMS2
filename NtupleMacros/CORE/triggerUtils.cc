//--------------------------------------------------
// Trigger utilities stolen from Derek and adapted
//-------------------------------------------------

//int nTriggeredJets(char*);
//LorentzVector singleJetBias( char*, int) ;
#include "CMS2.h"
#include "triggerUtils.h"

//--------------------------------------------------
// EG trigger selection from 5 July 2010
// data: Photon10 OR Electron 10 OR Photon 15
// MC  : E10 OR Photon15
// (Note: in data Photon10 is prescaled at some point)
// (Note: in MC Photon15 has a different L1 requirement, so we 
//        use Photon10 and we tighten the pt threshold)
//--------------------------------------------------
bool goodEGTrigger5July2010 (bool mc) {
  if (mc) {
 
   int e10 = nHLTObjects("HLT_Ele10_LW_L1R");
    if (e10 != 0) return true;

    int p10 = nHLTObjects("HLT_Photon10_L1R");
    for (int i=0; i<p10; i++) {
      LorentzVector v = p4HLTObject("HLT_Photon10_L1R", i);
      if (v.pt() > 15.) return true;
    }
 
  } else {  // data now
    int e10 = nHLTObjects("HLT_Ele10_LW_L1R");
    if (e10 != 0) return true;

    int p10 = nHLTObjects("HLT_Photon10_L1R");
    if (p10 != 0) return true;

    p10 = nHLTObjects("HLT_Photon10_Cleaned_L1R");
    if (p10 != 0) return true;

    int p15 = nHLTObjects("HLT_Photon15_L1R");
    if (p15 != 0) return true;

    p15 = nHLTObjects("HLT_Photon15_Cleaned_L1R");
    if (p15 != 0) return true;

 
  }
  return false;
}



//-----------------------------------------------------
// Returns the nth object that passes a given trigger
// (n starts from 0)
///----------------------------------------------------
LorentzVector p4HLTObject( char* arg, int objNumber){
 
  TString HLTTrigger( arg );
  int trigIndx = -1;
  vector<TString>::const_iterator begin_it = cms2.hlt_trigNames().begin();
  vector<TString>::const_iterator end_it = cms2.hlt_trigNames().end();
  vector<TString>::const_iterator found_it = find(begin_it, end_it, HLTTrigger );
  if( (found_it != end_it) ){
    trigIndx = found_it - begin_it;
    //cout << "p4HLTObject: Found Trigger: " << arg << endl;
  }
  else {
    cout << "p4HLTObject: Cannot find Trigger: " << arg << endl;
    gSystem->Exit(1);
  }

  int nobj = cms2.hlt_trigObjs_p4().at(trigIndx).size();
  if (nobj == 0 ) {
    cout << "ERROR: nobj == 0" << endl;
    gSystem->Exit(1);
  }

  if (objNumber > (nobj-1)) {
    cout << "ERROR: requested object number " << objNumber << " but we only have " << nobj <<endl;
    gSystem->Exit(1);
  }

  return cms2.hlt_trigObjs_p4().at(trigIndx).at(objNumber);

}
//--------------------------------------------------------
// Returns the number of objects passing a given trigger
// Returns zero if the trigger failed
// Returns -1 if the trigger passed but no onjects were found
//--------------------------------------------------------
int nHLTObjects( char* arg ){

  // put the trigger name into a string
  TString HLTTrigger( arg );

  // Did the trigger pass?
  if ( !(cms2.passHLTTrigger(HLTTrigger)) ) return 0;

  // The trigger passed, see how many associated objects there are
  int trigIndx = -1;
  vector<TString>::const_iterator begin_it = cms2.hlt_trigNames().begin();
  vector<TString>::const_iterator end_it = cms2.hlt_trigNames().end();
  vector<TString>::const_iterator found_it = find(begin_it, end_it, HLTTrigger );
  if( (found_it != end_it) ){
    trigIndx = found_it - begin_it;
    //cout << "nHLTObjects: Found Trigger: " << arg << endl;
  }
  else {
    cout << "nHLTObjects: Cannot find Trigger " << arg << endl;
    return 0;
  }

  int nobj = 0;
  for( unsigned int i=0; i < cms2.hlt_trigObjs_p4().at(trigIndx).size(); i++ ){
    nobj++;
    //cout << "\t" << i << ", (pt, eta, phi): " << hlt_trigObjs_p4().at(trigIndx).at(i).pt() << " "
    //              << hlt_trigObjs_p4().at(trigIndx).at(i).eta() << " " << hlt_trigObjs_p4().at(trigIndx).at(i).phi() << endl;
  }

  // cout << " Number of jets = " << njets << endl;

  if (nobj == 0) return -1;
  return nobj;
}
//---------------------------------------------
// Utility to print trigger names in the file
//---------------------------------------------
void PrintTriggers(){
  for( unsigned int i = 0; i < cms2.hlt_trigNames().size(); i++ ){
    cout << cms2.hlt_trigNames().at(i).Data() << endl;
  } 
  cout << endl;
}
