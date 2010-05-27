//--------------------------------------------------
// Trigger utilities stolen from Derek and adapted
//-------------------------------------------------

//int nTriggeredJets(char*);
//LorentzVector singleJetBias( char*, int) ;
#include "CMS2.h"
#include "triggerUtils.h"

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
