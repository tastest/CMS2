#include "CMS2.h"

int getDrellYanType();
bool isDYee();
bool isDYmm();
bool isDYtt();
bool isWe();
bool isWm();
bool isWt();
bool isWW();
bool isWZ();
bool isZZ();
int getZZType();
void dumpDocLines();

//--------------------------------------------------------
// Determines if the lepton in question is from W/Z
// and if its charge is correct
//
// Note that if we have
//     W->lepton->lepton gamma
// where the gamma is at large angles and it is the
// gamma that gives the lepton signature in the detector,
// then this returns "not from W/Z".  This is by design.
//
// Note W->tau->lepton is tagged as "from W"
//
// Inputs:  idx   = index in the els or mus block
//          id    = lepton ID (11 or 13 or -11 or -13)
//          v     = 4-vector of reco lepton
//
// Output:  2 = from W/Z incorrect charge
//          1 = from W/Z   correct charge
//          0 = not matched to a lepton (= fake)
//         -1 = lepton from b decay
//         -2 = lepton from c decay
//         -3 = lepton from some other source
//        
// Authors: Claudio in consultation with fkw 22-july-09    
//---------------------------------------------------------
int leptonIsFromW(int idx, int id, LorentzVector v);
bool trueGammaFromMuon(int electron);
// count genp leptons
//-------------------------------------------------- 
// Returns the number of e,mu, and tau in the doc lines 
//----------------------------------------------------- 
int leptonGenpCount(int& nele, int& nmuon, int& ntau);
int leptonGenpCount_lepTauDecays(int& nele, int& nmuon, int& ntau);
int genpDileptonType();
  //0 mumu; 1 emu; 2 ee
  
// -----------------------------------------
// MC helper functions for fakerate tests:
// -----------------------------------------
int elFakeMCCategory(int i_el);
int muFakeMCCategory(int i_mu);
bool idIsCharm(int id);
bool idIsBeauty(int id);
