#ifndef TTBARSELECTIONS_H
#define TTBARSELECTIONS_H

#include <vector>
#include "Math/LorentzVector.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

/******************************************************************************************/     
// good lepton (either mu or electron, no isolation cuts)
/******************************************************************************************/
bool isGoodLeptonNoIso(int id, int lepIdx);

/******************************************************************************************/     
// isolated lepton (either mu or electron)
/******************************************************************************************/
bool isGoodLeptonwIso(int id, int lepIdx);

/******************************************************************************************/     
// are the leptons in the hypothesis good (all cuts but isolation?)
/******************************************************************************************/
bool isGoodHypNoIso(int hypIdx);
  
/******************************************************************************************/     
// are the leptons in the hypothesis isolated?
/******************************************************************************************/     
bool isGoodHypwIso(int hypIdx);
  
/******************************************************************************************/     
// is it a good jet?
/******************************************************************************************/     
bool isGoodDilHypJet(LorentzVector jetp4, unsigned int& hypIdx, double ptCut, double absEtaCut, double dRCut, bool muJetClean);
             
/******************************************************************************************/     
//return the MET and the MET phi instead of a bool because the MT2 needs it
/******************************************************************************************/     
std::pair<float,float> getMet(string& algo);
  
/******************************************************************************************/     
//trigger requirement
/******************************************************************************************/         
bool passTriggersMu9orLisoE15(int dilType);
  
/******************************************************************************************/     
//hypothesis disabmiguation
/******************************************************************************************/     
int eventDilIndexByWeightTTDil08(const std::vector<unsigned int>& goodHyps, int& strasbourgDilType, bool printDebug, bool usePtOnlyForWeighting);
// a simpler version of the above maximising pt(ll) + pt(lt)
unsigned int eventDilIndexByWeightTTDil10(const std::vector<unsigned int> &hyp_index_selected);

/******************************************************************************************/
// MET cut - MET > 20 (OF) MET > 30 (SF)
/******************************************************************************************/
bool passMetAsIs_OF20_SF30(float met, int hyp_type);

#endif

