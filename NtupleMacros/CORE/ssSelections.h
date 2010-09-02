#ifndef SSSELECTIONS_H
#define SSSELECTIONS_H

#include <vector>
#include "Math/LorentzVector.h"


typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;
/******************************************************************************************/     
// good lepton (either mu or electron, no isolation cuts)
/******************************************************************************************/
bool isGoodLeptonNoIsoSS(int id, int lepIdx, bool applyAlignmentCorrection, bool removedEtaCutInEndcap=true);

/******************************************************************************************/     
// isolated lepton (either mu or electron)
/******************************************************************************************/
bool isGoodLeptonwIsoSS(int id, int lepIdx, bool applyAlignmentCorrection, bool removedEtaCutInEndcap=true);

/******************************************************************************************/     
// are the leptons in the hypothesis good (all cuts but isolation?)
/******************************************************************************************/
bool isGoodHypNoIsoSS(int hypIdx, bool applyAlignmentCorrection, bool removedEtaCutInEndcap=true);

/******************************************************************************************/     
// are the leptons in the hypothesis isolated?
/******************************************************************************************/     
bool isGoodHypwIsoSS(int hypIdx, bool applyAlignmentCorrection, bool removedEtaCutInEndcap=true);

/*****************************************************************************************/
//passes the EGamma triggers
/*****************************************************************************************/
bool passEGTriggerSS(unsigned int hypIdx, bool mc);

/*****************************************************************************************/
//passes the muon triggers
/*****************************************************************************************/
bool passMuTriggerSS(unsigned int hypIdx);

/*****************************************************************************************/
// extra Z veto
/*****************************************************************************************/
bool makesExtraZ(int hypIdx);

#endif

