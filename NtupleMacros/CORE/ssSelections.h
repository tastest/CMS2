#ifndef SSSELECTIONS_H
#define SSSELECTIONS_H

#include <vector>
#include "Math/LorentzVector.h"


typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;
/******************************************************************************************/     
// good lepton (either mu or electron, no isolation cuts)
/******************************************************************************************/
bool isGoodLeptonNoIsoSS(int id, int lepIdx, bool applyAlignmentCorrection=false, bool removedEtaCutInEndcap=false);

/******************************************************************************************/     
// isolated lepton (either mu or electron, no ID cuts)
/******************************************************************************************/
bool isGoodLeptonwIsoSS(int id, int lepIdx, bool applyAlignmentCorrection=false, bool removedEtaCutInEndcap=false);

/******************************************************************************************/     
// are the leptons in the hypothesis good (all cuts but isolation?)
/******************************************************************************************/
bool isGoodHypNoIsoSS(int hypIdx, bool applyAlignmentCorrection=false, bool removedEtaCutInEndcap=false);

/******************************************************************************************/     
// are the leptons in the hypothesis isolated? (no ID or other cuts applied)
/******************************************************************************************/     
bool isGoodHypwIsoSS(int hypIdx, bool applyAlignmentCorrection=false, bool removedEtaCutInEndcap=false);

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

/*****************************************************************************************/
// require electron GSF, CTF and SC charges agree
/*****************************************************************************************/
bool passThreeChargeRequirement(int elIdx);

#endif

