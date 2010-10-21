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
bool passMuTriggerSS(unsigned int hypIdx, bool mc);

/*****************************************************************************************/
// extra Z veto
/*****************************************************************************************/
bool makesExtraZ(int hypIdx, bool applyAlignmentCorrection=false, bool removedEtaCutInEndcap=false);

/*****************************************************************************************/
// require electron GSF, CTF and SC charges agree
/*****************************************************************************************/
bool passThreeChargeRequirement(int elIdx);

/******************************************************************************************/     
// is it a good jet?
/******************************************************************************************/     
bool isGoodJet(LorentzVector jetp4, double ptCut, double absEtaCut, double dRCut, bool muJetClean);



/******************************************************************************************/     
// good lepton (either mu or electron, no isolation cuts), d0 corrected for PV
/******************************************************************************************/
bool isGoodLeptonNoIsoSSd0PV(int id, int lepIdx, bool applyAlignmentCorrection=false, bool removedEtaCutInEndcap=false);

/******************************************************************************************/     
// are the leptons in the hypothesis good (all cuts but isolation?), d0 corrected for PV
/******************************************************************************************/
bool isGoodHypNoIsoSSd0PV(int hypIdx, bool applyAlignmentCorrection=false, bool removedEtaCutInEndcap=false);



/******************************************************************************************/     
// good lepton (either mu or electron, no isolation cuts), no d0 cut
/******************************************************************************************/
bool isGoodLeptonNoIsoSSnod0(int id, int lepIdx, bool applyAlignmentCorrection=false, bool removedEtaCutInEndcap=false);

/******************************************************************************************/     
// are the leptons in the hypothesis good (all cuts but isolation?), no d0 cut
/******************************************************************************************/
bool isGoodHypNoIsoSSnod0(int hypIdx, bool applyAlignmentCorrection=false, bool removedEtaCutInEndcap=false);


#endif

