#ifndef TTBARSELECTIONS_H
#define TTBARSELECTIONS_H

#include <vector>
#include "Math/LorentzVector.h"


typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;
/******************************************************************************************/     
// good lepton (either mu or electron, no isolation cuts)
/******************************************************************************************/
bool isGoodLeptonNoIso(int id, int lepIdx, bool applyAlignmentCorrection, bool removedEtaCutInEndcap=true);

/******************************************************************************************/     
// isolated lepton (either mu or electron)
/******************************************************************************************/
bool isGoodLeptonwIso(int id, int lepIdx, bool applyAlignmentCorrection, bool removedEtaCutInEndcap=true);

/******************************************************************************************/     
// are the leptons in the hypothesis good (all cuts but isolation?)
/******************************************************************************************/
bool isGoodHypNoIso(int hypIdx, bool applyAlignmentCorrection, bool removedEtaCutInEndcap=true);

/******************************************************************************************/     
// are the leptons in the hypothesis isolated?
/******************************************************************************************/     
bool isGoodHypwIso(int hypIdx, bool applyAlignmentCorrection, bool removedEtaCutInEndcap=true);

/******************************************************************************************/     
// is it a good jet?
/******************************************************************************************/     
bool isGoodDilHypJet(LorentzVector jetp4, unsigned int& hypIdx, double ptCut, double absEtaCut, double dRCut, bool muJetClean);

/******************************************************************************************/     
//return the MET and the MET phi, correcting for mus that are not corrected for by default
/******************************************************************************************/     
std::pair<float,float> getMet(const string algo, unsigned int hypIdx);

/*****************************************************************************************/
//hypothesis disambiguation. Returns the hypothesis that has the highest sum Pt
/*****************************************************************************************/
unsigned int selectHypByHighestSumPt(const vector<unsigned int> &v_goodHyps);

/******************************************************************************************/     
//corrects tcMET for mus that are not corrected for by default
/******************************************************************************************/     
//void correctTcMETForHypMus(unsigned int hypIdx, double& met, double& metPhi);

/*****************************************************************************************/
//passes the EGamma triggers
/*****************************************************************************************/
bool passEGTrigger(bool mc);

/*****************************************************************************************/
//get the impact parameter wrt the PV
/*****************************************************************************************/
double getd0wrtPV(LorentzVector p4, float d0);

/*****************************************************************************************/
//get the number of jets passing btag discriminator cuts
// takes as arguments a vector of indices, the jet algorithm
// and the btag discriminator
// the working points are hard coded
/*****************************************************************************************/
int getNbtags(const vector<unsigned int> v_jetsIdx, const string jetAlgo, const string bTagDiscriminator);

#endif

