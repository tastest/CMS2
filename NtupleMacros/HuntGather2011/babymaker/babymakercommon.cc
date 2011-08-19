#include "babymakercommon.h"
#include "CORE/CMS2.h"

#include "TMath.h"
#include "CORE/jetSelections.h"

float deltaPhi (float phi1, float phi2)
{
    float dphi = phi1-phi2;
    if (dphi < 0.) dphi = -1.*dphi;
    if (dphi > TMath::Pi()) dphi = 2.*TMath::Pi() - dphi;
    return dphi;
}

double dRbetweenVectors (const LorentzVector &vec1, const LorentzVector &vec2)
{
    float deta = vec1.eta()-vec2.eta();
    float dphi = deltaPhi(vec1.phi(), vec2.phi());
    return sqrt(dphi*dphi + deta*deta);
}

bool sortByPt (const LorentzVector &vec1, const LorentzVector &vec2)
{
    return vec1.pt() > vec2.pt();
}

bool sortByPFJetPt (const unsigned int &index1, const unsigned int &index2)
{
    return cms2.pfjets_p4()[index1].pt() > cms2.pfjets_p4()[index2].pt();
}
/*
   bool isGoodPFJet(unsigned int pfjeti) {
   bool pfjetisGood = false;

   float pfjet_pt_   = cms2.pfjets_p4()[pfjeti].pt();
   float pfjet_eta_  = cms2.pfjets_p4()[pfjeti].eta();
//float pfjet_phi_  = cms2.pfjets_p4()[pfjeti].phi();
float pfjet_chf_  = cms2.pfjets_chargedHadronE()[pfjeti]/cms2.pfjets_p4()[pfjeti].energy();
float pfjet_nhf_  = cms2.pfjets_neutralHadronE()[pfjeti]/cms2.pfjets_p4()[pfjeti].energy();
float pfjet_cef_  = cms2.pfjets_chargedEmE()[pfjeti]/cms2.pfjets_p4()[pfjeti].energy();
float pfjet_nef_  = cms2.pfjets_neutralEmE()[pfjeti]/cms2.pfjets_p4()[pfjeti].energy();

if (pfjet_pt_ > 25. && fabs(pfjet_eta_) < 3.0 && pfjet_nhf_ < 1. && pfjet_cef_ < 1. && pfjet_nef_ < 1.)
{
if (fabs(pfjet_eta_) > 2.4)
pfjetisGood = true;
else if (pfjet_chf_ > 0.)
pfjetisGood = true;
}

return pfjetisGood;
}
 */


bool isGoodPFJet(unsigned int pfjeti) {

    bool pfjetisGood = false;
    float pfjet_pt_   = cms2.pfjets_p4()[pfjeti].pt(); 
    float pfjet_eta_  = cms2.pfjets_p4()[pfjeti].eta();

    if (pfjet_pt_ > 25. && fabs(pfjet_eta_) < 3.0 )                                                                                  
    {                                         
        if( passesPFJetID(pfjeti) ) pfjetisGood = true;
    }

    return pfjetisGood;

}

