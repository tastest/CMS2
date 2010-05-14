
//
// electron selections
//

// CMS2 includes
#include "electronSelections.h"
#include "CMS2.h"

bool pass_electronSelection(const unsigned int index, const cuts_t selectionType)
{
    cuts_t cuts_passed = electronSelection(index);
    if ((cuts_passed & selectionType) == selectionType) return true;
    return false;
}

cuts_t electronSelection(const unsigned int index) 
{

    //
    // keep track of which cuts passed
    //

    cuts_t cuts_passed = 0;

    //
    // isolation
    //

    if (electronIsolation_rel(index, true) < 0.10) cuts_passed |= (1ll<<ELEISO_REL010);
    if (electronIsolation_rel(index, true) < 0.15) cuts_passed |= (1ll<<ELEISO_REL015);
    if (electronIsolation_rel(index, true) < 0.40) cuts_passed |= (1ll<<ELEISO_REL040);

    //
    // ip
    //

    if (fabs(cms2.els_d0corr()[index]) < 0.02) cuts_passed |= (1ll<<ELEIP_200);
    if (fabs(cms2.els_d0corr()[index]) < 0.04) cuts_passed |= (1ll<<ELEIP_400);

    //
    // id
    //

    if (electronId_cand(index, CAND_01)) cuts_passed |= (1ll<<ELEID_CAND01);
    if (electronId_cand(index, CAND_02)) cuts_passed |= (1ll<<ELEID_CAND02);
    if (electronId_extra(index)) cuts_passed |= (1ll<<ELEID_EXTRA);

    //
    // conversion rejection cuts
    //

    if (!isFromConversionPartnerTrack(index)) cuts_passed |= (1ll<<ELENOTCONV_DISTDCOT002);
    if (!isFromConversionHitPattern(index)) cuts_passed |= (1ll<<ELENOTCONV_HITPATTERN);

    //
    // fiduciality/other cuts
    //

    if ((cms2.els_type()[index] & (1ll<<ISECALDRIVEN))) cuts_passed |= (1ll<<ELESEED_ECAL);
    if (fabs(cms2.els_p4()[index].eta()) < 2.5) cuts_passed |= (1ll<<ELEETA_250);
    if (electronId_noMuon(index)) cuts_passed |= (1ll<<ELENOMUON_010);

    //
    // chargeflip
    //

    if (!isChargeFlip(index)) cuts_passed |= (1ll<<ELECHARGE_NOTFLIP);

    //
    // return which selections passed
    //

    return cuts_passed;

}

//
// if fbrem is low then cut on e/p_in
//
bool electronId_extra(const unsigned int index)
{
    if (cms2.els_fbrem()[index] < 0.2) {
        if (cms2.els_eOverPIn()[index] < 0.7 || cms2.els_eOverPIn()[index] > 1.5) return false;
    }

    return true;
}

//
// remove if close to a muon
//
bool electronId_noMuon(const unsigned int index)
{
    if ( cms2.els_closestMuon().at(index) != -1) return false;
    return true;
}

//
// candidate electron id function
//
bool electronId_cand(const unsigned int index, const cand_tightness tightness)
{

    std::vector<double> relisoThresholds;
    std::vector<double> dEtaInThresholds;
    std::vector<double> dPhiInThresholds;
    std::vector<double> hoeThresholds;
    std::vector<double> latThresholds;

    eidGetCand(tightness, dEtaInThresholds, dPhiInThresholds, hoeThresholds, latThresholds);

    //
    // apply cuts
    //
    if (fabs(cms2.els_etaSC()[index]) < 1.479) {
        if (fabs(cms2.els_dEtaIn()[index]) > dEtaInThresholds[0]) 	return false;
        if (fabs(cms2.els_dPhiIn()[index]) > dPhiInThresholds[0]) 	return false;
        if (cms2.els_hOverE()[index] > hoeThresholds[0]) 		    return false;
        if ((cms2.els_e2x5Max()[index]/cms2.els_e5x5()[index]) < latThresholds[0]) return false;
    }
    if (fabs(cms2.els_etaSC()[index]) > 1.479) {
        if (fabs(cms2.els_dEtaIn()[index]) > dEtaInThresholds[1]) 	return false;
        if (fabs(cms2.els_dPhiIn()[index]) > dPhiInThresholds[1]) 	return false;
        if (cms2.els_hOverE()[index] > hoeThresholds[1]) 		    return false;
        if (cms2.els_sigmaIEtaIEta()[index] > latThresholds[1])  return false;	
    }

    return true;

}

//
// class based electron id that we have used before
// this was optimised in the 2X series of CMSSW
// and should be considered depracated
//
bool electronId_classBasedLoose(const unsigned int index)
{
    if (cms2.els_egamma_looseId()[index]) return true;
    return false;
}

bool electronId_classBasedTight(const unsigned int index)
{
    if (cms2.els_egamma_tightId()[index]) return true;
    return false;
}

//
// class based id that is new/experimental
//

unsigned int electronId_CIC(const unsigned int index, const cic_tightness tightness)
{

    //
    // set the parameters for the chosen tightness
    //

    std::vector<double> cutdeta;
    std::vector<double> cutdphi;
    std::vector<double> cuteopin;
    std::vector<double> cutet;
    std::vector<double> cuthoe;
    std::vector<double> cutip;
    std::vector<double> cutisoecal;
    std::vector<double> cutisohcal;
    std::vector<double> cutisotk;
    std::vector<double> cutmishits;
    std::vector<double> cutsee;

    eidGetCIC(tightness, cutdeta, cutdphi, cuteopin, cutet, cuthoe, 
            cutip, cutisoecal, cutisohcal, cutisotk, cutmishits, cutsee);

    //
    // get the variables that are going to be cut on
    //

    double scTheta = (2*atan(exp(-1*cms2.els_etaSC()[index])));
    double scEt = cms2.els_eSC()[index]*sin(scTheta); 
    //double eta = cms2.els_p4()[index].Eta(); // not used
    // this is eSuperClusterOverP
    //double eOverP = cms2.els_eOverPIn()[index];

    //    double eSeed = cms2.els_eSeed()[index];
    //    double pin  = cms2.els_p4In()[index].R();   
    //    double pout = cms2.els_p4Out()[index].R(); 
    //    double eSeedOverPin = eSeed/pin; 
    //    double fBrem = (pin-pout)/pin;
    double eSeedOverPin = cms2.els_eSeedOverPIn()[index];
    double fBrem = cms2.els_fbrem()[index];


    double hOverE = cms2.els_hOverE()[index];
    double sigmaee = cms2.els_sigmaIEtaIEta()[index];
    //double e25Max = cms2.els_e2x5Max()[index];
    //double e15 = cms2.els_e1x5()[index];
    //double e55 = cms2.els_e5x5()[index];
    //double e25Maxoe55 = e25Max/e55 ;
    //double e15oe55 = e15/e55 ;
    double deltaPhiIn = cms2.els_dPhiIn()[index];
    double deltaEtaIn = cms2.els_dEtaIn()[index];

    double ip = 0;
    //int mishits = electron->gsfTrack()->trackerExpectedHitsInner().numberOfHits();
    int mishits = cms2.els_exp_innerlayers()[index];
    double tkIso = cms2.els_tkIso()[index];
    double ecalIso = cms2.els_ecalIso04()[index];
    double hcalIso = cms2.els_hcalIso04()[index];
    //double hcalIso1 = cms2.els_hcalDepth1TowerSumEt04()[index];
    //double hcalIso2 = cms2.els_hcalDepth2TowerSumEt04()[index];

    ip = fabs(cms2.els_d0corr()[index]);

    // should be sigmaIEtaIEta with respect to the SC...
    // (in the EE it is with respect to the seed BC...)
    //if (cms2.els_fiduciality()[index] & (1<<ISEB))
    //    sigmaee = cms2.els_sigmaIEtaIEtaSC()[index];

    // version for classify is 2
    unsigned int cat = classify(2, index);

    bool isEB = (cms2.els_fiduciality()[index] & (1<<ISEB));
    int eb;
    if (isEB) 
        eb = 0;
    else 
        eb = 1; 

    unsigned int result = 0;
    int bin = 0;

    if (scEt < 20.)
        bin = 2;
    else if (scEt > 30.)
        bin = 0;
    else
        bin = 1;

    if (fBrem > 0)
        eSeedOverPin = eSeedOverPin + fBrem;

    if (bin != 2) {     
        tkIso = tkIso*pow(40./scEt, 2); 
        ecalIso = ecalIso*pow(40./scEt, 2); 
        hcalIso = hcalIso*pow(40./scEt, 2); 
    }

    if ((tkIso < cutisotk[cat+3*eb+bin*6]) &&
            (ecalIso < cutisoecal[cat+3*eb+bin*6]) &&
            (hcalIso < cutisohcal[cat+3*eb+bin*6]))
        result |= (1<<ELEID_ISO);

    if (fBrem < -2)
        return result;

    if (hOverE < cuthoe[cat+3*eb+bin*6] &&
            sigmaee < cutsee[cat+3*eb+bin*6] &&
            fabs(deltaPhiIn) < cutdphi[cat+3*eb+bin*6] &&
            fabs(deltaEtaIn) < cutdeta[cat+3*eb+bin*6] &&
            eSeedOverPin > cuteopin[cat+3*eb+bin*6])
        result |= (1<<ELEID_ID);

    if (ip < cutip[cat+3*eb+bin*6])
        result |= (1<<ELEID_IP);

    if (mishits < cutmishits[cat+3*eb+bin*6])
        result |= (1<<ELEID_CONV);

    return result;

}

unsigned int classify(const unsigned int version, const unsigned int index) {

    // this is eSuperClusterOverP
    double eOverP = cms2.els_eOverPIn()[index];
    double fBrem = cms2.els_fbrem()[index];
    bool isEB = (cms2.els_fiduciality()[index] & (1<<ISEB));
    bool isEE = (cms2.els_fiduciality()[index] & (1<<ISEE));

    int cat = -1;
    if (version == 0 || version == 1) {
        if((isEB && fBrem<0.06) || (isEE && fBrem<0.1)) 
            cat=1;
        else if (eOverP < 1.2 && eOverP > 0.8) 
            cat=0;
        else 
            cat=2;
        return cat;
    } else {
        if (isEB) {       // BARREL
            if(fBrem < 0.12)
                cat=1;
            else if (eOverP < 1.2 && eOverP > 0.9) 
                cat=0;
            else 
                cat=2;
        } else {                     // ENDCAP
            if(fBrem < 0.2)
                cat=1;
            else if (eOverP < 1.22 && eOverP > 0.82) 
                cat=0;
            else 
                cat=2;
        }
        return cat;
    }
}

//
// VBTF stuff
//

unsigned int electronId_VBTF(const unsigned int index, const vbtf_tightness tightness)
{

    unsigned int answer = 0;

    std::vector<double> relisoThresholds;
    std::vector<double> dEtaInThresholds;
    std::vector<double> dPhiInThresholds;
    std::vector<double> hoeThresholds;
    std::vector<double> sigmaIEtaIEtaThresholds;

    eidGetVBTF(tightness, dEtaInThresholds, dPhiInThresholds, hoeThresholds, 
                sigmaIEtaIEtaThresholds, relisoThresholds);

    // barrel
    if (fabs(cms2.els_etaSC()[index]) < 1.479) {

        if (electronIsolation_rel(index, true) < relisoThresholds[0])
            answer |= (1<<ELEID_ISO);

        if (fabs(cms2.els_dEtaIn()[index]) < dEtaInThresholds[0] &&
                fabs(cms2.els_dPhiIn()[index]) < dPhiInThresholds[0] &&
                cms2.els_hOverE()[index] < hoeThresholds[0] &&
                cms2.els_sigmaIEtaIEta()[index] < sigmaIEtaIEtaThresholds[0])
            answer |= (1<<ELEID_ID);
    }

    // endcap
    if (fabs(cms2.els_etaSC()[index]) > 1.479) {
        if (electronIsolation_rel(index, true) < relisoThresholds[1])
            answer |= (1<<ELEID_ISO);

        if (fabs(cms2.els_dEtaIn()[index]) < dEtaInThresholds[1] &&
                fabs(cms2.els_dPhiIn()[index]) < dPhiInThresholds[1] &&
                cms2.els_hOverE()[index] < hoeThresholds[1] &&
                cms2.els_sigmaIEtaIEta()[index] < sigmaIEtaIEtaThresholds[1])
            answer |= (1<<ELEID_ID);
    }

    return answer;

}

//
// electron isolation definitions 
//
float electronIsolation_rel(const unsigned int index, bool use_calo_iso)
{
    float sum = cms2.els_tkIso().at(index);
    if (use_calo_iso) {
        if (fabs(cms2.els_etaSC().at(index)) > 1.479) sum += cms2.els_ecalIso().at(index);
        if (fabs(cms2.els_etaSC().at(index)) <= 1.479) sum += max(0., (cms2.els_ecalIso().at(index) -1.));
        sum += cms2.els_hcalIso().at(index);
    }
    double pt = cms2.els_p4().at(index).pt();
    return sum/max(pt, 20.);
}

//
//conversion rejection
//
bool isFromConversionHitPattern(const unsigned int index)
{
    if(cms2.els_exp_innerlayers().at(index) > 1) return true;
    return false;
}

bool isFromConversionPartnerTrack(const unsigned int index) {

    if(fabs(cms2.els_conv_dist().at(index)) < 0.02 &&
            fabs(cms2.els_conv_dcot().at(index)) < 0.02)
        return true;

    return false;

}

int getChargeUsingMajorityLogic(int elIdx, float minFracSharedHits) {

    if(cms2.els_sccharge()[elIdx]*cms2.els_trk_charge()[elIdx] > 0 || (cms2.els_trkidx()[elIdx] < 0 || cms2.els_trkshFrac()[elIdx] < minFracSharedHits))
        return cms2.els_sccharge()[elIdx];
    else 
        return  cms2.trks_charge().at(cms2.els_trkidx().at(elIdx));

}

bool isChargeFlip(int elIndex){
    //true if electron is likely to be a charge flip
    if ((cms2.els_trkidx().at(elIndex) >= 0) && (cms2.els_trk_charge().at(elIndex) != cms2.trks_charge().at(cms2.els_trkidx().at(elIndex))) ) return true;
    if ((cms2.els_trkidx().at(elIndex) < 0)  && (cms2.els_trk_charge().at(elIndex) != cms2.els_sccharge().at(elIndex))) return true;

    return false;
}
