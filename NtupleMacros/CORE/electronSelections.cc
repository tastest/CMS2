
//
// electron selections
//

// CMS2 includes
#include "electronSelections.h"
#include "CMS2.h"


bool pass_electronSelectionCompareMask(const cuts_t cuts_passed, const cuts_t selectionType)
{
    if ((cuts_passed & selectionType) == selectionType) return true;
    return false;
}

bool pass_electronSelection(const unsigned int index, const cuts_t selectionType, bool applyAlignmentCorrection)
{
    cuts_t cuts_passed = electronSelection(index, applyAlignmentCorrection);
    if ((cuts_passed & selectionType) == selectionType) return true;
    return false;
}

cuts_t electronSelection(const unsigned int index, bool applyAlignmentCorrection) 
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

    // 
    // "CAND" ID
    //
    // CAND01
    if (electronId_cand(index, CAND_01, applyAlignmentCorrection)) cuts_passed |= (1ll<<ELEID_CAND01);
    // CAND02
    if (electronId_cand(index, CAND_02, applyAlignmentCorrection)) cuts_passed |= (1ll<<ELEID_CAND02);
    // "EXTRA"
    if (electronId_extra(index)) cuts_passed |= (1ll<<ELEID_EXTRA);
    //
    // VBTF ID
    //
    electronIdComponent_t answer_vbtf = 0;
    // VBTF90 (optimised in 35X)
    answer_vbtf = electronId_VBTF(index, VBTF_35X_90, applyAlignmentCorrection);
    if ((answer_vbtf & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) cuts_passed |= (1ll<<ELEID_VBTF_35X_90);
    // VBTF70 (optimised in 35X)
    answer_vbtf = electronId_VBTF(index, VBTF_35X_80, applyAlignmentCorrection);
    if ((answer_vbtf & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) cuts_passed |= (1ll<<ELEID_VBTF_35X_80);
    // VBTF70 (optimised in 35X)
    answer_vbtf = electronId_VBTF(index, VBTF_35Xr2_70, applyAlignmentCorrection);
    if ((answer_vbtf & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) cuts_passed |= (1ll<<ELEID_VBTF_35X_70);
    //
    // CIC ID  
    //
    // MEDIUM (V03 optimisation)
    electronIdComponent_t answer_cic = electronId_CIC(index, 3, CIC_MEDIUM, applyAlignmentCorrection);
    if ((answer_cic & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) cuts_passed |= (1ll<<ELEID_CIC_V03_MEDIUM);

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
    // Pt
    //
    if( cms2.els_p4()[index].pt() > 10.0 ) cuts_passed |= (1ll<<ELEPT_010);

    //
    // Super Cluster Et
    //
    if( ( cms2.els_eSC()[index] / cosh(cms2.els_etaSC()[index]) ) > 10.0 ) cuts_passed |= (1ll<<ELESCET_010);

    //
    // chargeflip
    //

    if (!isChargeFlip(index)) cuts_passed |= (1ll<<ELECHARGE_NOTFLIP);

    //
    // spike rejection
    //
    if (!isSpikeElectron(index)) cuts_passed |= (1ll<<ELENOSPIKE_SWISS005);

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
bool electronId_cand(const unsigned int index, const cand_tightness tightness, bool applyAlignementCorrection)
{

    std::vector<double> relisoThresholds;
    std::vector<double> dEtaInThresholds;
    std::vector<double> dPhiInThresholds;
    std::vector<double> hoeThresholds;
    std::vector<double> latThresholds;

    eidGetCand(tightness, dEtaInThresholds, dPhiInThresholds, hoeThresholds, latThresholds);

    //
    // get corrected dEtaIn and dPhiIn
    //

    float dEtaIn = cms2.els_dEtaIn()[index];
    float dPhiIn = cms2.els_dPhiIn()[index];
    if (applyAlignementCorrection) electronCorrection_pos(index, dEtaIn, dPhiIn);

    //
    // apply cuts
    //
    if (fabs(cms2.els_etaSC()[index]) < 1.479) {
        if (fabs(dEtaIn) > dEtaInThresholds[0]) 	return false;
        if (fabs(dPhiIn) > dPhiInThresholds[0]) 	return false;
        if (cms2.els_hOverE()[index] > hoeThresholds[0]) 		    return false;
        if ((cms2.els_e2x5Max()[index]/cms2.els_e5x5()[index]) < latThresholds[0]) return false;
    }
    if (fabs(cms2.els_etaSC()[index]) > 1.479) {
        if (fabs(dEtaIn) > dEtaInThresholds[1]) 	return false;
        if (fabs(dPhiIn) > dPhiInThresholds[1]) 	return false;
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

electronIdComponent_t electronId_CIC(const unsigned int index, const unsigned int version, const cic_tightness tightness, bool applyAlignementCorrection)
{

    // check that a valid version number was supplied
    if (version != 2 && version != 3 && version != 4) {
        std::cout << "[electronId_CIC] Error! Version must be 2, 3 or 4 - fail" << std::endl;
        return 0;
    }

    //
    // get the variables that are going to be cut on
    //

    double scTheta = (2*atan(exp(-1*cms2.els_etaSC()[index])));
    double scEt = cms2.els_eSC()[index]*sin(scTheta); 
    double eSeedOverPin = cms2.els_eSeedOverPIn()[index];
    double fBrem = cms2.els_fbrem()[index];
    double hOverE = cms2.els_hOverE()[index];
    double sigmaee = cms2.els_sigmaIEtaIEta()[index];
    int mishits = cms2.els_exp_innerlayers()[index];
    double dist = cms2.els_conv_dist()[index];
    double dcot = cms2.els_conv_dcot()[index];
    double tkIso = cms2.els_tkIso()[index];
    double ecalIso = cms2.els_ecalIso04()[index];
    double hcalIso = cms2.els_hcalIso04()[index];
    double ip = fabs(cms2.els_d0corr()[index]);
    // should be sigmaIEtaIEta with respect to the SC...
    // (in the EE it is with respect to the seed BC...)
    if (cms2.els_fiduciality()[index] & (1<<ISEB))
        sigmaee = cms2.els_sigmaIEtaIEtaSC()[index];

    //
    // get corrected dEtaIn and dPhiIn
    //
    
    float deltaEtaIn = cms2.els_dEtaIn()[index];
    float deltaPhiIn = cms2.els_dPhiIn()[index];
    if (applyAlignementCorrection) electronCorrection_pos(index, deltaEtaIn, deltaPhiIn);

    // find the catagory for this electron
    unsigned int cat = classify(version, index);

    // determine if in EB or EE
    int eb;
    if (cms2.els_fiduciality()[index] & (1<<ISEB)) eb = 0;
    else eb = 1;

    //
    // Version V02
    //

    if (version == 2) {

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

       eidGetCIC_V02(tightness, cutdeta, cutdphi, cuteopin, cutet, cuthoe, 
                cutip, cutisoecal, cutisohcal, cutisotk, cutmishits, cutsee);


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

    //
    // version V03, V04 or V05
    //

    if (version == 3 || version == 4 || version == 5) {

        //
        // set the parameters for the chosen tightness
        //
        std::vector<double> cutdcotdist;
        std::vector<double> cutdetain;
        std::vector<double> cutdphiin;
        std::vector<double> cuteseedopcor;
        std::vector<double> cutet;
        std::vector<double> cutfmishits;
        std::vector<double> cuthoe;
        std::vector<double> cutip_gsf;
        std::vector<double> cutiso_sum;
        std::vector<double> cutiso_sumoet;
        std::vector<double> cutsee;

        // V03 uses Et binning
        bool wantBinning = true;

        if (version == 3) {
            eidGetCIC_V03(tightness, cutdcotdist, cutdetain, cutdphiin, cuteseedopcor, cutet,
                    cutfmishits, cuthoe, cutip_gsf, cutiso_sum, cutiso_sumoet, cutsee);
        }

        if (version == 4) {
            eidGetCIC_V04(tightness, cutdcotdist, cutdetain, cutdphiin, cuteseedopcor, cutet,
                    cutfmishits, cuthoe, cutip_gsf, cutiso_sum, cutiso_sumoet, cutsee);
            // V04 does not use Et binning
            wantBinning = false;
        }

        // this is certainly true for V03
        // but not sure the meaning of V04,05
        unsigned int result = 0;
        int bin = 0;
        if (wantBinning) {
            if (scEt < 20.)
                bin = 2;
            else if (scEt > 30.)
                bin = 0;
            else
                bin = 1;
        }

        if (fBrem > 0)
            eSeedOverPin = eSeedOverPin + fBrem;

        float iso_sum = tkIso + ecalIso + hcalIso;
        float iso_sum_corrected = iso_sum*pow(40./scEt, 2);

        if ((iso_sum < cutiso_sum[cat+bin*9]) &&
                (iso_sum_corrected < cutiso_sumoet[cat+bin*9]))
            result |= (1<<ELEID_ISO);

        if (fBrem > -2) {
            if ((hOverE < cuthoe[cat+bin*9]) and
                    (sigmaee < cutsee[cat+bin*9]) and
                    (fabs(deltaPhiIn) < cutdphiin[cat+bin*9]) and
                    (fabs(deltaEtaIn) < cutdetain[cat+bin*9]) and
                    (eSeedOverPin > cuteseedopcor[cat+bin*9]) and
                    (scEt > cutet[cat+bin*9]))
                result |= (1<<ELEID_ID);
        }

        if (ip < cutip_gsf[cat+bin*9])
            result |= (1<<ELEID_IP);

        float dcotdistcomb = ((0.4 - std::max(dist, dcot)) > 0?(0.4 - std::max(dist, dcot)):0);
        if ((mishits < cutfmishits[cat+bin*9]) and 
                (dcotdistcomb < cutdcotdist[cat+bin*9]))
            result |= (1<<ELEID_CONV);

        return result;
    }

    std::cout << "[electronId_CIC] Error! got to the end and didn't apply any cuts - fail" << std::endl;
    return 0;

}

unsigned int classify(const unsigned int version, const unsigned int index) {

    // variables used for classification
    double eta = fabs(cms2.els_etaSC()[index]);
    double eOverP = cms2.els_eOverPIn()[index];
    double fBrem = cms2.els_fbrem()[index];
    bool isEB = (cms2.els_fiduciality()[index] & (1<<ISEB));
    bool isEE = (cms2.els_fiduciality()[index] & (1<<ISEE));
    bool ecalDriven = (cms2.els_type()[index] & (1<<ISECALDRIVEN));
    bool trackerDriven = (cms2.els_type()[index] & (1<<ISTRACKERDRIVEN));

    int cat = -1;

    //
    // version V00 or V01
    //
    if (version == 0 || version == 1) {
        if((isEB && fBrem<0.06) || (isEE && fBrem<0.1)) 
            cat=1;
        else if (eOverP < 1.2 && eOverP > 0.8) 
            cat=0;
        else 
            cat=2;
        return cat;
    }

    //
    // version V02
    //
    if (version == 2) {
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

    //
    // version V03, V04 or V05
    // took this from V00-03-07-02 of 
    // ElectronIdentification/src/CutBasedElectronID.cc
    //

    if (version == 3 || version == 4 || version == 5) {

        // this is certainly true for V03 and V04
        // not sure about V05
        bool newCategories = false;

        if (isEB) {
            if ((fBrem >= 0.12) && (eOverP > 0.9) && (eOverP < 1.2))
                cat = 0;
            else if (((eta >  .445   && eta <  .45  ) ||
                        (eta >  .79    && eta <  .81  ) ||
                        (eta > 1.137   && eta < 1.157 ) ||
                        (eta > 1.47285 && eta < 1.4744)) && newCategories)
                cat = 6;
            else if (trackerDriven && !ecalDriven && newCategories)
                cat = 8;
            else if (fBrem < 0.12)
                cat = 1;
            else
                cat = 2;
        } else {
            if ((fBrem >= 0.2) && (eOverP > 0.82) && (eOverP < 1.22))
                cat = 3;
            else if (eta > 1.5 && eta <  1.58 && newCategories)
                cat = 7;
            else if (trackerDriven && !ecalDriven && newCategories)
                cat = 8;
            else if (fBrem < 0.2)
                cat = 4;
            else
                cat = 5;
        }

        return cat;
    }


    return cat;

}

//
// VBTF stuff
//

electronIdComponent_t electronId_VBTF(const unsigned int index, const vbtf_tightness tightness, bool applyAlignementCorrection)
{

    unsigned int answer = 0;

    std::vector<double> relisoThresholds;
    std::vector<double> dEtaInThresholds;
    std::vector<double> dPhiInThresholds;
    std::vector<double> hoeThresholds;
    std::vector<double> sigmaIEtaIEtaThresholds;

    eidGetVBTF(tightness, dEtaInThresholds, dPhiInThresholds, hoeThresholds, 
            sigmaIEtaIEtaThresholds, relisoThresholds);

    //
    // get corrected dEtaIn and dPhiIn
    //

    float dEtaIn = cms2.els_dEtaIn()[index];
    float dPhiIn = cms2.els_dPhiIn()[index];
    if (applyAlignementCorrection) electronCorrection_pos(index, dEtaIn, dPhiIn);

    // barrel
    if (fabs(cms2.els_etaSC()[index]) < 1.479) {

        if (electronIsolation_rel(index, true) < relisoThresholds[0])
            answer |= (1<<ELEID_ISO);

        if (fabs(dEtaIn) < dEtaInThresholds[0] &&
                fabs(dPhiIn) < dPhiInThresholds[0] &&
                cms2.els_hOverE()[index] < hoeThresholds[0] &&
                cms2.els_sigmaIEtaIEta()[index] < sigmaIEtaIEtaThresholds[0])
            answer |= (1<<ELEID_ID);
    }

    // endcap
    if (fabs(cms2.els_etaSC()[index]) > 1.479) {
        if (electronIsolation_rel(index, true) < relisoThresholds[1])
            answer |= (1<<ELEID_ISO);

        if (fabs(dEtaIn) < dEtaInThresholds[1] &&
                fabs(dPhiIn) < dPhiInThresholds[1] &&
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

//
// spike rejection for electrons
//

bool isSpikeElectron(const unsigned int index) {

    const int scidx = cms2.els_scindex()[index];
    bool isSpike = false;
    if (scidx != -1) {
        //subtract twice max since max is in both 1x3 and 3x1, and we want neither
        const float r4 = (cms2.scs_e1x3()[scidx] + cms2.scs_e3x1()[scidx] - 2*cms2.scs_eMax()[scidx])/cms2.scs_eMax()[scidx];
        if (r4 < 0.05) isSpike = true;
    }

    return isSpike;

}

//
// position correction for electrons
//

void electronCorrection_pos(const unsigned int index, float &dEtaIn, float &dPhiIn)
{

    //
    // uncorrected dEtaIn and dPhiIn
    //

    dEtaIn = cms2.els_dEtaIn()[index];
    dPhiIn = cms2.els_dPhiIn()[index];

    //
    // if configered not to apply correction
    // or in barrel or no valid super cluster
    // return uncorrected values
    //

    if (!(cms2.els_fiduciality()[index]) & (1<<ISEE)) return;
    if (cms2.els_scindex()[index] == -1) return;

    //
    // set up correction parameters for EE+ and EE-
    // RecoEgamma/EgammaTools/python/correctedElectronsProducer_cfi.py?revision=1.2
    //

    //                                      X',     Y',     Z'
    float scPositionCorrectionEEP[3] = {   0.52,   -0.81,  0.81};
    float scPositionCorrectionEEM[3] = {    -0.02,  -0.81,  -0.94};

    LorentzVector initial_pos = cms2.scs_pos_p4()[cms2.els_scindex()[index]];
    LorentzVector corrected_pos;

    //
    // work out corrected position
    //

    if (cms2.els_etaSC()[index] < 0) {
            corrected_pos = LorentzVector(  initial_pos.x() + scPositionCorrectionEEM[0],
                                            initial_pos.y() + scPositionCorrectionEEM[1],
                                            initial_pos.z() + scPositionCorrectionEEM[2], 0.0);
    }
    if (cms2.els_etaSC()[index] > 0) {
            corrected_pos = LorentzVector(  initial_pos.x() + scPositionCorrectionEEP[0],
                                            initial_pos.y() + scPositionCorrectionEEP[1],
                                            initial_pos.z() + scPositionCorrectionEEP[2], 0.0);
    }

    //
    // work out correction to dEtaIn and dPhiIn
    //

    float deta_sc = corrected_pos.Eta() - initial_pos.Eta();
    float dphi_sc = acos(cos(corrected_pos.Phi() - initial_pos.Phi()));
    dEtaIn = deta_sc + cms2.els_dEtaIn()[index];
    dPhiIn = acos(cos(dphi_sc + cms2.els_dPhiIn()[index]));

}

