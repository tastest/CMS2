
//
// electron selections
//

// CMS2 includes
#include "electronSelections.h"
#include "CMS2.h"

elecuts_t electronSelections_debug_;
elecuts_t electronId_debug_;

bool electronSelection_cand01(const unsigned int index)
{

    electronSelections_debug_ = 0;

    if (cms2.els_type()[index] & (1<<ISECALDRIVEN)) electronSelections_debug_ |= (1<<ELEPASS_TYPE);
    if (fabs(cms2.els_p4()[index].eta()) < 2.5)electronSelections_debug_ |= (1<<ELEPASS_FIDUCIAL);
    if (electronId_noMuon(index)) electronSelections_debug_ |= (1<<ELEPASS_NOMUON);
    if (electronId_cand01(index)) electronSelections_debug_ |= (1<<ELEPASS_ID);
    if (electronImpact_cand01(index)) electronSelections_debug_ |= (1<<ELEPASS_D0);
    if (electronIsolation_relsusy_cand1(index, true) < 0.10) electronSelections_debug_ |= (1<<ELEPASS_ISO);
    if (!isFromConversionPartnerTrack(index)) electronSelections_debug_ |= (1<<ELEPASS_NOTCONV);

    if ((electronSelections_debug_ & electronSelections_passall_) == electronSelections_passall_) return true;
    return false;
}

bool electronSelection_cand02(const unsigned int index)
{

    electronSelections_debug_ = 0;

    if (cms2.els_type()[index] & (1<<ISECALDRIVEN)) electronSelections_debug_ |= (1<<ELEPASS_TYPE);
    if (fabs(cms2.els_p4()[index].eta()) < 2.5) electronSelections_debug_ |= (1<<ELEPASS_FIDUCIAL);
    if (electronId_noMuon(index)) electronSelections_debug_ |= (1<<ELEPASS_NOMUON);
    if (electronId_cand02(index)) electronSelections_debug_ |= (1<<ELEPASS_ID);
    if (electronImpact_cand01(index)) electronSelections_debug_ |= (1<<ELEPASS_D0);
    if (electronIsolation_relsusy_cand1(index, true) < 0.10) electronSelections_debug_ |= (1<<ELEPASS_ISO);
    if (!isFromConversionPartnerTrack(index)) electronSelections_debug_ |= (1<<ELEPASS_NOTCONV);

    if ((electronSelections_debug_ & electronSelections_passall_) == electronSelections_passall_) return true;
    return false;

}

//but based on a simpler cut-based egamma ID
// with more better rejection: JuraTrackIso, susy-style (ped subtracted in EB)
// and conversion rejection
bool electronSelectionTTbar_cand01(const unsigned int index) 
{

    electronSelections_debug_ = 0;

    if (cms2.els_type()[index] & (1<<ISECALDRIVEN)) electronSelections_debug_ |= (1<<ELEPASS_TYPE);
    if (fabs(cms2.els_p4()[index].eta()) < 2.5) electronSelections_debug_ |= (1<<ELEPASS_FIDUCIAL);
    if (electronId_noMuon(index)) electronSelections_debug_ |= (1<<ELEPASS_NOMUON);
    if (electronId_cand01(index)) electronSelections_debug_ |= (1<<ELEPASS_ID); 
    if (electronImpactTTbar(index)) electronSelections_debug_ |= (1<<ELEPASS_D0);
    if (electronIsolation_relsusy_cand1(index, true) < 0.10) electronSelections_debug_ |= (1<<ELEPASS_ISO);
    if (!isFromConversionPartnerTrack(index)) electronSelections_debug_ |= (1<<ELEPASS_NOTCONV);

    if ((electronSelections_debug_ & electronSelections_passall_) == electronSelections_passall_) return true;
    return false;

}

bool electronImpactTTbar(const unsigned int index) 
{
    //
    // define thresholds for EB, EE
    //
    float d0Thresholds[2]               = {0.04, 0.04};

    //
    // apply cuts
    //
    if (fabs(cms2.els_etaSC()[index]) < 1.479) {
        if (fabs(cms2.els_d0corr()[index]) > d0Thresholds[0]) return false;
        return true;
    }
    if (fabs(cms2.els_etaSC()[index]) > 1.479) {
        if (fabs(cms2.els_d0corr()[index]) > d0Thresholds[1]) return false;
        return true;
    }
    return false;
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
bool electronId_cand01(const unsigned int index)
{

    electronId_debug_ = 0;

    //
    // define thresholds for EB, EE
    //
    float dEtaInThresholds[2]               = {0.007, 0.010};
    float dPhiInThresholds[2]               = {0.020, 0.025};
    float hoeThresholds[2]                  = {0.01, 0.01};
    float sigmaIEtaIEtaThresholds[2]        = {9999.99, 0.03};
    float e2x5Over5x5Thresholds[2]          = {0.90, 0.00};

    //
    // apply cuts
    //
    if (fabs(cms2.els_etaSC()[index]) < 1.479) {
        if (fabs(cms2.els_dEtaIn()[index]) < dEtaInThresholds[0]) 	electronId_debug_ |= (1<<ELEPASS_DETA);
        if (fabs(cms2.els_dPhiIn()[index]) < dPhiInThresholds[0]) 	electronId_debug_ |= (1<<ELEPASS_DPHI);
        if (cms2.els_hOverE()[index] < hoeThresholds[0]) 		electronId_debug_ |= (1<<ELEPASS_HOE);
        if ((cms2.els_e2x5Max()[index]/cms2.els_e5x5()[index]) > e2x5Over5x5Thresholds[0]) electronId_debug_ |= (1<<ELEPASS_LSHAPE);
    }
    if (fabs(cms2.els_etaSC()[index]) > 1.479) {
        if (fabs(cms2.els_dEtaIn()[index]) < dEtaInThresholds[1]) 	electronId_debug_ |= (1<<ELEPASS_DETA);
        if (fabs(cms2.els_dPhiIn()[index]) < dPhiInThresholds[1]) 	electronId_debug_ |= (1<<ELEPASS_DPHI);
        if (cms2.els_hOverE()[index] < hoeThresholds[1]) 		electronId_debug_ |= (1<<ELEPASS_HOE);
        if (cms2.els_sigmaIEtaIEta()[index] < sigmaIEtaIEtaThresholds[1])  electronId_debug_ |= (1<<ELEPASS_LSHAPE);	
    }

    if ((electronId_debug_ & electronSelections_passid_) == electronSelections_passid_) return true;
    return false;

}

bool electronId_cand02(const unsigned int index)
{

    electronId_debug_ = 0;

    //
    // define thresholds for EB, EE
    //
    float dEtaInThresholds[2]               = {0.005, 0.007};
    float dPhiInThresholds[2]               = {0.020, 0.025};
    float hoeThresholds[2]                  = {0.01, 0.01};
    float sigmaIEtaIEtaThresholds[2]        = {9999.99, 0.03};
    float e2x5Over5x5Thresholds[2]          = {0.94, 0.00};

    //
    // apply cuts
    //
    if (fabs(cms2.els_etaSC()[index]) < 1.479) {
        if (fabs(cms2.els_dEtaIn()[index]) < dEtaInThresholds[0])   electronId_debug_ |= (1<<ELEPASS_DETA);
        if (fabs(cms2.els_dPhiIn()[index]) < dPhiInThresholds[0])   electronId_debug_ |= (1<<ELEPASS_DPHI);
        if (cms2.els_hOverE()[index] < hoeThresholds[0])        electronId_debug_ |= (1<<ELEPASS_HOE);
        if ((cms2.els_e2x5Max()[index]/cms2.els_e5x5()[index]) > e2x5Over5x5Thresholds[0]) electronId_debug_ |= (1<<ELEPASS_LSHAPE);
    }
    if (fabs(cms2.els_etaSC()[index]) > 1.479) {
        if (fabs(cms2.els_dEtaIn()[index]) < dEtaInThresholds[1])   electronId_debug_ |= (1<<ELEPASS_DETA);
        if (fabs(cms2.els_dPhiIn()[index]) < dPhiInThresholds[1])   electronId_debug_ |= (1<<ELEPASS_DPHI);
        if (cms2.els_hOverE()[index] < hoeThresholds[1])        electronId_debug_ |= (1<<ELEPASS_HOE);
        if (cms2.els_sigmaIEtaIEta()[index] < sigmaIEtaIEtaThresholds[1])  electronId_debug_ |= (1<<ELEPASS_LSHAPE);
    }

    if ((electronId_debug_ & electronSelections_passid_) == electronSelections_passid_) return true;
    return false;

}

bool electronImpact_cand01(const unsigned int index)
{
    //
    // define thresholds for EB, EE
    //
    float d0Thresholds[2]               = {0.02, 0.02};

    //
    // apply cuts
    //
    if (fabs(cms2.els_etaSC()[index]) < 1.479) {
        if (fabs(cms2.els_d0corr()[index]) > d0Thresholds[0]) return false;
        return true;
    }
    if (fabs(cms2.els_etaSC()[index]) > 1.479) {
        if (fabs(cms2.els_d0corr()[index]) > d0Thresholds[1]) return false;
        return true;
    }

    return false;

}

//
// candidate electron isolation function
//
bool electronIsolation_cand01(const unsigned int index)
{

    //
    // define thresholds for EB, EE
    //
    //float tkThresholds[2]         =       {4.5, 6.0};
    float tkThresholds[2] 	= 	{2.5, 2.0};
    float ecalThresholds[2] = 	{2.5, 2.0};
    float hcalThresholds[2] = 	{1.0, 1.0};

    //
    // apply cuts
    //


    if (fabs(cms2.els_etaSC()[index]) < 1.479) {
        // if (cms2.els_tkIso()[index] > tkThresholds[0])    return false;
        if (cms2.els_tkJuraIso()[index] > tkThresholds[0]) 	return false;
        if (cms2.els_ecalIso()[index] 	> ecalThresholds[0]) 	return false;
        if (cms2.els_hcalIso()[index] 	> hcalThresholds[0]) 	return false;
        return true;
    }
    if (fabs(cms2.els_etaSC()[index]) > 1.479) {
        //if (cms2.els_tkIso()[index] > tkThresholds[1])      return false;
        if (cms2.els_tkJuraIso()[index] > tkThresholds[1])      return false;
        if (cms2.els_ecalIso()[index]   > ecalThresholds[1])    return false;
        if (cms2.els_hcalIso()[index]   > hcalThresholds[1])    return false;
        return true;
    }

    return false;
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

void electronId_classBasedIdAssign(std::vector<double> &cutarr, double cutvals[])
{
    cutarr.clear();
    for (unsigned int i = 0; i < 18; ++i) {
        cutarr.push_back(cutvals[i]);
    }
}

unsigned int electronId_CIC(const cic_tightness tightness, const unsigned int index)
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

    eidGet(tightness, cutdeta, cutdphi, cuteopin, cutet, cuthoe, 
            cutip, cutisoecal, cutisohcal, cutisotk, cutmishits, cutsee);

    //
    // get the variables that are going to be cut on
    //

    double scTheta = (2*atan(exp(cms2.els_etaSC()[index])));
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
        result |= (1<<ELEPASS_CIC_ISO);

    if (fBrem < -2)
        return result;

    if (hOverE < cuthoe[cat+3*eb+bin*6] &&
        sigmaee < cutsee[cat+3*eb+bin*6] &&
        fabs(deltaPhiIn) < cutdphi[cat+3*eb+bin*6] &&
        fabs(deltaEtaIn) < cutdeta[cat+3*eb+bin*6] &&
        eSeedOverPin > cuteopin[cat+3*eb+bin*6])
        result |= (1<<ELEPASS_CIC_ID);

    if (ip < cutip[cat+3*eb+bin*6])
        result |= (1<<ELEPASS_CIC_IP);

    if (mishits < cutmishits[cat+3*eb+bin*6])
        result |= (1<<ELEPASS_CIC_CONV);

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

//WP70 (70%)
//=======
//EB
//--
//track_iso  2.5 GeV
//ecal_iso   3.0 GeV
//hcal_iso   5.0 GeV
//sihih      0.01
//Dphi@vtx   0.02
//Deta@vtx   0.006
//H/E        0.02

//EE
//--
//track_iso  0.8 GeV
//ecal_iso   2.5 GeV
//hcal_iso   0.25 GeV
//sihih      0.03
//Dphi@vtx   0.02
//Deta@vtx   0.003
//H/E        0.0025

unsigned int electronId_VBTF70(const unsigned int index)
{

    unsigned int answer = 0;
    float tkThresholds[2]   =   {2.5, 0.8};
    float ecalThresholds[2] =   {3.0, 2.5};
    float hcalThresholds[2] =   {5.0, 0.25};
    float dEtaInThresholds[2]               = {0.006, 0.003};
    float dPhiInThresholds[2]               = {0.020, 0.020};
    float hoeThresholds[2]                  = {0.02, 0.0025};
    float sigmaIEtaIEtaThresholds[2]        = {0.01, 0.03};

    // barrel
    if (fabs(cms2.els_etaSC()[index]) < 1.479) {

        if (cms2.els_tkJuraIso()[index] < tkThresholds[0] &&
            cms2.els_ecalIso()[index]   < ecalThresholds[0] &&
            cms2.els_hcalIso()[index]   < hcalThresholds[0])
            answer |= (1<<1);

        if (fabs(cms2.els_dEtaIn()[index]) < dEtaInThresholds[0] &&
            fabs(cms2.els_dPhiIn()[index]) < dPhiInThresholds[0] &&
            cms2.els_hOverE()[index] < hoeThresholds[0] &&
            cms2.els_sigmaIEtaIEta()[index] < sigmaIEtaIEtaThresholds[0])
            answer |= (1<<0);
    }

    // endcap
    if (fabs(cms2.els_etaSC()[index]) > 1.479) {
        if (cms2.els_tkJuraIso()[index] < tkThresholds[1] &&
            cms2.els_ecalIso()[index]   < ecalThresholds[1] &&
            cms2.els_hcalIso()[index]   < hcalThresholds[1])
            answer |= (1<<1);

        if (fabs(cms2.els_dEtaIn()[index]) < dEtaInThresholds[1] &&
            fabs(cms2.els_dPhiIn()[index]) < dPhiInThresholds[1] &&
            cms2.els_hOverE()[index] < hoeThresholds[1] &&
            cms2.els_sigmaIEtaIEta()[index] < sigmaIEtaIEtaThresholds[1])
            answer |= (1<<0);
    }

    return answer;

}



//
//
//




//
// electron isolation definitions that we have used before
//

float electronIsolation_rel(const unsigned int index, bool use_calo_iso)
{
    float sum = cms2.els_tkIso().at(index);
    if (use_calo_iso) {
        sum += cms2.els_ecalIso()[index];
        sum += cms2.els_hcalIso()[index];
    }
    float pt = cms2.els_p4().at(index).pt();
    return pt / (pt + sum + 1e-5);
}

float electronIsolation_relsusy(const unsigned int index, bool use_calo_iso)
{
    float sum = cms2.els_tkIso().at(index);
    if (use_calo_iso) {
        sum += max(0., (cms2.els_ecalIso().at(index) -2.));
        sum += cms2.els_hcalIso().at(index);
    }
    double pt = cms2.els_p4().at(index).pt();
    return sum/max(pt, 20.);
}

float electronIsolation_relsusy_cand0(const unsigned int index, bool use_calo_iso)
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

float electronIsolation_relsusy_cand1(const unsigned int index, bool use_calo_iso)
{
    float sum = cms2.els_tkJuraIso().at(index);
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
