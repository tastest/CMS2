
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

unsigned int electronId_classBasedId(const unsigned int tightness, const unsigned int index)
{

    if (tightness != 0 && tightness != 1) {
        std::cout << "[electronId_classBasedExperimental] ERROR - invalid tightness" << std::endl;
        return 0;
    }

    //
    // set the parameters for the chosen tightness
    //

    std::vector<double> cutdeta;
    std::vector<double> cutdphi;
    std::vector<double> cuteopin;
    std::vector<double> cuthoe;
    std::vector<double> cutip;
    std::vector<double> cutisoecal;
    std::vector<double> cutisohcal;
    std::vector<double> cutisotk;
    std::vector<double> cutmishits;
    std::vector<double> cutsee;

    // loose
    if (tightness == 0) {

        double cutdeta_tmp[18] = { 
            9.58e-03, 4.06e-03, 1.22e-02, 1.37e-02, 8.37e-03, 1.27e-02,    
            1.10e-02, 3.36e-03, 9.77e-03, 1.50e-02, 6.75e-03, 1.09e-02,    
            1.40e-02, 5.08e-03, 1.09e-02, 1.46e-02, 5.06e-03, 1.27e-02 };
        double cutdphi_tmp[18] = {
            3.72e-02, 1.14e-01, 1.18e-01, 4.88e-02, 1.17e-01, 1.19e-01,
            6.06e-02, 5.48e-02, 1.17e-01, 7.00e-02, 3.55e-02, 1.17e-01,
            8.80e-02, 4.50e-02, 1.18e-01, 9.19e-02, 2.36e-02, 5.15e-02 };
        double cuteopin_tmp[18] = {
            8.78e-01, 8.02e-01, 8.14e-01, 9.42e-01, 7.35e-01, 7.74e-01,
            8.29e-01, 9.09e-01, 8.29e-01, 8.13e-01, 8.60e-01, 8.97e-01,
            8.17e-01, 8.31e-01, 8.18e-01, 8.61e-01, 7.87e-01, 7.89e-01 };
        double cuthoe_tmp[18] = {
            8.87e-02, 9.34e-02, 9.49e-02, 9.86e-02, 4.31e-02, 8.78e-02,
            9.70e-02, 5.09e-02, 9.80e-02, 9.91e-02, 3.21e-02, 9.28e-02,
            6.63e-02, 7.17e-02, 9.66e-02, 7.58e-02, 1.49e-02, 1.31e-02 };
        double cutip_tmp[18] = {
            2.46e-02, 7.60e-02, 9.66e-02, 8.85e-02, 4.41e-01, 2.05e-01,
            2.92e-02, 2.93e-02, 6.19e-02, 2.51e-02, 1.59e-01, 8.15e-02,
            7.29e+00, 1.06e-02, 5.76e+00, 6.89e+00, 1.27e+00, 5.89e+00 };
        double cutisoecal_tmp[18] = {
            3.34e+01, 2.81e+01, 7.32e+00, 2.74e+01, 7.33e+00, 2.17e+01,
            9.38e+01, 1.02e+02, 1.21e+01, 2.60e+01, 8.91e+00, 1.00e+01,
            1.61e+01, 3.13e+01, 1.69e+01, 1.54e+01, 1.33e+01, 3.77e+01 };
        double cutisohcal_tmp[18] = {
            1.35e+01, 9.93e+00, 7.56e+00, 1.48e+01, 8.10e+00, 1.08e+01,
            4.27e+01, 2.01e+01, 9.11e+00, 1.04e+01, 6.89e+00, 5.59e+00,
            8.53e+00, 9.59e+00, 2.42e+01, 2.78e+00, 8.67e+00, 2.88e-01 };
        double cutisotk_tmp[18] = {
            2.43e+01, 8.45e+00, 1.44e+01, 2.78e+01, 6.02e+00, 1.05e+01,
            1.41e+01, 1.02e+01, 1.45e+01, 1.91e+01, 6.10e+00, 1.41e+01,
            8.59e+00, 8.33e+00, 8.30e+00, 8.93e+00, 8.60e+00, 1.60e+01 };
        double cutmishits_tmp[18] = {
            5.50e+00, 1.50e+00, 5.50e+00, 2.50e+00, 2.50e+00, 2.50e+00,
            3.50e+00, 5.50e+00, 5.00e-01, 1.50e+00, 2.50e+00, 5.00e-01,
            1.50e+00, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01 };
        double cutsee_tmp[18] = {
            1.72e-02, 1.15e-02, 1.43e-02, 3.44e-02, 2.95e-02, 3.04e-02,
            1.45e-02, 1.08e-02, 1.28e-02, 3.47e-02, 3.07e-02, 3.16e-02,
            1.80e-02, 1.10e-02, 1.32e-02, 3.49e-02, 3.10e-02, 3.27e-02 };

        electronId_classBasedExperimentalAssign(cutdeta, cutdeta_tmp);
        electronId_classBasedExperimentalAssign(cutdphi, cutdphi_tmp);
        electronId_classBasedExperimentalAssign(cuteopin, cuteopin_tmp);
        electronId_classBasedExperimentalAssign(cuthoe, cuthoe_tmp);
        electronId_classBasedExperimentalAssign(cutip, cutip_tmp);
        electronId_classBasedExperimentalAssign(cutisoecal, cutisoecal_tmp);
        electronId_classBasedExperimentalAssign(cutisohcal, cutisohcal_tmp);
        electronId_classBasedExperimentalAssign(cutisotk, cutisotk_tmp);
        electronId_classBasedExperimentalAssign(cutmishits, cutmishits_tmp);
        electronId_classBasedExperimentalAssign(cutsee, cutsee_tmp);
    }

    if (tightness == 1) {
        double cutdeta_tmp[18] = {
            9.15e-03, 3.02e-03, 6.10e-03, 1.35e-02, 5.65e-03, 7.93e-03,
            1.02e-02, 2.66e-03, 1.06e-02, 9.03e-03, 7.66e-03, 7.23e-03,
            1.16e-02, 2.03e-03, 6.59e-03, 1.48e-02, 5.55e-03, 1.28e-02};
        double cutdphi_tmp[18] = {
            3.69e-02, 3.07e-02, 1.17e-01, 4.75e-02, 2.16e-02, 1.17e-01,
            3.72e-02, 2.46e-02, 4.26e-02, 6.12e-02, 1.42e-02, 3.90e-02,
            7.37e-02, 5.66e-02, 3.59e-02, 1.87e-02, 1.20e-02, 3.58e-02};
        double cuteopin_tmp[18] = {
            8.78e-01, 8.59e-01, 8.74e-01, 9.44e-01, 7.37e-01, 7.73e-01,
            8.60e-01, 9.67e-01, 9.17e-01, 8.12e-01, 9.15e-01, 1.01e+00,
            8.47e-01, 9.53e-01, 9.79e-01, 8.41e-01, 7.71e-01, 1.09e+00};
        double cuthoe_tmp[18] = {
            8.71e-02, 2.89e-02, 7.83e-02, 9.46e-02, 2.45e-02, 3.63e-02,
            6.71e-02, 4.80e-02, 6.14e-02, 9.24e-02, 1.58e-02, 4.90e-02,
            3.82e-02, 9.15e-02, 4.51e-02, 4.52e-02, 1.96e-03, 4.30e-03};
        double cutip_tmp[18] = {
            2.39e-02, 2.70e-02, 7.68e-02, 2.31e-02, 1.78e-01, 9.57e-02,
            1.02e-02, 1.68e-02, 4.30e-02, 1.66e-02, 5.94e-02, 3.08e-02,
            2.10e+00, 5.27e-03, 3.17e+00, 4.91e+00, 7.69e-01, 5.90e+00};
        double cutisoecal_tmp[18] = {
            2.00e+01, 2.72e+01, 4.48e+00, 1.35e+01, 4.56e+00, 3.19e+00,
            1.22e+01, 1.31e+01, 7.42e+00, 7.67e+00, 4.12e+00, 4.85e+00,
            1.01e+01, 1.24e+01, 1.11e+01, 1.10e+01, 1.06e+01, 1.34e+01};
        double cutisohcal_tmp[18] = {
            1.09e+01, 7.01e+00, 8.75e+00, 3.51e+00, 7.75e+00, 1.62e+00,
            1.16e+01, 9.90e+00, 4.97e+00, 5.33e+00, 3.18e+00, 2.32e+00,
            1.64e-01, 5.46e+00, 1.20e+01, 6.04e-03, 4.10e+00, 6.28e-04};
        double cutisotk_tmp[18] = {
            6.53e+00, 4.60e+00, 6.00e+00, 8.63e+00, 3.11e+00, 7.77e+00,
            5.42e+00, 4.81e+00, 4.06e+00, 6.47e+00, 2.80e+00, 3.45e+00,
            5.29e+00, 5.18e+00, 1.54e+01, 5.38e+00, 4.47e+00, 3.47e-02};
        double cutmishits_tmp[18] = {
            5.50e+00, 1.50e+00, 5.00e-01, 1.50e+00, 2.50e+00, 5.00e-01,
            3.50e+00, 5.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01,
            5.00e-01, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01};
        double cutsee_tmp[18] = {
            1.31e-02, 1.06e-02, 1.15e-02, 3.06e-02, 2.80e-02, 2.93e-02,
            1.31e-02, 1.06e-02, 1.15e-02, 3.17e-02, 2.90e-02, 2.89e-02,
            1.42e-02, 1.06e-02, 1.03e-02, 3.50e-02, 2.96e-02, 3.33e-02};

        electronId_classBasedExperimentalAssign(cutdeta, cutdeta_tmp);
        electronId_classBasedExperimentalAssign(cutdphi, cutdphi_tmp);
        electronId_classBasedExperimentalAssign(cuteopin, cuteopin_tmp);
        electronId_classBasedExperimentalAssign(cuthoe, cuthoe_tmp);
        electronId_classBasedExperimentalAssign(cutip, cutip_tmp);
        electronId_classBasedExperimentalAssign(cutisoecal, cutisoecal_tmp);
        electronId_classBasedExperimentalAssign(cutisohcal, cutisohcal_tmp);
        electronId_classBasedExperimentalAssign(cutisotk, cutisotk_tmp);
        electronId_classBasedExperimentalAssign(cutmishits, cutmishits_tmp);
        electronId_classBasedExperimentalAssign(cutsee, cutsee_tmp);

    }

    double scTheta = (2*atan(exp(cms2.els_etaSC()[index])));
    double scEt = cms2.els_eSC()[index]*sin(scTheta); 
    //double eta = cms2.els_p4()[index].Eta(); // not used
    // this is eSuperClusterOverP
    double eOverP = cms2.els_eOverPIn()[index];

    //    double eSeed = cms2.els_eSeed()[index];
    //    double pin  = cms2.els_p4In()[index].R();   
    //    double pout = cms2.els_p4Out()[index].R(); 
    //    double eSeedOverPin = eSeed/pin; 
    //    double fBrem = (pin-pout)/pin;
    double eSeedOverPin = cms2.els_eSeedOverPIn()[index];
    double fBrem = cms2.els_fbrem()[index];


    double hOverE = cms2.els_hOverE()[index];
    double sigmaee = cms2.els_sigmaIEtaIEta()[index];
    double e25Max = cms2.els_e2x5Max()[index];
    double e15 = cms2.els_e1x5()[index];
    double e55 = cms2.els_e5x5()[index];
    double e25Maxoe55 = e25Max/e55 ;
    double e15oe55 = e15/e55 ;
    double deltaPhiIn = cms2.els_dPhiIn()[index];
    double deltaEtaIn = cms2.els_dEtaIn()[index];

    double ip = 0;
    //int mishits = electron->gsfTrack()->trackerExpectedHitsInner().numberOfHits();
    int mishits = cms2.els_exp_innerlayers()[index];
    double tkIso = cms2.els_tkIso()[index];
    double ecalIso = cms2.els_ecalIso04()[index];
    double hcalIso = cms2.els_hcalIso04()[index];
    double hcalIso1 = cms2.els_hcalDepth1TowerSumEt04()[index];
    double hcalIso2 = cms2.els_hcalDepth2TowerSumEt04()[index];

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

    if ((tkIso > cutisotk[cat+3*eb+bin*6]) ||
            (ecalIso > cutisoecal[cat+3*eb+bin*6]) ||
            (hcalIso > cutisohcal[cat+3*eb+bin*6]))
        result = 0;
    else
        result = 2;

    if (fBrem < -2)
        return result;

    //std::cout << "hoe" << hOverE << std::endl;
    if (hOverE > cuthoe[cat+3*eb+bin*6])
        return result;

    //std::cout << "see" << sigmaee << std::endl;
    if (sigmaee > cutsee[cat+3*eb+bin*6])
        return result;

    //std::cout << "dphiin" << fabs(deltaPhiIn) << std::endl;
    if (fabs(deltaPhiIn) > cutdphi[cat+3*eb+bin*6])
        return result;  

    //std::cout << "detain" << fabs(deltaEtaIn) << std::endl;
    if (fabs(deltaEtaIn) > cutdeta[cat+3*eb+bin*6])
        return result;

    //std::cout << "eseedopin " << eSeedOverPin << std::endl;
    if (eSeedOverPin < cuteopin[cat+3*eb+bin*6])
        return result;

    //std::cout << "ip" << ip << std::endl;
    if (ip > cutip[cat+3*eb+bin*6])
        return result;

    if (mishits > cutmishits[cat+3*eb+bin*6])
        return result;

    result = result + 1;
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
