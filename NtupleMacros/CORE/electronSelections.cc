
//
// electron selections
//

// CMS2 includes
#include "electronSelections.h"
#include "CMS2.h"

//
// enums and typedefs
//

enum EgammaFiduciality {
	ISEB,
	ISEBEEGAP,
	ISEE,
	ISEEGAP,
	ISEBETAGAP,
	ISEBPHIGAP,
	ISEEDEEGAP,
	ISEERINGGAP,
	ISGAP
};

// seeding type used and corrections applied

enum EgammaElectronType {
	ISECALENERGYCORRECTED,	// if false, the electron "ecalEnergy" is just the supercluster energy 
	ISMOMENTUMCORRECTED,  	// has E-p combination been applied
	ISECALDRIVEN,
	ISTRACKERDRIVEN
};

//
// candidate electron id function
//
bool electronId_cand01(const unsigned int index)
{

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
        if (abs(cms2.els_etaSC()[index]) < 1.5) {
	//if (cms2.els_fiduciality()[index] & (1<<ISEB)) {
		if (fabs(cms2.els_dEtaIn()[index] > dEtaInThresholds[0])) 	return false;
		if (fabs(cms2.els_dPhiIn()[index] > dPhiInThresholds[0])) 	return false;
		if (cms2.els_hOverE()[index] > hoeThresholds[0]) 		return false;
		if ((cms2.els_e2x5Max()[index]/cms2.els_e5x5()[index]) < e2x5Over5x5Thresholds[0]) return false;
		return true;
	}
        if (abs(cms2.els_etaSC()[index]) > 1.5) {
	//if (cms2.els_fiduciality()[index] & (1<<ISEE)) {
		if (fabs(cms2.els_dEtaIn()[index] > dEtaInThresholds[0])) 	return false;
		if (fabs(cms2.els_dPhiIn()[index] > dPhiInThresholds[0])) 	return false;
		if (cms2.els_hOverE()[index] > hoeThresholds[0]) 		return false;
		if (cms2.els_sigmaIEtaIEta()[index] > sigmaIEtaIEtaThresholds[1]) 	return false;
		return true;
	}

	return false;
}

bool electronImpact_cand01(const unsigned int index)
{
        //
        // define thresholds for EB, EE
        //
        float d0Thresholds[2]               = {0.25, 0.35};

        //
        // apply cut
        //
        if (abs(cms2.els_etaSC()[index]) < 1.5) {
        //if (cms2.els_fiduciality()[index] & (1<<ISEB)) {
		if (cms2.els_d0corr()[index] > d0Thresholds[0]) return false;
                return true;
        }
        if (abs(cms2.els_etaSC()[index]) > 1.5) {
        //if (cms2.els_fiduciality()[index] & (1<<ISEE)) {
                if (cms2.els_d0corr()[index] > d0Thresholds[1]) return false;
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
        float tkThresholds[2]         =       {4.5, 6.0};
	//float tkThresholds[2] 	= 	{2.5, 2.0};
	float ecalThresholds[2] = 	{2.5, 2.0};
	float hcalThresholds[2] = 	{1.0, 1.0};

	//
	// apply cuts
	//


	if (abs(cms2.els_etaSC()[index]) < 1.5) {
	//if (cms2.els_fiduciality()[index] & (1<<ISEB)) {
                if (cms2.els_tkIso()[index] > tkThresholds[0])    return false;
		//if (cms2.els_tkJuraIso()[index] > tkThresholds[0]) 	return false;
		if (cms2.els_ecalIso()[index] 	> ecalThresholds[0]) 	return false;
		if (cms2.els_ecalIso()[index] 	> hcalThresholds[0]) 	return false;
		return true;
	}
        if (abs(cms2.els_etaSC()[index]) > 1.5) {
	//if (cms2.els_fiduciality()[index] & (1<<ISEE)) {
                if (cms2.els_tkIso()[index] > tkThresholds[0])      return false;
		//if (cms2.els_tkJuraIso()[index] > tkThresholds[0])      return false;
		if (cms2.els_ecalIso()[index]   > ecalThresholds[0])    return false;
		if (cms2.els_ecalIso()[index]   > hcalThresholds[0])    return false;
		return true;
	}

	return false;
}

//
// class based electron id that we have used before
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

//
//conversion rejection
//
bool isFromConversionHitPattern(const unsigned int index)
{
	//if(cms2.els_exp_innerlayers().at(index) > 1) return true;
	return false;
}

bool isFromConversionPartnerTrack(const unsigned int index) {

  if(fabs(cms2.els_conv_dist().at(index)) < 0.02 &&
     fabs(cms2.els_conv_dcot().at(index)) < 0.02)
    return true;

  return false;

}
