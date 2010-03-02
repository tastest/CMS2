

#ifndef ELECTRONSELECTIONS_H
#define ELECTRONSELECTIONS_H


//
// electron selections 
//

typedef ULong64_t elecuts_t;

//
// This is for debugging and should not be used
// by anyone who does not know what it is
//
elecuts_t electronSelections_debug_;
elecuts_t electronId_debug_;

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
    ISECALENERGYCORRECTED,  // if false, the electron "ecalEnergy" is just the supercluster energy 
    ISMOMENTUMCORRECTED,    // has E-p combination been applied
    ISECALDRIVEN,
    ISTRACKERDRIVEN
};

enum ElectronSelection {
    ELEPASS_PT10,
    ELEPASS_PT20,
    ELEPASS_PT10NOT20,
    ELEPASS_DPHI,
    ELEPASS_DETA,
    ELEPASS_HOE,
    ELEPASS_LSHAPE,
    ELEPASS_ISO,
    ELEPASS_EXTRA,
    ELEPASS_D0,
	ELEPASS_ID,
    ELEPASS_NOTCONV,
    ELEPASS_NOMUON,
    ELEPASS_TYPE,
    ELEPASS_FIDUCIAL,
    ELEPASS_FULLSELECTION
};


//
// These are the components of the 
// electronSelection_candXY functions.
//
static const elecuts_t electronSelections_passall_ =   (1<<ELEPASS_TYPE)       |
                                (1<<ELEPASS_FIDUCIAL)   |
                                (1<<ELEPASS_NOMUON)     |
                                (1<<ELEPASS_ID)         |
                                (1<<ELEPASS_D0)         |
                                (1<<ELEPASS_ISO)        |
                                (1<<ELEPASS_NOTCONV);

//
// pass the candXY cut based ID part
//

static const elecuts_t electronSelections_passid_ =  (1<<ELEPASS_DETA) |
                                            (1<<ELEPASS_DPHI) |
                                            (1<<ELEPASS_HOE) |
                                            (1<<ELEPASS_LSHAPE);

//
// combined selection functions
//
bool electronSelection_cand01(const unsigned int index);
bool electronSelection_cand02(const unsigned int index);

//
// candidate electron id function
// see http://www.t2.ucsd.edu/tastwiki/bin/view/CMS/ElectronID
//
bool electronId_cand01(const unsigned int index);
bool electronId_cand02(const unsigned int index);

bool electronId_extra(const unsigned int index);
bool electronImpact_cand01(const unsigned int index);

//
// remove electrons that are overlapping with a muon
//
bool electronId_noMuon(const unsigned int index);

//
// candidate electron isolation function
// see http://www.t2.ucsd.edu/tastwiki/bin/view/CMS/EfficiencyMeasurement
//
bool electronIsolation_cand01(const unsigned int index);

//
// class based electron id that we have used before
//
bool electronId_classBasedLoose(const unsigned int index);
bool electronId_classBasedTight(const unsigned int index);

//
// electron isolation definitions that we have used before
//
float electronIsolation_rel(const unsigned int index, bool use_calo_iso);
float electronIsolation_relsusy(const unsigned int index, bool use_calo_iso);
float electronIsolation_relsusy_cand0(const unsigned int index, bool use_calo_iso);
float electronIsolation_relsusy_cand1(const unsigned int index, bool use_calo_iso);

//
// conversion rejection
//
bool isFromConversionHitPattern(const unsigned int index);
bool isFromConversionPartnerTrack(const unsigned int index);

//electron charge using the majority logic of the egamma group
int getChargeUsingMajorityLogic(int elIdx, float minFracSharedHits = 0.45);

//charge flip rejection
bool isChargeFlip(int elIndex);

#endif

