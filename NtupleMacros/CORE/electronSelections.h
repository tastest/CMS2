
#ifndef ELECTRONSELECTIONS_H
#define ELECTRONSELECTIONS_H

#include <stdint.h>
#include <vector>

#include "electronSelectionsParameters.h"

//
// Configurable components
// of the selection
//

enum EleSelectionType {
    //
    // iso cuts
    //
    ELEISO_REL010,
    ELEISO_REL015,
    //
    // ip cuts
    //
    ELEIP_200,
    ELEIP_400,
    //
    // id cuts
    //
    ELEID_CAND01,
    ELEID_CAND02,
    ELEID_EXTRA,
    //
    // conv rej cuts
    //
    ELENOTCONV_DISTDCOT002,
    ELENOTCONV_HITPATTERN,
    //
    // eta cuts
    //
    ELEETA_250,
    //
    // no muon cuts
    //
    ELENOMUON_010,
    //
    // ecal driven requirement
    //
    ELESEED_ECAL,
};

//
// The standard TTBar selection
//
//---------------------------------------------------------
static const unsigned int electronSelection_ttbar =
                    (1<<ELEISO_REL010) |
                    (1<<ELEIP_400) |
                    (1<<ELEID_CAND01) |
                    (1<<ELENOTCONV_DISTDCOT002) |
                    (1<<ELEETA_250) |
                    (1<<ELENOMUON_010) |
                    (1<<ELESEED_ECAL);
//---------------------------------------------------------

//
// The standard WW selection
//
//---------------------------------------------------------
static const unsigned int electronSelection_ww =
                    (1<<ELEISO_REL010) |
                    (1<<ELEIP_400) |
                    (1<<ELEID_CAND01) |
                    (1<<ELENOTCONV_DISTDCOT002) |
                    (1<<ELEETA_250) |
                    (1<<ELENOMUON_010) |
                    (1<<ELESEED_ECAL);
//---------------------------------------------------------

//
// The standard SS susy selection
//
//---------------------------------------------------------
static const unsigned int electronSelection_ss =
                    (1<<ELEISO_REL010) |
                    (1<<ELEIP_400) |
                    (1<<ELEID_CAND01) |
                    (1<<ELENOTCONV_DISTDCOT002) |
                    (1<<ELEETA_250) |
                    (1<<ELENOMUON_010) |
                    (1<<ELESEED_ECAL);
//---------------------------------------------------------

//
// The standard OS susy selection
//
//---------------------------------------------------------
static const unsigned int electronSelection_os =
                    (1<<ELEISO_REL010) |
                    (1<<ELEIP_400) |
                    (1<<ELEID_CAND01) |
                    (1<<ELENOTCONV_DISTDCOT002) |
                    (1<<ELEETA_250) |
                    (1<<ELENOMUON_010) |
                    (1<<ELESEED_ECAL);
//---------------------------------------------------------


//
//
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

enum EgammaElectronType {
    ISECALENERGYCORRECTED,  // if false, the electron "ecalEnergy" is just the supercluster energy 
    ISMOMENTUMCORRECTED,    // has E-p combination been applied
    ISECALDRIVEN,
    ISTRACKERDRIVEN
};

enum ElectronIDComponent {
    ELEID_ID,
    ELEID_ISO,
    ELEID_CONV,
    ELEID_IP,
};

// 
// master selection function
//
bool pass_electronSelection(const unsigned int index, const unsigned int selectionType);
unsigned int electronSelection(const unsigned int index, const unsigned int selectionType);

//
// "cand" electron id
//
bool electronId_cand(const unsigned int index, const cand_tightness tightness);
bool electronId_extra(const unsigned int index);

//
// class based electron id that we have used before
// e.g. Jim Branson V01
//
bool electronId_classBasedLoose(const unsigned int index);
bool electronId_classBasedTight(const unsigned int index);

//
// "VBTF" id
//
unsigned int electronId_VBTF(const unsigned int index, const vbtf_tightness tightness);

//
// "CIC" id
//
unsigned int electronId_CIC(const unsigned int index, const cic_tightness tightness);
unsigned int classify(const unsigned int version, const unsigned int index);

//
// relative isolation 
// - standard track isolation from CMSSW
// --- CMSSW >= 3_5_X track jurassic strip half width 0.015
// --- CMSSW < 3_5_X no jurassic iso for tracks
// - ecal pedestal of 1 GeV subtracted in EB as Max(0, ecalIso - 1.0)
// - hcal iso as usual
float electronIsolation_rel(const unsigned int index, bool use_calo_iso);

//
// remove electrons that are overlapping with a muon
//
bool electronId_noMuon(const unsigned int index);

//
// conversion rejection
//
bool isFromConversionHitPattern(const unsigned int index);
bool isFromConversionPartnerTrack(const unsigned int index);

//
//electron charge using the majority logic of the egamma group
//
int getChargeUsingMajorityLogic(int elIdx, float minFracSharedHits = 0.45);

//
//charge flip rejection
//
bool isChargeFlip(int elIndex);

#endif

