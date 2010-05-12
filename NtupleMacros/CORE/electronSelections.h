
#ifndef ELECTRONSELECTIONS_H
#define ELECTRONSELECTIONS_H

#include <stdint.h>
#include <vector>

#include "electronSelectionsParameters.h"

//
// This is the menu of electron selections
//

enum EleSelectionType {
    //
    // iso cuts
    //
    // rel iso (fixed below 20 GeV) < 0.10
    // 0.3 cone size for all, 1 GeV pedestal sub in EB
    ELEISO_REL010,
    // rel iso (fixed below 20 GeV) < 0.15
    // 0.3 cone size for all, 1 GeV pedestal sub in EB
    ELEISO_REL015,
    // rel iso (fixed below 20 GeV) < 0.40
    // 0.3 cone size for all, 1 GeV pedestal sub in EB
    ELEISO_REL040,
    //
    // ip cuts
    //
    // d0 corrected for beamspot < 0.02
    ELEIP_200,
    // d0 corrected for beamspot < 0.04
    ELEIP_400,
    //
    // id cuts
    //
    // pass "CAND01" electron ID
    ELEID_CAND01,
    // pass "CAND02" electron ID
    ELEID_CAND02,
    // pass "EXTRA" electron ID
    ELEID_EXTRA,
    //
    // conv rej cuts
    //
    // dist < 0.02 && dcot(theta) < 0.02
    ELENOTCONV_DISTDCOT002,
    // < 2 missing hits
    ELENOTCONV_HITPATTERN,
    //
    // eta cuts
    //
    // |eta| < 2.50 where eta is the ecal eta
    ELEETA_250,
    //
    // no muon cuts
    //
    // no muon within dR < 0.1
    ELENOMUON_010,
    //
    // ecal driven requirement
    //
    // seed must have been found by at least the ecal algo
    ELESEED_ECAL,
    //
    // charge flipping
    //
    // is not a charge flip
    ELECHARGE_NOTFLIP,

};

// ===================== FAKERATES ========================
//
// The Fakeable Objects
//
// These cuts are applied in all
// fakeable object denominators
//---------------------------------------------------------
static const unsigned int electronSelectionFO_baseline =
        (1<<ELEETA_250) |               // |eta| < 2.5
        (1<<ELENOMUON_010) |            // no mu (dR < 0.1)
        (1<<ELESEED_ECAL) |             // not only tracker seeded
        (1<<ELENOTCONV_DISTDCOT002);    // dist/dcot conv. rejection
//---------------------------------------------------------
//
// TTBar fakeable object definition v1
// extrapolating in isolation and id
//---------------------------------------------------------
static const unsigned int electronSelectionFO_el_ttbar_v1 =
        electronSelectionFO_baseline |
        (1<<ELEISO_REL040);
//---------------------------------------------------------
//
// TTBar fakeable object definition v2
// extrapolating in id
//---------------------------------------------------------
static const unsigned int electronSelectionFO_el_ttbar_v2 =
        electronSelectionFO_baseline |
        (1<<ELEISO_REL010);
//---------------------------------------------------------
//
// TTBar fakeable object definition v3
// extrapolating in iso
//---------------------------------------------------------
static const unsigned int electronSelectionFO_el_ttbar_v3 =
        electronSelectionFO_baseline |
        (1<<ELEID_CAND01);
//---------------------------------------------------------
//  
// V1 fakeable object for cand01
// extrapolating in iso and id
// NOTE - this is the same as el_ttbar_v1
//---------------------------------------------------------
static const unsigned int electronSelectionFO_el_v1_cand01 =
        electronSelectionFO_baseline |
        (1<<ELEISO_REL040);
//---------------------------------------------------------
//  
// V1 fakeable object for cand02 + extra
// extrapolating in iso and id
// NOTE - this is the same as el_ttbar_v1
//---------------------------------------------------------
static const unsigned int electronSelectionFO_el_v1_cand02 =
        electronSelectionFO_baseline |
        (1<<ELEISO_REL040);
//---------------------------------------------------------
//  
// V1 fakeable object for cand02 + extra + flipveto
// extrapolating in iso and id with charge flip veto
//---------------------------------------------------------
static const unsigned int electronSelectionFO_el_v1_cand02flip =
        electronSelectionFO_baseline |
        (1<<ELEISO_REL040) |
        (1<<ELECHARGE_NOTFLIP);
//---------------------------------------------------------
//  
// V2 fakeable object for cand01
// extrapolating in id
//---------------------------------------------------------
static const unsigned int electronSelectionFO_el_v2_cand01 =
        electronSelectionFO_baseline |
        (1<<ELEISO_REL010);
//---------------------------------------------------------
//  
// V2 fakeable object for cand02 + extra
// extrapolating in id
// NOTE: this is the same as above
//---------------------------------------------------------
static const unsigned int electronSelectionFO_el_v2_cand02 =
        electronSelectionFO_baseline |
        (1<<ELEISO_REL010);
//---------------------------------------------------------
//  
// V2 fakeable object for cand02 + extra + flipveto
// extrapolating in id with charge flip veto
//---------------------------------------------------------
static const unsigned int electronSelectionFO_el_v2_cand02flip =
        electronSelectionFO_baseline |
        (1<<ELEISO_REL010) |
        (1<<ELECHARGE_NOTFLIP);
//---------------------------------------------------------
//  
// V3 fakeable object for cand01
// extrapolating in isolation
//---------------------------------------------------------
static const unsigned int electronSelectionFO_el_v3_cand01 =
        electronSelectionFO_baseline |
        (1<<ELEISO_REL040) |
        (1<<ELEID_CAND01);
//---------------------------------------------------------
//  
// V3 fakeable object for cand02 + extra
// extrapolating in isolation
//---------------------------------------------------------
static const unsigned int electronSelectionFO_el_v3_cand02 =
        electronSelectionFO_baseline |
        (1<<ELEISO_REL040) |
        (1<<ELEID_CAND02) |
        (1<<ELEID_EXTRA);
//---------------------------------------------------------
//  
// V3 fakeable object for cand02 + extra + flipveto
// extrapolating in isolation with charge flip veto
//---------------------------------------------------------
static const unsigned int electronSelectionFO_el_v3_cand02flip =
        electronSelectionFO_baseline |
        (1<<ELEISO_REL040) |
        (1<<ELEID_CAND02) |
        (1<<ELEID_EXTRA) |
        (1<<ELECHARGE_NOTFLIP);
//---------------------------------------------------------




//========================================================
// ANALYSIS SELECTIONS
//========================================================



// ===================== TTBar ============================
//
// The standard TTBar selection without isolation
//---------------------------------------------------------
static const unsigned int electronSelection_ttbar =
                    (1<<ELEIP_400) |
                    (1<<ELEID_CAND01) |
                    (1<<ELENOTCONV_DISTDCOT002) |
                    (1<<ELEETA_250) |
                    (1<<ELENOMUON_010) |
                    (1<<ELESEED_ECAL);
//---------------------------------------------------------
//
// The standard TTBar isolation cut
//---------------------------------------------------------
static const unsigned int electronSelection_ttbar_iso =
                    (1<<ELEISO_REL010);
//---------------------------------------------------------



// ======================== WW ============================
//
// The standard WW selection
//---------------------------------------------------------
static const unsigned int electronSelection_ww =
                    (1<<ELEISO_REL010) |
                    (1<<ELEID_CAND02) |
                    (1<<ELEID_EXTRA) |
                    (1<<ELENOTCONV_DISTDCOT002) |
                    (1<<ELEETA_250) |
                    (1<<ELENOMUON_010) |
                    (1<<ELESEED_ECAL);
//---------------------------------------------------------
//
// The standard WW selection with no isolation
//---------------------------------------------------------
static const unsigned int electronSelection_ww_noiso =
                    (1<<ELEID_CAND02) |
                    (1<<ELEID_EXTRA) |
                    (1<<ELENOTCONV_DISTDCOT002) |
                    (1<<ELEETA_250) |
                    (1<<ELENOMUON_010) |
                    (1<<ELESEED_ECAL);
//---------------------------------------------------------
//
// The standard WW impact parameter cut
//---------------------------------------------------------
static const unsigned int electronSelection_ww_ip =
                    (1<<ELEIP_200);
//---------------------------------------------------------
//
// The standard WW isolation cut
//---------------------------------------------------------
static const unsigned int electronSelection_ww_iso =
                    (1<<ELEISO_REL010);
//---------------------------------------------------------



// ======================== OS ============================
//
// The standard OS selection
//---------------------------------------------------------
static const unsigned int electronSelection_os =
                    (1<<ELEISO_REL010) |
                    (1<<ELEIP_200) |
                    (1<<ELEID_CAND01) |
                    (1<<ELENOTCONV_DISTDCOT002) |
                    (1<<ELEETA_250) |
                    (1<<ELENOMUON_010) |
                    (1<<ELESEED_ECAL);
//---------------------------------------------------------



// ======================== SS ============================
//
// The standard SS selection
//---------------------------------------------------------

//---------------------------------------------------------





//
// other enums for internal use
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
unsigned int electronSelection(const unsigned int index);

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

