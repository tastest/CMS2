
#ifndef ELECTRONSELECTIONS_H
#define ELECTRONSELECTIONS_H

#include <stdint.h>
#include <vector>
#include "TMath.h"

#include "electronSelectionsParameters.h"

typedef ULong64_t   uint64;
typedef uint64      cuts_t;
typedef uint64      electronIdComponent_t;
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
    // VBTF90 electron ID (35X)
    ELEID_VBTF_35X_90,
    // CIC_MEDIUM electron ID (V03)
    ELEID_CIC_V03_MEDIUM,
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
    // Pt
    //
    ELEPT_010,
    //
    // Super Cluster Et
    //
    ELESCET_010,
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
    //
    // spike rejection
    //
    // swiss variable should be > 0.05 to not be a spike
    ELENOSPIKE_SWISS005,

};

// ===================== FAKERATES ========================
//
// The Fakeable Objects
//
// These cuts are applied in all
// fakeable object denominators
//---------------------------------------------------------
static const cuts_t electronSelectionFO_baseline =
        (1ll<<ELEETA_250) |               // |eta| < 2.5
        (1ll<<ELENOMUON_010) |            // no mu (dR < 0.1)
        (1ll<<ELESEED_ECAL) |             // not only tracker seeded
        (1ll<<ELENOTCONV_DISTDCOT002);    // dist/dcot conv. rejection
//---------------------------------------------------------
//
// TTBar fakeable object definition v1
// extrapolating in isolation and id
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_ttbar_v1 =
        electronSelectionFO_baseline |
        (1ll<<ELEISO_REL040);
//---------------------------------------------------------
//
// TTBar fakeable object definition v2
// extrapolating in id
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_ttbar_v2 =
        electronSelectionFO_baseline |
        (1ll<<ELEISO_REL010);
//---------------------------------------------------------
//
// TTBar fakeable object definition v3
// extrapolating in iso
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_ttbar_v3 =
        electronSelectionFO_baseline |
        (1ll<<ELEID_CAND01) |
        (1ll<<ELEISO_REL040);
//---------------------------------------------------------
//  
// V1 fakeable object for cand01
// extrapolating in iso and id
// NOTE - this is the same as el_ttbar_v1
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_v1_cand01 =
        electronSelectionFO_baseline |
        (1ll<<ELEISO_REL040);
//---------------------------------------------------------
//  
// V1 fakeable object for cand02 + extra
// extrapolating in iso and id
// NOTE - this is the same as el_ttbar_v1
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_v1_cand02 =
        electronSelectionFO_baseline |
        (1ll<<ELEISO_REL040);
//---------------------------------------------------------
//  
// V1 fakeable object for cand02 + extra + flipveto
// extrapolating in iso and id with charge flip veto
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_v1_cand02flip =
        electronSelectionFO_baseline |
        (1ll<<ELEISO_REL040) |
        (1ll<<ELECHARGE_NOTFLIP);
//---------------------------------------------------------
//  
// V2 fakeable object for cand01
// extrapolating in id
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_v2_cand01 =
        electronSelectionFO_baseline |
        (1ll<<ELEISO_REL010);
//---------------------------------------------------------
//  
// V2 fakeable object for cand02 + extra
// extrapolating in id
// NOTE: this is the same as above
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_v2_cand02 =
        electronSelectionFO_baseline |
        (1ll<<ELEISO_REL010);
//---------------------------------------------------------
//  
// V2 fakeable object for cand02 + extra + flipveto
// extrapolating in id with charge flip veto
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_v2_cand02flip =
        electronSelectionFO_baseline |
        (1ll<<ELEISO_REL010) |
        (1ll<<ELECHARGE_NOTFLIP);
//---------------------------------------------------------
//  
// V3 fakeable object for cand01
// extrapolating in isolation
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_v3_cand01 =
        electronSelectionFO_baseline |
        (1ll<<ELEISO_REL040) |
        (1ll<<ELEID_CAND01);
//---------------------------------------------------------
//  
// V3 fakeable object for cand02 + extra
// extrapolating in isolation
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_v3_cand02 =
        electronSelectionFO_baseline |
        (1ll<<ELEISO_REL040) |
        (1ll<<ELEID_CAND02) |
        (1ll<<ELEID_EXTRA);
//---------------------------------------------------------
//  
// V3 fakeable object for cand02 + extra + flipveto
// extrapolating in isolation with charge flip veto
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_v3_cand02flip =
        electronSelectionFO_baseline |
        (1ll<<ELEISO_REL040) |
        (1ll<<ELEID_CAND02) |
        (1ll<<ELEID_EXTRA) |
        (1ll<<ELECHARGE_NOTFLIP);
//---------------------------------------------------------




//========================================================
// ANALYSIS SELECTIONS
//========================================================



// ===================== TTBar ============================
//
// The standard TTBar selection without isolation
//---------------------------------------------------------
static const cuts_t electronSelection_ttbar_noiso =
                    (1ll<<ELEIP_400) |
                    (1ll<<ELEID_CAND01) |
                    (1ll<<ELENOTCONV_DISTDCOT002) |
                    (1ll<<ELEETA_250) |
                    (1ll<<ELENOMUON_010) |
                    (1ll<<ELESEED_ECAL);
//---------------------------------------------------------
//
// The standard TTBar isolation cut
//---------------------------------------------------------
static const cuts_t electronSelection_ttbar_iso =
                    (1ll<<ELEISO_REL010);
//---------------------------------------------------------
//
// The standard TTBar selection with isolation
//---------------------------------------------------------
static const cuts_t electronSelection_ttbar =
        electronSelection_ttbar_noiso |
        electronSelection_ttbar_iso;
//---------------------------------------------------------

//---------------------------------------------------------
// TTBar selection with VBTF90 ID
//---------------------------------------------------------
static const cuts_t electronSelection_ttbarV1 =
  (1ll<<ELEID_VBTF_35X_90) |
  (1ll<<ELEIP_400) |										
  (1ll<<ELEISO_REL015) |
  (1ll<<ELENOMUON_010) |
  (1ll<<ELENOTCONV_HITPATTERN) |
  (1ll<<ELENOTCONV_DISTDCOT002) |
  (1ll<<ELESCET_010) |
  (1ll<<ELEPT_010) |
  (1ll<<ELEETA_250) |
  (1ll<<ELESEED_ECAL);
//---------------------------------------------------------
// TTBarV1 fakeable object definition v1
// extrapolating in isolation and id
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_ttbarV1_v1 =
  electronSelectionFO_baseline |
  (1ll<<ELENOTCONV_HITPATTERN) |
  (1ll<<ELESCET_010) |
  (1ll<<ELEPT_010) |

  (1ll<<ELEISO_REL040);
//---------------------------------------------------------
// TTBarV1 fakeable object definition v2
// extrapolating in id
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_ttbarV1_v2 =
  electronSelectionFO_baseline |
  (1ll<<ELENOTCONV_HITPATTERN) |
  (1ll<<ELESCET_010) |
  (1ll<<ELEPT_010) |

  (1ll<<ELEISO_REL015);
//---------------------------------------------------------
// TTBarV1 fakeable object definition v3
// extrapolating in iso
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_ttbarV1_v3 =
  electronSelectionFO_baseline |
  (1ll<<ELENOTCONV_HITPATTERN) |
  (1ll<<ELESCET_010) |
  (1ll<<ELEPT_010) |

  (1ll<<ELEID_VBTF_35X_90) |
  (1ll<<ELEISO_REL040);
//---------------------------------------------------------


//
// ======================== WW ============================
//
// The standard WW selection with no isolation
//
//---------------------------------------------------------
static const cuts_t electronSelection_ww_noiso =
                    (1ll<<ELEID_CAND02) |
                    (1ll<<ELEID_EXTRA) |
                    (1ll<<ELENOTCONV_DISTDCOT002) |
                    (1ll<<ELEETA_250) |
                    (1ll<<ELENOMUON_010) |
                    (1ll<<ELESEED_ECAL);
//---------------------------------------------------------
//
// The standard WW impact parameter cut
//
//---------------------------------------------------------
static const cuts_t electronSelection_ww_ip =
                    (1ll<<ELEIP_200);
//---------------------------------------------------------
//
// The standard WW isolation cut
//
//---------------------------------------------------------
static const cuts_t electronSelection_ww_iso =
                    (1ll<<ELEISO_REL010);
//---------------------------------------------------------
// FIXME: for fakerates
static const cuts_t electronSelection_cand02 = 
                    electronSelection_ww_noiso |
                    electronSelection_ww_ip |
                    electronSelection_ww_iso;
//---------------------------------------------------------



// ======================== OS ============================
//
// The standard OS selection
//---------------------------------------------------------
static const cuts_t electronSelection_os =
                    (1ll<<ELEISO_REL010) |
                    (1ll<<ELEIP_200) |
                    (1ll<<ELEID_CAND01) |
                    (1ll<<ELENOTCONV_DISTDCOT002) |
                    (1ll<<ELEETA_250) |
                    (1ll<<ELENOMUON_010) |
                    (1ll<<ELESEED_ECAL);
//---------------------------------------------------------
// FIXME: for fake rates
static const cuts_t electronSelection_cand01 = electronSelection_os;
//---------------------------------------------------------



// ======================== SS ============================
//
// The standard SS selection
//---------------------------------------------------------
static const cuts_t electronSelection_ss =
                    (1ll<<ELEISO_REL010) |
                    (1ll<<ELEIP_200) |
                    (1ll<<ELEID_CAND02) |
                    (1ll<<ELEID_EXTRA) |
                    (1ll<<ELENOTCONV_DISTDCOT002) |
                    (1ll<<ELEETA_250) |
                    (1ll<<ELENOMUON_010) |
                    (1ll<<ELESEED_ECAL) |
                    (1ll<<ELECHARGE_NOTFLIP);
//---------------------------------------------------------
// FIXME: for fake rates
static const cuts_t electronSelection_cand02flip = electronSelection_ss;
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
bool pass_electronSelectionCompareMask(const cuts_t cuts_passed, const cuts_t selectionType);
bool pass_electronSelection(const unsigned int index, const cuts_t selectionType);
cuts_t electronSelection(const unsigned int index);

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
electronIdComponent_t electronId_VBTF(const unsigned int index, const vbtf_tightness tightness);

//
// "CIC" id
//
electronIdComponent_t electronId_CIC(const unsigned int index, const unsigned int version, const cic_tightness tightness);
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

//
// spike rejection for electrons
//
bool isSpikeElectron(const unsigned int index);

#endif

