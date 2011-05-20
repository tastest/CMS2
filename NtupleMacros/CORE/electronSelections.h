#ifndef ELECTRONSELECTIONS_H
#define ELECTRONSELECTIONS_H

// C++
#include <stdint.h>
#include <vector>

// Root
#include "TMath.h"
#include "CMS2.h"

// Header
#include "electronSelectionsParameters.h"

//
typedef ULong64_t   uint64;
typedef uint64      cuts_t;
typedef uint64      electronIdComponent_t;

/////////////////////////////////////////////
// This is the menu of electron selections //
/////////////////////////////////////////////

enum EleSelectionType {

///////////////
// Isolation //
///////////////
 
	 ELEISO_REL010,             // rel iso (fixed below 20 GeV) < 0.10, 0.3 cone size for all, 1 GeV pedestal sub in EB
	 ELEISO_REL015,             // rel iso (fixed below 20 GeV) < 0.15, 0.3 cone size for all, 1 GeV pedestal sub in EB
	 ELEISO_REL040,             // rel iso (fixed below 20 GeV) < 0.40, 0.3 cone size for all, 1 GeV pedestal sub in EB
	 ELEISO_REL100,             // rel iso (fixed below 20 GeV) < 1.00, 0.3 cone size for all, 1 GeV pedestal sub in EB
	 ELEISO_REL010_WW,          // rel iso (fixed below 20 GeV) < 0.10, 0.3 cone size for all, 1 GeV pedestal sub in EB/EE
	 ELEISO_REL040_WW,          // rel iso (fixed below 20 GeV) < 0.40, 0.3 cone size for all, 1 GeV pedestal sub in EB/EE
	 ELEISO_REL100_WW,          // rel iso (fixed below 20 GeV) < 1.00, 0.3 cone size for all, 1 GeV pedestal sub in EB/EE
	 ELEISO_SMURFV4,            // non-truncated relative iso with cut [0.05,0.07,0.10] for pT [10,15,20]
	 ELEISO_RELNT010,           // non-truncated relative iso < 0.10, 0.3 cone size for all, 1 GeV pedestal sub in EB
	 ELEISO_RELNT015,           // non-truncated relative iso < 0.15, 0.3 cone size for all, 1 GeV pedestal sub in EB
	 ELEISO_RELNT040,           // non-truncated relative iso < 0.40, 0.3 cone size for all, 1 GeV pedestal sub in EB
	 ELEISO_TRK_RELNT020,       // non-truncated relative Tracker iso < 0.20, 0.3 cone size for all
	 ELEISO_ECAL_RELNT020,      // non-truncated relative ECAL    iso < 0.20, 0.3 cone size for all, 1 GeV pedestal sub in EB
	 ELEISO_ECAL_RELNT020_NPS,  // non-truncated relative ECAL    iso < 0.20, 0.3 cone size for all, no pedestal sub in EB
	 ELEISO_HCAL_RELNT020,      // non-truncated relative HCAL    iso < 0.20, 0.3 cone size for all
	 ELEISO_ECAL_REL020,        // truncated relative ECAL    iso < 0.20, 0.3 cone size for all, 1 GeV pedestal sub in EB
	 ELEISO_HCAL_REL020,        // truncated relative HCAL    iso < 0.20, 0.3 cone size for all
	 ELEISO_FASTJET_REL005,     // truncated reliso < 0.05, 0.3 cone size, 1 GeV pedestal subtraction, fastjet-corrected
	 ELEISO_FASTJET_REL010,     // truncated reliso < 0.05, 0.3 cone size, 1 GeV pedestal subtraction, fastjet-corrected
	 ELEISO_FASTJET_REL015,     // truncated reliso < 0.05, 0.3 cone size, 1 GeV pedestal subtraction, fastjet-corrected

//////////////////////
// Impact Parameter //
//////////////////////

  ELEIP_200,         // d0 corrected for beamspot < 0.02
  ELEIP_400,         // d0 corrected for beamspot < 0.04
  ELEIP_PV_200,      // d0 corrected for primary vertex < 0.02
  ELEIP_PV_wwV1,     // d0 (PV) < 0.02 and dz (PV) < 1.0
  ELEIP_PV_SMURFV3,  // d0 (PV) < 0.02 and dz (PV) < 0.2, using first DA PV
  ELEIP_PV_SMURFV4,  // d0 (PV) < 0.02 and dz (PV) < 0.1, using first DA PV
  ELEIP_PV_OSV2,     // d0 (PV) < 0.04 and dz (PV) < 1, using first DA PV
  ELEIP_PV_OSV2_FO,  // d0 (PV) < 0.2 and dz (PV) < 1, using first DA PV
  ELEIP_SS200,       // 2011 SS numerator d0 cut
  ELEIP_SS2000,      // 2011 SS potential denominator d0 cut


/////////////////////////////
// Electron Identification //
/////////////////////////////

  ELEID_SMURFV1_EXTRA,                // pass smurf v1 electron ID
  ELEID_SMURFV2_EXTRA,                // pass smurf v2 electron ID
  ELEID_SMURFV3_EXTRA,                // pass smurf v3 electron ID
  ELEID_SMURFV1SS_EXTRA,              // electron ID with VBTF80
  ELEID_CAND01,                       // pass "CAND01" electron ID
  ELEID_CAND02,                       // pass "CAND02" electron ID
	ELEID_EXTRA,                        // pass "EXTRA" electron ID
	ELEID_VBTF_35X_95,                  // VBTF95 electron ID (35X)
	ELEID_VBTF_35X_90,                  // VBTF90 electron ID (35X)
	ELEID_VBTF_35X_80,                  // VBTF80 electron ID (35X)
  ELEID_VBTF_35X_70,                  // VBTF70 electron ID (35X)
  ELEID_VBTF_80_NOHOEEND,             // VBTF80 electron ID no HoE in endcap
	ELEID_VBTF_85_NOHOEEND,             // VBTF85 electron ID no HoE in endcap
	ELEID_VBTF_85,                      // VBTF85 electron ID
	ELEID_VBTF_70_NOHOEEND,             // VBTF70 electron ID no HoE in endcap
  ELEID_VBTF_90_HLT,                  // VBTF90 electron ID with HoE and dPhiIn cuts tuned to represent HLT requirements for CaloIdL_TrkIdVL
  ELEID_VBTF_90_HLT_CALOIDT_TRKIDVL,  // VBTF90 electron ID with HoE and dPhiIn cuts tuned to represent HLT requirements for CaloIdT_TrkIdVL
	ELEID_CIC_V03_MEDIUM,               // CIC_MEDIUM electron ID (V03)

//////////////////////////
// Conversion Rejection //
//////////////////////////

  ELENOTCONV_MIT,                // mit conversion rejection v11 
  ELENOTCONV_DISTDCOT002,        // dist < 0.02 && dcot(theta) < 0.02
  ELENOTCONV_HITPATTERN,         // < 2 missing hits
  ELENOTCONV_HITPATTERN_0MHITS,  // < 1 missing hits

//////////////////////
// Basic Selections //
//////////////////////

  ELEETA_250,               // |eta| < 2.50
  ELEETA_240,               // |eta| < 2.40 

	ELEPT_010,                // Pt > 10
	ELEPT_015,                // Pt > 15
	ELEPT_020,                // Pt > 20

  ELESCET_010,              // SC Et > 10
  ELESCET_015,              // SC Et > 20

  ELENOMUON_010,            // no muon within dR < 0.1

  ELESEED_ECAL,             // seed must have been found by at least the ecal algo

  ELECHARGE_NOTFLIP,        // Not a charge flip
  ELECHARGE_NOTFLIP3AGREE,  // Not a charge flip and CTF, GSF, and Pixel-SuperCluster charges agree

 
  //*************************************************************************************//
  // This is the last enumerator element                                                 //
  // The electron selection bitmask will not work with more than 63 element in this enum //
  // DO NOT ADD ANY ELEMENTS AFTER THIS!                                                 //
  //*************************************************************************************//

  ELENOSPIKE_SWISS005, // Spike Rejection, swiss variable should be > 0.05 to not be a spike

};

// Assuming the constants in EleSelectionType have default numeric values ( 0, 1, 2, ... N ),
// the last enum constant will have integer value N
// For the bitmasks to work, N must be <= 63
static bool shown = false;
inline void checkElectronSelections(void){
  int n    = (int) EleSelectionType(ELENOSPIKE_SWISS005);
  int nMax = (int) 8*sizeof(1ll) - 1;
  if( n > nMax ){
    cout << endl << "ERROR at line " << __LINE__ << " in " << __FILE__ << ":" << endl;
    cout << "enum \"EleSelectionType\" has " << n << " elements but cannot have more than " << nMax << " elements... Exiting." << endl << endl;
    exit(1);
  }
  else{
    if( !shown ){
      cout << endl << "There are " << ( nMax - n ) << " available selectors left in enum EleSelectionType" << endl;
      cout << "\t( " << __FILE__ << " )" << endl << endl;
    }
    shown = true;
  }
}

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
// TTBarV1 selection with VBTF90 ID
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
// TTBarV1 selection with VBTF90 ID for pass5 TTbar selections
// Only difference here is that SC Et goes from 10 to 15
// When pass_electronSelection is called, you need to remove the 
// alignment correction and NOT remove the dEta cut in the endcaps
// if running on 38X data
//---------------------------------------------------------
static const cuts_t electronSelection_ttbarV1_pass5 =
	 (1ll<<ELEID_VBTF_35X_90) |
	 (1ll<<ELEIP_400) |										
	 (1ll<<ELEISO_REL015) |
	 (1ll<<ELENOMUON_010) |
	 (1ll<<ELENOTCONV_HITPATTERN) |
	 (1ll<<ELENOTCONV_DISTDCOT002) |
	 (1ll<<ELESCET_015) |
	 (1ll<<ELEPT_015) |
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

//-------------------------------------------------------
//    Same as the above 3 selections
//    but with pt > 15, SCET > 15
//    This is for the ttbar dilepton pass5
//-------------------------------------------------------
static const cuts_t electronSelectionFO_el_ttbarV1_v1_pass5 =
	 electronSelectionFO_baseline |
	 (1ll<<ELENOTCONV_HITPATTERN) |
	 (1ll<<ELESCET_015) |
	 (1ll<<ELEPT_015) |
	 (1ll<<ELEISO_REL040);

//just defining a new FO object for pass6 for consistency with the muon FO ttbar naming scheme
//nothing has changed from pass5 -> pass6 however.
static const cuts_t electronSelectionFO_el_ttbarV1_v1_pass6 = electronSelectionFO_el_ttbarV1_v1_pass5;
//---------------------------------------------------------
// TTBarV1 fakeable object definition v2
// extrapolating in id
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_ttbarV1_v2_pass5 =
	 electronSelectionFO_baseline |
	 (1ll<<ELENOTCONV_HITPATTERN) |
	 (1ll<<ELESCET_015) |
	 (1ll<<ELEPT_015) |
	 (1ll<<ELEISO_REL015);

//just defining a new FO object for pass6 for consistency with the muon FO ttbar naming scheme
//nothing has changed from pass5 -> pass6 however.
static const cuts_t electronSelectionFO_el_ttbarV1_v2_pass6 =  electronSelectionFO_el_ttbarV1_v2_pass5;
//---------------------------------------------------------
// TTBarV1 fakeable object definition v3
// extrapolating in iso
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_ttbarV1_v3_pass5 =
	 electronSelectionFO_baseline |
	 (1ll<<ELENOTCONV_HITPATTERN) |
	 (1ll<<ELESCET_015) |
	 (1ll<<ELEPT_015) |
	 (1ll<<ELEID_VBTF_35X_90) |
	 (1ll<<ELEISO_REL040);

//just defining a new FO object for pass6 for consistency with the muon FO ttbar naming scheme
//nothing has changed from pass5 -> pass6 however.
static const cuts_t electronSelectionFO_el_ttbarV1_v3_pass6 = electronSelectionFO_el_ttbarV1_v3_pass5;
//---------------------------------------------------------


//---------------------------------------------------------
// TTBarV2 selection with VBTF90 ID
//---------------------------------------------------------
static const cuts_t electronSelection_ttbarV2 =
	 (1ll<<ELEID_VBTF_35X_90) |
	 (1ll<<ELEIP_400) |
	 (1ll<<ELEISO_REL015) |
	 (1ll<<ELENOMUON_010) |
	 (1ll<<ELENOTCONV_HITPATTERN) |
	 (1ll<<ELENOTCONV_DISTDCOT002) |
	 (1ll<<ELESCET_010) |
	 (1ll<<ELEPT_010) |
	 (1ll<<ELEETA_250) |
	 (1ll<<ELESEED_ECAL) |
	 (1ll<<ELENOSPIKE_SWISS005); 

//------------------------------------------------------- FOs and selections for OS -----------------------------------//
//---------------------------------------------- Same as TTbar, but different SC Et cut ----------------//

//---------------------------------------------------------
// TTbarV1 fakeable object definition v1
// extrapolating in isolation and id - > For SS
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_OSV1_v1 =
	 electronSelectionFO_baseline |
	 (1ll<<ELENOTCONV_HITPATTERN) |
	 (1ll<<ELESCET_010) |
	 (1ll<<ELEPT_010) |
	 (1ll<<ELEISO_REL040);
//---------------------------------------------------------
// TTBarV1 fakeable object definition v2
// extrapolating in id
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_OSV1_v2 =
	 electronSelectionFO_baseline |
	 (1ll<<ELENOTCONV_HITPATTERN) |
	 (1ll<<ELESCET_010) |
	 (1ll<<ELEPT_010) |
	 (1ll<<ELEISO_REL015);
//---------------------------------------------------------
// TTBarV1 fakeable object definition v3
// extrapolating in iso
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_OSV1_v3 =
	 electronSelectionFO_baseline |
	 (1ll<<ELENOTCONV_HITPATTERN) |
	 (1ll<<ELESCET_010) |
	 (1ll<<ELEPT_010) |
	 (1ll<<ELEID_VBTF_35X_90) |
	 (1ll<<ELEISO_REL040);
//---------------------------------------------------------


//Same as  electronSelection_ttbarV1_pass5, but uses pt 10, sc 10 cut instead of 15
//---------------------------------------------------------

static const cuts_t electronSelection_el_OSV1_noiso = 
   (1ll<<ELEID_VBTF_35X_90) |
   (1ll<<ELEIP_400) |    
   (1ll<<ELENOMUON_010) |
   (1ll<<ELENOTCONV_HITPATTERN) |
   (1ll<<ELENOTCONV_DISTDCOT002) |
   (1ll<<ELESCET_010) |
   (1ll<<ELEPT_010) |
   (1ll<<ELEETA_250) |
   (1ll<<ELESEED_ECAL);

static const cuts_t electronSelection_el_OSV1 = 
   (1ll<<ELEISO_REL015) | electronSelection_el_OSV1_noiso;

//---------------------------------------------------------
// baseline selection for OS 2011 analysis
//---------------------------------------------------------

static const cuts_t electronSelection_el_OSV2_noiso = 
  (1ll<<ELEID_VBTF_90_HLT_CALOIDT_TRKIDVL) | // VBTF90, tightened to match CaloIdT+TrkIdVL
  (1ll<<ELEIP_PV_OSV2)                     | // d0(PV) < 0.04 cm, dz(PV) < 1.0 cm
  (1ll<<ELENOMUON_010)                     | // no muon dR < 0.1
  (1ll<<ELENOTCONV_HITPATTERN)             | // <=1 missing hits
  (1ll<<ELENOTCONV_MIT)                    | // MIT conversion rejection
  (1ll<<ELEPT_010)                         | // electron p_T > 10 GeV
  (1ll<<ELEETA_250);                         // |eta| < 2.5

static const cuts_t electronSelection_el_OSV2_iso = 
  (1ll<<ELEISO_ECAL_RELNT020_NPS)          | // ecal/pt < 0.2 (matches HLT requirement)
  (1ll<<ELEISO_REL015);                      // reliso < 0.15, truncated, 1 GeV EB PS

static const cuts_t electronSelection_el_OSV2 = 
  electronSelection_el_OSV2_iso | electronSelection_el_OSV2_noiso;

static const cuts_t electronSelection_el_OSV2_FO = 
  (1ll<<ELEID_VBTF_90_HLT_CALOIDT_TRKIDVL) | // VBTF90, tightened to match CaloIdT+TrkIdVL
  (1ll<<ELEIP_PV_OSV2_FO)                  | // d0(PV) < 0.2 cm, dz(PV) < 1.0 cm
  (1ll<<ELENOMUON_010)                     | // no muon dR < 0.1
  (1ll<<ELENOTCONV_HITPATTERN)             | // <=1 missing hits
  (1ll<<ELENOTCONV_MIT)                    | // MIT conversion rejection
  (1ll<<ELEPT_010)                         | // electron p_T > 10 GeV
  (1ll<<ELEISO_ECAL_RELNT020_NPS)          | // ecal/pt < 0.2 (matches HLT requirement)
  (1ll<<ELEISO_REL040)                     | // reliso < 0.4, truncated, 1 GeV EB PS
  (1ll<<ELEETA_250);                         // |eta| < 2.5


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


//--------beginning of WW V0 cuts--------------------------

//---------------------------------------------------------
// WWV0 base cut
//---------------------------------------------------------
static const cuts_t electronSelection_wwV0_base  = 
	 (1ll<<ELEETA_250) |
	 (1ll<<ELENOTCONV_DISTDCOT002);

//---------------------------------------------------------
// WWV0 impact parameter cut
//---------------------------------------------------------
static const cuts_t electronSelection_wwV0_ip  = 
	 (1ll<<ELEIP_PV_200);

//---------------------------------------------------------
// WWV0 id cut
//---------------------------------------------------------
static const cuts_t electronSelection_wwV0_id  = 
	 (1ll<<ELEID_VBTF_35X_80) ;

//---------------------------------------------------------
// WWV0 isolation cut
//---------------------------------------------------------
static const cuts_t electronSelection_wwV0_iso  = 
	 (1ll<<ELEISO_REL010_WW); 

//---------------------------------------------------------
// WWV0 selection
//--------------------------------------------------------
static const cuts_t electronSelection_wwV0  = 
	 electronSelection_wwV0_base |
	 electronSelection_wwV0_ip   |
	 electronSelection_wwV0_id   |
	 electronSelection_wwV0_iso;

//---------------------------------------------------------
// WWV0 fakeable object baseline definition
//---------------------------------------------------------
static const cuts_t electronSelectionFO_wwV0_baseline =
	 electronSelection_wwV0_base;

//---------------------------------------------------------
// WWV0 fakeable object definition v1
// extrapolating in isolation and id
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_wwV0_v1 =
	 electronSelectionFO_wwV0_baseline |
	 (1ll<<ELEISO_REL040_WW); 

//---------------------------------------------------------
// WWV0 fakeable object definition v2
// extrapolating in id
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_wwV0_v2 =
	 electronSelectionFO_wwV0_baseline |
	 electronSelection_wwV0_iso;

//---------------------------------------------------------
// WWV0 fakeable object definition v3
// extrapolating in iso
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_wwV0_v3 =
	 electronSelectionFO_wwV0_baseline |
	 electronSelection_wwV0_id |
	 (1ll<<ELEISO_REL100_WW); 

//---------------------------------------------------------
// WWV0 fakeable object definition v4
// extrapolating in iso
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_wwV0_v4 =
	 electronSelectionFO_wwV0_baseline |
	 (1ll<<ELEID_VBTF_35X_90) |
	 (1ll<<ELEISO_REL100_WW); 

//--------end of WW V0 cuts--------------------------------

//--------beginning of WW V0b cuts--------------------------
// At the time of this writing the only difference between V0
// and V0b is that WW V0b additionally requires zero missing
// hits conversion rejection in the electron selection. Conversion
// rejection is in the baseline selection, so the following
// is really an exact copy of the above, except that the
// baseline includes ELENOTCONV_HITPATTERN_0MHITS, and V0 ->
// V0b everywhere. However, I will write the selection out in
// its full glory, rather than just set things equal to the
// above, in recognition that we may choose to evolve V0b into
// a real second generation selection that is more than just
// a copy of V0 plus the 0mhits conversion rejection ;)

//---------------------------------------------------------
// WWV0b base cut
//---------------------------------------------------------
static const cuts_t electronSelection_wwV0b_base  = 
	 (1ll<<ELEETA_250) |
	 (1ll<<ELENOTCONV_DISTDCOT002) |
	 (1ll<<ELENOTCONV_HITPATTERN_0MHITS);

//---------------------------------------------------------
// WWV0b impact parameter cut
//---------------------------------------------------------
static const cuts_t electronSelection_wwV0b_ip  = 
	 (1ll<<ELEIP_PV_200);

//---------------------------------------------------------
// WWV0b id cut
//---------------------------------------------------------
static const cuts_t electronSelection_wwV0b_id  = 
	 (1ll<<ELEID_VBTF_35X_80) ;

//---------------------------------------------------------
// WWV0b isolation cut
//---------------------------------------------------------
static const cuts_t electronSelection_wwV0b_iso  = 
	 (1ll<<ELEISO_REL010); 

//---------------------------------------------------------
// WWV0b selection
//--------------------------------------------------------
static const cuts_t electronSelection_wwV0b  = 
	 electronSelection_wwV0b_base |
	 electronSelection_wwV0b_ip   |
	 electronSelection_wwV0b_id   |
	 electronSelection_wwV0b_iso;

//---------------------------------------------------------
// WWV0b fakeable object baseline definition
//---------------------------------------------------------
static const cuts_t electronSelectionFO_wwV0b_baseline =
	 electronSelection_wwV0b_base;

//---------------------------------------------------------
// WWV0b fakeable object definition v1
// extrapolating in isolation and id
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_wwV0b_v1 =
	 electronSelectionFO_wwV0b_baseline |
	 (1ll<<ELEISO_REL040); 

//---------------------------------------------------------
// WWV0b fakeable object definition v2
// extrapolating in id
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_wwV0b_v2 =
	 electronSelectionFO_wwV0b_baseline |
	 electronSelection_wwV0b_iso;

//---------------------------------------------------------
// WWV0b fakeable object definition v3
// extrapolating in iso
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_wwV0b_v3 =
	 electronSelectionFO_wwV0b_baseline |
	 electronSelection_wwV0b_id |
	 (1ll<<ELEISO_REL100); 

//---------------------------------------------------------
// WWV0b fakeable object definition v4
// extrapolating in iso
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_wwV0b_v4 =
	 electronSelectionFO_wwV0b_baseline |
	 (1ll<<ELEID_VBTF_35X_90) |
	 (1ll<<ELEISO_REL100); 

//--------end of WW V0b cuts--------------------------------

//---------------------------------------------------------
// WWV1 base cut
//---------------------------------------------------------
static const cuts_t electronSelection_wwV1_base  = 
         (1ll<<ELEETA_250) ;

//---------------------------------------------------------
// WWV1 convrej
//---------------------------------------------------------
static const cuts_t electronSelection_wwV1_convrej  = 
	 (1ll<<ELENOTCONV_DISTDCOT002) | 
       //(1ll<<ELENOTCONV_HITPATTERN39X_0MHITS);
	 (1ll<<ELENOTCONV_HITPATTERN_0MHITS);

//---------------------------------------------------------
// WWV1 impact parameter cut
//---------------------------------------------------------
static const cuts_t electronSelection_wwV1_ip  = 
	 (1ll<<ELEIP_PV_wwV1);

//---------------------------------------------------------
// WWV1 id cut
//---------------------------------------------------------
static const cuts_t electronSelection_wwV1_id  = 
	 (1ll<<ELEID_VBTF_35X_80) ;

//---------------------------------------------------------
// WWV1 isolation cut
//---------------------------------------------------------
static const cuts_t electronSelection_wwV1_iso  = 
	 (1ll<<ELEISO_REL010); 

//---------------------------------------------------------
// WWV1 selection
//--------------------------------------------------------
static const cuts_t electronSelection_wwV1  = 
	 electronSelection_wwV1_base |
	 electronSelection_wwV1_convrej |
	 electronSelection_wwV1_ip   |
	 electronSelection_wwV1_id   |
	 electronSelection_wwV1_iso;

//---------------------------------------------------------
// WWV1 fakeable object baseline definition
//---------------------------------------------------------
static const cuts_t electronSelectionFO_wwV1_baseline =
	 electronSelection_wwV1_base | 
         electronSelection_wwV1_convrej;


//---------------------------------------------------------
// WWV1 fakeable object definition v1
// extrapolating in isolation and id
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_wwV1_v1 =
	 electronSelectionFO_wwV1_baseline |
	 (1ll<<ELEISO_REL040); 

//---------------------------------------------------------
// WWV1 fakeable object definition v2
// extrapolating in id
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_wwV1_v2 =
	 electronSelectionFO_wwV1_baseline |
	 electronSelection_wwV1_iso;

//---------------------------------------------------------
// WWV1 fakeable object definition v3
// extrapolating in iso
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_wwV1_v3 =
         electronSelectionFO_wwV1_baseline |
	 electronSelection_wwV1_id |
	 (1ll<<ELEISO_REL100); 

//---------------------------------------------------------
// WWV1 fakeable object definition v4
// extrapolating in iso
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_wwV1_v4 =
	 electronSelectionFO_wwV1_baseline |
	 (1ll<<ELEID_VBTF_35X_90) |
	 (1ll<<ELEISO_REL100); 

//--------end of WW V1 cuts--------------------------------

//--------SMURF V1 cuts--------------------------------
/*
static const cuts_t electronSelection_smurfV1_baseline  = 
	 electronSelection_wwV1_base |
	 electronSelection_wwV1_convrej |
	 electronSelection_wwV1_ip;
static const cuts_t electronSelection_smurfV1_iso  = 
         (1ll<<ELEISO_SMURFV1);
static const cuts_t electronSelection_smurfV1_id  = 
	 electronSelection_wwV1_id |
         (1ll<<ELEID_SMURFV1_EXTRA);
static const cuts_t electronSelection_smurfV1  = 
         electronSelection_smurfV1_baseline |
         electronSelection_smurfV1_iso |
         electronSelection_smurfV1_id;
*/
//--------end of SMURF V1 cuts--------------------------------

//--------SMURF V2 cuts--------------------------------
static const cuts_t electronSelection_smurfV2_baseline  = 
	 electronSelection_wwV1_base |
	 electronSelection_wwV1_convrej |
	 electronSelection_wwV1_ip;
static const cuts_t electronSelection_smurfV2_iso  = 
         (1ll<<ELEISO_RELNT010);
static const cuts_t electronSelection_smurfV2_id  = 
	 electronSelection_wwV1_id |
         (1ll<<ELEID_SMURFV2_EXTRA);
static const cuts_t electronSelection_smurfV2  = 
         electronSelection_smurfV2_baseline |
         electronSelection_smurfV2_iso |
         electronSelection_smurfV2_id;
//--------end of SMURF V2 cuts--------------------------------

//--------SMURF V3 cuts--------------------------------
static const cuts_t electronSelection_smurfV3_ip  = 
         (1ll<<ELEIP_PV_SMURFV3);
static const cuts_t electronSelection_smurfV3_baseline  = 
	 electronSelection_wwV1_base |
	 electronSelection_smurfV3_ip;
static const cuts_t electronSelection_smurfV3_convrej  = 
	 (1ll<<ELENOTCONV_HITPATTERN_0MHITS) |
         (1ll<<ELENOTCONV_MIT);
static const cuts_t electronSelection_smurfV3_iso  = 
         (1ll<<ELEISO_RELNT010);
static const cuts_t electronSelection_smurfV3_id  = 
	 (1ll<<ELEID_VBTF_80_NOHOEEND) |
         (1ll<<ELEID_SMURFV3_EXTRA);

static const cuts_t electronSelection_smurfV1ss_id  =
         (1ll<<ELEID_VBTF_80_NOHOEEND) |
         (1ll<<ELEID_SMURFV1SS_EXTRA);

static const cuts_t electronSelection_smurfV3  = 
         electronSelection_smurfV3_baseline |
         electronSelection_smurfV3_convrej |
         electronSelection_smurfV3_iso |
         electronSelection_smurfV3_id;
//--------end of SMURF V3 cuts--------------------------------

//--------SMURF V4 cuts--------------------------------
static const cuts_t electronSelection_smurfV4_ip  = 
         (1ll<<ELEIP_PV_SMURFV4);
static const cuts_t electronSelection_smurfV4_baseline  = 
	 electronSelection_wwV1_base |
	 electronSelection_smurfV4_ip;
static const cuts_t electronSelection_smurfV4_convrej  = 
	 (1ll<<ELENOTCONV_HITPATTERN_0MHITS) |
         (1ll<<ELENOTCONV_MIT);
static const cuts_t electronSelection_smurfV4_iso  = 
         (1ll<<ELEISO_SMURFV4);
static const cuts_t electronSelection_smurfV4_id  = 
	 (1ll<<ELEID_VBTF_80_NOHOEEND) |
         (1ll<<ELEID_SMURFV3_EXTRA);
static const cuts_t electronSelection_smurfV4  = 
         electronSelection_smurfV4_baseline |
         electronSelection_smurfV4_convrej |
         electronSelection_smurfV4_iso |
         electronSelection_smurfV4_id;
//--------end of SMURF V4 cuts--------------------------------

//--------SMURF FakableObject cuts--------------------------------
static const cuts_t electronSelectionFO_el_smurf_base =
  (1ll<<ELEETA_250) |
  electronSelection_smurfV3_convrej;
//---------------------------------------------------------
// Fakeable object definition (option V3)
// extrapolating in isolation as much as the trigger allows
// *! USE WITH CARE !*
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_smurf_v3 =
  electronSelectionFO_el_smurf_base |
  (1ll<<ELEID_VBTF_80_NOHOEEND);
//---------------------------------------------------------
// Fakeable object definition (option V1)
// extrapolating in isolation as much as the trigger allows
// and in partial id
// *! USE WITH CARE !*
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_smurf_v1 =
  electronSelectionFO_el_smurf_base |
  (1ll<<ELEID_VBTF_90_HLT);
//---------------------------------------------------------
// Fakeable object definition (option V4)
// extrapolating in partial id and partial isolation
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_smurf_v4 =
  electronSelectionFO_el_smurf_base |
  (1ll<<ELEID_VBTF_90_HLT) |
  (1ll<<ELEISO_TRK_RELNT020) |
  (1ll<<ELEISO_ECAL_RELNT020) |
  (1ll<<ELEISO_HCAL_RELNT020);
//---------------------------------------------------------
// Fakeable object definition (option V2)
// extrapolating in partial id
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_smurf_v2 =
  electronSelectionFO_el_smurf_base |
  (1ll<<ELEID_VBTF_90_HLT) |
  (1ll<<ELEISO_RELNT010); 
//--------end of SMURF FakableObject cuts------------------


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
// SS NoIso selections

static const cuts_t electronSelection_ss_NoIso =
	 (1ll<<ELEID_VBTF_35X_80) |
	 (1ll<<ELEIP_200) |
	 (1ll<<ELENOMUON_010) |
	 (1ll<<ELENOTCONV_HITPATTERN) |
	 (1ll<<ELENOTCONV_DISTDCOT002) |
	 (1ll<<ELESCET_010) |
	 (1ll<<ELEPT_010) |
	 (1ll<<ELEETA_250) |
	 (1ll<<ELECHARGE_NOTFLIP3AGREE) |
	 (1ll<<ELESEED_ECAL);

// second version of SS electron ID defined on SS analysis twiki
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/SameSignDilepton2010
static const cuts_t electronSelection_ssV2_NoIso =
	 (1ll<<ELEID_VBTF_35X_80) |
	 (1ll<<ELEIP_200) |
	 (1ll<<ELENOMUON_010) |
	 (1ll<<ELENOTCONV_HITPATTERN) |
	 (1ll<<ELENOTCONV_DISTDCOT002) |
	 (1ll<<ELESCET_010) |
	 (1ll<<ELEPT_010) |
	 (1ll<<ELEETA_240) |
	 (1ll<<ELESEED_ECAL) |
	 (1ll<<ELECHARGE_NOTFLIP3AGREE) |
	 (1ll<<ELENOTCONV_HITPATTERN_0MHITS);


// SS v2 electron ID, cut on d0 corrected for PV
// instead of d0 corrected for BS
static const cuts_t electronSelection_ssV2d0PV_NoIso =
	 (1ll<<ELEID_VBTF_35X_80) |
	 (1ll<<ELEIP_PV_200) |
	 (1ll<<ELENOMUON_010) |
	 (1ll<<ELENOTCONV_HITPATTERN) |
	 (1ll<<ELENOTCONV_DISTDCOT002) |
	 (1ll<<ELESCET_010) |
	 (1ll<<ELEPT_010) |
	 (1ll<<ELEETA_240) |
	 (1ll<<ELESEED_ECAL) |
	 (1ll<<ELECHARGE_NOTFLIP3AGREE) |
	 (1ll<<ELENOTCONV_HITPATTERN_0MHITS);

// SS v2 electron ID with no d0 cut
static const cuts_t electronSelection_ssV2nod0_NoIso =
	 (1ll<<ELEID_VBTF_35X_80) |
	 (1ll<<ELEIP_PV_200) |
	 (1ll<<ELENOMUON_010) |
	 (1ll<<ELENOTCONV_HITPATTERN) |
	 (1ll<<ELENOTCONV_DISTDCOT002) |
	 (1ll<<ELESCET_010) |
	 (1ll<<ELEPT_010) |
	 (1ll<<ELEETA_240) |
	 (1ll<<ELESEED_ECAL) |
	 (1ll<<ELECHARGE_NOTFLIP3AGREE) |
	 (1ll<<ELENOTCONV_HITPATTERN_0MHITS);

// SS v2 electron ID with no triple charge consistency requirement
static const cuts_t electronSelection_ssV2noTripleChargeReq_NoIso =
	 (1ll<<ELEID_VBTF_35X_80) |
	 (1ll<<ELEIP_200) |
	 (1ll<<ELENOMUON_010) |
	 (1ll<<ELENOTCONV_HITPATTERN) |
	 (1ll<<ELENOTCONV_DISTDCOT002) |
	 (1ll<<ELESCET_010) |
	 (1ll<<ELEPT_010) |
	 (1ll<<ELEETA_240) |
	 (1ll<<ELESEED_ECAL) |
	 (1ll<<ELENOTCONV_HITPATTERN_0MHITS);

// SS Iso selections

/////////////////////////////////////
// 2011 SS Selections              //
/////////////////////////////////////

//baseline SS FO selections
static const cuts_t electronSelectionFO_SS_baseline =
           (1ll<<ELEETA_240) |
           (1ll<<ELENOMUON_010);

// Analysis Selection (fake rate numerator)
static const cuts_t electronSelection_ssV3_noIso = 
           electronSelectionFO_SS_baseline    |
           electronSelection_smurfV3_convrej  |
           electronSelection_smurfV1ss_id     |
           (1ll<<ELEIP_SS200)                 |
           (1ll<<ELECHARGE_NOTFLIP3AGREE);

static const cuts_t electronSelection_ssV3_iso =
                 (1ll<<ELEISO_RELNT015) |
                 (1ll<<ELEISO_ECAL_RELNT020_NPS);

static const cuts_t electronSelection_ssV3 = 
                       electronSelection_ssV3_noIso |
                       electronSelection_ssV3_iso;

// Loose "Fakeable Object" Selection (fake rate denominators)

static const cuts_t electronSelectionFOV3_ssVBTF80_v1 =       // V1 - relaxed Id & Isolation
                 electronSelectionFO_SS_baseline    |
                 electronSelection_smurfV3_convrej  |                
                 (1ll<<ELECHARGE_NOTFLIP3AGREE)     |
                 (1ll<<ELEISO_ECAL_RELNT020_NPS)         |
                 (1ll<<ELEISO_HCAL_RELNT020);

static const cuts_t electronSelectionFOV3_ssVBTF80_v2 =       // V2 - relaxed Id
                 electronSelectionFO_SS_baseline    |
                 electronSelection_smurfV3_convrej  |                
                 (1ll<<ELECHARGE_NOTFLIP3AGREE)     |
                 (1ll<<ELEISO_RELNT015)             |
                 (1ll<<ELEISO_ECAL_RELNT020_NPS);

static const cuts_t electronSelectionFOV3_ssVBTF80_v3 =       // V3 - relaxed isolation (relaxed all the way; we store the relIso and can cut on it separately in the babies or elsewhere)
                 electronSelectionFO_SS_baseline    |
                 electronSelection_smurfV3_convrej  |
                 electronSelection_smurfV1ss_id     |              
                 (1ll<<ELECHARGE_NOTFLIP3AGREE)     |
                 (1ll<<ELEISO_ECAL_RELNT020_NPS)    |
                 (1ll<<ELEISO_HCAL_RELNT020);

/////////////////////////////////////
// End 2011 SS Selections          //
/////////////////////////////////////




static const cuts_t electronSelection_ss_Iso =
	 (1ll<<ELEISO_REL010);

// SS Final selections

static const cuts_t electronSelection_ss =
	 electronSelection_ss_NoIso |
	 electronSelection_ss_Iso;               

// second version of SS electron selections defined on SS analysis twiki
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/SameSignDilepton2010
static const cuts_t electronSelection_ssV2 =
	 electronSelection_ssV2_NoIso |
	 electronSelection_ss_Iso;               

// SS v2 electron ID, cut on d0 corrected for PV
// instead of d0 corrected for BS
static const cuts_t electronSelection_ssV2d0PV =
	 electronSelection_ssV2d0PV_NoIso |
	 electronSelection_ss_Iso;               

// SS v2 electron ID with no d0 cut
static const cuts_t electronSelection_ssV2nod0 =
	 electronSelection_ssV2nod0_NoIso |
	 electronSelection_ss_Iso;               

// SS v2 electron ID with no triple charge consistency requirement
static const cuts_t electronSelection_ssV2noTripleChargeReq =
	 electronSelection_ssV2noTripleChargeReq_NoIso |
	 electronSelection_ss_Iso;               

// SS Flip Veto 
const cuts_t electronSelection_ss_Flip =
       (1ll<<ELECHARGE_NOTFLIP3AGREE);


//---------------------------------------------------------
//  SS With VBTFID V1, V2, V3
//---------------------------------------------------------

static const cuts_t electronSelectionFO_ssVBTF80_v1 =
	 electronSelectionFO_baseline |
	 (1ll<<ELENOTCONV_HITPATTERN) |
	 (1ll<<ELESCET_010) |
	 (1ll<<ELEPT_010) |
	 (1ll<<ELEISO_REL040);

static const cuts_t electronSelectionFO_ssVBTF80_v2 =
	 electronSelectionFO_baseline |
	 (1ll<<ELENOTCONV_HITPATTERN) |
	 (1ll<<ELESCET_010) |
	 (1ll<<ELEPT_010) |
	 (1ll<<ELEISO_REL010);

static const cuts_t electronSelectionFO_ssVBTF80_v3 =
	 electronSelectionFO_baseline |
	 (1ll<<ELENOTCONV_HITPATTERN) |
	 (1ll<<ELESCET_010) |
	 (1ll<<ELEPT_010) |
	 (1ll<<ELEID_VBTF_35X_80) |
	 (1ll<<ELEISO_REL040);

//---------------------------------------------------------
//  SS V2 With VBTFID V1, V2, V3
//---------------------------------------------------------

static const cuts_t electronSelectionFOV2_ssVBTF80_v1 =
	 electronSelectionFO_baseline |
	 (1ll<<ELEETA_240) |
	 (1ll<<ELESCET_010) |
	 (1ll<<ELEPT_010) |
	 (1ll<<ELEISO_REL040) |
	 (1ll<<ELENOTCONV_HITPATTERN) |
	 (1ll<<ELENOTCONV_HITPATTERN_0MHITS) |
	 (1ll<<ELECHARGE_NOTFLIP3AGREE);

static const cuts_t electronSelectionFOV2_ssVBTF80_v2 =
	 electronSelectionFO_baseline |
	 (1ll<<ELEETA_240) |
	 (1ll<<ELENOTCONV_HITPATTERN) |
	 (1ll<<ELESCET_010) |
	 (1ll<<ELEPT_010) |
	 (1ll<<ELENOTCONV_HITPATTERN_0MHITS) |
	 (1ll<<ELECHARGE_NOTFLIP3AGREE) |
	 (1ll<<ELEISO_REL010);


static const cuts_t electronSelectionFOV2_ssVBTF80_v3 =
	 electronSelectionFO_baseline |
	 (1ll<<ELEETA_240) |
	 (1ll<<ELENOTCONV_HITPATTERN) |
	 (1ll<<ELESCET_010) |
	 (1ll<<ELEPT_010) |
	 (1ll<<ELENOTCONV_HITPATTERN_0MHITS) |
	 (1ll<<ELECHARGE_NOTFLIP3AGREE) |
	 (1ll<<ELEID_VBTF_35X_80) |
	 (1ll<<ELEISO_REL040);


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
bool pass_electronSelection(const unsigned int index, const cuts_t selectionType, bool applyAlignmentCorrection = false, bool removedEtaCutInEndcap = false, int vertex_index = -1);
cuts_t electronSelection(const unsigned int index, bool applyAlignmentCorrection = false, bool removedEtaCutInEndcap = false, int vertex_index = -1);

//
// "smurf" electron id
// WARNING!!! this is not the full smurf selection, just the additional ID on top of VBTF80 for low pt guys
//
bool electronId_smurf_v1(const unsigned int index);
bool electronId_smurf_v2(const unsigned int index);
bool electronId_smurf_v3(const unsigned int index);
bool electronId_smurf_v1ss(const unsigned int index);

//
// "cand" electron id
//
bool electronId_cand(const unsigned int index, const cand_tightness tightness, bool applyAlignementCorrection = false, bool removedEtaCutInEndcap = false);
bool electronId_extra(const unsigned int index);


//
// "VBTF" id
//
electronIdComponent_t electronId_VBTF(const unsigned int index, const vbtf_tightness tightness,  bool applyAlignementCorrection = false, bool removedEtaCutInEndcap = false);

//
// "CIC" id
//
electronIdComponent_t electronId_CIC(const unsigned int index, const unsigned int version, const cic_tightness tightness, bool applyAlignementCorrection = false, bool removedEtaCutInEndcap =false);

unsigned int eidClassify(const unsigned int version, const unsigned int index);
bool eidComputeCut(double x, double et, double cut_min, double cut_max, bool gtn=false);

//
// relative isolation 
// - standard track isolation from CMSSW
// --- CMSSW >= 3_5_X track jurassic strip half width 0.015
// --- CMSSW < 3_5_X no jurassic iso for tracks
// - ecal pedestal of 1 GeV subtracted in EB as Max(0, ecalIso - 1.0)
// - hcal iso as usual
float electronIsolation_rel(const unsigned int index, bool use_calo_iso);
//non-truncated relative iso
//float electronIsolation_rel_v1Original(const unsigned int index, bool use_calo_iso);
float electronIsolation_rel_v1(const unsigned int, bool);
float electronIsolation_ECAL_rel_v1(const unsigned int, bool useEBps = true);
float electronIsolation_HCAL_rel_v1(const unsigned int);
float electronIsolation_ECAL_rel(const unsigned int);
float electronIsolation_HCAL_rel(const unsigned int);

float electronIsolation_rel_FastJet(const unsigned int, bool);
float electronIsolation_rel_v1_FastJet(const unsigned int, bool);
float el_fastjet_rel_offset(const unsigned int);

// the difference from above is that the pedestal sub is applied on both EB/EE
float electronIsolation_rel_ww(const unsigned int index, bool use_calo_iso);

float electronIsoValuePF(const unsigned int iel, unsigned int idavtx, float coner=0.4, float minptn=1.0, float dzcut=0.1, 
			 float footprintdr=0.07, float gammastripveto=0.025, float elestripveto=0.025);


//
// remove electrons that are overlapping with a muon
//
bool electronId_noMuon(const unsigned int index);

//
// conversion rejection
//
bool isFromConversionHitPattern(const unsigned int index);
bool isFromConversionPartnerTrack(const unsigned int index);
bool isFromConversionMIT(const unsigned int index);

//
//electron charge using the majority logic of the egamma group
//
int getChargeUsingMajorityLogic(int elIdx, float minFracSharedHits = 0.45);

//
//charge flip rejection
//
bool isChargeFlip(int elIndex);
bool isChargeFlip3agree(int elIndex); 
//
// spike rejection for electrons
//
bool isSpikeElectron(const unsigned int index);

//
// position correction for electrons
//

void electronCorrection_pos(const unsigned int index, float &dEtaIn, float &dPhiIn);

//
// d0 corrected by the primary vertex
//
double electron_d0PV(unsigned int index);
double electron_d0PV_wwV1(unsigned int index);
double electron_d0PV_mindz(unsigned int index);
double electron_dzPV_wwV1(unsigned int index);
double electron_d0PV_smurfV3(unsigned int index);
double electron_dzPV_smurfV3(unsigned int index);

double dzPV(const LorentzVector& vtx, const LorentzVector& p4, const LorentzVector& pv);

#endif

