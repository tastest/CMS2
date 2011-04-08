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
	 // rel iso (fixed below 20 GeV) < 1.0
	 // 0.3 cone size for all, 1 GeV pedestal sub in EB
	 ELEISO_REL100,
	 // rel iso (fixed below 20 GeV) < 0.10
	 // 0.3 cone size for all, 1 GeV pedestal sub in EB/EE
	 ELEISO_REL010_WW,
	 // rel iso (fixed below 20 GeV) < 0.40
	 // 0.3 cone size for all, 1 GeV pedestal sub in EB/EE
	 ELEISO_REL040_WW,
	 // rel iso (fixed below 20 GeV) < 1.00
	 // 0.3 cone size for all, 1 GeV pedestal sub in EB/EE
	 ELEISO_REL100_WW,
	 // non-truncated relative iso with cut [0.05,0.07,0.10] for pT [10,15,20]
	 ELEISO_SMURFV1,

     ELEISO_RELNT010,       // non-truncated relative iso < 0.10, 0.3 cone size for all, 1 GeV pedestal sub in EB
     ELEISO_RELNT015,       // non-truncated relative iso < 0.15, 0.3 cone size for all, 1 GeV pedestal sub in EB
     ELEISO_RELNT040,       // non-truncated relative iso < 0.40, 0.3 cone size for all, 1 GeV pedestal sub in EB
     ELEISO_TRK_RELNT020,   // non-truncated relative Tracker iso < 0.20, 0.3 cone size for all
     ELEISO_ECAL_RELNT020,  // non-truncated relative ECAL    iso < 0.20, 0.3 cone size for all, 1 GeV pedestal sub in EB
     ELEISO_HCAL_RELNT020,  // non-truncated relative HCAL    iso < 0.20, 0.3 cone size for all


	 //
	 // ip cuts
	 //
	 // d0 corrected for beamspot < 0.02
	 ELEIP_200,
	 // d0 corrected for beamspot < 0.04
	 ELEIP_400,
	 // d0 corrected for primary vertex < 0.02
	 ELEIP_PV_200,
	 // d0 (PV) < 0.02 and dz (PV) < 1.0
	 ELEIP_PV_wwV1,


	 //
	 // id cuts
	 //
	 // pass smurf v1 electron ID
	 ELEID_SMURFV1_EXTRA,
	 // pass smurf v2 electron ID
	 ELEID_SMURFV2_EXTRA,
	 // pass smurf v2 electron ID
	 ELEID_SMURFV3_EXTRA,
	 // pass "CAND01" electron ID
	 ELEID_CAND01,
	 // pass "CAND02" electron ID
	 ELEID_CAND02,
	 // pass "EXTRA" electron ID
	 ELEID_EXTRA,
	 // VBTF95 electron ID (35X)
	 ELEID_VBTF_35X_95,
	 // VBTF90 electron ID (35X)
	 ELEID_VBTF_35X_90,
	 // VBTF80 electron ID (35X)
	 ELEID_VBTF_35X_80,
	 // VBTF70 electron ID (35X)
	 ELEID_VBTF_35X_70,
	 // VBTF80 electron ID no HoE in endcap
	 ELEID_VBTF_80_NOHOEEND,
	 // VBTF70 electron ID no HoE in endcap
	 ELEID_VBTF_70_NOHOEEND,
	 // VBTF90 electron ID with HoE and dPhiIn cuts tuned to 
	 // represent HLT requirements for CaloIdL_TrkIdVL
	 ELEID_VBTF_90_HLT,
	 // CIC_MEDIUM electron ID (V03)
	 ELEID_CIC_V03_MEDIUM,
	 //
	 // conv rej cuts
	 //
	 // mit conversion rejection v11 
	 //(maxhitsbeforevtx, minprob, minlxy, allowckfmatch, requirearbitratedmerged) = (0,   1e-6,   2.0,   true,  false)
	 ELENOTCONV_MIT,
	 // dist < 0.02 && dcot(theta) < 0.02
	 ELENOTCONV_DISTDCOT002,
	 // < 2 missing hits
	 ELENOTCONV_HITPATTERN,
	 // < 1 missing hits
	 ELENOTCONV_HITPATTERN_0MHITS,
	 //ELENOTCONV_HITPATTERN39X_0MHITS,
	 //
	 // eta cuts
	 //
	 // |eta| < 2.50 where eta is the ecal eta
	 ELEETA_250,
	 ELEETA_240,
	 //
	 // Pt
	 //
	 ELEPT_010,
	 ELEPT_015,
	 ELEPT_020,
	 //
	 // Super Cluster Et
	 //
	 ELESCET_010,
	 ELESCET_015,
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
	 ELECHARGE_NOTFLIP3AGREE,
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
static const cuts_t electronSelection_smurfV3_baseline  = 
	 electronSelection_wwV1_base |
	 electronSelection_wwV1_ip;
static const cuts_t electronSelection_smurfV3_convrej  = 
	 (1ll<<ELENOTCONV_HITPATTERN_0MHITS) |
         (1ll<<ELENOTCONV_MIT);
static const cuts_t electronSelection_smurfV3_iso  = 
         (1ll<<ELEISO_RELNT010);
static const cuts_t electronSelection_smurfV3_id  = 
	 (1ll<<ELEID_VBTF_80_NOHOEEND) |
         (1ll<<ELEID_SMURFV3_EXTRA);
static const cuts_t electronSelection_smurfV3  = 
         electronSelection_smurfV3_baseline |
         electronSelection_smurfV3_convrej |
         electronSelection_smurfV3_iso |
         electronSelection_smurfV3_id;
//--------end of SMURF V3 cuts--------------------------------

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
           (1ll<<ELEETA_250) |
           (1ll<<ELENOMUON_010);

// Analysis Selection (fake rate numerator)
static const cuts_t electronSelection_ssV3 = 
           electronSelection_smurfV3_baseline |
           electronSelection_smurfV3_convrej  |
           electronSelection_smurfV3_id       |
           (1ll<<ELEISO_RELNT015)             |
           (1ll<<ELECHARGE_NOTFLIP3AGREE)     |
           (1ll<<ELENOMUON_010);

// Loose "Fakeable Object" Selection (fake rate denominators)

static const cuts_t electronSelectionFOV3_ssVBTF80_v1 =       // V1 - relaxed Id & Isolation (same story here for the relIso cut)
                 electronSelectionFO_SS_baseline    |
                 electronSelection_smurfV3_convrej  |                
                 (1ll<<ELECHARGE_NOTFLIP3AGREE);

static const cuts_t electronSelectionFOV3_ssVBTF80_v2 =       // V2 - relaxed Id
                 electronSelectionFO_SS_baseline    |
                 electronSelection_smurfV3_convrej  |                
                 (1ll<<ELECHARGE_NOTFLIP3AGREE)     |
                 (1ll<<ELEISO_RELNT015);

static const cuts_t electronSelectionFOV3_ssVBTF80_v3 =       // V3 - relaxed isolation (relaxed all the way; we store the relIso and can cut on it separately in the babies or elsewhere)
                 electronSelectionFO_SS_baseline    |
                 electronSelection_smurfV3_convrej  |
                 electronSelection_smurfV3_id       |              
                 (1ll<<ELECHARGE_NOTFLIP3AGREE);

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
bool pass_electronSelection(const unsigned int index, const cuts_t selectionType, bool applyAlignmentCorrection = false, bool removedEtaCutInEndcap = false);
cuts_t electronSelection(const unsigned int index, bool applyAlignmentCorrection = false, bool removedEtaCutInEndcap = false);

//
// "smurf" electron id
// WARNING!!! this is not the full smurf selection, just the additional ID on top of VBTF80 for low pt guys
//
bool electronId_smurf_v1(const unsigned int index);
bool electronId_smurf_v2(const unsigned int index);
bool electronId_smurf_v3(const unsigned int index);

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
float electronIsolation_ECAL_rel_v1(const unsigned int);
float electronIsolation_HCAL_rel_v1(const unsigned int);

// the difference from above is that the pedestal sub is applied on both EB/EE
float electronIsolation_rel_ww(const unsigned int index, bool use_calo_iso);

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

#endif

