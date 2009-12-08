#include "Cuts.h"
//Sorry for having to make this file. The Powers That Be say it shall be so.

//this is used to turn off *ALL W AND Z CUTS* and jet cuts
cuts_t non_lepcuts = ( CUT_BIT(DILEP_PT) |
					   CUT_BIT(DILEP_PT_1020) |
					   CUT_BIT(DILEP_PT_1030) |
					   CUT_BIT(DILEP_PT_2025) |
					   CUT_BIT(DILEP_PT_2030) |
					   CUT_BIT(DILEP_OS) |
					   CUT_BIT(DILEP_MASS) |   
					   CUT_BIT(DILEP_MASS_WIDE) |
					   CUT_BIT(TCMET) |
					   CUT_BIT(TCMET_20) |
					   CUT_BIT(TCMET_25) |
					   CUT_BIT(TCMET_30) |
					   CUT_BIT(TCMET_35) |
					   CUT_BIT(TMASS_100) |
					   CUT_BIT(TMASS_120) |
					   CUT_BIT(JET_VETO_20) |
					   CUT_BIT(JET_VETO_30) 
					   );

//lep cuts
cuts_t baselepcuts  = ( CUT_BIT(LEP_PT_20) | CUT_BIT(LEP_GOOD) | CUT_BIT(LEP_ISO) ); //FOR BASE LEP COUNTING
cuts_t baselepnopt  = ( CUT_BIT(LEP_GOOD) | CUT_BIT(LEP_ISO) );


//w cuts
cuts_t w_basecuts   = ( baselepcuts | CUT_BIT(TCMET_20) | CUT_BIT(TMASS_100) | CUT_BIT(JET_VETO_20) ); //define this as 'baseline w selection'
cuts_t w_basenopt   = ( baselepnopt | CUT_BIT(TCMET_20) | CUT_BIT(TMASS_100) | CUT_BIT(JET_VETO_20) ); //don't use this for selection

cuts_t w_cuts_pt25  = ( baselepnopt | CUT_BIT(TCMET_20) | CUT_BIT(TMASS_100) | CUT_BIT(JET_VETO_20) | CUT_BIT(LEP_PT_25) ); //tighter pt
cuts_t w_cuts_pt30  = ( baselepnopt | CUT_BIT(TCMET_20) | CUT_BIT(TMASS_100) | CUT_BIT(JET_VETO_20) | CUT_BIT(LEP_PT_30) ); //tightest pt

cuts_t w_cuts_met25 = ( baselepcuts | CUT_BIT(TCMET_25) | CUT_BIT(TMASS_100) | CUT_BIT(JET_VETO_20) ); //baseline but change met 
cuts_t w_cuts_met30 = ( baselepcuts | CUT_BIT(TCMET_30) | CUT_BIT(TMASS_100) | CUT_BIT(JET_VETO_20) ); //baseline but change met 
cuts_t w_cuts_met35 = ( baselepcuts | CUT_BIT(TCMET_35) | CUT_BIT(TMASS_100) | CUT_BIT(JET_VETO_20) ); //baseline but change met 

cuts_t w_cuts_jet20 = ( baselepcuts | CUT_BIT(TCMET_20) | CUT_BIT(TMASS_100) | CUT_BIT(JET_VETO_20) ); //same as w_basecuts for now
cuts_t w_cuts_jet30 = ( baselepcuts | CUT_BIT(TCMET_20) | CUT_BIT(TMASS_100) | CUT_BIT(JET_VETO_30) ); 
cuts_t w_cuts_nojet = ( baselepcuts | CUT_BIT(TCMET_20) | CUT_BIT(TMASS_100) );


//z cuts
cuts_t z_basecuts    = ( CUT_BIT(DILEP_PT) | CUT_BIT(DILEP_MASS) | CUT_BIT(DILEP_OS) );
cuts_t z_basenopt    = ( CUT_BIT(DILEP_MASS) | CUT_BIT(DILEP_OS) );
//the LEP_PT_X cut is the pt of the LOOSE leg so that after turning off all Z cuts, can select leps with this
cuts_t z_cuts_pt1020 = ( baselepnopt | z_basenopt | CUT_BIT(LEP_PT_10) | CUT_BIT(DILEP_PT_1020) | CUT_BIT(JET_VETO_20) ); 
cuts_t z_cuts_pt1030 = ( baselepnopt | z_basenopt | CUT_BIT(LEP_PT_10) | CUT_BIT(DILEP_PT_1030) | CUT_BIT(JET_VETO_20) );
cuts_t z_cuts_pt2025 = ( baselepnopt | z_basenopt | CUT_BIT(LEP_PT_20) | CUT_BIT(DILEP_PT_2025) | CUT_BIT(JET_VETO_20) );
cuts_t z_cuts_pt2030 = ( baselepnopt | z_basenopt | CUT_BIT(LEP_PT_20) | CUT_BIT(DILEP_PT_2030) | CUT_BIT(JET_VETO_20) );

cuts_t z_cuts_jet20  = ( baselepcuts | z_basecuts | CUT_BIT(JET_VETO_20) );
cuts_t z_cuts_jet30  = ( baselepcuts | z_basecuts | CUT_BIT(JET_VETO_30) );
cuts_t z_cuts_nojet  = ( baselepcuts | z_basecuts );
