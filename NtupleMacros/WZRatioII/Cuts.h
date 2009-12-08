#ifndef CUTS_H
#define CUTS_H

#include "Tools/LooperTypes.h"
#include "Tools/LooperBase.h"

enum {
  //lepton selection cuts
  LEP_PT,          //kept for backward compatibility with ABCDLooper--don't use in RLooper nor Cuts.cc
  LEP_PT_10,       //10 -- loose--use for counting, but not W selection
  LEP_PT_20,       //20 -- nominal
  LEP_PT_25,       //25 -- tighter
  LEP_PT_30,       //30 -- tightest
  LEP_ISO,         //0.1 -- nominal cut on susy_rel_iso
  LEP_ISO_TIGHT,   //0.08 -- tighter
  LEP_GOOD,
  LEP_GOOD_NOD0,
  LEP_D0,
  //Z cuts
  DILEP_PT,      //20,20 -- nominal nominal
  DILEP_PT_1020, //10,20 -- loose nominal
  DILEP_PT_1030, //10,30 -- loose tightest
  DILEP_PT_2025, //20,25 -- nominal tighter
  DILEP_PT_2030, //20,30 -- nominal tightest
  DILEP_OS,
  DILEP_MASS,      //76-106
  DILEP_MASS_WIDE, //71-111
  //W cuts
  TCMET,       //kept for backward compatibility with ABCDLooper--don't use in RLooper nor Cuts.cc
  TCMET_20,
  TCMET_25,
  TCMET_30,
  TCMET_35,
  TMASS_100, //transverse mass for W selection
  TMASS_120,
  //jet cuts
  JET_VETO_20, //if njet_20 > 0, this cut fails
  JET_VETO_30,
};

//this one is special--used to switch off W specific cuts (anything not checked in LepSelect) for use as baselepcuts in RLooper
// w_cuts has lep_cuts + met, tmass. non_baselepcuts has all W+jet cuts.
// ie, baselepcuts = (w_cuts & ~non_baselepcuts);
extern cuts_t non_lepcuts;

extern cuts_t baselepcuts ;// = ( CUT_BIT(LEP_PT) | CUT_BIT(LEP_GOOD) | CUT_BIT(LEP_ISO) ); //FOR BASE LEP COUNTING
//extern cuts_t w_evt_sel   ;// = ( CUT_BIT(TCMET) | CUT_BIT(TMASS_100) | CUT_BIT(JET_VETO_20) );
extern cuts_t w_basecuts;
extern cuts_t w_cuts_pt25;
extern cuts_t w_cuts_pt30;
extern cuts_t w_cuts_jet20;// = ( CUT_BIT(LEP_PT) | CUT_BIT(LEP_GOOD) | w_evt_sel | CUT_BIT(LEP_ISO) | CUT_BIT(JET_VETO_20) ); //ALL W CUTS
extern cuts_t w_cuts_jet30;// = ( CUT_BIT(LEP_PT) | CUT_BIT(LEP_GOOD) | w_evt_sel | CUT_BIT(LEP_ISO) | CUT_BIT(JET_VETO_30) ); //ALL W CUTS
extern cuts_t w_cuts_nojet;// = ( CUT_BIT(LEP_PT) | CUT_BIT(LEP_GOOD) | w_evt_sel | CUT_BIT(LEP_ISO) ); 
extern cuts_t w_cuts_met25;
extern cuts_t w_cuts_met30;
extern cuts_t w_cuts_met35;
				   
extern cuts_t z_basecuts  ;// = ( CUT_BIT(DILEP_PT) | CUT_BIT(DILEP_MASS) | CUT_BIT(DILEP_OS) );
extern cuts_t z_cuts_pt1020;
extern cuts_t z_cuts_pt1030;
extern cuts_t z_cuts_pt2025;
extern cuts_t z_cuts_pt2030;
extern cuts_t z_cuts_jet20;// = ( CUT_BIT(LEP_PT) | CUT_BIT(LEP_GOOD) | z_basecuts | CUT_BIT(LEP_ISO) | CUT_BIT(JET_VETO_20) ); //ALL Z CUTS
extern cuts_t z_cuts_jet30;// = ( CUT_BIT(LEP_PT) | CUT_BIT(LEP_GOOD) | z_basecuts | CUT_BIT(LEP_ISO) | CUT_BIT(JET_VETO_30) ); //ALL Z CUTS
extern cuts_t z_cuts_nojet;// = ( CUT_BIT(LEP_PT) | CUT_BIT(LEP_GOOD) | z_basecuts | CUT_BIT(LEP_ISO) );


#endif

