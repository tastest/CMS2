
#ifndef CUTS_H
#define CUTS_H

enum {
	//EL_GOOD,
	//MU_GOOD,
	DILEP_PT,
	DILEP_OS,
	DILEP_MASS,
	LEP_PT,
	LEP_GOOD,
	LEP_ISO,
	LEP_GOOD_NOD0,
	LEP_D0,
	JET_VETO_20, //if njet_20 > 0, this cut fails
	JET_VETO_30,
	TCMET,
	TMASS, //transverse mass for W selection
	ZMASS,
	//DUMMY,
};

// baseline cuts
// placeholder 
// this does not do anything
//const static cuts_t event_cuts = 
//(CUT_BIT(DUMMY)) ;



#endif

