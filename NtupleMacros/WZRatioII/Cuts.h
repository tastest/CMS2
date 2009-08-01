
#ifndef CUTS_H
#define CUTS_H

enum {
	EVT_LEP,
	EVT_DILEP,
};

// baseline cuts
// placeholder 
const static cuts_t event_cuts = 
  (CUT_BIT(EVT_LEP) | CUT_BIT(EVT_DILEP) ) ;



#endif

