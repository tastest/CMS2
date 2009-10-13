
#ifndef CUTS_H
#define CUTS_H

// CMS2 related 
#include "CORE/CMS2.h"

enum {
	EVT_NELE_ONE,
	ELE_ISO_15,
	ELE_ISO_10,
	EVT_TCMET_30,
	EVT_TCMET_20,
	EVT_JPT_25,
	EVT_JPT_PHIMAX_100,
	EVT_JPT_PHIMAX_110,
	ELE_NOCONV,
};

// baseline cuts
const static cuts_t event_cuts = 
  (CUT_BIT(EVT_NELE_ONE)  ) ;

/*
bool tasElectron_v1(int eleIndex)
{

	// ecal driven electron
	if (! cms2.els_type()[0] & (1<<ISECALDRIVEN)) return false;

	// eb
        if (cms2.els_fiduciality()[eleIndex] & (1<<ISEB)) {
		if (cms2.els_dEtaIn()[eleIndex] > 0.007) return false;
		if (cms2.els_dPhiIn()[eleIndex] > 0.02) return false;
		if (cms2.els_hOverE()[eleIndex] > 0.01) return false;
		if (cms2.els_e1x5()[eleIndex]/cms2.els_e5x5()[eleIndex] < 0.9) return false;
	}

	// ee
	if (cms2.els_fiduciality()[eleIndex] & (1<<ISEE)) {
                if (cms2.els_dEtaIn()[eleIndex] > 0.01) return false;
                if (cms2.els_dPhiIn()[eleIndex] > 0.025) return false;
                if (cms2.els_hOverE()[eleIndex] > 0.01) return false;
                if (cms2.els_sigmaIEtaIEta()[eleIndex] < 0.9) return false;
	}

	return true;
}
*/

#endif

