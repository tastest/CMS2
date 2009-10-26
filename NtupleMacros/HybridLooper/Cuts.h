
#ifndef CUTS_H
#define CUTS_H

// CMS2 related 
#include "CORE/CMS2.h"

// electron fiduciality flags
#include "../../NtupleMaker/interface/EgammaFiduciality.h"

enum {
	EVT_NELE_ONE,
	ELE_PT_20,
	ELE_PT_30,
	ELE_ISO_15,
	ELE_ISO_10,
	ELE_ISO_V0,
	ELE_ISO_V1,
	ELE_ISO_V2,
	EVT_TCMET_30,
	EVT_TCMET_20,
	EVT_JPT_25,
	EVT_JPT_PHIMAX_100,
	EVT_JPT_PHIMAX_110,
	EVT_JPT_PHIMAX_130,
	ELE_NOCONV,
	ELE_TAS_V1,
};

// baseline cuts
const static cuts_t event_cuts = 
(CUT_BIT(EVT_NELE_ONE)  ) ;

namespace ele {

	static bool tasElectron_v1(int eleIndex)
	{

		// ecal driven electron
		if ( (cms2.els_type()[eleIndex] & (1<<ISECALDRIVEN)) == 0 ) return false;

		// eb
		if ((cms2.els_fiduciality()[eleIndex] & (1<<ISEB))) {

			if (fabs(cms2.els_dEtaIn()[eleIndex]) > 0.007) return false;
			if (fabs(cms2.els_dPhiIn()[eleIndex]) > 0.02) return false;
			if (cms2.els_hOverE()[eleIndex] > 0.01) return false;
			if (cms2.els_e2x5Max()[eleIndex]/cms2.els_e5x5()[eleIndex] < 0.9) return false;

			return true;
		}

		// ee
		else if ((cms2.els_fiduciality()[eleIndex] & (1<<ISEE))) {
			if (fabs(cms2.els_dEtaIn()[eleIndex]) > 0.01) return false;
			if (fabs(cms2.els_dPhiIn()[eleIndex]) > 0.025) return false;
			if (cms2.els_hOverE()[eleIndex] > 0.01) return false;
			if (cms2.els_sigmaIEtaIEta()[eleIndex] > 0.03) return false;
			return true;
		}

		return false;
	}

}

#endif

