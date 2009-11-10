
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

	ELEID_TAS_ECAL,
	ELEID_TAS_DETAIN,
	ELEID_TAS_DPHIIN,
	ELEID_TAS_HOE,
	ELEID_TAS_LSHAPE,

	CONTROL_STUDYID,
	CONTROL_STUDYW,
};

// baseline cuts
const static cuts_t event_cuts = 
	(CUT_BIT(EVT_NELE_ONE)  ) ;

const static cuts_t eleid_tasElectron_v1 = 
	(CUT_BIT(ELEID_TAS_ECAL)) |
	(CUT_BIT(ELEID_TAS_DETAIN)) |
	(CUT_BIT(ELEID_TAS_DPHIIN)) |
	(CUT_BIT(ELEID_TAS_HOE)) |
	(CUT_BIT(ELEID_TAS_LSHAPE));

namespace ele {

	static cuts_t tasElectron_v1(int eleIndex)
	{

		cuts_t ret = 0;
		// ecal driven electron
		if ( (cms2.els_type()[eleIndex] & (1<<ISECALDRIVEN))) ret |= CUT_BIT(ELEID_TAS_ECAL);

		// eb
		if ((cms2.els_fiduciality()[eleIndex] & (1<<ISEB)) == (1<<ISEB)) {

			if (fabs(cms2.els_dEtaIn()[eleIndex]) < 0.007) ret |= CUT_BIT(ELEID_TAS_DETAIN);
			if (fabs(cms2.els_dPhiIn()[eleIndex]) < 0.02) ret |= CUT_BIT(ELEID_TAS_DPHIIN);
			if (cms2.els_hOverE()[eleIndex] < 0.01) ret |= CUT_BIT(ELEID_TAS_HOE);
			if (cms2.els_e2x5Max()[eleIndex]/cms2.els_e5x5()[eleIndex] > 0.9) ret |= CUT_BIT(ELEID_TAS_LSHAPE);
			return ret;
		}

		// ee
		else if ((cms2.els_fiduciality()[eleIndex] & (1<<ISEE)) == (1<<ISEE)) {

			if (fabs(cms2.els_dEtaIn()[eleIndex]) < 0.01) ret |= CUT_BIT(ELEID_TAS_DETAIN);
			if (fabs(cms2.els_dPhiIn()[eleIndex]) < 0.025) ret |= CUT_BIT(ELEID_TAS_DPHIIN);
			if (cms2.els_hOverE()[eleIndex] < 0.01) ret |= CUT_BIT(ELEID_TAS_HOE);
			if (cms2.els_sigmaIEtaIEta()[eleIndex] < 0.03) ret |= CUT_BIT(ELEID_TAS_LSHAPE);
			return ret;
		}

		return ret;
	}

}

#endif

