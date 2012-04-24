#ifndef MVAMET_H
#define MVAMET_H

#include "../Tools/MVAMET/MVAMet.h" 
#include "jetcorr/FactorizedJetCorrector.h"

float mvaMET( int hyp_index, double deltaZCut = 0.2, MVAMet* mvaMet = 0, 
			  FactorizedJetCorrector *jet_corrector = 0 );

#endif

