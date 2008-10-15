#include <unistd.h>
#include "selections.h"
#include "CMS2.h"
#include "WWLooper.h"

WWTopEstimateLooper::WWTopEstimateLooper (Sample s, const double *eff,
					  cuts_t c, cuts_t muc, const char *fname) 
     : WWLooperBase(s, c, fname),
       mutag_cuts(muc)
{
     memcpy(mutag_eff, eff, sizeof(mutag_eff));
     memset(top_prediction, 0, sizeof(top_prediction));
     memset(top_prediction_w2, 0, sizeof(top_prediction_w2));
}

void WWTopEstimateLooper::Dilep (int i_hyp)
{
     const enum DileptonHypType myType = hyp_typeToHypType(cms2.hyp_type()[i_hyp]);
     const double weight = Weight(i_hyp);
     cuts_t cuts_passed = DilepSelect(i_hyp);
     // need to invert the muon veto
     if ((cuts_passed & cuts) != cuts) 
	  return;
     if ((cuts_passed & mutag_cuts) != mutag_cuts)
	  return;
     double vetoed_top_pred = 1 / mutag_eff[myType]; // (1 - mutag_eff[myType]) / mutag_eff[myType];
     top_prediction[myType] += weight * vetoed_top_pred;
     top_prediction_w2[myType] += weight * weight * vetoed_top_pred * vetoed_top_pred;
     top_prediction[DILEPTON_ALL] += weight * vetoed_top_pred;
     top_prediction_w2[DILEPTON_ALL] += weight * weight * vetoed_top_pred * vetoed_top_pred;
}
