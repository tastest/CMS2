#include "selections.h"
#include "CMS2.h"
#include "WWLooper.h"
#include "fakerates.h"

WWFakeRateLooper::WWFakeRateLooper (Sample s, cuts_t cuts, const char *fname) 
     : WWLooperBase(s, cuts, fname)
{
     
}

cuts_t WWFakeRateLooper::DilepSelect (int i_hyp)
{
     // this doesn't really do anything special for now
     return WWLooperBase::DilepSelect(i_hyp);
}

void WWFakeRateLooper::Dilep (int i_hyp)
{
     cuts_t cuts_passed = DilepSelect(i_hyp);
     
     if ((cuts_passed & cuts) != cuts)
	  return;

     FillHists(i_hyp);
}

double WWFakeRateLooper::Weight (int i_hyp)
{ 
     double weight = WWLooperBase::Weight(i_hyp);
     double fr = 0;
     // this doesn't work for ee, since it assumes only one of
     // the hyp leptons is an electron
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11) {
	  fr = elFakeProb(cms2.hyp_lt_index()[i_hyp]);
     } else if (abs(cms2.hyp_ll_id()[i_hyp]) == 11) {
	  fr = elFakeProb(cms2.hyp_ll_index()[i_hyp]);
     }
     
     return weight * fr / (1 - fr); 
}
