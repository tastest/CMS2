#include "../CORE/selections.h"
#include "../CORE/CMS2.h"
#include "WWLooper.h"

WWResultsLooper::WWResultsLooper (Sample s, cuts_t c, const char *fname) 
     : WWLooperBase(s, c, fname)
{
     
}

void WWResultsLooper::Dilep (int i_hyp)
{
     cuts_t cuts_passed = DilepSelect(i_hyp);
//      printf("cuts passed: %x, need: %x\n", cuts_passed, ww_old_baseline_cuts);

     // eventually, we want to do the N - 1 plots 
     // ...

     FillHists(i_hyp);

     // for now, we shoot first, ask questions later
     if ((cuts_passed & cuts) != cuts)
	  return;

     // if we get an extra muon, print it out
     int tag_muon_idx = tagMuonIdx(i_hyp);
     enum DileptonHypType myType = hyp_typeToHypType(cms2.hyp_type()[i_hyp]);
     if (tag_muon_idx != -1 && myType == DILEPTON_EMU)
	  printf("tag muon: %d (out of %u), pt = %.1f (%.1f), iso = %f\n",
		 tag_muon_idx, cms2.mus_p4().size(), cms2.mus_p4()[tag_muon_idx].pt(),
		 tagMuonPt(i_hyp), tagMuonRelIso(i_hyp));
}
