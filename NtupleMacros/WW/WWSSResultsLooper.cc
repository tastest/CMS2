#include "selections.h"
#include "CMS2.h"
#include "WWLooper.h"

WWSSResultsLooper::WWSSResultsLooper (Sample s, cuts_t cuts, const char *log) 
     : WWLooperBase(s, cuts, log)
{
     
}

void WWSSResultsLooper::Dilep (int i_hyp)
{
     cuts_t cuts_passed = DilepSelect(i_hyp);
     
     if ((cuts_passed & cuts) != cuts)
	  return;
     
     FillHists(i_hyp);
}
