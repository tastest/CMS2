#include "../CORE/selections.h"
#include "../CORE/CMS2.h"
#include "WWLooper.h"

WWDumboLooper::WWDumboLooper (Sample s, cuts_t c, const char *fname) 
     : WWLooperBase(s, c, fname)
{
     
}

void WWDumboLooper::Dilep (int i_hyp)
{
     cuts_t cuts_passed = DilepSelect(i_hyp);
     
     if ((cuts_passed & cuts) != cuts)
	  return;
     
     FillHists(i_hyp);
}
