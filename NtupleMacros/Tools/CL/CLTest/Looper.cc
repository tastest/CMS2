#include <math.h>
#include "TVector3.h"
#include "CORE/electronSelections.h"
#include "CORE/utilities.h"
#include "CORE/CMS2.h"
#include "Tools/tools.h"
#include "Looper.h"

Looper::Looper (Sample s, cuts_t c, const char *fname) 
     : LooperBase(s, c, fname)
{

}

void Looper::BookHistos ()
{
     hdilMass_		= new TH1F("m_ee", "m_{ee};m_{ee}", 150, 0, 150);
}


cuts_t Looper::DilepSelect (int i_hyp)
{
     //------------------------------------------------------------
     // Example dilepton cuts; edit for your application
     //------------------------------------------------------------

     // cuts are failed until proven otherwise
     cuts_t ret = 0;

     // sign cuts
     if ( cms2.hyp_lt_id()[i_hyp] * cms2.hyp_ll_id()[i_hyp] < 0 ) 
	  ret |= (CUT_BIT(CUT_OPP_SIGN));
     if ( cms2.hyp_lt_id()[i_hyp] * cms2.hyp_ll_id()[i_hyp] > 0 ) 
	  ret |= (CUT_BIT(CUT_SAME_SIGN));
     // goodness cuts
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11 && electronSelection_cand01(cms2.hyp_lt_index()[i_hyp]))
	  ret |= CUT_BIT(CUT_LT_GOOD);
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 11 && electronSelection_cand01(cms2.hyp_ll_index()[i_hyp]))
	  ret |= CUT_BIT(CUT_LL_GOOD);

     return ret;
}
void Looper::FillDilepHistos (int i_hyp)
{
     // these are the cuts that the candidate passes:
     cuts_t cuts_passed = DilepSelect(i_hyp);

     // this is how to test that the candidate passes the cuts (which
     // we specified in the constructor when we made the looper)
     // (*note: the parentheses are important*):
     if ((cuts_passed & cuts_) == cuts_) {
	  if (cms2.hyp_p4()[i_hyp].M2() > 0)
	       hdilMass_->Fill(cms2.hyp_p4()[i_hyp].M());
     }
}

void Looper::End ()
{
     //------------------------------------------------------------
     //Example status message at the end of a looper; edit for your
     //application
     //------------------------------------------------------------
     hdilMass_->Draw();
}
