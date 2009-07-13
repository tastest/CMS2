#include <math.h>
#include <sstream>
#include "TVector3.h"
#include "CORE/selections.h"
#include "CORE/utilities.h"
#include "CORE/CMS2.h"
#include "Tools/tools.h"
#include "Tools/fakerates.h"
#include "Looper.h"

FakeRateLooper::FakeRateLooper (Sample s, cuts_t c, const char *fname) 
     : Looper(s, c, fname)
{
     // zero out the candidate counters (don't comment this out)
     memset(cands_passing_syst_hi, 0, sizeof(cands_passing_syst_hi));
     memset(cands_passing_syst_lo, 0, sizeof(cands_passing_syst_lo));
}

void FakeRateLooper::BookHistos 	()
{
     fake_syst = new TH2F("fake_syst", "fake syst uncertainty;#eta;pt", 
			  fakeRate().GetNbinsX(), fakeRate().GetXaxis()->GetXbins()->GetArray(),
			  fakeRate().GetNbinsY(), fakeRate().GetYaxis()->GetXbins()->GetArray());
     Looper::BookHistos();
}

cuts_t FakeRateLooper::DilepSelect (int i_hyp)
{
     // this doesn't do anything special
     return Looper::DilepSelect(i_hyp);
}

void FakeRateLooper::FillDilepHistos (int i_hyp)
{
     const enum DileptonHypType myType = hyp_typeToHypType(cms2.hyp_type()[i_hyp]);
     cuts_t cuts_passed = DilepSelect(i_hyp);
     
     if ((cuts_passed & cuts_) != cuts_)
	  return;

     // the logic is as follows:
     //
     // first, we assume the lt to be real (by testing whether it
     // passes the reco+truth tags) and calculate the probability that the
     // ll is fake.
     //
     // then we switch it around and assume the ll to be real (by
     // testing whether it passes the reco+truth tags) and calculate the probability that the
     // lt is fake.
     //
     // this neglects double fakes (but our W+jets shouldn't contain a
     // lot of those anyway)
     //
     // it also neglects the horrible things that will happen if the
     // trigger cuts are tighter than the denominator cuts 

     double weight = Weight(i_hyp);
     
     double fr_ll = 0;
     if ((cuts_passed & (CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_LT_ISO))) == 
	 (CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_LT_ISO))) { // muon truth tag is already in global cuts
// 	  const double eta = cms2.hyp_ll_p4()[i_hyp].eta();
// 	  const double pt = cms2.hyp_ll_p4()[i_hyp].pt();
	  switch (abs(cms2.hyp_ll_id()[i_hyp])) {
	  case 11:
	       if (isFakeable(cms2.hyp_ll_index()[i_hyp]) &&
		   not isNumeratorElectron(cms2.hyp_ll_index()[i_hyp]))
		    fr_ll = elFakeProb(cms2.hyp_ll_index()[i_hyp], 0);
	       break;
 	  case 13:
// 	       if (isFakeableMuon(cms2.hyp_ll_index()[i_hyp]) &&
// 		   not isNumeratorMuon(cms2.hyp_ll_index()[i_hyp]))
// 		    fr_ll = muFakeProb(cms2.hyp_ll_index()[i_hyp], 0);
 	       break;
	  default:
	       assert(0);
	  }
     }

     double fr_lt = 0;
     if ((cuts_passed & (CUT_BIT(CUT_LL_GOOD) | CUT_BIT(CUT_LL_ISO))) == 
	 (CUT_BIT(CUT_LL_GOOD) | CUT_BIT(CUT_LL_ISO))) {
// 	  const double eta = cms2.hyp_lt_p4()[i_hyp].eta();
// 	  const double pt = cms2.hyp_lt_p4()[i_hyp].pt();
	  switch (abs(cms2.hyp_lt_id()[i_hyp])) {
	  case 11:
	       if (isFakeable(cms2.hyp_lt_index()[i_hyp]) &&
		   not isNumeratorElectron(cms2.hyp_lt_index()[i_hyp]))
		    fr_lt = elFakeProb(cms2.hyp_lt_index()[i_hyp], 0);
	       break;
	  case 13:
// 	       if (isFakeableMuon(cms2.hyp_lt_index()[i_hyp]) &&
// 		   not isNumeratorMuon(cms2.hyp_lt_index()[i_hyp]))
// 		    fr_lt = muFakeProb(cms2.hyp_lt_index()[i_hyp], 0);
	       break;
	  default:
	       assert(0);
	  }
     }

     weight *= fr_lt / (1 - fr_lt) + fr_ll / (1 - fr_ll);

     // just to make sure: only one of fr_lt and fr_ll should ever be non-zero
     printf("fr_lt = %g\tfr_ll = %g\n", fr_lt, fr_ll);
     assert(not ((fr_lt != 0) && (fr_ll != 0)));

     cands_passing_[myType] += weight;
     cands_passing_w2_[myType] += weight * weight;
     cands_count_[myType]++;
     cands_passing_[DILEPTON_ALL] += weight;
     cands_passing_w2_[DILEPTON_ALL] += weight * weight;
     cands_count_[DILEPTON_ALL]++;
}
