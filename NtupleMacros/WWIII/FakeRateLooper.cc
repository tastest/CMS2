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

     const double weight_hi = Weight(i_hyp, 1);
     const double weight_lo = Weight(i_hyp, -1);
     cands_passing_syst_hi[myType] += weight_hi;
     cands_passing_syst_lo[myType] += weight_lo;
     cands_passing_syst_hi[DILEPTON_ALL] += weight_hi;
     cands_passing_syst_lo[DILEPTON_ALL] += weight_lo;
     const double weight = Looper::Weight(i_hyp);
     // this doesn't work for ee, since it assumes only one of
     // the hyp leptons is an electron
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11) {
	  const double err = elFakeProb(cms2.hyp_lt_index()[i_hyp], 1) - 
	       elFakeProb(cms2.hyp_lt_index()[i_hyp], 0);
	  const double eta = cms2.els_p4()[cms2.hyp_lt_index()[i_hyp]].eta();
	  const double pt = cms2.els_p4()[cms2.hyp_lt_index()[i_hyp]].pt();
	  fake_syst->Fill(eta, pt, weight * err);
     } else if (abs(cms2.hyp_ll_id()[i_hyp]) == 11) {
	  const double err = elFakeProb(cms2.hyp_ll_index()[i_hyp], 1) - 
	       elFakeProb(cms2.hyp_ll_index()[i_hyp], 0);
	  const double eta = cms2.els_p4()[cms2.hyp_ll_index()[i_hyp]].eta();
	  const double pt = cms2.els_p4()[cms2.hyp_ll_index()[i_hyp]].pt();
	  fake_syst->Fill(eta, pt, weight * err);
     }
     Looper::FillDilepHistos(i_hyp);
}

double FakeRateLooper::Weight (int i_hyp)
{ 
     return Weight(i_hyp, 0);
}

double FakeRateLooper::Weight (int i_hyp, int n_sig_syst)
{
     double weight = Looper::Weight(i_hyp);
     double fr = 0;
     // this doesn't work for ee, since it assumes only one of
     // the hyp leptons is an electron
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11) {
	  fr = elFakeProb(cms2.hyp_lt_index()[i_hyp], n_sig_syst);
     } else if (abs(cms2.hyp_ll_id()[i_hyp]) == 11) {
	  fr = elFakeProb(cms2.hyp_ll_index()[i_hyp], n_sig_syst);
     }
     
     return weight * fr / (1 - fr); 
}

double FakeRateLooper::FakeSyst (enum DileptonHypType i) const
{
     switch (i) {
     case DILEPTON_EE: case DILEPTON_EMU:
	  return 0;
     default:
	  break;
     }
     double err2 = 0;
     for (int i = 0; i <= fake_syst->GetNbinsX() + 1; ++i) {
	  for (int j = 0; j <= fake_syst->GetNbinsY() + 1; ++j) {
	       double err = fake_syst->GetBinContent(i, j);
	       err2 += err * err;
	  }
     }
     return sqrt(err2);
}
