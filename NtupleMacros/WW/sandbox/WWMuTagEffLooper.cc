#include <unistd.h>
#include "../CORE/selections.h"
#include "../CORE/CMS2.h"
#include "WWLooper.h"

WWMuTagEffLooper::WWMuTagEffLooper (Sample s, cuts_t c, cuts_t muc, const char *fname) 
     : WWLooperBase(s, c, fname),
       mutag_cuts(muc)
{
     memset(cand_twojet, 0, sizeof(cand_twojet));
     memset(cand_twojet_w2, 0, sizeof(cand_twojet_w2));
     memset(cand_twojet_mutagged, 0, sizeof(cand_twojet_mutagged));
     memset(cand_twojet_mutagged_w2, 0, sizeof(cand_twojet_mutagged_w2));
}

void WWMuTagEffLooper::Dilep (int i_hyp)
{
     const enum DileptonHypType myType = hyp_typeToHypType(cms2.hyp_type()[i_hyp]);
     const double weight = Weight(i_hyp);
     cuts_t cuts_passed = DilepSelect(i_hyp);
     if ((cuts_passed & cuts) != cuts) 
	  return;
     if (cms2.hyp_njets()[i_hyp] != 2)
	  return;
     cand_twojet[myType] += weight;
     cand_twojet_w2[myType] += weight * weight;
     cand_twojet[DILEPTON_ALL] += weight;
     cand_twojet_w2[DILEPTON_ALL] += weight * weight;
     if ((cuts_passed & mutag_cuts) != mutag_cuts) // no extra muon
	  return;
     cand_twojet_mutagged[myType] += weight;
     cand_twojet_mutagged_w2[myType] += weight * weight;
     cand_twojet_mutagged[DILEPTON_ALL] += weight;
     cand_twojet_mutagged_w2[DILEPTON_ALL] += weight * weight;
}

void WWMuTagEffLooper::End ()
{
     for (int i = 0; i < 4; ++i) 
	  mutag_eff[i] = cand_twojet_mutagged[i] / cand_twojet[i];
     int ret = fprintf(logfile, 
		       "Sample %10s: (ee emu mumu all)\n"
		       " Njet == 2            %10.1f +- %10.1f %10.1f +- %10.1f %10.1f +- %10.1f %10.1f +- %10.1f\n"   
		       " Njet == 2, mu-tagged %10.1f +- %10.1f %10.1f +- %10.1f %10.1f +- %10.1f %10.1f +- %10.1f\n",   
		       sample.name.c_str(),
		       CandTwoJet(DILEPTON_EE)  , CandTwoJetRMS(DILEPTON_EE),  
		       CandTwoJet(DILEPTON_EMU) , CandTwoJetRMS(DILEPTON_EMU),  
		       CandTwoJet(DILEPTON_MUMU), CandTwoJetRMS(DILEPTON_MUMU), 
		       CandTwoJet(DILEPTON_ALL) , CandTwoJetRMS(DILEPTON_ALL),
		       CandTwoJetMuTagged(DILEPTON_EE)  , CandTwoJetMuTaggedRMS(DILEPTON_EE),  
		       CandTwoJetMuTagged(DILEPTON_EMU) , CandTwoJetMuTaggedRMS(DILEPTON_EMU),  
		       CandTwoJetMuTagged(DILEPTON_MUMU), CandTwoJetMuTaggedRMS(DILEPTON_MUMU), 
		       CandTwoJetMuTagged(DILEPTON_ALL) , CandTwoJetMuTaggedRMS(DILEPTON_ALL));
     if (ret < 0)
	  perror("writing to log file");
     if (logfile != stdout)
	  fclose(logfile);
}

WWMuTagEffLooper WWMuTagEffLooper::operator + (const WWMuTagEffLooper &other) const
{
     WWMuTagEffLooper ret = *this;
     for (int i = 0; i < 4; ++i) {
	  ret.cand_twojet[i]		 = cand_twojet[i]             + other.cand_twojet[i]            ;
	  ret.cand_twojet_mutagged[i]	 = cand_twojet_mutagged[i]    + other.cand_twojet_mutagged[i]   ;
	  ret.cand_twojet_w2[i]		 = cand_twojet_w2[i]          + other.cand_twojet_w2[i]         ;
	  ret.cand_twojet_mutagged_w2[i] = cand_twojet_mutagged_w2[i] + other.cand_twojet_mutagged_w2[i];
	  ret.mutag_eff[i] = ret.cand_twojet_mutagged[i] / ret.cand_twojet[i];
     }
     return ret;
}
