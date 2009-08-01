#include <math.h>
#include "TVector3.h"
#include "CORE/utilities.h"
#include "CORE/CMS2.h"
#include "Tools/tools.h"
#include "Looper.h"

#include "ElectronId.h"

Looper::Looper (Sample s, cuts_t c, const char *fname) 
     : LooperBase(s, c, fname)
{
     // zero out the candidate counters (don't comment this out)
     memset(cands_passing_	, 0, sizeof(cands_passing_       ));
     memset(cands_passing_w2_	, 0, sizeof(cands_passing_w2_    ));
     memset(cands_count_		, 0, sizeof(cands_count_         ));
}

void Looper::FormatHist(TH1* hist)
{
	hist->SetFillColor(sample_.histo_color);
}

void Looper::BookHistos ()
{

     // book histograms the manual way:

	for (unsigned int i = 0; i < 2; ++i)
	{

		std::string det = "eb";
		if (i == 1) det = "ee";

     		h1_pt_[i] = new TH1F(Form("%s_%s_%s", SampleName().c_str(), "h1_pt", det.c_str()), 
			"h1_pt", 100, 0.0, 100.0);
     		FormatHist(h1_pt_[i]);
	}

	

}


bool Looper::FilterEvent()
{ 
  return false; 
}

cuts_t Looper::EventSelect ()
{
     cuts_t ret = 0;

	// only one electron
	if (cms2.evt_nels() == 1)
		ret |= CUT_BIT(EVT_NELE_ONE);

	return ret;

}

void Looper::FillEventHistos ()
{

	// get the event weight
	float weight = cms2.evt_scale1fb() * sample_.kFactor;

	// did the event pass the baseline cuts
	cuts_t cuts_passed = EventSelect();
	if ((cuts_passed & cuts_) == cuts_) {

                // determine what detector the electron is in
                unsigned int det = 0;
                if (fabs(cms2.els_etaSC()[0]) > 1.5) det = 1;
		h1_pt_[det]->Fill(cms2.els_p4()[0].Pt(), weight);

	} // end event level cuts passed

}

void Looper::End ()
{
     int ret = fprintf(logfile_, 
		       "Sample %10s: Total candidate count (ee em mm all): %8u %8u %8u %8u."
		       " Total weight %10.1f +- %10.1f %10.1f +- %10.1f %10.1f +- %10.1f %10.1f +- %10.1f\n",   
		       sample_.name.c_str(),
		       CandsCount(DILEPTON_EE), CandsCount(DILEPTON_EMU), CandsCount(DILEPTON_MUMU), CandsCount(DILEPTON_ALL), 
		       CandsPassing(DILEPTON_EE)  , RMS(DILEPTON_EE),  
		       CandsPassing(DILEPTON_EMU) , RMS(DILEPTON_EMU),  
		       CandsPassing(DILEPTON_MUMU), RMS(DILEPTON_MUMU), 
		       CandsPassing(DILEPTON_ALL) , RMS(DILEPTON_ALL));
     if (ret < 0)
	  perror("writing to log file");

}

