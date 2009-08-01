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

	// single lepton histograms (two + 1types)
	for (unsigned int i = 0; i < 3; ++i)
	{
		std::string hyp = "e";
		if (i == 1) hyp = "m";
		if (i == 2) hyp = "all";

     		h1_lep_pt_[i] = new TH1F(Form("%s_%s_%s", SampleName().c_str(), "lep_pt", hyp.c_str()), 
			"lep_pt", 100, 0.0, 100.0);
     		FormatHist(h1_lep_pt_[i]);

                h1_lep_met_[i] = new TH1F(Form("%s_%s_%s", SampleName().c_str(), "lep_met", hyp.c_str()),
                        "lep_met", 100, 0.0, 100.0);
                FormatHist(h1_lep_met_[i]);


	}

	// di-lepton histograms (three + 1 types)
        for (unsigned int i = 0; i < 4; ++i)
        {
                h1_dilep_0_pt_[i] = new TH1F(Form("%s_%s_%s", SampleName().c_str(), "dilep_0_pt", dilepton_hypo_names[i]),
                        "dilep_0_pt", 100, 0.0, 100.0);
                FormatHist(h1_dilep_0_pt_[i]);

                h1_dilep_1_pt_[i] = new TH1F(Form("%s_%s_%s", SampleName().c_str(), "dilep_1_pt", dilepton_hypo_names[i]),
                        "dilep_1_pt", 100, 0.0, 100.0);
                FormatHist(h1_dilep_1_pt_[i]);

        }
	

}


bool Looper::FilterEvent()
{ 
  return false; 
}

cuts_t Looper::EventSelect ()
{
     cuts_t ret = 0;

	// is this a single lepton event
	// or a dilepton event?
	// also which flavour


	if (cms2.evt_nels() > 1 || cms2.mus_p4().size() > 1)
		ret |= CUT_BIT(EVT_DILEP);

	if ((cms2.evt_nels() == 0 && cms2.mus_p4().size() == 1) ||
		(cms2.evt_nels() == 1 && cms2.mus_p4().size() == 0))
		ret |= CUT_BIT(EVT_LEP);

	return ret;

}

void Looper::FillEventHistos ()
{

	// need to determine if this is a di-lepton
	// or a single lepton event

	cuts_t cuts_passed = EventSelect();
	if ((cuts_passed & CUT_BIT(EVT_LEP)) == CUT_BIT(EVT_LEP)) WEvent();
        if ((cuts_passed & CUT_BIT(EVT_DILEP)) == CUT_BIT(EVT_LEP)) ZEvent();

}

void Looper::WEvent ()
{

        // get the event weight
        float weight = cms2.evt_scale1fb() * sample_.kFactor; 

	// histogram indices are e, m, all (0, 1, 2)
	unsigned int hyp = 1;
	if (cms2.mus_p4().size() == 0) {
		hyp = 0;
		h1_lep_pt_[hyp]->Fill(cms2.els_p4()[0].Pt(), weight);
                h1_lep_pt_[2]->Fill(cms2.els_p4()[0].Pt(), weight);
		
		h1_lep_met_[hyp]->Fill(cms2.evt_tcmet(), weight);
                h1_lep_met_[2]->Fill(cms2.evt_tcmet(), weight);

	}

        if (cms2.evt_nels() == 0) {
                hyp = 1;
                h1_lep_pt_[hyp]->Fill(cms2.mus_p4()[0].Pt(), weight);
                h1_lep_pt_[2]->Fill(cms2.mus_p4()[0].Pt(), weight);

                h1_lep_met_[hyp]->Fill(cms2.evt_tcmet(), weight);
                h1_lep_met_[2]->Fill(cms2.evt_tcmet(), weight);
        }



}

void Looper::ZEvent ()
{

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

