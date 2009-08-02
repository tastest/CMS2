# include <math.h>
#include "TVector3.h"
#include "CORE/utilities.h"
#include "CORE/CMS2.h"
#include "Tools/tools.h"
#include "Looper.h"


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

	// single lepton histograms (two + 1 types)
	for (unsigned int i = 0; i < 3; ++i)
	{
		std::string hyp = "e";
		if (i == 1) hyp = "m";
		if (i == 2) hyp = "all";

     		h1_lep_pt_[i] = new TH1F(
			Form("%s_%s_%s", SampleName().c_str(), "lep_pt", hyp.c_str()), 
			"lep_tkIso", 100, 0.0, 100.0);
     		FormatHist(h1_lep_pt_[i]);

                h1_lep_met_[i] = new TH1F(
			Form("%s_%s_%s", SampleName().c_str(), "lep_met", hyp.c_str()),
                        "lep_met", 100, 0.0, 100.0);
                FormatHist(h1_lep_met_[i]);

                h1_lep_met_dphi_[i] = new TH1F(
                        Form("%s_%s_%s", SampleName().c_str(), "lep_met_dphi", hyp.c_str()),
                        "lep_met_dphi", 100, 0, 2 * 3.14159);
                FormatHist(h1_lep_met_dphi_[i]);

                h1_lep_tkIso_[i] = new TH1F(
			Form("%s_%s_%s", SampleName().c_str(), "lep_tkIso", hyp.c_str()),
                        "lep_tkIso", 100, 0.0, 100.0);
                FormatHist(h1_lep_tkIso_[i]);


	}

	// di-lepton histograms (three + 1 types)
        for (unsigned int i = 0; i < 4; ++i)
        {
                h1_dilep_0_pt_[i] = new TH1F(
			Form("%s_%s_%s", SampleName().c_str(), "dilep_0_pt", dilepton_hypo_names[i]),
                        "dilep_0_pt", 100, 0.0, 100.0);
                FormatHist(h1_dilep_0_pt_[i]);

                h1_dilep_1_pt_[i] = new TH1F(
			Form("%s_%s_%s", SampleName().c_str(), "dilep_1_pt", dilepton_hypo_names[i]),
                        "dilep_1_pt", 100, 0.0, 100.0);
                FormatHist(h1_dilep_1_pt_[i]);

                h1_dilep_mass_[i] = new TH1F(
                        Form("%s_%s_%s", SampleName().c_str(), "dilep_mass", dilepton_hypo_names[i]),
                        "dilep_mass", 200, 0.0, 200.0);
                FormatHist(h1_dilep_mass_[i]);

                h1_dilep_met_[i] = new TH1F(
                        Form("%s_%s_%s", SampleName().c_str(), "dilep_met", dilepton_hypo_names[i]),
                        "lep_met", 100, 0.0, 100.0);
                FormatHist(h1_dilep_met_[i]);


        }
	
	// event level histograms
	h1_dilep_nhyp_ = new TH1F(
		Form("%s_%s_%s", SampleName().c_str(), "dilep_nhyp", "all"),
		"dilep_nhyp", 10, -0.5, 9.5);
	FormatHist(h1_dilep_nhyp_);

}

cuts_t Looper::DilepSelect (int i_hyp)
{
	cuts_t ret = 0;

	float ptcut = 20.0;
	if (cms2.hyp_lt_p4()[i_hyp].pt() > ptcut && cms2.hyp_ll_p4()[i_hyp].pt() > ptcut)
		ret |= CUT_BIT(DILEP_PT);

	return ret;
}

cuts_t Looper::LepSelect(int lep_type, int i)
{
	cuts_t ret = 0;
	return ret;
}

void Looper::FillEventHistos ()
{

	// need to determine if this is a di-lepton
	// or a single lepton event

        if (cms2.evt_nels() > 1 || cms2.mus_p4().size() > 1)
		ZEvent();
        else if ((cms2.evt_nels() == 0 && cms2.mus_p4().size() == 1) ||
                (cms2.evt_nels() == 1 && cms2.mus_p4().size() == 0)) 
		WEvent();

}

void Looper::WEvent ()
{

        // get the event weight
        float weight = cms2.evt_scale1fb() * sample_.kFactor; 

	// histogram indices are e, m, all (0, 1, 2)
	unsigned int hyp = 1;
	if (cms2.mus_p4().size() == 0) hyp = 0;

	if (hyp == 0) {	
		h1_lep_pt_[hyp]->Fill(cms2.els_p4()[0].pt(), weight);
        	h1_lep_pt_[2]->Fill(cms2.els_p4()[0].pt(), weight);
		
		h1_lep_met_[hyp]->Fill(cms2.evt_tcmet(), weight);
	        h1_lep_met_[2]->Fill(cms2.evt_tcmet(), weight);

		float dphi = acos(cos(cms2.evt_tcmetPhi() - cms2.els_p4()[0].Phi() ));
	        h1_lep_met_dphi_[hyp]->Fill(dphi, weight);
        	h1_lep_met_dphi_[2]->Fill(dphi, weight);

	        h1_lep_tkIso_[hyp]->Fill(cms2.els_tkIso()[0], weight);
		h1_lep_tkIso_[2]->Fill(cms2.els_tkIso()[0], weight);
	}
	if (hyp == 1) {
                h1_lep_pt_[hyp]->Fill(cms2.mus_p4()[0].pt(), weight);
                h1_lep_pt_[2]->Fill(cms2.mus_p4()[0].pt(), weight);
                
                h1_lep_met_[hyp]->Fill(cms2.evt_tcmet(), weight);
                h1_lep_met_[2]->Fill(cms2.evt_tcmet(), weight);

                float dphi = acos(cos(cms2.evt_tcmetPhi() - cms2.mus_p4()[0].Phi() ));
                h1_lep_met_dphi_[hyp]->Fill(dphi, weight);
                h1_lep_met_dphi_[2]->Fill(dphi, weight);

                h1_lep_tkIso_[hyp]->Fill(cms2.mus_iso03_sumPt()[0], weight);
                h1_lep_tkIso_[2]->Fill(cms2.mus_iso03_sumPt()[0], weight);

	}

}

void Looper::ZEvent ()
{

        // get the event weight
        float weight = cms2.evt_scale1fb() * sample_.kFactor;
	h1_dilep_nhyp_->Fill(cms2.hyp_p4().size(), weight);

	// define the cuts to be used
	cuts_t cuts = CUT_BIT(DILEP_PT);

	for (unsigned int i = 0; i < cms2.hyp_p4().size();  ++i)
	{
		// get what type of di-lepton hypothesis this is
		const enum DileptonHypType myType = hyp_typeToHypType(cms2.hyp_type()[i]);

		// does this hypothesis pass the required cuts?
		cuts_t cuts_passed = DilepSelect(i);
		if ((cuts_passed & cuts) == cuts) {

			// fill histograms for the 0th lepton
			h1_dilep_0_pt_[myType]->Fill(cms2.hyp_lt_p4()[i].pt(), weight);
			h1_dilep_0_pt_[DILEPTON_ALL]->Fill(cms2.hyp_lt_p4()[i].pt(), weight);
		
			// fill histograms for the 1th lepton
        	        h1_dilep_1_pt_[myType]->Fill(cms2.hyp_ll_p4()[i].pt(), weight);
                	h1_dilep_1_pt_[DILEPTON_ALL]->Fill(cms2.hyp_ll_p4()[i].pt(), weight);

			// fill hypothesis level histograms
			h1_dilep_mass_[myType]->Fill(cms2.hyp_p4()[i].mass(), weight);
			h1_dilep_mass_[DILEPTON_ALL]->Fill(cms2.hyp_p4()[i].mass(), weight);

                        h1_dilep_met_[myType]->Fill(cms2.evt_tcmet(), weight);
                        h1_dilep_met_[DILEPTON_ALL]->Fill(cms2.evt_tcmet(), weight);

		}

	}


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

