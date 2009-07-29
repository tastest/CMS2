#include <math.h>
#include "TVector3.h"
//#include "CORE/selections.h"
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

     // book histograms the manual way:

	for (unsigned int i = 0; i < 2; ++i)
	{

		std::string det = "eb";
		if (i == 1) det = "ee";

     		h1_pt_[i] = new TH1F(Form("%s_%s_%s_", SampleName().c_str(), "h1_pt", det.c_str()), 
			"h1_pt", 100, 0.0, 100.0);
     		FormatHist(h1_pt_[i]);

     		h1_eta_[i] = new TH1F(Form("%s_%s_%s_", SampleName().c_str(), "h1_eta", det.c_str()), 
			"h1_eta", 60, -3.0, 3.0);
     		FormatHist(h1_eta_[i]);

     		h1_ecalIso03_[i] = new TH1F(Form("%s_%s_%s_", SampleName().c_str(), "h1_ecalIso03", det.c_str()), 
			"h1_ecalIso03", 100, 0.0, 1.0);
     		FormatHist(h1_ecalIso03_[i]);

     		h1_hcalIso03_[i] = new TH1F(Form("%s_%s_%s_", SampleName().c_str(), "h1_hcalIso03", det.c_str()), 
			"h1_hcalIso03", 100, 0.0, 0.5);
     		FormatHist(h1_hcalIso03_[i]);

    		h1_tkIso03_[i] = new TH1F(Form("%s_%s_%s_", SampleName().c_str(), "h1_tkIso03", det.c_str()), 
			"h1_tkIso03", 100, 0.0, 0.5);
     		FormatHist(h1_tkIso03_[i]);

                h1_esJuraIso03_[i] = new TH1F(Form("%s_%s_%s_", SampleName().c_str(), "h1_esJuraIso03", det.c_str()), 
			"h1_esJuraIso03", 100, 0.0, 1.0);
                FormatHist(h1_esJuraIso03_[i]);

     		h1_wwIso_[i] = new TH1F(Form("%s_%s_%s_", SampleName().c_str(), "h1_wwIso", det.c_str()), 
			"h1_wwIso", 100, 0.0, 1.0);
     		FormatHist(h1_wwIso_[i]);


		//
		// electron ID related
		//

		h1_dEtaIn_[i] = new TH1F(Form("%s_%s_%s_", SampleName().c_str(), "h1_dEtaIn", det.c_str()),
			"h1_dEtaIn", 100, 0.0, 0.1);
                FormatHist(h1_dEtaIn_[i]);

                h1_dPhiIn_[i] = new TH1F(Form("%s_%s_%s_", SampleName().c_str(), "h1_dPhiIn", det.c_str()),
                        "h1_dPhiIn", 100, 0.0, 0.1);
                FormatHist(h1_dPhiIn_[i]);

                h1_hoe_[i] = new TH1F(Form("%s_%s_%s_", SampleName().c_str(), "h1_hoe", det.c_str()),
                        "h1_hoe", 100, 0.0, 0.2);
                FormatHist(h1_hoe_[i]);

                h1_sigmaIEtaIEta_[i] = new TH1F(Form("%s_%s_%s_", SampleName().c_str(), "h1_sigmaIEtaIEta", det.c_str()),
                        "h1_sigmaIEtaIEta", 100, 0.0, 0.1);
                FormatHist(h1_sigmaIEtaIEta_[i]);


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

cuts_t Looper::DilepSelect (int i_hyp)
{
     cuts_t ret = 0;
     return ret;
}

cuts_t Looper::TrilepSelect (int i_hyp)
{
     cuts_t ret = 0;
     return ret;
}

cuts_t Looper::QuadlepSelect (int i_hyp)
{
     cuts_t ret = 0;
     return ret;
}

void Looper::FillEventHistos ()
{

	// get the event weight
	float weight = cms2.evt_scale1fb() * sample_.kFactor;

	// did the event pass the baseline cuts
	cuts_t cuts_passed = EventSelect();
	if ((cuts_passed & cuts_) == cuts_) {

		// 20 GeV electrons in the barrel and avoid -ve endcap
		if (cms2.els_p4()[0].Pt() < 20.0 || cms2.els_etaSC()[0] < -1.5) return;

		// determine what detector the electron is in
		unsigned int det = 0;
		if (cms2.els_etaSC()[0] > 1.5) det = 1;

		// 3_1_X
		//float ecalIso = cms2.els_tkIso03()[0];
		//float hcalIso = cms2.els_tkIso03()[0];
		//float tkIso = cms2.els_tkIso03()[0];
		// 2_2_X
                float ecalIso = cms2.els_ecalIso()[0];
                float hcalIso = cms2.els_hcalIso()[0];
                float tkIso = cms2.els_tkIso()[0];
                float isoSum = ecalIso + hcalIso + tkIso;

		// electron id related
		//
		if (isoSum / cms2.els_p4()[0].Pt() < 0.15) {

			h1_dEtaIn_[det]->Fill(fabs(cms2.els_dEtaIn()[0]), weight);
			h1_dPhiIn_[det]->Fill(fabs(cms2.els_dPhiIn()[0]), weight);
			h1_hoe_[det]->Fill(cms2.els_hOverE()[0], weight);
			h1_sigmaIEtaIEta_[det]->Fill(cms2.els_sigmaIEtaIEta()[0], weight);
		}

		// iso related
		//

		// min matteo tight
		if (cms2.els_tightId22XMinMatteo()[0] == 1) {

	                h1_pt_[det]->Fill(cms2.els_p4()[0].Pt(), weight);
			h1_eta_[det]->Fill(cms2.els_p4()[0].Eta(), weight);

			h1_ecalIso03_[det]->Fill(ecalIso/cms2.els_p4()[0].Pt(), weight);
	                h1_hcalIso03_[det]->Fill(hcalIso/cms2.els_p4()[0].Pt(), weight);
	                h1_tkIso03_[det]->Fill(tkIso/cms2.els_p4()[0].Pt(), weight);
			// for 3_1_X
        	        //h1_esJuraIso03_[det]->Fill(cms2.els_esJuraIso03()[0]/cms2.els_p4()[0].Pt(), weight);
			h1_wwIso_[det]->Fill(isoSum / cms2.els_p4()[0].Pt(), weight);

		}

	} // end event level cuts passed

}

void Looper::FillDilepHistos (int i_hyp)
{
}

void Looper::FillTrilepHistos (int i_hyp)
{
}

void Looper::FillQuadlepHistos (int i_hyp)
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

