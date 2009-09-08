#include <math.h>
#include "TVector3.h"
//#include "CORE/selections.h"
#include "CORE/utilities.h"
#include "CORE/CMS2.h"
#include "Tools/tools.h"
#include "Looper.h"

// CMS2 related
#include "../../NtupleMaker/interface/EgammaFiduciality.h"

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

		h1_sigmaIPhiIPhi_[i] = new TH1F(Form("%s_%s_%s_", SampleName().c_str(), "h1_sigmaIPhiIPhi", det.c_str()),
				"h1_sigmaIPhiIPhi", 100, 0.0, 0.1);
		FormatHist(h1_sigmaIPhiIPhi_[i]);

		h1_E2x5Norm5x5_[i] = new TH1F(Form("%s_%s_%s_", SampleName().c_str(), "h1_E2x5Norm5x5", det.c_str()),
				"h1_E2x5Norm5x5", 100, 0.0, 1.0);
		FormatHist(h1_E2x5Norm5x5_[i]);

		h1_E1x5Norm5x5_[i] = new TH1F(Form("%s_%s_%s_", SampleName().c_str(), "h1_E1x5Norm5x5", det.c_str()),
				"h1_E1x5Norm5x5", 100, 0.0, 1.0);
		FormatHist(h1_E1x5Norm5x5_[i]);

		h1_eopIn_[i] = new TH1F(Form("%s_%s_%s_", SampleName().c_str(), "h1_eopIn", det.c_str()),
				"h1_eopIn", 100, 0.0, 5.0);
		FormatHist(h1_eopIn_[i]);




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

	// get the event weight (for 1 pb^{-1})
	float weight = cms2.evt_scale1fb() * sample_.kFactor;
	weight /= 1000;	

	// did the event pass the baseline cuts
	cuts_t cuts_passed = EventSelect();
	if ((cuts_passed & cuts_) == cuts_) {

		for (size_t i = 0; i < cms2.evt_nels(); ++i)
		{

					// matched to an mc electron	
			if (		abs(cms2.els_mc_id()[i]) == 11
					// 20 GeV Pt
					&& cms2.els_p4()[i].Pt() < 20.0 
					// Is a regular 'Egamma' electron - not PFlow!
					&& (cms2.els_type()[i] & (1<<ISECALDRIVEN)));

			// determine what detector the electron is in
			int det = -1;
			if (cms2.els_fiduciality()[i] & (1<<ISEB)) det = 0;
			else if (cms2.els_fiduciality()[i] & (1<<ISEE)) det = 1;

			// check that electron is in either EE or EB
			if (det == -1) {
				std::cout << "Not fiducial: " << cms2.els_etaSC()[i] << std::endl;
				continue;
			}

			float ecalIso = cms2.els_ecalIso()[i];
			float hcalIso = cms2.els_hcalIso()[i];
			float tkIso = cms2.els_tkIso()[i];
			float isoSum = ecalIso + hcalIso + tkIso;

			// electron id related
			//
			if (isoSum / cms2.els_p4()[i].Pt() < 0.15) {

				h1_dEtaIn_[det]->Fill(fabs(cms2.els_dEtaIn()[i]), weight);
				h1_dPhiIn_[det]->Fill(fabs(cms2.els_dPhiIn()[i]), weight);
				h1_hoe_[det]->Fill(cms2.els_hOverE()[i], weight);
				h1_sigmaIEtaIEta_[det]->Fill(cms2.els_sigmaIEtaIEta()[i], weight);
				h1_sigmaIPhiIPhi_[det]->Fill(cms2.els_sigmaIPhiIPhi()[i], weight);
				h1_E2x5Norm5x5_[det]->Fill(cms2.els_e2x5Max()[i]/cms2.els_e5x5()[i], weight);
				h1_E1x5Norm5x5_[det]->Fill(cms2.els_e1x5()[i]/cms2.els_e5x5()[i], weight);
				h1_eopIn_[det]->Fill(cms2.els_eOverPIn()[i], weight);
			}

			// iso related
			//

			// catagory based tight
			if (cms2.els_egamma_tightId()[i] == 1) {

				h1_pt_[det]->Fill(cms2.els_p4()[i].Pt(), weight);
				h1_eta_[det]->Fill(cms2.els_p4()[i].Eta(), weight);

				h1_ecalIso03_[det]->Fill(ecalIso/cms2.els_p4()[i].Pt(), weight);
				h1_hcalIso03_[det]->Fill(hcalIso/cms2.els_p4()[i].Pt(), weight);
				h1_tkIso03_[det]->Fill(tkIso/cms2.els_p4()[i].Pt(), weight);
				h1_wwIso_[det]->Fill(isoSum / cms2.els_p4()[i].Pt(), weight);
			}

		} // end loop on electrons

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

