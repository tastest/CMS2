#include <math.h>
#include "TVector3.h"
//#include "CORE/selections.h"
#include "CORE/utilities.h"
#include "CORE/CMS2.h"
#include "Tools/tools.h"
#include "Looper.h"
#include "TMath.h"

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

void Looper::FormatHist(TH1F** hist, std::string name, Int_t n, Float_t min, Float_t max)
{
	// loop on EB, EE
	for (unsigned int i = 0; i < 2; ++i)
	{
		std::string det = "eb";
		if (i == 1) det = "ee";
		hist[i] = new TH1F(Form("%s_%s_%s", SampleName().c_str(), name.c_str(), det.c_str()), 
				name.c_str(), n, min, max);
		hist[i]->SetFillColor(sample_.histo_color);
		hist[i]->GetXaxis()->SetTitle(name.c_str());
	}
}

void Looper::FormatEffHist(EffMulti** hist, bool lessThan, float thresholdEB, float thresholdEE, std::string name)
{
	// loop on EB, EE
	for (unsigned int i = 0; i < 2; ++i)
	{
		std::string det = "eb";
		if (i == 1) det = "ee";
		hist[i] = new EffMulti(lessThan, thresholdEB, thresholdEE, SampleName(), name, det);
	}
}

void Looper::BookHistos ()
{

	// General
	//
	FormatHist(h1_pt_, "pt", 100, 0, 100);
	FormatHist(h1_eta_, "eta", 100, -3, 3);
	FormatHist(h1_phi_, "phi", 100, -4, 4);
	FormatHist(h1_wwIsoAll_, "wwIsoAll", 100, 0.0, 1.0);

	// for "W eff studies"
	FormatHist(h1_weff_pt_, "weff_pt", 100, 0, 100);
	FormatHist(h1_weff_iso_, "weff_iso", 100, 0, 1);
	FormatHist(h1_weff_tcmet_, "weff_tcmet", 100, 0, 100);
	FormatHist(h1_weff_jptpt_, "weff_jptpt", 100, 0, 100);
	FormatHist(h1_weff_tcmet_after_iso_, "weff_tcmet_after_iso", 100, 0, 100);
        FormatHist(h1_weff_jptpt_after_iso_, "weff_jptpt_after_iso", 100, 0, 100);
        FormatHist(h1_weff_tcmet_after_iso_jetveto_, "weff_tcmet_after_iso_jetveto", 100, 0, 100);
	
	// Isolation related
	//
	FormatHist(h1_ecalIso03_, "ecalIso03", 100, 0, 1);
	FormatHist(h1_hcalIso03_, "hcalIso03", 100, 0, 1);
	FormatHist(h1_tkIso03_, "tkIso03", 100, 0, 1);
	FormatHist(h1_wwIso_, "wwIso", 100, 0, 1);

	// electron ID related
	//
	FormatHist(h1_dEtaIn_, "dEtaIn", 100, 0.0, 0.1);
	FormatHist(h1_dPhiIn_, "dPhiIn", 100, 0.0, 0.1);
	FormatHist(h1_hoe_, "hoe", 100, 0.0, 0.2);
	FormatHist(h1_sigmaIEtaIEta_, "sigmaIEtaIEta", 100, 0.0, 0.1);
	FormatHist(h1_sigmaIPhiIPhi_, "sigmaIPhiIPhi", 100, 0.0, 0.1);
	FormatHist(h1_E2x5Norm5x5_, "E2x5Norm5x5", 100, 0.0, 1.0);
	FormatHist(h1_E1x5Norm5x5_, "E1x5Norm5x5", 100, 0.0, 1.0);
	FormatHist(h1_eopIn_, "eopIn", 100, 0.0, 5.0);
	FormatHist(h1_d0corr_, "d0corr", 100, -0.2, 0.2);
	FormatHist(h1_closestMuon_, "closestMuon", 100, -1, 5);

	// The "Egamma robust tight V1 (2_2_X tune)"
	//
	FormatEffHist(em_dEtaIn_, true, 0.0040, 0.0066, "dEtaIn");
	FormatEffHist(em_dPhiIn_, true, 0.025, 0.020, "dPhiIn");
	FormatEffHist(em_hoe_, true, 0.01, 0.01, "hoe");
	FormatEffHist(em_sieie_, true, 0.0099, 0.028, "sieie");

	FormatEffHist(em_classBasedTight_, false, 0, 0, "classBasedTight");	
	FormatEffHist(em_robustTight_, false, 0, 0, "robustTight");

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

// to decdie if to fill the EB histo (zero in the array)
// or the EE histo (one in the array)
int Looper::getSubdet(int eleIndex)
{
	// determine what detector the electron is in
	if (cms2.els_fiduciality()[eleIndex] & (1<<ISEB)) return 0;
	else if (cms2.els_fiduciality()[eleIndex] & (1<<ISEE)) return 1;
	return -1;
}

void Looper::wEfficiency()
{

        // get the event weight (for 1 pb^{-1})
        float weight = cms2.evt_scale1fb() * sample_.kFactor;
        weight /= 1000;

        // did the event pass the baseline cuts
	// e.g. require one electron exactly
        cuts_t cuts_passed = EventSelect();
        if ((cuts_passed & cuts_) == cuts_) {

		// construct the isolation variable
                float ecalIso = cms2.els_ecalIso()[0];
                float hcalIso = cms2.els_hcalIso()[0];
                float tkIso = cms2.els_tkIso()[0];
                float isoSum = ecalIso + hcalIso + tkIso;

		// require egamma electron
		if (! cms2.els_type()[0] & (1<<ISECALDRIVEN)) return;
		int det = getSubdet(0);

                h1_weff_pt_[det]->Fill(cms2.els_p4()[0].Pt(), weight);

		// 20 GeV Pt cut
		if (cms2.els_p4()[0].Pt() < 20.0) return;

             	h1_weff_iso_[det]->Fill(isoSum / cms2.els_p4()[0].Pt(), weight);
                h1_weff_tcmet_[det]->Fill(cms2.evt_tcmet(), weight);

		float leadingJPT = 0.0;
		for (size_t j = 0; j < cms2.jpts_p4().size(); ++j)
		{
         	 	if ( TMath::Abs(cms2.jpts_p4()[j].eta()) > 2.5 ) continue;
        	  	if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.els_p4()[0], cms2.jpts_p4()[j])) < 0.4) continue;
			if (cms2.jpts_p4()[j].Et() > leadingJPT) leadingJPT = cms2.jpts_p4()[j].Et();
		}

		h1_weff_jptpt_[det]->Fill(leadingJPT, weight);

		// isolation cut
                if (isoSum / cms2.els_p4()[0].Pt() > 0.15) return;

		h1_weff_tcmet_after_iso_[det]->Fill(cms2.evt_tcmet(), weight);
                h1_weff_jptpt_after_iso_[det]->Fill(leadingJPT, weight);

		// jet veto
		if (leadingJPT > 30.0) return;
	
		h1_weff_tcmet_after_iso_jetveto_[det]->Fill(cms2.evt_tcmet(), weight);

	}

}

void Looper::FillEventHistos ()
{

	// do the W efficiency studies
	wEfficiency();


	// get the event weight (for 1 pb^{-1})
	float weight = cms2.evt_scale1fb() * sample_.kFactor;
	weight /= 1000;	

	// did the event pass the baseline cuts
	cuts_t cuts_passed = EventSelect();
	if ((cuts_passed & cuts_) == cuts_) {

		for (size_t i = 0; i < cms2.evt_nels(); ++i)
		{

			// matched to an mc electron if a signal process
			// a bit of a fudge but do something better later
			if (	
					(sample_.name == "wenu" && abs(cms2.els_mc3_id()[i]) != 11)
					// 20 GeV Pt
					//|| cms2.els_p4()[i].Pt() < 20.0
					// Is a regular 'Egamma' electron - not PFlow!
					|| (! (cms2.els_type()[i] & (1<<ISECALDRIVEN))) )
				continue;


			// determine what detector the electron is in
			int det = getSubdet(i);

			// construct isolation variable
			float ecalIso = cms2.els_ecalIso()[i];
			float hcalIso = cms2.els_hcalIso()[i];
			float tkIso = cms2.els_tkIso()[i];
			float isoSum = ecalIso + hcalIso + tkIso;

			// electron id related
			//
			h1_pt_[det]->Fill(cms2.els_p4()[i].Pt(), weight);
			h1_eta_[det]->Fill(cms2.els_etaSC()[i], weight);
			h1_phi_[det]->Fill(cms2.els_phiSC()[i], weight);
			h1_wwIsoAll_[det]->Fill(isoSum / cms2.els_p4()[i].Pt(), weight);

			if (isoSum / cms2.els_p4()[i].Pt() < 0.15) {

				if (cms2.els_p4()[i].Pt() > 20.0) {

					h1_dEtaIn_[det]->Fill(fabs(cms2.els_dEtaIn()[i]), weight);
					h1_dPhiIn_[det]->Fill(fabs(cms2.els_dPhiIn()[i]), weight);
					h1_hoe_[det]->Fill(cms2.els_hOverE()[i], weight);
					h1_sigmaIEtaIEta_[det]->Fill(cms2.els_sigmaIEtaIEta()[i], weight);
					h1_sigmaIPhiIPhi_[det]->Fill(cms2.els_sigmaIPhiIPhi()[i], weight);
					h1_E2x5Norm5x5_[det]->Fill(cms2.els_e2x5Max()[i]/cms2.els_e5x5()[i], weight);
					h1_E1x5Norm5x5_[det]->Fill(cms2.els_e1x5()[i]/cms2.els_e5x5()[i], weight);
					h1_eopIn_[det]->Fill(cms2.els_eOverPIn()[i], weight);
					h1_d0corr_[det]->Fill(cms2.els_d0corr()[i], weight);
					h1_closestMuon_[det]->Fill(cms2.els_closestMuon()[i], weight);
				}

				// Efficiency histograms for the Egamma robust tight V1 (2_2_X) tune...
				// ... note that all histograms except pt have a pt cut of 20.0
				em_dEtaIn_[det]->Fill(fabs(cms2.els_dEtaIn()[i]),
						cms2.els_p4()[i].Pt(), cms2.els_etaSC()[i], cms2.els_phiSC()[i], weight);

				em_dPhiIn_[det]->Fill(fabs(cms2.els_dPhiIn()[i]),
						cms2.els_p4()[i].Pt(), cms2.els_etaSC()[i], cms2.els_phiSC()[i], weight);

				em_hoe_[det]->Fill(fabs(cms2.els_hOverE()[i]),
						cms2.els_p4()[i].Pt(), cms2.els_etaSC()[i], cms2.els_phiSC()[i], weight);

				em_sieie_[det]->Fill(fabs(cms2.els_sigmaIEtaIEta()[i]),
						cms2.els_p4()[i].Pt(), cms2.els_etaSC()[i], cms2.els_phiSC()[i], weight);

				em_classBasedTight_[det]->Fill(cms2.els_egamma_tightId()[i],
						cms2.els_p4()[i].Pt(), cms2.els_etaSC()[i], cms2.els_phiSC()[i], weight);

				em_robustTight_[det]->Fill(cms2.els_egamma_robustTightId()[i],
						cms2.els_p4()[i].Pt(), cms2.els_etaSC()[i], cms2.els_phiSC()[i], weight);



			}

			// iso related
			//

			// catagory based tight
			if (cms2.els_egamma_tightId()[i] == 1) {

				// only for 20 GeV pT
				if (cms2.els_p4()[i].Pt() < 20.0) continue;

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
