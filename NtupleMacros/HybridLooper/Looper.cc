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
        memset(AN2009_98Events_passing_, 0, sizeof(AN2009_98Events_passing_));
        memset(AN2009_98Events_passing_w2_, 0, sizeof(AN2009_98Events_passing_w2_));
        memset(AN2009_98Events_count_, 0, sizeof(AN2009_98Events_count_));

	memset(wEvents_passing_, 0, sizeof(wEvents_passing_));
	memset(wEvents_passing_w2_, 0, sizeof(wEvents_passing_w2_));
	memset(wEvents_count_, 0, sizeof(wEvents_count_));

	memset(cands_passing_	, 0, sizeof(cands_passing_       ));
	memset(cands_passing_w2_	, 0, sizeof(cands_passing_w2_    ));
	memset(cands_count_		, 0, sizeof(cands_count_         ));
}

void Looper::Format2DHist(TH2F** hist, std::string name, Int_t nx, Float_t minx, Float_t maxx,
				Int_t ny, Float_t miny, Float_t maxy)
{       
        // loop on EB, EE
        for (unsigned int i = 0; i < 2; ++i)
        {
                std::string det = "eb";
                if (i == 1) det = "ee";

                hist[i] = new TH2F(Form("%s_%s_%s", SampleName().c_str(), name.c_str(), det.c_str()),
                                name.c_str(), nx, minx, maxx, ny, miny, maxy);
        }
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

	// for "W eff studies"
	FormatHist(h1_weff_pt_, "weff_pt", 100, 0, 100);
	FormatHist(h1_weff_iso_, "weff_iso", 100, 0, 1);
	FormatHist(h1_weff_tcmet_, "weff_tcmet", 100, 0, 100);
	FormatHist(h1_weff_jptpt_, "weff_jptpt", 100, 0, 100);
	FormatHist(h1_weff_tcmet_after_iso_, "weff_tcmet_after_iso", 100, 0, 100);
	FormatHist(h1_weff_jptpt_after_iso_, "weff_jptpt_after_iso", 100, 0, 100);
	FormatHist(h1_weff_leadjptphi_after_iso_, "weff_leadjptphi_after_iso", 100, 0, 180);
	FormatHist(h1_weff_jptphimax_after_iso_, "weff_jptphimax_after_iso", 100, 0, 180);
        FormatHist(h1_weff_d0corr_after_iso_, "weff_d0corr_after_iso", 100, 0, 0.2);
        FormatHist(h1_weff_d0corr_after_iso_jpt_, "weff_d0corr_after_iso_jpt", 100, 0, 0.2);

	FormatHist(h1_weff_tcmet_after_iso_jpt_, "weff_tcmet_after_iso_jpt", 100, 0, 100);
	FormatHist(h1_weff_leadjptphi_after_iso_jpt_, "weff_leadjptphi_after_iso_jpt", 100, 0, 180);
	FormatHist(h1_weff_leadjptphi_after_iso_jpt_tcmet_, "weff_leadjptphi_after_iso_jpt_tcmet", 100, 0, 180);

	FormatHist(h1_weffs_sigmaIEtaIEta_, "weffs_sigmaIEtaIEta", 100, 0.0, 0.06);
	FormatHist(h1_weffbg_sigmaIEtaIEta_, "weffbg_sigmaIEtaIEta", 100, 0.0, 0.06);

        FormatHist(h1_weff_tcmet_after_iso_jpt_conv_, "weff_tcmet_after_iso_jpt_conv", 100, 0, 100);

	// electron ID related
	//
	FormatHist(h1_dEtaIn_, "dEtaIn", 100, 0.0, 0.04);
	FormatHist(h1_dPhiIn_, "dPhiIn", 100, 0.0, 0.15);
	FormatHist(h1_dPhiInSigned_, "dPhiInSigned", 200, -0.15, 0.15);
	FormatHist(h1_hoe_, "hoe", 100, 0.0, 0.2);
	FormatHist(h1_sigmaIEtaIEta_, "sigmaIEtaIEta", 100, 0.0, 0.06);
	FormatHist(h1_sigmaIPhiIPhi_, "sigmaIPhiIPhi", 100, 0.0, 0.04);
	FormatHist(h1_E2x5Norm5x5_, "E2x5Norm5x5", 100, 0.4, 1.0);
	FormatHist(h1_E1x5Norm5x5_, "E1x5Norm5x5", 100, 0.0, 1.0);
	FormatHist(h1_eopIn_, "eopIn", 100, 0.0, 5.0);
	FormatHist(h1_d0corr_, "d0corr", 100, 0.0, 0.2);
	FormatHist(h1_closestMuon_, "closestMuon", 100, -1, 5);

        // Isolation
        //
        FormatHist(h1_wwIsoAll_, "wwIsoAll", 100, 0.0, 1.0);
        FormatHist(h1_tkIso03All_, "tkIso03All", 100, 0.0, 10);
        FormatHist(h1_ecalIso03All_, "ecalIso03All", 100, 0.0, 10);
        FormatHist(h1_hcalIso03All_, "hcalIso03All", 100, 0.0, 10);

	Format2DHist(h2_tkIso03All_, "tkIso03All2D", 30, 0.0, 150.0, 30, 0.0, 15.0);
        Format2DHist(h2_ecalIso03All_, "ecalIso03All2D", 30, 0.0, 150.0, 30, 0.0, 15.0);
        Format2DHist(h2_hcalIso03All_, "hcalIso03All2D", 30, 0.0, 150.0, 30, 0.0, 15.0);

	// N-1
        FormatHist(h1_tkIso03AllNM1_, "tkIso03AllNM1", 100, 0.0, 10);
        FormatHist(h1_ecalIso03AllNM1_, "ecalIso03AllNM1", 100, 0.0, 10);
        FormatHist(h1_hcalIso03AllNM1_, "hcalIso03AllNM1", 100, 0.0, 10);


	// The "Egamma robust tight V1 (2_2_X tune)"
	//
	FormatEffHist(em_dEtaIn_, true, 0.0040, 0.0066, "dEtaIn");
	FormatEffHist(em_dPhiIn_, true, 0.025, 0.020, "dPhiIn");
	FormatEffHist(em_hoe_, true, 0.01, 0.01, "hoe");
	FormatEffHist(em_sieie_, true, 0.0099, 0.028, "sieie");

	FormatEffHist(em_classBasedTight_, false, 0, 0, "classBasedTight");	
	FormatEffHist(em_robustTight_, false, 0, 0, "robustTight");

	FormatEffHist(em_eopInLT30_, true, 3.0, 3.0, "eopInLT30");
	FormatEffHist(em_eopInGT05_, true, 3.0, 3.0, "eopInGT05");

	// AN2009-98 related
	//
        FormatHist(h1_AN2009_098_pt2_, "AN2009_098_pt2", 100, 0, 100);
        FormatHist(h1_AN2009_098_eta1_, "AN2009_098_eta1", 100, -3, 3);

        FormatHist(h1_AN2009_098_ecalIso_, "AN2009_098_ecalIso", 100, 0, 10.0);
        FormatHist(h1_AN2009_098_hcalIso_, "AN2009_098_hcalIso", 100, 0, 10.0);
        FormatHist(h1_AN2009_098_tkIso_, "AN2009_098_tkIso", 100, 0, 10.0);

        FormatHist(h1_AN2009_098_tcmet_, "AN2009_098_tcmet", 100, 0, 100);
        FormatHist(h1_AN2009_098_tcmet_after_selection_, "AN2009_098_tcmet_after_selection", 100, 0, 100);

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

// to decdie if to fill the EB histo (zero in the array)
// or the EE histo (one in the array)
int Looper::getSubdet(int eleIndex)
{
	// determine what detector the electron is in
	if (cms2.els_fiduciality()[eleIndex] & (1<<ISEB)) return 0;
	else if (cms2.els_fiduciality()[eleIndex] & (1<<ISEE)) return 1;
	return -1;
}

void Looper::AN2009_98()
{

        // get the event weight (for 1 pb^{-1})
        float weight = cms2.evt_scale1fb() * sample_.kFactor;
        weight /= 1000;

	// find candidate electron
	// assumes electrons are sorted by pT descending
	int eleIndex = 0;
	int eleSecondIndex = 0;
	bool foundFirst = false;
	bool foundSecond = false;
	for (size_t i = 0; i < cms2.evt_nels(); ++i)
	{
		// no particle flow
	        if (! cms2.els_type()[i] & (1<<ISECALDRIVEN)) continue;

		if (foundFirst && !foundSecond) {
			eleSecondIndex = i;
			foundSecond = true;
			break;
		}
                if (!foundFirst) {
                        eleIndex = i;
                        foundFirst = true;
                }
	}

	// get isolation
        const float &ecalIso = cms2.els_ecalIso()[0];
        const float &hcalIso = cms2.els_hcalIso()[0];
        const float &tkIso = cms2.els_tkIso()[0];

	// get subdetector (for histogramming)
        int det = getSubdet(eleIndex);

	//
	// Apply cuts here
	//

	// require highest pT electron to have pT > 30.0 GeV
	//
	if (cms2.els_p4()[eleIndex].Pt() < 30.0) return;
	//

	if (foundSecond) h1_AN2009_098_pt2_[det]->Fill(cms2.els_p4()[eleSecondIndex].Pt(), weight);
	else h1_AN2009_098_pt2_[det]->Fill(0.0, weight);

        h1_AN2009_098_eta1_[det]->Fill(cms2.els_etaSC()[eleIndex], weight);

	// don't allow events with a second electron above 20.0 GeV
	//
	if (foundSecond && cms2.els_p4()[eleSecondIndex].Pt() > 20.0) return;
	//

	// impose fiducial cuts in Eta
	// veto barrel-endcap gap
	if (fabs(cms2.els_etaSC()[eleIndex]) > 1.4442 && fabs(cms2.els_etaSC()[eleIndex]) < 1.560) return;
	// veto high eta endcap
	if (fabs(cms2.els_etaSC()[eleIndex]) > 2.500) return;
	//

        h1_AN2009_098_tkIso_[det]->Fill(tkIso, weight);
        h1_AN2009_098_ecalIso_[det]->Fill(ecalIso, weight);
        h1_AN2009_098_hcalIso_[det]->Fill(hcalIso, weight);
	h1_AN2009_098_tcmet_[det]->Fill(cms2.evt_tcmet(), weight);

	// isolations cuts
	//
        // some candidate thresholds to look at N-1 with
        float tkThresholds[2] = {2.2, 1.1};
        float ecalThresholds[2] = {4.2, 3.4};
        float hcalThresholds[2] = {2.0, 1.3};

	if (tkIso > tkThresholds[det]) return;
	if (ecalIso > ecalThresholds[det]) return;
	if (hcalIso > hcalThresholds[det]) return;

	else std::cout << "Error! Not in barrel or endcap - something is wrong" << std::endl;
	//

	h1_AN2009_098_tcmet_after_selection_[det]->Fill(cms2.evt_tcmet(), weight);

	if (cms2.evt_tcmet() < 30.0) return;

        AN2009_98Events_passing_[det] += weight;
        AN2009_98Events_passing_w2_[det] += weight;
        AN2009_98Events_count_[det] ++;
        AN2009_98Events_passing_[2] += weight;
        AN2009_98Events_passing_w2_[2] += weight;
        AN2009_98Events_count_[2] ++;

}


void Looper::wEfficiency()
{

	// get the event weight (for 1 pb^{-1})
	float weight = cms2.evt_scale1fb() * sample_.kFactor;
	weight /= 1000;

	//
	// apply preselection cuts are always the same
	//

	// one egamma electron
	if (cms2.evt_nels() != 1) return;
	// 
	// electron is egamma type
	if (! cms2.els_type()[0] & (1<<ISECALDRIVEN)) return;
	//
	// get the subdestector (EE or EB)
	int det = getSubdet(0);
	//
	// 20 GeV Pt cut
	h1_weff_pt_[det]->Fill(cms2.els_p4()[0].Pt(), weight);
	if (cms2.els_p4()[0].Pt() < 20.0) return;

	//
	// construct variables that are not already in ntuple
	//

	// isolation
	const float &ecalIso = cms2.els_ecalIso()[0];
	const float &hcalIso = cms2.els_hcalIso()[0];
	const float &tkIso = cms2.els_tkIso()[0];
	float isoSum = ecalIso + hcalIso + tkIso;

	// leading pT JPT jet
	// that is dR > 0.4 from the nearest electron
	float leadingJPT = 0.0;
	float mostBackToBackJPT = 0.0;
	int leadingJPTIndex = 0;
	int mostBackToBackJPTIndex = 0;
	for (size_t j = 0; j < cms2.jpts_p4().size(); ++j)
	{
		if ( TMath::Abs(cms2.jpts_p4()[j].eta()) > 2.5 ) continue;
		if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.els_p4()[0], cms2.jpts_p4()[j])) < 0.4) continue;
		if (cms2.jpts_p4()[j].Et() > leadingJPT) {
			leadingJPT = cms2.jpts_p4()[j].Et();
			leadingJPTIndex = j;
		}
		float dPhi = acos(cos(cms2.els_p4()[0].Phi() - cms2.jpts_p4()[j].Phi()));
		if (dPhi > mostBackToBackJPT) {
			mostBackToBackJPT = dPhi;
			mostBackToBackJPTIndex = j;
		}

	}

	// distance in degrees of the leading JPT from the electron
	float leadingJPTAngle = 0.0;
	float mostBackToBackJPTAngle = 0.0;
	if (leadingJPT > 0.0) {
		// angle between the electron and the highest pT JPT
		leadingJPTAngle = (180.0 * acos(cos(cms2.els_p4()[0].Phi()
						- cms2.jpts_p4()[leadingJPTIndex].Phi())) / 3.14159265358979312);
		// angle between the electron and the JPT that is most back to back to it
		mostBackToBackJPTAngle = (180.0 * acos(cos(cms2.els_p4()[0].Phi()
						- cms2.jpts_p4()[mostBackToBackJPTIndex].Phi())) / 3.14159265358979312);	
	}

	// 
	// set up the selections to apply
	//

	float isolationThreshold = 0;
	if ((cuts_ & (CUT_BIT(ELE_ISO_15))) == (CUT_BIT(ELE_ISO_15))) isolationThreshold = 0.15;
	if ((cuts_ & (CUT_BIT(ELE_ISO_10))) == (CUT_BIT(ELE_ISO_10))) isolationThreshold = 0.10;
	float jptThreshold = 999999.99;
	if ((cuts_ & (CUT_BIT(EVT_JPT_25))) == (CUT_BIT(EVT_JPT_25))) jptThreshold = 25.0;
	float tcMetThreshold = 0;
	if ((cuts_ & (CUT_BIT(EVT_TCMET_30))) == (CUT_BIT(EVT_TCMET_30))) tcMetThreshold = 30.0;
	if ((cuts_ & (CUT_BIT(EVT_TCMET_20))) == (CUT_BIT(EVT_TCMET_20))) tcMetThreshold = 20.0;
	float jptPhiThreshold = 180.0;
        if ((cuts_ & (CUT_BIT(EVT_JPT_PHIMAX_100))) == (CUT_BIT(EVT_JPT_PHIMAX_100))) jptPhiThreshold = 100.0;
        if ((cuts_ & (CUT_BIT(EVT_JPT_PHIMAX_100))) == (CUT_BIT(EVT_JPT_PHIMAX_100))) jptPhiThreshold = 110.0;
	bool rejectConversions = false;
	if ((cuts_ & (CUT_BIT(ELE_NOCONV))) == (CUT_BIT(ELE_NOCONV))) rejectConversions = true;

	//
	// plots of key quantities before selection applied
	h1_weff_iso_[det]->Fill(isoSum / cms2.els_p4()[0].Pt(), weight);
	h1_weff_tcmet_[det]->Fill(cms2.evt_tcmet(), weight);
	h1_weff_jptpt_[det]->Fill(leadingJPT, weight);

	//
	// apply selections
	//

	// isolation cut
	if (isoSum / cms2.els_p4()[0].Pt() > isolationThreshold) return;
	//

	// plots of tcmet and leading jpt pt
	h1_weff_tcmet_after_iso_[det]->Fill(cms2.evt_tcmet(), weight);
	h1_weff_jptpt_after_iso_[det]->Fill(leadingJPT, weight);
	// angle between leading JPT and electron
	h1_weff_leadjptphi_after_iso_[det]->Fill(leadingJPTAngle, weight); 
	// angle between electron and JPT that is most back to back to it
	h1_weff_jptphimax_after_iso_[det]->Fill(mostBackToBackJPTAngle, weight);
	// d0 corrected
	h1_weff_d0corr_after_iso_[det]->Fill(fabs(cms2.els_d0corr()[0]), weight);

	// leading JPT cut
	if (leadingJPT > jptThreshold) return;
	//

        // JPT phimax cut
        if (mostBackToBackJPTAngle > jptPhiThreshold) return;
        //


	// plot of tcmet after isolation and jpt veto applied
	h1_weff_tcmet_after_iso_jpt_[det]->Fill(cms2.evt_tcmet(), weight);
	// plot of angle between leading JPT and electron
	h1_weff_leadjptphi_after_iso_jpt_[det]->Fill(leadingJPTAngle, weight);
        // d0 corrected
        h1_weff_d0corr_after_iso_jpt_[det]->Fill(fabs(cms2.els_d0corr()[0]), weight);

	// conversion rejection cut
	if (rejectConversions && isconversionElectron09(0)) return;
	//

	h1_weff_tcmet_after_iso_jpt_conv_[det]->Fill(cms2.evt_tcmet(), weight);

	// BACKGROUND
	// plots of electron ID quantities
	//
        if (cms2.evt_tcmet() < 15.0) {
		h1_weffbg_sigmaIEtaIEta_[det]->Fill(cms2.els_sigmaIEtaIEta()[0], weight); 
	}

	// tcMet cut
	if (cms2.evt_tcmet() > tcMetThreshold) {
	//

		// plot of angle between leading JPT and electron
		h1_weff_leadjptphi_after_iso_jpt_tcmet_[det]->Fill(leadingJPTAngle, weight);

		// SIGNAL
		// plots of electron ID quantities
		//
		h1_weffs_sigmaIEtaIEta_[det]->Fill(cms2.els_sigmaIEtaIEta()[0], weight);	

		// add to count of events passing all cuts
		wEvents_passing_[det] += weight;
		wEvents_passing_w2_[det] += weight;
		wEvents_count_[det] ++;
		wEvents_passing_[2] += weight;
		wEvents_passing_w2_[2] += weight;
		wEvents_count_[2] ++;

	}

}

void Looper::FillEventHistos ()
{

	// do the W efficiency studies
	wEfficiency();

	// do the AN2009-98 studies
	AN2009_98();

	// get the event weight (for 1 pb^{-1})
	float weight = cms2.evt_scale1fb() * sample_.kFactor;
	weight /= 1000;	

	// did the event pass the baseline cuts
	//cuts_t cuts_passed = EventSelect();
	//if ((cuts_passed & cuts_) == cuts_) {
	if (cms2.evt_nels() == 1) {


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

			// for isolation optimisation
			//
			// some candidate thresholds to look at N-1 with
			float tkThresholdsNM1[2] = {4.5, 6.0};
			float ecalThresholdsNM1[2] = {2.5, 2.0};
			float hcalThresholdsNM1[2] = {1.0, 1.0};
			//
			if (cms2.els_p4()[i].Pt() > 20.0) {
				h1_wwIsoAll_[det]->Fill(isoSum / cms2.els_p4()[i].Pt(), weight);
				h1_ecalIso03All_[det]->Fill(ecalIso, weight);
				h1_hcalIso03All_[det]->Fill(hcalIso, weight);
				h1_tkIso03All_[det]->Fill(tkIso, weight);

				h2_tkIso03All_[det]->Fill(cms2.els_p4()[i].Pt(), tkIso, weight);
                                h2_ecalIso03All_[det]->Fill(cms2.els_p4()[i].Pt(), ecalIso, weight);
                                h2_hcalIso03All_[det]->Fill(cms2.els_p4()[i].Pt(), hcalIso, weight);

				// N-1
				if (ecalIso < ecalThresholdsNM1[det] && hcalIso < hcalThresholdsNM1[det])
					h1_tkIso03AllNM1_[det]->Fill(tkIso, weight);
                                if (tkIso < tkThresholdsNM1[det] && hcalIso < hcalThresholdsNM1[det])
					h1_ecalIso03AllNM1_[det]->Fill(ecalIso, weight);
                                if (tkIso < tkThresholdsNM1[det] && ecalIso < ecalThresholdsNM1[det])
					h1_hcalIso03AllNM1_[det]->Fill(hcalIso, weight);

			}

			if (isoSum / cms2.els_p4()[i].Pt() < 0.15) {

				if (cms2.els_p4()[i].Pt() > 20.0) {

					h1_dEtaIn_[det]->Fill(fabs(cms2.els_dEtaIn()[i]), weight);
					h1_dPhiIn_[det]->Fill(fabs(cms2.els_dPhiIn()[i]), weight);
					h1_dPhiInSigned_[det]->Fill(cms2.els_dPhiIn()[i]*cms2.els_charge()[i], weight);
					h1_hoe_[det]->Fill(cms2.els_hOverE()[i], weight);
					h1_sigmaIEtaIEta_[det]->Fill(cms2.els_sigmaIEtaIEta()[i], weight);
					h1_sigmaIPhiIPhi_[det]->Fill(cms2.els_sigmaIPhiIPhi()[i], weight);
					h1_E2x5Norm5x5_[det]->Fill(cms2.els_e2x5Max()[i]/cms2.els_e5x5()[i], weight);
					h1_E1x5Norm5x5_[det]->Fill(cms2.els_e1x5()[i]/cms2.els_e5x5()[i], weight);
					h1_eopIn_[det]->Fill(cms2.els_eOverPIn()[i], weight);
					h1_d0corr_[det]->Fill(fabs(cms2.els_d0corr()[i]), weight);
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

				em_eopInGT05_[det]->Fill(cms2.els_eOverPIn()[i],
						cms2.els_p4()[i].Pt(), cms2.els_etaSC()[i], cms2.els_phiSC()[i], weight);

				em_eopInLT30_[det]->Fill(cms2.els_eOverPIn()[i],
						cms2.els_p4()[i].Pt(), cms2.els_etaSC()[i], cms2.els_phiSC()[i], weight);

				//em_tasElectronV1_[det]->Fill(tasElectron_v1(i),
				//                cms2.els_p4()[i].Pt(), cms2.els_etaSC()[i], cms2.els_phiSC()[i], weight);


			}

		} // end loop on electrons

	} // end event level cuts passed

}

void Looper::End ()
{

	int ret = fprintf(logfile_,
			"Sample %10s: Total candidate count (EB EE ALL): %8u %8u %8u"
			" Total weight %10.1f +- %10.1f %10.1f +- %10.1f %10.1f +- %10.1f \n",
			sample_.name.c_str(),
			CandsCountW(0), CandsCountW(1), CandsCountW(2),
			CandsPassingW(0)  , RMSW(0),
			CandsPassingW(1) , RMSW(1),
			CandsPassingW(2) , RMSW(2));

	if (ret < 0)
		perror("HybridLooper: writing w study to log file");

        ret = fprintf(logfile_,
                        "Sample %10s: Total candidate count (EB EE ALL): %8u %8u %8u"
                        " Total weight %10.1f +- %10.1f %10.1f +- %10.1f %10.1f +- %10.1f \n",
                        sample_.name.c_str(),
                        CandsCountAN2009_98(0), CandsCountAN2009_98(1), CandsCountAN2009_98(2),
                        CandsPassingAN2009_98(0)  , RMSAN2009_98(0),
                        CandsPassingAN2009_98(1) , RMSAN2009_98(1),
                        CandsPassingAN2009_98(2) , RMSAN2009_98(2));

        if (ret < 0)
                perror("HybridLooper: writing AN2009_98 study to log file");


	ret = fprintf(logfile_, 
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
