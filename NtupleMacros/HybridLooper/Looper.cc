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
	FormatHist(h1_weff_leastemjpt_after_iso_, "weff_leastemjpt_after_iso", 100, 0, 1.0);

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
	FormatHist(h1_E2x5Norm5x5_, "E2x5Norm5x5", 80, 0.6, 1.0);
	FormatHist(h1_E1x5Norm5x5_, "E1x5Norm5x5", 100, 0.0, 1.0);
	FormatHist(h1_eopIn_, "eopIn", 100, 0.0, 5.0);
	FormatHist(h1_d0corr_, "d0corr", 100, 0.0, 0.2);
	FormatHist(h1_closestMuon_, "closestMuon", 100, -1, 5);

	// efficiencies in pt for the tas ID, N-1 etc.
        FormatHist(h1_dEtaInTasV1NM1_, "dEtaInTasV1NM1", 100, 0.0, 0.04);	
        FormatEffHist(em_dEtaInTasV1NM1_, true, 0.007, 0.010, "dEtaInTasV1NM1");

        FormatHist(h1_dPhiInTasV1NM1_, "dPhiInTasV1NM1", 100, 0.0, 0.15);       
        FormatEffHist(em_dPhiInTasV1NM1_, true, 0.020, 0.025, "dPhiInTasV1NM1");

        FormatHist(h1_hoeTasV1NM1_, "hoeTasV1NM1", 100, 0.0, 0.2);
        FormatEffHist(em_hoeTasV1NM1_, true, 0.01, 0.01, "hoeTasV1NM1");

        FormatHist(h1_sigmaIEtaIEtaTasV1NM1_, "sigmaIEtaIEtaTasV1NM1", 100, 0.0, 0.06);
        FormatEffHist(em_sigmaIEtaIEtaTasV1NM1_, true, 0.0, 0.03, "sigmaIEtaIEtaTasV1NM1");

        FormatHist(h1_E2x5Norm5x5TasV1NM1_, "E2x5Norm5x5TasV1NM1", 80, 0.6, 1.0);
        FormatEffHist(em_E2x5Norm5x5TasV1NM1_, false, 0.9, 0.0, "E2x5Norm5x5TasV1NM1");


        // Isolation
        //
        FormatHist(h1_wwIsoAll_, "wwIsoAll", 100, 0.0, 1.0);
        FormatHist(h1_tkIso03All_, "tkIso03All", 150, 0.0, 15);
        FormatHist(h1_ecalIso03All_, "ecalIso03All", 150, 0.0, 15);
        FormatHist(h1_hcalIso03All_, "hcalIso03All", 150, 0.0, 15);
	FormatHist(h1_caloIso03All_, "caloIso03All", 150, 0.0, 15);

	Format2DHist(h2_tkIso03All_, "tkIso03All2D", 30, 0.0, 150.0, 30, 0.0, 15.0);
        Format2DHist(h2_ecalIso03All_, "ecalIso03All2D", 30, 0.0, 150.0, 30, 0.0, 15.0);
        Format2DHist(h2_hcalIso03All_, "hcalIso03All2D", 30, 0.0, 150.0, 30, 0.0, 15.0);
        Format2DHist(h2_caloIso03All_, "caloIso03All2D", 30, 0.0, 150.0, 30, 0.0, 15.0);

	// N-1
        FormatHist(h1_tkIso03AllNM1_, "tkIso03AllNM1", 150, 0.0, 15);
        FormatHist(h1_tkIso03AllIDNM1_, "tkIso03AllIDNM1", 150, 0.0, 15);
        FormatHist(h1_tkIso03AllConvNM1_, "tkIso03AllConvNM1", 150, 0.0, 15);
        FormatHist(h1_tkIso03AllConvIDNM1_, "tkIso03AllConvIDNM1", 150, 0.0, 15);

        FormatHist(h1_ecalIso03AllNM1_, "ecalIso03AllNM1", 150, 0.0, 15);
        FormatHist(h1_hcalIso03AllNM1_, "hcalIso03AllNM1", 150, 0.0, 15);
	FormatHist(h1_tkIso03AllReJura01In015NM1_, "tkIso03AllReJura01In015NM1", 150, 0.0, 15);
        FormatHist(h1_tkIso03AllReJura01In015IDNM1_, "tkIso03AllReJura01In015IDNM1", 150, 0.0, 15);
        FormatHist(h1_tkIso03AllReJura01In015ConvNM1_, "tkIso03AllReJura01In015ConvNM1", 150, 0.0, 15);
        FormatHist(h1_tkIso03AllReJura01In015ConvIDNM1_, "tkIso03AllReJura01In015ConvIDNM1", 150, 0.0, 15);
	FormatHist(h1_tkIso03AllReShCutNM1_, "tkIso03AllReShCutNM1", 150, 0.0, 15);

	// track iso studies
        FormatHist(h1_tkIso03Alld0corr_, "tkIso03Alld0corr", 100, 0.0, 0.2);
        FormatHist(h1_tkIso03AllRe_, "tkIso03AllRe", 150, 0.0, 15.0);
        FormatHist(h1_tkIso03AllReShVeto_, "tkIso03AllReShVeto", 220, -1.1, 1.1);
        FormatHist(h1_tkIso03AllReRel_, "tkIso03AllReRel", 150, 0.0, 1.5);
        FormatHist(h1_tkIso03AllReJura01_, "tkIso03AllReJura01", 150, 0.0, 15.0);
        FormatHist(h1_tkIso03AllReJura02_, "tkIso03AllReJura02", 150, 0.0, 15.0);
        FormatHist(h1_tkIso03AllReJura03_, "tkIso03AllReJura03", 150, 0.0, 15.0);
	FormatHist(h1_tkIso03AllReJura01In015_, "tkIso03AllReJura01In015", 150, 0.0, 15.0);
        FormatHist(h1_tkIso03AllRedEta_, "tkIso03AllRedEta", 80, -0.4, 0.4);
        FormatHist(h1_tkIso03AllRedPhi_, "tkIso03AllRedPhi", 80, -0.4, 0.4);
        Format2DHist(h2_tkIso03AllRedR2D_, "tkIso03AllRedR2D", 80, -0.4, 0.4, 50, -0.1, 0.4);


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

	FormatEffHist(em_tasElectronV1_, false, 0, 0, "tasElectronV1");

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

float Looper::recomputeTrackIsolation(int eleIndex, float strip, float dRIn, float dROut, float &shCutSum)
{
        
        float isoSum = 0.0;
	shCutSum = 0.0;
        for (size_t i = 0; i < cms2.trks_trk_p4().size(); ++i)
        {

                float dEta = cms2.trks_trk_p4()[i].Eta() - cms2.els_p4()[eleIndex].Eta();
                float dPhi = acos(cos(cms2.trks_trk_p4()[i].Phi() - cms2.els_p4()[eleIndex].Phi()));
                float dR = sqrt(dEta*dEta + dPhi*dPhi);
                const float &pT = cms2.trks_trk_p4()[i].Pt();

                if (pT < 0.7) continue;
                if (fabs(cms2.trks_z0()[i] - cms2.els_z0()[eleIndex]) > 0.2) continue;
                if (dR < dRIn) continue;
		if (dR > dROut) continue;

		// jurassic isolation
		if (fabs(dEta) > strip) isoSum += pT;

		// hit sharing
		// if the ctf track is the one that shares the most hits
		// only add it to the sum if the shared fraction is < 0.45.
		if (cms2.els_trkidx()[eleIndex] == i && cms2.els_trkshFrac()[eleIndex] < 0.45) shCutSum += pT;
		else shCutSum += pT;

	}

	return isoSum;

}

void Looper::trackIsolationStudy(int eleIndex, int det)
{

        // get the event weight (for 1 pb^{-1})
        float weight = cms2.evt_scale1fb() * sample_.kFactor;
        weight /= 1000;

	float isoSum = 0.0;
	float isoSumJura01 = 0.0;
        float isoSumJura02 = 0.0;
        float isoSumJura03 = 0.0;
	float isoSumJura01In015 = 0.0;

	for (size_t i = 0; i < cms2.trks_trk_p4().size(); ++i)
	{

		float dEta = cms2.trks_trk_p4()[i].Eta() - cms2.els_p4()[eleIndex].Eta();
		float dPhi = acos(cos(cms2.trks_trk_p4()[i].Phi() - cms2.els_p4()[eleIndex].Phi()));
		float dR = sqrt(dEta*dEta + dPhi*dPhi);
		const float &pT = cms2.trks_trk_p4()[i].Pt();

		if (pT < 0.7) continue;
		if (dR > 0.3) continue;
		if (fabs(cms2.trks_z0()[i] - cms2.els_z0()[eleIndex]) > 0.2) continue;

		// use weight of 1.0 when filling the 2d hist
                h2_tkIso03AllRedR2D_[det]->Fill(dEta, dPhi, 1.0);

		if (dR > 0.015 && fabs(dEta) > 0.01) {
			isoSumJura01In015 += pT;
			// if the ctf track is the one with the highest number of 
			// shared hits with the electron, then plot the number of shared hits.
			if (cms2.els_trkidx()[eleIndex] == i) 
				h1_tkIso03AllReShVeto_[det]->Fill(cms2.els_trkshFrac()[eleIndex], weight);
			else 
				h1_tkIso03AllReShVeto_[det]->Fill(-1.0, weight);
		}

                if (dR > 0.015 && fabs(dEta) < 0.01) {
                        // if the ctf track is the one with the highest number of 
                        // shared hits with the electron, then plot the number of shared hits.
                        if (cms2.els_trkidx()[eleIndex] == i)
                                h1_tkIso03AllReShVeto_[det]->Fill(cms2.els_trkshFrac()[eleIndex], weight);
                        else
                                h1_tkIso03AllReShVeto_[det]->Fill(-1.0, weight);
                }


		// standard inner veto
                if (dR < 0.04) continue;

                h1_tkIso03AllRedEta_[det]->Fill(dEta, weight);
                h1_tkIso03AllRedPhi_[det]->Fill(dPhi, weight);

		isoSum += pT;

                if (fabs(dEta) > 0.01)
                        isoSumJura01 += pT;
		if (fabs(dEta) > 0.02)
			isoSumJura02 += pT;
                if (fabs(dEta) > 0.03)
                        isoSumJura03 += pT;

		h1_tkIso03AllRedEta_[det]->Fill(dEta, weight);
		h1_tkIso03AllRedPhi_[det]->Fill(dPhi, weight);
		h2_tkIso03AllRedR2D_[det]->Fill(dEta, dPhi, weight);

	}

	h1_tkIso03AllRe_[det]->Fill(isoSum, weight);
	h1_tkIso03AllReRel_[det]->Fill(isoSum/cms2.els_p4()[eleIndex].Pt(), weight);
        h1_tkIso03AllReJura01_[det]->Fill(isoSumJura01, weight);
        h1_tkIso03AllReJura02_[det]->Fill(isoSumJura02, weight);
        h1_tkIso03AllReJura03_[det]->Fill(isoSumJura03, weight);
	h1_tkIso03AllReJura01In015_[det]->Fill(isoSumJura01In015, weight);


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

	if (!foundFirst) return;

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
        float tkThresholds[2] = {2.2, 1.1};
        float ecalThresholds[2] = {4.2, 3.4};
        float hcalThresholds[2] = {2.0, 1.3};

	if (tkIso > tkThresholds[det]) return;
	if (ecalIso > ecalThresholds[det]) return;
	if (hcalIso > hcalThresholds[det]) return;

	if (det != 0 && det != 1) std::cout << "Error! Not in barrel or endcap - something is wrong" << std::endl;
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

	// must have found first electron
	if (!foundFirst) return;

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

	//
	// get the subdestector (EE or EB)
	int det = getSubdet(eleIndex);

	//
	// Pt cut
        float ptThreshold = 0.0;
        if ((cuts_ & (CUT_BIT(ELE_PT_20))) == (CUT_BIT(ELE_PT_20))) ptThreshold = 20.0;
        if ((cuts_ & (CUT_BIT(ELE_PT_30))) == (CUT_BIT(ELE_PT_30))) ptThreshold = 30.0;

	h1_weff_pt_[det]->Fill(cms2.els_p4()[eleIndex].Pt(), weight);
	if (cms2.els_p4()[eleIndex].Pt() < ptThreshold) return;

	//
	// construct variables that are not already in ntuple
	//

	// isolation
	const float &ecalIso = cms2.els_ecalIso()[eleIndex];
	const float &hcalIso = cms2.els_hcalIso()[eleIndex];
	const float &tkIso = cms2.els_tkIso()[eleIndex];
	float isoSum = ecalIso + hcalIso + tkIso;

	// leading pT JPT jet
	// that is dR > 0.4 from the nearest electron
	float leadingJPT = 0.0;
	float mostBackToBackJPT = 0.0;
	float leastEMJPT = 999.9;
	int leastEMJPTIndex = 0;
	int leadingJPTIndex = 0;
	int mostBackToBackJPTIndex = 0;
	for (size_t j = 0; j < cms2.jpts_p4().size(); ++j)
	{
		if ( TMath::Abs(cms2.jpts_p4()[j].eta()) > 2.5 ) continue;
		if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.els_p4()[eleIndex], cms2.jpts_p4()[j])) < 0.4) continue;

		// find leading pT JPT
		if (cms2.jpts_p4()[j].Et() > leadingJPT) {
			leadingJPT = cms2.jpts_p4()[j].Et();
			leadingJPTIndex = j;
		}

		// find most back to back JPT from the electron
		float dPhi = acos(cos(cms2.els_p4()[eleIndex].Phi() - cms2.jpts_p4()[j].Phi()));
		if (dPhi > mostBackToBackJPT) {
			mostBackToBackJPT = dPhi;
			mostBackToBackJPTIndex = j;
		}

		// find the JPT with the lowest EM fraction
		if (cms2.jpts_emFrac()[j] < leastEMJPT) {
			leastEMJPT = cms2.jpts_emFrac()[j];
			leastEMJPTIndex = j;
		}


	}

	// distance in degrees of the leading JPT from the electron
	float leadingJPTAngle = 0.0;
	float mostBackToBackJPTAngle = 0.0;
	if (leadingJPT > 0.0) {
		// angle between the electron and the highest pT JPT
		leadingJPTAngle = (180.0 * acos(cos(cms2.els_p4()[eleIndex].Phi()
						- cms2.jpts_p4()[leadingJPTIndex].Phi())) / 3.14159265358979312);
		// angle between the electron and the JPT that is most back to back to it
		mostBackToBackJPTAngle = (180.0 * acos(cos(cms2.els_p4()[eleIndex].Phi()
						- cms2.jpts_p4()[mostBackToBackJPTIndex].Phi())) / 3.14159265358979312);	
	}

	// 
	// set up the selections to apply
	//

	float isolationThreshold = 999;
	bool applyIsoV0 = false;
	bool applyIsoV1 = false;
	bool applyIsoV2 = false;
	if ((cuts_ & (CUT_BIT(ELE_ISO_15))) == (CUT_BIT(ELE_ISO_15))) isolationThreshold = 0.15;
	if ((cuts_ & (CUT_BIT(ELE_ISO_10))) == (CUT_BIT(ELE_ISO_10))) isolationThreshold = 0.10;
	if ((cuts_ & (CUT_BIT(ELE_ISO_V0))) == (CUT_BIT(ELE_ISO_V0))) applyIsoV0 = true;
        if ((cuts_ & (CUT_BIT(ELE_ISO_V1))) == (CUT_BIT(ELE_ISO_V1))) applyIsoV1 = true;
        if ((cuts_ & (CUT_BIT(ELE_ISO_V2))) == (CUT_BIT(ELE_ISO_V2))) applyIsoV2 = true;
	float jptThreshold = 999;
	if ((cuts_ & (CUT_BIT(EVT_JPT_25))) == (CUT_BIT(EVT_JPT_25))) jptThreshold = 25.0;
	float tcMetThreshold = 0;
	if ((cuts_ & (CUT_BIT(EVT_TCMET_30))) == (CUT_BIT(EVT_TCMET_30))) tcMetThreshold = 30.0;
	if ((cuts_ & (CUT_BIT(EVT_TCMET_20))) == (CUT_BIT(EVT_TCMET_20))) tcMetThreshold = 20.0;
	float jptPhiThreshold = 180.0;
        if ((cuts_ & (CUT_BIT(EVT_JPT_PHIMAX_100))) == (CUT_BIT(EVT_JPT_PHIMAX_100))) jptPhiThreshold = 100.0;
        if ((cuts_ & (CUT_BIT(EVT_JPT_PHIMAX_110))) == (CUT_BIT(EVT_JPT_PHIMAX_110))) jptPhiThreshold = 110.0;
        if ((cuts_ & (CUT_BIT(EVT_JPT_PHIMAX_130))) == (CUT_BIT(EVT_JPT_PHIMAX_130))) jptPhiThreshold = 130.0;
	bool rejectConversions = false;
	if ((cuts_ & (CUT_BIT(ELE_NOCONV))) == (CUT_BIT(ELE_NOCONV))) rejectConversions = true;

	bool use_tasElectron_v1 = false;
	if ((cuts_ & (CUT_BIT(ELE_TAS_V1))) == (CUT_BIT(ELE_TAS_V1))) use_tasElectron_v1 = true;

	//
	// plots of key quantities before selection applied
	h1_weff_iso_[det]->Fill(isoSum / cms2.els_p4()[eleIndex].Pt(), weight);
	h1_weff_tcmet_[det]->Fill(cms2.evt_tcmet(), weight);
	h1_weff_jptpt_[det]->Fill(leadingJPT, weight);

	//
	// apply selections
	//

	// isolation cut
	float dummy = 0.0;
	if (isoSum / cms2.els_p4()[0].Pt() > isolationThreshold && isolationThreshold != 999) return;
	if (applyIsoV0) {
	        float tkThresholds[2] = {4.5, 6.0};
        	float ecalThresholds[2] = {2.5, 2.0};
	        float hcalThresholds[2] = {1.0, 1.0};
	        if (tkIso > tkThresholds[det]) return;
	        if (ecalIso > ecalThresholds[det]) return;
	        if (hcalIso > hcalThresholds[det]) return;
	}
        if (applyIsoV1) {
		float tkIsoJura01In015 = recomputeTrackIsolation(eleIndex, 0.01, 0.015, 0.3, dummy);
                float tkThresholds[2] = {2.5, 2.0};
                float ecalThresholds[2] = {2.5, 2.0};
                float hcalThresholds[2] = {1.0, 1.0};
                if (tkIsoJura01In015 > tkThresholds[det]) return;
                if (ecalIso > ecalThresholds[det]) return;
                if (hcalIso > hcalThresholds[det]) return;
        }
        if (applyIsoV2) {
                float tkIsoJura01In015 = recomputeTrackIsolation(eleIndex, 0.01, 0.015, 0.3, dummy);
                float tkThresholds[2] = {2.5, 2.0};
                float caloThresholds[2] = {3.0, 2.5};
                if (tkIsoJura01In015 > tkThresholds[det]) return;
                if ((ecalIso+hcalIso) > caloThresholds[det]) return;
        }    
	//

	// plots of tcmet and leading jpt pt
	h1_weff_tcmet_after_iso_[det]->Fill(cms2.evt_tcmet(), weight);
	h1_weff_jptpt_after_iso_[det]->Fill(leadingJPT, weight);
	// angle between leading JPT and electron
	h1_weff_leadjptphi_after_iso_[det]->Fill(leadingJPTAngle, weight); 
	// angle between electron and JPT that is most back to back to it
	h1_weff_jptphimax_after_iso_[det]->Fill(mostBackToBackJPTAngle, weight);
	// d0 corrected
	h1_weff_d0corr_after_iso_[det]->Fill(fabs(cms2.els_d0corr()[eleIndex]), weight);

	// emFraction of JPT with lowest emFraction
	h1_weff_leastemjpt_after_iso_[det]->Fill(leastEMJPT, weight);


	// leading JPT cut
	if (leadingJPT > jptThreshold && jptThreshold != 999) return;
	//

        // JPT phimax cut
        if (mostBackToBackJPTAngle > jptPhiThreshold && jptPhiThreshold != 180.0) return;
        //

        // tas electron 
	cuts_t eleIdResult = ele::tasElectron_v1(eleIndex);
        if (use_tasElectron_v1 && !((eleIdResult & eleid_tasElectron_v1) == eleid_tasElectron_v1)) return;

	// plot of tcmet after isolation and jpt veto applied
	h1_weff_tcmet_after_iso_jpt_[det]->Fill(cms2.evt_tcmet(), weight);
	// plot of angle between leading JPT and electron
	h1_weff_leadjptphi_after_iso_jpt_[det]->Fill(leadingJPTAngle, weight);
        // d0 corrected
        h1_weff_d0corr_after_iso_jpt_[det]->Fill(fabs(cms2.els_d0corr()[eleIndex]), weight);

	// conversion rejection cut
	if (rejectConversions && isconversionElectron09(eleIndex)) return;
	//

	h1_weff_tcmet_after_iso_jpt_conv_[det]->Fill(cms2.evt_tcmet(), weight);

	// BACKGROUND
	// plots of electron ID quantities
	//
        if (cms2.evt_tcmet() < 15.0) {
		h1_weffbg_sigmaIEtaIEta_[det]->Fill(cms2.els_sigmaIEtaIEta()[eleIndex], weight); 
	}

	// tcMet cut
	if (cms2.evt_tcmet() > tcMetThreshold) {
	//

		// plot of angle between leading JPT and electron
		h1_weff_leadjptphi_after_iso_jpt_tcmet_[det]->Fill(leadingJPTAngle, weight);

		// SIGNAL
		// plots of electron ID quantities
		//
		h1_weffs_sigmaIEtaIEta_[det]->Fill(cms2.els_sigmaIEtaIEta()[eleIndex], weight);	

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
	//AN2009_98();

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
				h1_caloIso03All_[det]->Fill(ecalIso + hcalIso, weight);

				h2_tkIso03All_[det]->Fill(cms2.els_p4()[i].Pt(), tkIso, weight);
                                h2_ecalIso03All_[det]->Fill(cms2.els_p4()[i].Pt(), ecalIso, weight);
                                h2_hcalIso03All_[det]->Fill(cms2.els_p4()[i].Pt(), hcalIso, weight);
				h2_caloIso03All_[det]->Fill(cms2.els_p4()[i].Pt(), ecalIso + hcalIso, weight);

				// N-1
				if (ecalIso < ecalThresholdsNM1[det] && hcalIso < hcalThresholdsNM1[det]) {
					float shCutSum = 0.0;
                                        float tkIsoJura01In015 = recomputeTrackIsolation(i, 0.01, 0.015, 0.3, shCutSum);
					h1_tkIso03AllNM1_[det]->Fill(tkIso, weight);
					h1_tkIso03AllReJura01In015NM1_[det]->Fill(tkIsoJura01In015, weight);
					h1_tkIso03AllReShCutNM1_[det]->Fill(shCutSum, weight);				
	
					if (cms2.els_egamma_tightId()[i]) {
						h1_tkIso03AllIDNM1_[det]->Fill(tkIso, weight);
						h1_tkIso03AllReJura01In015IDNM1_[det]->Fill(tkIsoJura01In015, weight);	
					}
                                        if (!isconversionElectron09(i)) {
                                                h1_tkIso03AllConvNM1_[det]->Fill(tkIso, weight);
                                                h1_tkIso03AllReJura01In015ConvNM1_[det]->Fill(tkIsoJura01In015, weight);
					}                  
                                        if (!isconversionElectron09(i) && cms2.els_egamma_tightId()[i]) {
                                                h1_tkIso03AllConvIDNM1_[det]->Fill(tkIso, weight);
                                                h1_tkIso03AllReJura01In015ConvIDNM1_[det]->Fill(tkIsoJura01In015, weight);
					}
				}
                                if (tkIso < tkThresholdsNM1[det] && hcalIso < hcalThresholdsNM1[det])
					h1_ecalIso03AllNM1_[det]->Fill(ecalIso, weight);
                                if (tkIso < tkThresholdsNM1[det] && ecalIso < ecalThresholdsNM1[det])
					h1_hcalIso03AllNM1_[det]->Fill(hcalIso, weight);

				//trackIsolationStudy(i, det);

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
						cms2.els_p4()[i].Pt(), cms2.els_etaSC()[i], cms2.els_phiSC()[i], 1);

				em_dPhiIn_[det]->Fill(fabs(cms2.els_dPhiIn()[i]),
						cms2.els_p4()[i].Pt(), cms2.els_etaSC()[i], cms2.els_phiSC()[i], 1);

				em_hoe_[det]->Fill(fabs(cms2.els_hOverE()[i]),
						cms2.els_p4()[i].Pt(), cms2.els_etaSC()[i], cms2.els_phiSC()[i], 1);

				em_sieie_[det]->Fill(fabs(cms2.els_sigmaIEtaIEta()[i]),
						cms2.els_p4()[i].Pt(), cms2.els_etaSC()[i], cms2.els_phiSC()[i], 1);

				em_classBasedTight_[det]->Fill(cms2.els_egamma_tightId()[i],
						cms2.els_p4()[i].Pt(), cms2.els_etaSC()[i], cms2.els_phiSC()[i], 1);

				em_robustTight_[det]->Fill(cms2.els_egamma_robustTightId()[i],
						cms2.els_p4()[i].Pt(), cms2.els_etaSC()[i], cms2.els_phiSC()[i], 1);

				em_eopInGT05_[det]->Fill(cms2.els_eOverPIn()[i],
						cms2.els_p4()[i].Pt(), cms2.els_etaSC()[i], cms2.els_phiSC()[i], 1);

				em_eopInLT30_[det]->Fill(cms2.els_eOverPIn()[i],
						cms2.els_p4()[i].Pt(), cms2.els_etaSC()[i], cms2.els_phiSC()[i], 1);

				cuts_t eleIdResult = ele::tasElectron_v1(i);
				
				em_tasElectronV1_[det]->Fill(CheckCuts(eleid_tasElectron_v1, eleIdResult),
				                cms2.els_p4()[i].Pt(), cms2.els_etaSC()[i], cms2.els_phiSC()[i], 1);

				// N-1 for dEtaIn
				if (CheckCutsNM1(eleid_tasElectron_v1, CUT_BIT(ELEID_TAS_DETAIN), eleIdResult)) {
					em_dEtaInTasV1NM1_[det]->Fill(fabs(cms2.els_dEtaIn()[i]),
                                                cms2.els_p4()[i].Pt(), cms2.els_etaSC()[i], cms2.els_phiSC()[i], 1);			
					h1_dEtaInTasV1NM1_[det]->Fill(fabs(cms2.els_dEtaIn()[i]), weight);
				}

                                // N-1 for dPhiIn
                                if (CheckCutsNM1(eleid_tasElectron_v1, CUT_BIT(ELEID_TAS_DPHIIN), eleIdResult)) {
                                        em_dPhiInTasV1NM1_[det]->Fill(fabs(cms2.els_dPhiIn()[i]),
                                                cms2.els_p4()[i].Pt(), cms2.els_etaSC()[i], cms2.els_phiSC()[i], 1);
                                        h1_dPhiInTasV1NM1_[det]->Fill(fabs(cms2.els_dPhiIn()[i]), weight);
                                }

                                // N-1 for hoe
                                if (CheckCutsNM1(eleid_tasElectron_v1, CUT_BIT(ELEID_TAS_HOE), eleIdResult)) {
                                        em_hoeTasV1NM1_[det]->Fill(cms2.els_hOverE()[i],
                                                cms2.els_p4()[i].Pt(), cms2.els_etaSC()[i], cms2.els_phiSC()[i], 1);
                                        h1_hoeTasV1NM1_[det]->Fill(cms2.els_hOverE()[i], weight);
                                }

                                // N-1 for sigmaIEtaIEta
                                if (CheckCutsNM1(eleid_tasElectron_v1, CUT_BIT(ELEID_TAS_LSHAPE), eleIdResult)) {
                                        em_sigmaIEtaIEtaTasV1NM1_[det]->Fill(cms2.els_sigmaIEtaIEta()[i],
                                                cms2.els_p4()[i].Pt(), cms2.els_etaSC()[i], cms2.els_phiSC()[i], 1);
                                        h1_sigmaIEtaIEtaTasV1NM1_[det]->Fill(cms2.els_sigmaIEtaIEta()[i], weight);
                                }

                                // N-1 for E2x5Max/E5x5
                                if (CheckCutsNM1(eleid_tasElectron_v1, CUT_BIT(ELEID_TAS_LSHAPE), eleIdResult)) {
                                        em_E2x5Norm5x5TasV1NM1_[det]->Fill(cms2.els_e2x5Max()[i]/cms2.els_e5x5()[i],
                                                cms2.els_p4()[i].Pt(), cms2.els_etaSC()[i], cms2.els_phiSC()[i], 1);
                                        h1_E2x5Norm5x5TasV1NM1_[det]->Fill(cms2.els_e2x5Max()[i]/cms2.els_e5x5()[i], weight);
                                }


			}

		} // end loop on electrons

	} // end event level cuts passed

}

bool Looper::CheckCutsNM1(cuts_t apply, cuts_t remove, cuts_t passed)
{
	if ((passed & (apply & (~remove))) == (apply & (~remove))) return true;
	return false;
}

bool Looper::CheckCuts(cuts_t apply, cuts_t passed)
{
	if ((apply & passed) == apply) return true;
	return false;
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
