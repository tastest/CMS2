
//
// ttbar -> ll
// Dave "the one but not the only" Evans 
//

#include "MyScanChain.h"

// ROOT includes
#include "TChain.h"
#include "TChainElement.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"

#include "Math/LorentzVector.h"

// CMS2 includes
#include "CMS2.h"
#include "../../CORE/electronSelections.h"
#include "../../CORE/selections.h"
#include "../../Tools/DileptonHypType.h"

//
// Namespaces
//
using namespace tas;

//
//
//

enum ele_selection {
	PASS_DPHI,
	PASS_DETA,
	PASS_HOE,
	PASS_LSHAPE,
	PASS_D0,
	PASS_ISO,
	PASS_DETA_CAND02,
	PASS_LSHAPE_CAND02,
	PASS_EXTRA,
	PASS_EXTRA_V2,
	PASS_CONV,
	PASS_NOMUON,
};

//
// for hyps
//

double dRbetweenVectors(const LorentzVector &vec1,
		const LorentzVector &vec2 ){

	double dphi = std::min(::fabs(vec1.Phi() - vec2.Phi()), 2 * M_PI - fabs(vec1.Phi() - vec2.Phi()));
	double deta = vec1.Eta() - vec2.Eta();
	return sqrt(dphi*dphi + deta*deta);
} 

enum DileptonHypType hyp_typeToHypType (int hyp_type)
{
	switch (hyp_type) {
		case 0:
			return DILEPTON_MUMU;
		case 1: case 2:
			return DILEPTON_EMU;
		case 3:
			return DILEPTON_EE;
		default:
			assert(hyp_type < 4);
	}
	return DILEPTON_ALL;
}

void MyScanChain::Fill(TH1F** hist, const unsigned int hyp, const float &val, const float &weight)
{
	hist[hyp]->Fill(val, weight);
	hist[DILEPTON_ALL]->Fill(val, weight);
}

void MyScanChain::Fill2D(TH2F** hist, const unsigned int hyp, const float &valx, const float &valy, const float &weight)
{   
	hist[hyp]->Fill(valx, valy, weight); 
	hist[DILEPTON_ALL]->Fill(valx, valy, weight);
}


void MyScanChain::FormatHist(TH1F** hist, std::string sampleName, std::string name, int n, float min, float max)
{       
	// loop on EB, EE
	for (unsigned int i = 0; i < 4; ++i)
	{
		std::string str = dilepton_hypo_names[i];
		std::string title = name + "_" + str;
		hist[i] = new TH1F(Form("%s_%s_%s", sampleName.c_str(), name.c_str(), str.c_str()),
				title.c_str(), n, min, max);
		hist[i]->GetXaxis()->SetTitle(name.c_str());
		hist[i]->Sumw2();
	}
}    

void MyScanChain::FormatHist2D(TH2F** hist, std::string sampleName, std::string name, int nx, float minx, float maxx, int ny, float miny, float maxy)
{       
	// loop on EB, EE 
	for (unsigned int i = 0; i < 4; ++i)
	{   
		std::string str = dilepton_hypo_names[i];
		std::string title = name + "_" + str;
		hist[i] = new TH2F(Form("%s_%s_%s", sampleName.c_str(), name.c_str(), str.c_str()),
				title.c_str(), nx, minx, maxx, ny, miny, maxy);
		hist[i]->GetXaxis()->SetTitle(name.c_str());
		hist[i]->Sumw2();
	}
}

void MyScanChain::FormatAllEleIdHistograms(std::string sampleName)
{

	FormatHist(h1_hyp_lt_eb_pt_, sampleName, "hyp_lt_eb_pt", 20, 0.0, 100.0);
	FormatHist(h1_hyp_lt_ee_pt_, sampleName, "hyp_lt_ee_pt", 20, 0.0, 100.0);

	FormatHist(h1_hyp_lt_eb_hoe_, sampleName, "hyp_lt_eb_hoe", 100, 0.0, 0.1);
	FormatHist(h1_hyp_lt_ee_hoe_, sampleName, "hyp_lt_ee_hoe", 100, 0.0, 0.1);

	FormatHist(h1_hyp_lt_eb_sigmaIEtaIEta_, sampleName, "hyp_lt_eb_sigmaIEtaIEta", 50, 0.0, 0.05);
	FormatHist(h1_hyp_lt_ee_sigmaIEtaIEta_, sampleName, "hyp_lt_ee_sigmaIEtaIEta", 50, 0.0, 0.05);

	FormatHist(h1_hyp_lt_eb_dEtaIn_, sampleName, "hyp_lt_eb_dEtaIn", 100, 0.0, 0.05);
	FormatHist(h1_hyp_lt_ee_dEtaIn_, sampleName, "hyp_lt_ee_dEtaIn", 100, 0.0, 0.05);

	FormatHist(h1_hyp_lt_eb_dPhiIn_, sampleName, "hyp_lt_eb_dPhiIn", 100, 0.0, 0.1);
	FormatHist(h1_hyp_lt_ee_dPhiIn_, sampleName, "hyp_lt_ee_dPhiIn", 100, 0.0, 0.1);

	FormatHist(h1_hyp_lt_eb_d0_, sampleName, "hyp_lt_eb_d0", 100, 0.0, 0.05);
	FormatHist(h1_hyp_lt_ee_d0_, sampleName, "hyp_lt_ee_d0", 100, 0.0, 0.05);

	FormatHist(h1_hyp_lt_eb_E2x5MaxOver5x5_, sampleName, "hyp_lt_eb_E2x5MaxOver5x5", 50, 0.6, 1.1);
	FormatHist(h1_hyp_lt_ee_E2x5MaxOver5x5_, sampleName, "hyp_lt_ee_E2x5MaxOver5x5", 50, 0.6, 1.1);

	// NM1
	FormatHist(h1_hyp_lt_eb_nm1_hoe_, sampleName, "hyp_lt_eb_nm1_hoe", 100, 0.0, 0.1);
	FormatHist(h1_hyp_lt_ee_nm1_hoe_, sampleName, "hyp_lt_ee_nm1_hoe", 100, 0.0, 0.1);

	FormatHist(h1_hyp_lt_eb_nm1_sigmaIEtaIEta_, sampleName, "hyp_lt_eb_nm1_sigmaIEtaIEta", 50, 0.0, 0.05);
	FormatHist(h1_hyp_lt_ee_nm1_sigmaIEtaIEta_, sampleName, "hyp_lt_ee_nm1_sigmaIEtaIEta", 50, 0.0, 0.05);

	FormatHist(h1_hyp_lt_eb_nm1_dEtaIn_, sampleName, "hyp_lt_eb_nm1_dEtaIn", 100, 0.0, 0.05);
	FormatHist(h1_hyp_lt_ee_nm1_dEtaIn_, sampleName, "hyp_lt_ee_nm1_dEtaIn", 100, 0.0, 0.05);

	FormatHist(h1_hyp_lt_eb_nm1_dPhiIn_, sampleName, "hyp_lt_eb_nm1_dPhiIn", 100, 0.0, 0.1);
	FormatHist(h1_hyp_lt_ee_nm1_dPhiIn_, sampleName, "hyp_lt_ee_nm1_dPhiIn", 100, 0.0, 0.1);

	FormatHist(h1_hyp_lt_eb_nm1_d0_, sampleName, "hyp_lt_eb_nm1_d0", 100, 0.0, 0.05);
	FormatHist(h1_hyp_lt_ee_nm1_d0_, sampleName, "hyp_lt_ee_nm1_d0", 100, 0.0, 0.05);

	FormatHist(h1_hyp_lt_eb_nm1_E2x5MaxOver5x5_, sampleName, "hyp_lt_eb_nm1_E2x5MaxOver5x5", 50, 0.6, 1.1);
	FormatHist(h1_hyp_lt_ee_nm1_E2x5MaxOver5x5_, sampleName, "hyp_lt_ee_nm1_E2x5MaxOver5x5", 50, 0.6, 1.1);

	FormatHist2D(h1_hyp_lt_eb_nm1_lateral_, sampleName, "hyp_lt_eb_nm1_lateral", 110, 0, 1.1, 110, 0, 1.1);
	FormatHist2D(h1_hyp_lt_ee_nm1_lateral_, sampleName, "hyp_lt_ee_nm1_lateral", 110, 0, 1.1, 110, 0, 1.1);

	// old and new pt spectra for second pass at id
	//
	FormatHist(h1_hyp_lt_eb_pt_cand01_, sampleName, "hyp_lt_eb_pt_cand01", 20, 0.0, 100.0);
	FormatHist(h1_hyp_lt_ee_pt_cand01_, sampleName, "hyp_lt_ee_pt_cand01", 20, 0.0, 100.0);

	FormatHist(h1_hyp_lt_eb_pt_cand02_, sampleName, "hyp_lt_eb_pt_cand02", 20, 0.0, 100.0);
	FormatHist(h1_hyp_lt_ee_pt_cand02_, sampleName, "hyp_lt_ee_pt_cand02", 20, 0.0, 100.0);

	FormatHist(h1_hyp_lt_eb_pt_cand02_extra_, sampleName, "hyp_lt_eb_pt_cand02_extra", 20, 0.0, 100.0);
	FormatHist(h1_hyp_lt_ee_pt_cand02_extra_, sampleName, "hyp_lt_ee_pt_cand02_extra", 20, 0.0, 100.0);

	FormatHist(h1_hyp_lt_eb_pt_cand02_extra_v2_, sampleName, "hyp_lt_eb_pt_cand02_extra_v2", 20, 0.0, 100.0);
	FormatHist(h1_hyp_lt_ee_pt_cand02_extra_v2_, sampleName, "hyp_lt_ee_pt_cand02_extra_v2", 20, 0.0, 100.0);

	FormatHist(h1_hyp_lt_eb_eta_cand02_extra_v2_, sampleName, "hyp_lt_eb_eta_cand02_extra_v2", 120, -3, 3);
	FormatHist(h1_hyp_lt_ee_eta_cand02_extra_v2_, sampleName, "hyp_lt_ee_eta_cand02_extra_v2", 120, -3, 3);


	//
	//

	FormatHist(h1_hyp_lt_eb_ecalIso_, sampleName, "hyp_lt_eb_ecalIso", 100, 0, 10);
	FormatHist(h1_hyp_lt_eb_hcalIso_, sampleName, "hyp_lt_eb_hcalIso", 100, 0, 25);
	FormatHist(h1_hyp_lt_eb_tkIso_, sampleName, "hyp_lt_eb_tkIso", 100, 0, 25);

	FormatHist(h1_hyp_lt_ee_ecalIso_, sampleName, "hyp_lt_ee_ecalIso", 100, 0, 10);
	FormatHist(h1_hyp_lt_ee_hcalIso_, sampleName, "hyp_lt_ee_hcalIso", 100, 0, 25);
	FormatHist(h1_hyp_lt_ee_tkIso_, sampleName, "hyp_lt_ee_tkIso", 100, 0, 25);

	FormatHist(h1_hyp_lt_ee_relsusy_, sampleName, "hyp_lt_ee_relsusy", 100, 0, 1);
	FormatHist(h1_hyp_lt_eb_relsusy_, sampleName, "hyp_lt_eb_relsusy", 100, 0, 1);

	FormatHist(h1_hyp_lt_ee_afterid_relsusy_, sampleName, "hyp_lt_ee_afterid_relsusy", 100, 0, 1);
	FormatHist(h1_hyp_lt_eb_afterid_relsusy_, sampleName, "hyp_lt_eb_afterid_relsusy", 100, 0, 1);

	FormatHist(h1_hyp_lt_ee_afterid_fbrem_, sampleName, "hyp_lt_ee_afterid_fbrem", 100, 0, 1);
	FormatHist(h1_hyp_lt_eb_afterid_fbrem_, sampleName, "hyp_lt_eb_afterid_fbrem", 100, 0, 1);

	FormatHist(h1_hyp_lt_ee_afterid_eopin_, sampleName, "hyp_lt_ee_afterid_eopin", 100, 0, 5);
	FormatHist(h1_hyp_lt_eb_afterid_eopin_, sampleName, "hyp_lt_eb_afterid_eopin", 100, 0, 5);

	FormatHist(h1_hyp_lt_ee_afterid_relsusy_lowfbrem_, sampleName, "hyp_lt_ee_afterid_relsusy_lowfbrem", 100, 0, 1);
	FormatHist(h1_hyp_lt_eb_afterid_relsusy_lowfbrem_, sampleName, "hyp_lt_eb_afterid_relsusy_lowfbrem", 100, 0, 1);
	FormatHist(h1_hyp_lt_ee_afterid_relsusy_highfbrem_, sampleName, "hyp_lt_ee_afterid_relsusy_highfbrem", 100, 0, 1);
	FormatHist(h1_hyp_lt_eb_afterid_relsusy_highfbrem_, sampleName, "hyp_lt_eb_afterid_relsusy_highfbrem", 100, 0, 1);

	FormatHist(h1_hyp_lt_ee_afterid_eopin_lowfbrem_, sampleName, "hyp_lt_ee_afterid_eopin_lowfbrem", 100, 0, 5);
	FormatHist(h1_hyp_lt_eb_afterid_eopin_lowfbrem_, sampleName, "hyp_lt_eb_afterid_eopin_lowfbrem", 100, 0, 5);
	FormatHist(h1_hyp_lt_ee_afterid_eopin_highfbrem_, sampleName, "hyp_lt_ee_afterid_eopin_highfbrem", 100, 0, 5);
	FormatHist(h1_hyp_lt_eb_afterid_eopin_highfbrem_, sampleName, "hyp_lt_eb_afterid_eopin_highfbrem", 100, 0, 5);

	// preshower after id
	FormatHist(h1_hyp_lt_ee_afterid_preshowerEnergy_lowfbrem_, sampleName, "hyp_lt_ee_afterid_preshowerEnergy_lowfbrem", 100, 0, 2);
	FormatHist(h1_hyp_lt_eb_afterid_preshowerEnergy_lowfbrem_, sampleName, "hyp_lt_eb_afterid_preshowerEnergy_lowfbrem", 100, 0, 2);
	FormatHist(h1_hyp_lt_ee_afterid_preshowerEnergy_highfbrem_, sampleName, "hyp_lt_ee_afterid_preshowerEnergy_highfbrem", 100, 0, 2);
	FormatHist(h1_hyp_lt_eb_afterid_preshowerEnergy_highfbrem_, sampleName, "hyp_lt_eb_afterid_preshowerEnergy_highfbrem", 100, 0, 2);


	// dPhiIn after id (except dPhiIn)
	FormatHist(h1_hyp_lt_ee_afterid_dPhiIn_lowfbrem_, sampleName, "hyp_lt_ee_afterid_dPhiIn_lowfbrem", 100, 0, 0.1);
	FormatHist(h1_hyp_lt_eb_afterid_dPhiIn_lowfbrem_, sampleName, "hyp_lt_eb_afterid_dPhiIn_lowfbrem", 100, 0, 0.1);
	FormatHist(h1_hyp_lt_ee_afterid_dPhiIn_highfbrem_, sampleName, "hyp_lt_ee_afterid_dPhiIn_highfbrem", 100, 0, 0.1);
	FormatHist(h1_hyp_lt_eb_afterid_dPhiIn_highfbrem_, sampleName, "hyp_lt_eb_afterid_dPhiIn_highfbrem", 100, 0, 0.1);

	// dEtaIn after id (except dEtaIn)
	FormatHist(h1_hyp_lt_ee_afterid_dEtaIn_lowfbrem_, sampleName, "hyp_lt_ee_afterid_dEtaIn_lowfbrem", 40, 0, 0.04);
	FormatHist(h1_hyp_lt_eb_afterid_dEtaIn_lowfbrem_, sampleName, "hyp_lt_eb_afterid_dEtaIn_lowfbrem", 40, 0, 0.04);
	FormatHist(h1_hyp_lt_ee_afterid_dEtaIn_highfbrem_, sampleName, "hyp_lt_ee_afterid_dEtaIn_highfbrem", 40, 0, 0.04);
	FormatHist(h1_hyp_lt_eb_afterid_dEtaIn_highfbrem_, sampleName, "hyp_lt_eb_afterid_dEtaIn_highfbrem", 40, 0, 0.04);

	// closest mu
	FormatHist(h1_hyp_lt_ee_afterid_closestMu_, sampleName, "hyp_lt_ee_closestMu", 60, 0.0, 3); 
	FormatHist(h1_hyp_lt_eb_afterid_closestMu_, sampleName, "hyp_lt_eb_closestMu", 60, 0.0, 3); 


	FormatHist(h1_hyp_lt_eb_pt_idnew_, sampleName, "hyp_lt_eb_pt_idnew", 20, 0.0, 100.0);
	FormatHist(h1_hyp_lt_ee_pt_idnew_, sampleName, "hyp_lt_ee_pt_idnew", 20, 0.0, 100.0);
	FormatHist(h1_hyp_lt_eb_pt_idold_, sampleName, "hyp_lt_eb_pt_idold", 20, 0.0, 100.0);
	FormatHist(h1_hyp_lt_ee_pt_idold_, sampleName, "hyp_lt_ee_pt_idold", 20, 0.0, 100.0);

	FormatHist(h1_hyp_lt_eb_pt_isonew_cand1_, sampleName, "hyp_lt_eb_pt_isonew_cand1", 20, 0.0, 100.0);
	FormatHist(h1_hyp_lt_ee_pt_isonew_cand1_, sampleName, "hyp_lt_ee_pt_isonew_cand1", 20, 0.0, 100.0);
	FormatHist(h1_hyp_lt_eb_pt_isonew_, sampleName, "hyp_lt_eb_pt_isonew", 20, 0.0, 100.0);
	FormatHist(h1_hyp_lt_ee_pt_isonew_, sampleName, "hyp_lt_ee_pt_isonew", 20, 0.0, 100.0);
	FormatHist(h1_hyp_lt_eb_pt_isoold_, sampleName, "hyp_lt_eb_pt_isoold", 20, 0.0, 100.0);
	FormatHist(h1_hyp_lt_ee_pt_isoold_, sampleName, "hyp_lt_ee_pt_isoold", 20, 0.0, 100.0);
	FormatHist(h1_hyp_lt_eb_pt_conv_, sampleName, "hyp_lt_eb_pt_conv", 20, 0.0, 100.0);
	FormatHist(h1_hyp_lt_ee_pt_conv_, sampleName, "hyp_lt_ee_pt_conv", 20, 0.0, 100.0);

	FormatHist(h1_hyp_lt_eb_pt_id1_iso1_conv_, sampleName, "hyp_lt_eb_pt_id1_iso1_conv", 20, 0.0, 100.0);
	FormatHist(h1_hyp_lt_ee_pt_id1_iso1_conv_, sampleName, "hyp_lt_ee_pt_id1_iso1_conv", 20, 0.0, 100.0);
}

void MyScanChain::FillAllEleIdHistogramsNoHyp(const float &weight, const TString &sampleName)
{

	for (size_t i = 0; i < cms2.els_p4().size(); ++i) {
		if (!electron20Eta2p4(i)) continue;
		if (!(cms2.els_type()[i] & (1<<ISECALDRIVEN))) continue;
		FillAllEleIdHistograms(i, weight, sampleName);
	}
}

void MyScanChain::FillAllEleIdHistogramsHyp(const unsigned int h, const float &weight, const TString &sampleName)
{

	DileptonHypType hypType = hyp_typeToHypType(cms2.hyp_type()[h]);

	// apply part of electron denominator first here
	if (hypType == DILEPTON_EMU && sampleName == "ttbar") {
		if(abs(cms2.hyp_ll_id()[h]) == 13) {
			if (trueMuonFromW_WJets(cms2.hyp_ll_index()[h]) && (cms2.els_type()[cms2.hyp_lt_index()[h]] & (1<<ISECALDRIVEN)) && fabs(cms2.els_p4().at(cms2.hyp_lt_index()[h]).eta()) < 2.4)
				FillAllEleIdHistograms(cms2.hyp_lt_index()[h], weight, sampleName);
		}
		if(abs(cms2.hyp_lt_id()[h]) == 13) {
			if (trueMuonFromW_WJets(cms2.hyp_lt_index()[h]) && (cms2.els_type()[cms2.hyp_ll_index()[h]] & (1<<ISECALDRIVEN)) && fabs(cms2.els_p4().at(cms2.hyp_ll_index()[h]).eta()) < 2.4)
				FillAllEleIdHistograms(cms2.hyp_ll_index()[h], weight, sampleName);
		}
	}
	if (hypType == DILEPTON_EMU && sampleName == "wm") {
		if(abs(cms2.hyp_ll_id()[h]) == 13) {
			if (trueMuonFromW_WJets(cms2.hyp_ll_index()[h]) && (cms2.els_type()[cms2.hyp_lt_index()[h]] & (1<<ISECALDRIVEN)) && fabs(cms2.els_p4().at(cms2.hyp_lt_index()[h]).eta()) < 2.4)
				//electron20Eta2p4(cms2.hyp_lt_index()[h])) 
				FillAllEleIdHistograms(cms2.hyp_lt_index()[h], weight, sampleName);
		}
		if(abs(cms2.hyp_lt_id()[h]) == 13) {
			if (trueMuonFromW_WJets(cms2.hyp_lt_index()[h]) && (cms2.els_type()[cms2.hyp_ll_index()[h]] & (1<<ISECALDRIVEN)) && fabs(cms2.els_p4().at(cms2.hyp_ll_index()[h]).eta()) < 2.4)
				//electron20Eta2p4(cms2.hyp_ll_index()[h]))
				FillAllEleIdHistograms(cms2.hyp_ll_index()[h], weight, sampleName);
		}

	}

}

void MyScanChain::FillAllEleIdHistograms(const unsigned int index, const float &weight, const TString &sampleName)
{

	// apply truth match behavior if ttbar
	if (sampleName == "ttbar") {
		if(!((abs(cms2.els_mc_id()[index]) == 11) && abs(cms2.els_mc_motherid()[index]) == 24) ) return;
	}
	// apply truth match behavior if wjets
	if (sampleName == "wjets") {
		if(((abs(cms2.els_mc_id()[index]) == 11) && abs(cms2.els_mc_motherid()[index]) == 24) ) return;
	}
	// apply truth match behavior if dyee
	if (sampleName == "dyee") {
		if(!((abs(cms2.els_mc_id()[index]) == 11) && abs(cms2.els_mc_motherid()[index]) == 23) ) return;
	}  
	// if particle gun
	if (sampleName == "elegun") {
		if(!((abs(cms2.els_mc_id()[index]) == 11))) return;
	} 

	//
	// define thresholds for EB, EE
	//
	float dEtaInThresholds[2]               = {0.007, 0.010};
	float dPhiInThresholds[2]               = {0.020, 0.025};
	float hoeThresholds[2]                  = {0.01, 0.01};
	float sigmaIEtaIEtaThresholds[2]        = {9999.99, 0.03};
	float e2x5Over5x5Thresholds[2]          = {0.90, 0.00};
	float d0Thresholds[2]               = {0.02, 0.02};

	float e2x5Over5x5Thresholds_cand02[2]          = {0.94, 0.00};
	float dEtaInThresholds_cand02[2]               = {0.005, 0.007};

	int ele_result = 0;
	int ele_passall = (1<<PASS_DETA) | (1<<PASS_DPHI) | (1<<PASS_HOE) | (1<<PASS_LSHAPE) | (1<<PASS_D0) | (1<<PASS_NOMUON) | (1<<PASS_CONV);
	int ele_passall_id_cand02 = (1<<PASS_DETA_CAND02) | (1<<PASS_DPHI) | (1<<PASS_HOE) | (1<<PASS_LSHAPE_CAND02) | (1<<PASS_D0) | (1<<PASS_NOMUON) | (1<<PASS_CONV);

	int ele_passall_id_and_iso_cand01 =  (1<<PASS_DETA) | (1<<PASS_DPHI) | (1<<PASS_HOE) | (1<<PASS_LSHAPE) | (1<<PASS_D0) | (1<<PASS_NOMUON) | (1<<PASS_CONV) | (1<<PASS_ISO);
	int ele_passall_id_and_iso_cand02 =  (1<<PASS_DETA_CAND02) | (1<<PASS_DPHI) | (1<<PASS_HOE) | (1<<PASS_LSHAPE_CAND02) | (1<<PASS_D0) | (1<<PASS_NOMUON) | (1<<PASS_CONV) | (1<<PASS_ISO);
	int ele_passall_id_and_iso_cand02_extra = ele_passall_id_and_iso_cand02 | (1<<PASS_EXTRA);
	int ele_passall_id_and_iso_cand02_extra_v2 = ele_passall_id_and_iso_cand02 | (1<<PASS_EXTRA_V2);


	//
	// apply cuts
	//
	float iso_relsusy = -1;
	iso_relsusy = electronIsolation_relsusy_cand1(index, true);

	if (electronId_noMuon(index)) ele_result |= (1<<PASS_NOMUON);
	if (!isFromConversionPartnerTrack(index)) ele_result |= (1<<PASS_CONV);

	if (fabs(cms2.els_etaSC()[index]) < 1.479) {
		if (fabs(cms2.els_dEtaIn()[index]) < dEtaInThresholds[0]) 	ele_result |= (1<<PASS_DETA);
		if (fabs(cms2.els_dEtaIn()[index]) < dEtaInThresholds_cand02[0])   ele_result |= (1<<PASS_DETA_CAND02);
		if (fabs(cms2.els_dPhiIn()[index]) < dPhiInThresholds[0])	ele_result |= (1<<PASS_DPHI);
		if (cms2.els_hOverE()[index] < hoeThresholds[0]) 			ele_result |= (1<<PASS_HOE);
		if ((cms2.els_e2x5Max()[index]/cms2.els_e5x5()[index]) > e2x5Over5x5Thresholds[0]) ele_result |= (1<<PASS_LSHAPE);
		if ((cms2.els_e2x5Max()[index]/cms2.els_e5x5()[index]) > e2x5Over5x5Thresholds_cand02[0]) ele_result |= (1<<PASS_LSHAPE_CAND02);
		if (cms2.els_d0corr()[index] < d0Thresholds[0]) ele_result |= (1<<PASS_D0);
		if (iso_relsusy < 0.10) ele_result |= (1<<PASS_ISO);

		if (cms2.els_fbrem()[index] > 0.20) {
			ele_result |= (1<<PASS_EXTRA);
			ele_result |= (1<<PASS_EXTRA_V2);
		}
		if (cms2.els_fbrem()[index] < 0.20) {
			if (cms2.els_eOverPIn()[index] > 0.7 && cms2.els_eOverPIn()[index] < 1.5) ele_result |= (1<<PASS_EXTRA);
			if (cms2.els_eOverPIn()[index] > 0.9 && cms2.els_eOverPIn()[index] < 1.5 && fabs(cms2.els_dEtaIn()[index]) < 0.03) ele_result |= (1<<PASS_EXTRA_V2);
		}

	}
	if (fabs(cms2.els_etaSC()[index]) > 1.479) {
		if (fabs(cms2.els_dEtaIn()[index]) < dEtaInThresholds[1]) 	ele_result |= (1<<PASS_DETA);
		if (fabs(cms2.els_dEtaIn()[index]) < dEtaInThresholds_cand02[1])   ele_result |= (1<<PASS_DETA_CAND02);
		if (fabs(cms2.els_dPhiIn()[index]) < dPhiInThresholds[1]) 	ele_result |= (1<<PASS_DPHI);
		if (cms2.els_hOverE()[index] < hoeThresholds[1]) 			ele_result |= (1<<PASS_HOE);
		if (cms2.els_sigmaIEtaIEta()[index] < sigmaIEtaIEtaThresholds[1]) 	ele_result |= (1<<PASS_LSHAPE) | (1<<PASS_LSHAPE_CAND02);
		if (cms2.els_d0corr()[index] < d0Thresholds[1]) ele_result |= (1<<PASS_D0);
		if (iso_relsusy < 0.10) ele_result |= (1<<PASS_ISO);

		if (cms2.els_fbrem()[index] > 0.20) {
			ele_result |= (1<<PASS_EXTRA);
			ele_result |= (1<<PASS_EXTRA_V2);
		}
		//if (cms2.els_fbrem()[index] < 0.15) {
		//    if (cms2.els_eOverPIn()[index] > 0.7 && cms2.els_eOverPIn()[index] < 1.5) ele_result |= (1<<PASS_EXTRA);
		//}


	}

	//
	// fill denominator histograms
	//

	DileptonHypType hypType = DILEPTON_EE;

	// find sc index
	int scidx = -1;
	for (size_t s = 0; s < cms2.evt_nscs(); ++s) {
		if (cms2.scs_elsidx()[s] == index) {
			scidx = s;
			break;
		}
	}

	// find closest mu
	int muidx = -1;
	float closestMu = 999.99;
	for (size_t m = 0; m < cms2.mus_p4().size(); ++m) {
		double dR = dRbetweenVectors(cms2.mus_p4()[m], cms2.els_p4()[index]); 
		if (dR < closestMu) {
			muidx = m;
			closestMu = dR;
		}
	}

	float E2x5MaxOver5x5 = cms2.els_e2x5Max()[index] / cms2.els_e5x5()[index];
	float E1x5Over5x5 = cms2.els_e1x5()[index] / cms2.els_e5x5()[index];

	if (fabs(cms2.els_etaSC()[index]) > 1.479) {
		Fill(h1_hyp_lt_ee_pt_, hypType, cms2.els_p4()[index].Pt(), weight);
		Fill(h1_hyp_lt_ee_hoe_, hypType, cms2.els_hOverE()[index], weight);
		Fill(h1_hyp_lt_ee_d0_, hypType, fabs(cms2.els_d0corr()[index]), weight);
		Fill(h1_hyp_lt_ee_dPhiIn_, hypType, fabs(cms2.els_dPhiIn()[index]), weight);
		Fill(h1_hyp_lt_ee_dEtaIn_, hypType, fabs(cms2.els_dEtaIn()[index]), weight);
		Fill(h1_hyp_lt_ee_sigmaIEtaIEta_, hypType, cms2.els_sigmaIEtaIEta()[index], weight);
		Fill(h1_hyp_lt_ee_E2x5MaxOver5x5_, hypType, E2x5MaxOver5x5, weight);
		Fill(h1_hyp_lt_ee_ecalIso_, hypType, cms2.els_ecalIso()[index], weight);
		Fill(h1_hyp_lt_ee_hcalIso_, hypType, cms2.els_hcalIso()[index], weight);
		Fill(h1_hyp_lt_ee_tkIso_, hypType, cms2.els_tkIso()[index], weight);
		Fill(h1_hyp_lt_ee_relsusy_, hypType, iso_relsusy, weight);

		if ((ele_result & (ele_passall & ~(1<<PASS_DETA))) == (ele_passall & ~(1<<PASS_DETA))) Fill(h1_hyp_lt_ee_nm1_dEtaIn_, hypType, fabs(cms2.els_dEtaIn()[index]), weight);
		if ((ele_result & (ele_passall & ~(1<<PASS_DPHI))) == (ele_passall & ~(1<<PASS_DPHI))) Fill(h1_hyp_lt_ee_nm1_dPhiIn_, hypType, fabs(cms2.els_dPhiIn()[index]), weight);
		if ((ele_result & (ele_passall & ~(1<<PASS_HOE))) == (ele_passall & ~(1<<PASS_HOE))) Fill(h1_hyp_lt_ee_nm1_hoe_, hypType, fabs(cms2.els_hOverE()[index]), weight);
		if ((ele_result & (ele_passall & ~(1<<PASS_LSHAPE))) == (ele_passall & ~(1<<PASS_LSHAPE))) {
			Fill(h1_hyp_lt_ee_nm1_sigmaIEtaIEta_, hypType, fabs(cms2.els_sigmaIEtaIEta()[index]), weight);
			Fill2D(h1_hyp_lt_ee_nm1_lateral_, hypType, E2x5MaxOver5x5, E1x5Over5x5, weight);
		}
		if ((ele_result & (ele_passall & ~(1<<PASS_D0))) == (ele_passall & ~(1<<PASS_D0))) Fill(h1_hyp_lt_ee_nm1_d0_, hypType, fabs(cms2.els_d0corr()[index]), weight);

		if ((ele_result &ele_passall_id_cand02) == (ele_passall_id_cand02)) {
			Fill(h1_hyp_lt_ee_afterid_relsusy_, hypType, iso_relsusy, weight);
			Fill(h1_hyp_lt_ee_afterid_fbrem_, hypType, cms2.els_fbrem()[index], weight);
			Fill(h1_hyp_lt_ee_afterid_eopin_, hypType, cms2.els_eOverPIn()[index], weight);
			if (cms2.els_fbrem()[index] > 0.2) {
				Fill(h1_hyp_lt_ee_afterid_relsusy_highfbrem_, hypType, iso_relsusy, weight);
				Fill(h1_hyp_lt_ee_afterid_eopin_highfbrem_, hypType, cms2.els_eOverPIn()[index], weight);
				Fill(h1_hyp_lt_ee_afterid_dPhiIn_highfbrem_, hypType, fabs(cms2.els_dPhiIn()[index]), weight);				
				Fill(h1_hyp_lt_ee_afterid_dEtaIn_highfbrem_, hypType, fabs(cms2.els_dEtaIn()[index]), weight);
				Fill(h1_hyp_lt_ee_afterid_preshowerEnergy_highfbrem_, hypType, cms2.scs_preshowerEnergy()[scidx]/cms2.scs_rawEnergy()[scidx], weight);
			}
			else {
				Fill(h1_hyp_lt_ee_afterid_relsusy_lowfbrem_, hypType, iso_relsusy, weight);
				Fill(h1_hyp_lt_ee_afterid_eopin_lowfbrem_, hypType, cms2.els_eOverPIn()[index], weight);
				Fill(h1_hyp_lt_ee_afterid_dPhiIn_lowfbrem_, hypType, fabs(cms2.els_dPhiIn()[index]), weight); 
				Fill(h1_hyp_lt_ee_afterid_dEtaIn_lowfbrem_, hypType, fabs(cms2.els_dEtaIn()[index]), weight);
				Fill(h1_hyp_lt_ee_afterid_preshowerEnergy_lowfbrem_, hypType, cms2.scs_preshowerEnergy()[scidx]/cms2.scs_rawEnergy()[scidx], weight);
			}
		}

		if ((ele_result & ele_passall_id_and_iso_cand01) == (ele_passall_id_and_iso_cand01)) Fill(h1_hyp_lt_ee_pt_cand01_, hypType, cms2.els_p4()[index].Pt(), weight);
		if ((ele_result & ele_passall_id_and_iso_cand02) == (ele_passall_id_and_iso_cand02)) Fill(h1_hyp_lt_ee_pt_cand02_, hypType, cms2.els_p4()[index].Pt(), weight);
		if ((ele_result & ele_passall_id_and_iso_cand02_extra) == (ele_passall_id_and_iso_cand02_extra)) Fill(h1_hyp_lt_ee_pt_cand02_extra_, hypType, cms2.els_p4()[index].Pt(), weight);
		if ((ele_result & ele_passall_id_and_iso_cand02_extra_v2) == (ele_passall_id_and_iso_cand02_extra_v2)) {
			Fill(h1_hyp_lt_ee_pt_cand02_extra_v2_, hypType, cms2.els_p4()[index].Pt(), weight);
			Fill(h1_hyp_lt_ee_eta_cand02_extra_v2_, hypType, cms2.els_etaSC()[index], weight);
			std::cout << cms2.evt_run() << ", " << cms2.evt_event() << ", " << index << std::endl;
			Fill(h1_hyp_lt_ee_afterid_closestMu_, hypType, closestMu, weight);
		}
	}

	if (fabs(cms2.els_etaSC()[index]) < 1.479) {
		Fill(h1_hyp_lt_eb_pt_, hypType, cms2.els_p4()[index].Pt(), weight);
		Fill(h1_hyp_lt_eb_hoe_, hypType, cms2.els_hOverE()[index], weight);
		Fill(h1_hyp_lt_eb_d0_, hypType, fabs(cms2.els_d0corr()[index]), weight);
		Fill(h1_hyp_lt_eb_dPhiIn_, hypType, fabs(cms2.els_dPhiIn()[index]), weight);
		Fill(h1_hyp_lt_eb_dEtaIn_, hypType, fabs(cms2.els_dEtaIn()[index]), weight);
		Fill(h1_hyp_lt_eb_sigmaIEtaIEta_, hypType, cms2.els_sigmaIEtaIEta()[index], weight);
		Fill(h1_hyp_lt_eb_E2x5MaxOver5x5_, hypType, E2x5MaxOver5x5, weight);
		Fill(h1_hyp_lt_eb_ecalIso_, hypType, cms2.els_ecalIso()[index], weight);
		Fill(h1_hyp_lt_eb_hcalIso_, hypType, cms2.els_hcalIso()[index], weight);
		Fill(h1_hyp_lt_eb_tkIso_, hypType, cms2.els_tkIso()[index], weight);
		Fill(h1_hyp_lt_eb_relsusy_, hypType, iso_relsusy, weight);

		if ((ele_result & (ele_passall & ~(1<<PASS_DETA))) == (ele_passall & ~(1<<PASS_DETA))) Fill(h1_hyp_lt_eb_nm1_dEtaIn_, hypType, fabs(cms2.els_dEtaIn()[index]), weight);
		if ((ele_result & (ele_passall & ~(1<<PASS_DPHI))) == (ele_passall & ~(1<<PASS_DPHI))) Fill(h1_hyp_lt_eb_nm1_dPhiIn_, hypType, fabs(cms2.els_dPhiIn()[index]), weight);
		if ((ele_result & (ele_passall & ~(1<<PASS_HOE))) == (ele_passall & ~(1<<PASS_HOE))) Fill(h1_hyp_lt_eb_nm1_hoe_, hypType, fabs(cms2.els_hOverE()[index]), weight);
		if ((ele_result & (ele_passall & ~(1<<PASS_LSHAPE))) == (ele_passall & ~(1<<PASS_LSHAPE))) {
			Fill(h1_hyp_lt_eb_nm1_E2x5MaxOver5x5_, hypType, E2x5MaxOver5x5, weight);
			Fill2D(h1_hyp_lt_eb_nm1_lateral_, hypType, E2x5MaxOver5x5, E1x5Over5x5, weight);
		}
		if ((ele_result & (ele_passall & ~(1<<PASS_D0))) == (ele_passall & ~(1<<PASS_D0))) Fill(h1_hyp_lt_eb_nm1_d0_, hypType, fabs(cms2.els_d0corr()[index]), weight);

		if ((ele_result & ele_passall_id_cand02) == (ele_passall_id_cand02)) {
			Fill(h1_hyp_lt_eb_afterid_relsusy_, hypType, iso_relsusy, weight);
			Fill(h1_hyp_lt_eb_afterid_fbrem_, hypType, cms2.els_fbrem()[index], weight);
			Fill(h1_hyp_lt_eb_afterid_eopin_, hypType, cms2.els_eOverPIn()[index], weight);
			if (cms2.els_fbrem()[index] > 0.2) {
				Fill(h1_hyp_lt_eb_afterid_relsusy_highfbrem_, hypType, iso_relsusy, weight);
				Fill(h1_hyp_lt_eb_afterid_eopin_highfbrem_, hypType, cms2.els_eOverPIn()[index], weight);
				Fill(h1_hyp_lt_eb_afterid_dPhiIn_highfbrem_, hypType, fabs(cms2.els_dPhiIn()[index]), weight);
				Fill(h1_hyp_lt_eb_afterid_dEtaIn_highfbrem_, hypType, fabs(cms2.els_dEtaIn()[index]), weight);
				Fill(h1_hyp_lt_eb_afterid_preshowerEnergy_highfbrem_, hypType, cms2.scs_preshowerEnergy()[scidx]/cms2.scs_rawEnergy()[scidx], weight);
			}
			else {
				Fill(h1_hyp_lt_eb_afterid_relsusy_lowfbrem_, hypType, iso_relsusy, weight);
				Fill(h1_hyp_lt_eb_afterid_eopin_lowfbrem_, hypType, cms2.els_eOverPIn()[index], weight);
				Fill(h1_hyp_lt_eb_afterid_dPhiIn_lowfbrem_, hypType, fabs(cms2.els_dPhiIn()[index]), weight); 
				Fill(h1_hyp_lt_eb_afterid_dEtaIn_lowfbrem_, hypType, fabs(cms2.els_dEtaIn()[index]), weight); 
				Fill(h1_hyp_lt_eb_afterid_preshowerEnergy_lowfbrem_, hypType, cms2.scs_preshowerEnergy()[scidx]/cms2.scs_rawEnergy()[scidx], weight);
			}
		}

		if ((ele_result & ele_passall_id_and_iso_cand01) == (ele_passall_id_and_iso_cand01)) Fill(h1_hyp_lt_eb_pt_cand01_, hypType, cms2.els_p4()[index].Pt(), weight);
		if ((ele_result & ele_passall_id_and_iso_cand02) == (ele_passall_id_and_iso_cand02)) Fill(h1_hyp_lt_eb_pt_cand02_, hypType, cms2.els_p4()[index].Pt(), weight);
		if ((ele_result & ele_passall_id_and_iso_cand02_extra) == (ele_passall_id_and_iso_cand02_extra)) Fill(h1_hyp_lt_eb_pt_cand02_extra_, hypType, cms2.els_p4()[index].Pt(), weight);
		if ((ele_result & ele_passall_id_and_iso_cand02_extra_v2) == (ele_passall_id_and_iso_cand02_extra_v2)) {
			Fill(h1_hyp_lt_eb_pt_cand02_extra_v2_, hypType, cms2.els_p4()[index].Pt(), weight);
			Fill(h1_hyp_lt_eb_eta_cand02_extra_v2_, hypType, cms2.els_etaSC()[index], weight);
			Fill(h1_hyp_lt_eb_afterid_closestMu_, hypType, closestMu, weight);
		}
	}

	// find out what passed
	bool ltPassOld = false;
	bool ltPassNew = false;
	bool relSusyIso = false;
	bool relSusyIso_cand0 = false;
	bool relSusyIso_cand1 = false;
	bool isConv = false;
	if (cms2.els_egamma_looseId().at(index)) ltPassOld = true;
	if (electronId_cand01(index)) ltPassNew = true;
	if (electronIsolation_relsusy(index, true) < 0.10) relSusyIso = true;
	if (electronIsolation_relsusy_cand0(index, true) < 0.10) relSusyIso_cand0 = true;
	if (electronIsolation_relsusy_cand1(index, true) < 0.10) relSusyIso_cand1 = true;
	if (isFromConversionPartnerTrack(index)) isConv = true;

	//
	// fill numerator histograms
	//
	if (abs(cms2.els_p4()[index].eta()) > 1.5) {
		if (ltPassOld) Fill(h1_hyp_lt_ee_pt_idold_, hypType, cms2.els_p4()[index].Pt(), weight);
		if (ltPassNew) Fill(h1_hyp_lt_ee_pt_idnew_, hypType, cms2.els_p4()[index].Pt(), weight);
		if (relSusyIso) Fill(h1_hyp_lt_ee_pt_isoold_, hypType, cms2.els_p4()[index].Pt(), weight);
		if (relSusyIso_cand0) Fill(h1_hyp_lt_ee_pt_isonew_, hypType, cms2.els_p4()[index].Pt(), weight);
		if (relSusyIso_cand1) Fill(h1_hyp_lt_ee_pt_isonew_cand1_, hypType, cms2.els_p4()[index].Pt(), weight);
		if (relSusyIso_cand1 && ltPassNew && !isConv) Fill(h1_hyp_lt_ee_pt_id1_iso1_conv_, hypType, cms2.els_p4()[index].Pt(), weight);
		if (!isConv) Fill(h1_hyp_lt_ee_pt_conv_, hypType, cms2.els_p4()[index].Pt(), weight);
	}
	if (abs(cms2.els_p4()[index].eta()) < 1.5)  {
		if (ltPassOld) Fill(h1_hyp_lt_eb_pt_idold_, hypType, cms2.els_p4()[index].Pt(), weight);
		if (ltPassNew) Fill(h1_hyp_lt_eb_pt_idnew_, hypType, cms2.els_p4()[index].Pt(), weight);
		if (relSusyIso) Fill(h1_hyp_lt_eb_pt_isoold_, hypType, cms2.els_p4()[index].Pt(), weight);
		if (relSusyIso_cand0) Fill(h1_hyp_lt_eb_pt_isonew_, hypType, cms2.els_p4()[index].Pt(), weight);
		if (relSusyIso_cand1) Fill(h1_hyp_lt_eb_pt_isonew_cand1_, hypType, cms2.els_p4()[index].Pt(), weight);
		if (relSusyIso_cand1 && ltPassNew && !isConv) Fill(h1_hyp_lt_eb_pt_id1_iso1_conv_, hypType, cms2.els_p4()[index].Pt(), weight);
		if (!isConv) Fill(h1_hyp_lt_eb_pt_conv_, hypType, cms2.els_p4()[index].Pt(), weight);
	}

}

//
// Main function
//
int MyScanChain::ScanChain(bool isData, std::string sampleName, TChain *chain, int nEvents, std::string skimFilePrefix) {

	//
	//
	//
	TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
	if (rootdir == 0){
		std::cout<<"Head directory root: not found. Try Rint: ..."<<std::endl;
		rootdir = gROOT->GetDirectory("Rint:");
		if (rootdir){
			std::cout<<"OK: Got Rint:"<<std::endl;
		} else {
			std::cout<<"ERROR: no root: or Rint: found. Histograms will likely be lost"<<std::endl;
		}
	} 

	//
	// format histograms
	//
	FormatAllEleIdHistograms(sampleName);

	// file loop
	//

	unsigned int nEventsChain=0;
	if(nEvents == -1) nEvents = chain->GetEntries();
	nEventsChain = nEvents;
	unsigned int nEventsTotal = 0;
	int i_permille_old = 0;

	TObjArray *listOfFiles = chain->GetListOfFiles();
	TIter fileIter(listOfFiles);
	TFile *currentFile = 0;
	while (TChainElement *currentFile = (TChainElement*)fileIter.Next()) {
		TFile *f = TFile::Open(currentFile->GetTitle());
		TTree *tree = (TTree*)f->Get("Events");
		cms2.Init(tree);

		//Event Loop
		ULong64_t nEvents = tree->GetEntries();
		for(ULong64_t event = 0; event < nEvents; ++event) {
			cms2.GetEntry(event);
			++nEventsTotal;

			// Progress feedback to the user
			int i_permille = (int)floor(1000 * nEventsTotal / float(nEventsChain));
			if (i_permille != i_permille_old) {
				// xterm magic from L. Vacavant and A. Cerri
				if (isatty(1)) {
					printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
							"\033[0m\033[32m <---\033[0m\015", i_permille/10.);
					fflush(stdout);
				}
				i_permille_old = i_permille;
			}

			// work out event weight
			float weight = cms2.evt_scale1fb()*0.01;

			//
			// Fill event level electron histograms
			//
			//FillAllEleIdHistogramsNoHyp(weight, sampleName);

			//
			// loop on hypothesis
			//
			for (size_t h = 0; h < cms2.hyp_type().size(); ++h) {

				//
				// fill basic electron ID histograms
				//
				FillAllEleIdHistogramsHyp(h, weight, sampleName);

			} // end loop on hypothesis

		}

	} // end loop on files

} // end loop on events

if ( nEventsChain != nEventsTotal ) {
	std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
}

rootdir = gROOT->GetDirectory("root:");
if (rootdir) rootdir->cd();
else{
	std::cout<<"Cant find root: . Current dir is "<<gDirectory->GetName()<<std::endl;
	rootdir = gROOT->GetDirectory("Rint:");
	if (rootdir){
		std::cout<<"OK, got Rint: "<<std::endl;
		rootdir->cd();
	} else {
		std::cout<<"Cant find Rint: either . Current dir is "<<gDirectory->GetName()<<std::endl;
	}
}

return 0;
}

