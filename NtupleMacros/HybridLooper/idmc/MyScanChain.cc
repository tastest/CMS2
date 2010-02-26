
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

static const char det_names[][128] = { "EB", "EE"};

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

void printCuts(elecuts_t result_electronSelections_cand01)
{
    if(result_electronSelections_cand01 & (1<<ELEPASS_ISO)) std::cout << "pass iso\n";
    else std::cout << "fail iso\n";

    if(result_electronSelections_cand01 & (1<<ELEPASS_ID)) std::cout << "pass id\n";
    else std::cout << "fail id\n";

    if(result_electronSelections_cand01 & (1<<ELEPASS_D0)) std::cout << "pass d0\n";
    else std::cout << "fail d0\n";

    if(result_electronSelections_cand01 & (1<<ELEPASS_NOTCONV)) std::cout << "pass conv\n";
    else std::cout << "fail conv\n";

    if(result_electronSelections_cand01 & (1<<ELEPASS_TYPE)) std::cout << "pass type\n";
    else std::cout << "fail type\n";

    if(result_electronSelections_cand01 & (1<<ELEPASS_FIDUCIAL)) std::cout << "pass eta\n";
    else std::cout << "fail eta\n";

    if(result_electronSelections_cand01 & (1<<ELEPASS_NOMUON)) std::cout << "pass muon\n";
    else std::cout << "fail muon\n";
}

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

    for (unsigned int i = 0; i < 2; ++i) {
        std::string detname = det_names[i];
        FormatHist(h1_hyp_ltid_sigmaIEtaIEta_[i], sampleName, "h1_hyp_ltid_sigmaIEtaIEta_" + detname, 50, 0.0, 0.05);
    }


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
	// fill denominator histograms
	//

	DileptonHypType hypType = DILEPTON_EE;

	bool pass_electronSelection_cand01 = electronSelection_cand01(index);
    elecuts_t result_electronSelections_cand01 = electronSelections_debug_;

    bool pass_electronId_cand01 = electronId_cand01(index);
    elecuts_t result_electronId_cand01 = electronId_debug_;

    bool checkAll = false;
    if ((result_electronSelections_cand01 & electronSelections_passall_) == result_electronSelections_cand01) checkAll = true;

    bool checkId = false;
    if ((result_electronId_cand01 & electronSelections_passid_) == result_electronId_cand01) checkId = true;

    unsigned int det = 0;
    if (fabs(cms2.els_etaSC()[index]) > 1.479) det = 1;


    Fill(h1_hyp_ltid_sigmaIEtaIEta_[det], hypType, cms2.els_sigmaIEtaIEta()[index], weight);


    std::cout << "id " << checkId << " all " << checkAll << std::endl;

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

		} // end loop on events

	} // end loop on files


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

