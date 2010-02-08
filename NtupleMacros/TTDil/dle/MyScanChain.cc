
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
// for jets
//

static const char jetbin_names[][128] = { "0j", "1j", "2j"};

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

	FormatHist(h1_hyp_lt_eb_dPhiIn_, sampleName, "hyp_lt_eb_dPhiIn", 100, 0.0, 0.05);
	FormatHist(h1_hyp_lt_ee_dPhiIn_, sampleName, "hyp_lt_ee_dPhiIn", 100, 0.0, 0.05);

	FormatHist(h1_hyp_lt_eb_d0_, sampleName, "hyp_lt_eb_d0", 100, 0.0, 0.05);
	FormatHist(h1_hyp_lt_ee_d0_, sampleName, "hyp_lt_ee_d0", 100, 0.0, 0.05);

	FormatHist(h1_hyp_lt_eb_E2x5MaxOver5x5_, sampleName, "hyp_lt_eb_E2x5MaxOver5x5", 110, 0.0, 1.1);
	FormatHist(h1_hyp_lt_ee_E2x5MaxOver5x5_, sampleName, "hyp_lt_ee_E2x5MaxOver5x5", 110, 0.0, 1.1);

	FormatHist(h1_hyp_lt_eb_ecalIso_, sampleName, "hyp_lt_eb_ecalIso", 100, 0, 10);
	FormatHist(h1_hyp_lt_eb_hcalIso_, sampleName, "hyp_lt_eb_hcalIso", 100, 0, 25);
	FormatHist(h1_hyp_lt_eb_tkIso_, sampleName, "hyp_lt_eb_tkIso", 100, 0, 25);

	FormatHist(h1_hyp_lt_ee_ecalIso_, sampleName, "hyp_lt_ee_ecalIso", 100, 0, 10);
	FormatHist(h1_hyp_lt_ee_hcalIso_, sampleName, "hyp_lt_ee_hcalIso", 100, 0, 25);
	FormatHist(h1_hyp_lt_ee_tkIso_, sampleName, "hyp_lt_ee_tkIso", 100, 0, 25);

	FormatHist(h1_hyp_lt_ee_relsusy_, sampleName, "hyp_lt_ee_relsusy", 100, 0, 1);
    FormatHist(h1_hyp_lt_eb_relsusy_, sampleName, "hyp_lt_eb_relsusy", 100, 0, 1);

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

void MyScanChain::FormatAllAnaHistograms(std::string sampleName)
{
	FormatHist(h1_hyp_njets_, sampleName, "hyp_njets", 10, -0.5, 9.5);
}

void MyScanChain::FormatAllDYEstHistograms(std::string sampleName)
{

	for (unsigned int j = 0; j < 3; ++j) {
		std::string jetbin = jetbin_names[j];
    	FormatHist(h1_dyest_mll_met_[j], sampleName, "dyest_mll_met_" + jetbin, 40, 0.0, 200.0);
        FormatHist(h1_dyest_mll_nomet_[j], sampleName, "dyest_mll_nomet_" + jetbin, 40, 0.0, 200.0);
		FormatHist(h1_dyest_met_in_[j], sampleName, "dyest_met_in_" + jetbin, 40, 0.0, 200.0);
        FormatHist(h1_dyest_met_out_[j], sampleName, "dyest_met_out_" + jetbin, 40, 0.0, 200.0);
	}

}

void MyScanChain::FillAllDYEstHistograms(const unsigned int h, const float &weight, const unsigned int njet)
{

	// which hypothesis type
    DileptonHypType hypType = hyp_typeToHypType(cms2.hyp_type()[h]); 

	// which jet bin to fill
	unsigned int jetbin = 0;
	if (njet == 1) jetbin = 1;
	if (njet >= 2) jetbin = 2;

	// fill the mass histogram
	float mass = cms2.hyp_p4()[h].mass();
	Fill(h1_dyest_mll_nomet_[jetbin], hypType, mass, weight);
    if (passMet_OF20_SF30(h, false)) Fill(h1_dyest_mll_met_[jetbin], hypType, mass, weight);

	// fill the met histograms for "in" and "out" regions
	float mymet = met_pat_metCor_hyp(h);
	if (inZmassWindow(mass)) {
		Fill(h1_dyest_met_in_[jetbin], hypType, mymet, weight);
	}
	else Fill(h1_dyest_met_out_[jetbin], hypType, mymet, weight);

}

void MyScanChain::FillAllEleIdHistograms(const unsigned int h, const float &weight, const TString &sampleName)
{

	DileptonHypType hypType = hyp_typeToHypType(cms2.hyp_type()[h]);

	// apply truth match behavior if ttbar
	if (sampleName == "ttbar") {
		if (abs(cms2.hyp_lt_id()[h]) == 11 && 
			!((abs(cms2.hyp_lt_mc_id()[h]) == 11) && abs(cms2.hyp_lt_mc_motherid()[h]) == 24) ) return;
	}
	// apply truth match behavior if wjets
    if (sampleName == "wjets") {
    	if (abs(cms2.hyp_lt_id()[h]) == 11 &&
            ((abs(cms2.hyp_lt_mc_id()[h]) == 11) && abs(cms2.hyp_lt_mc_motherid()[h]) == 24) ) return;
    }
	// apply truth match behavior if dyee
    if (sampleName == "dyee") {
        if (abs(cms2.hyp_lt_id()[h]) == 11 &&
            !((abs(cms2.hyp_lt_mc_id()[h]) == 11) && abs(cms2.hyp_lt_mc_motherid()[h]) == 23) ) return;
    }


	// apply part of electron denominator first here
	if (abs(cms2.hyp_lt_id()[h]) == 11) {
		if (!electron20Eta2p4(cms2.hyp_lt_index()[h])) return;
		if (!(cms2.els_type()[cms2.hyp_lt_index()[h]] & (1<<ISECALDRIVEN))) return;
    }
	if (abs(cms2.hyp_ll_id()[h]) == 11) {
		if (!electron20Eta2p4(cms2.hyp_ll_index()[h])) return;
        if (!(cms2.els_type()[cms2.hyp_ll_index()[h]] & (1<<ISECALDRIVEN))) return;
	}

	//
	// fill denominator histograms
	//

	float iso_relsusy = -1;
    if (abs(cms2.hyp_lt_id()[h]) == 11)	
		iso_relsusy = electronIsolation_relsusy_cand1(cms2.hyp_lt_index()[h], true);

	if (abs(cms2.hyp_lt_p4()[h].eta()) > 1.5 && abs(cms2.hyp_lt_id()[h]) == 11) {
		Fill(h1_hyp_lt_ee_pt_, hypType, cms2.hyp_lt_p4()[h].Pt(), weight);
		Fill(h1_hyp_lt_ee_hoe_, hypType, cms2.els_hOverE()[cms2.hyp_lt_index()[h]], weight);
		Fill(h1_hyp_lt_ee_d0_, hypType, fabs(cms2.els_d0corr()[cms2.hyp_lt_index()[h]]), weight);
		Fill(h1_hyp_lt_ee_dPhiIn_, hypType, fabs(cms2.els_dPhiIn()[cms2.hyp_lt_index()[h]]), weight);
		Fill(h1_hyp_lt_ee_dEtaIn_, hypType, fabs(cms2.els_dEtaIn()[cms2.hyp_lt_index()[h]]), weight);
		Fill(h1_hyp_lt_ee_sigmaIEtaIEta_, hypType, cms2.els_sigmaIEtaIEta()[cms2.hyp_lt_index()[h]], weight);
		float E2x5MaxOver5x5 = cms2.els_e2x5Max()[cms2.hyp_lt_index()[h]] / cms2.els_e5x5()[cms2.hyp_lt_index()[h]];
		Fill(h1_hyp_lt_ee_E2x5MaxOver5x5_, hypType, E2x5MaxOver5x5, weight);
		Fill(h1_hyp_lt_ee_ecalIso_, hypType, cms2.els_ecalIso()[h], weight);
		Fill(h1_hyp_lt_ee_hcalIso_, hypType, cms2.els_hcalIso()[h], weight);
		Fill(h1_hyp_lt_ee_tkIso_, hypType, cms2.els_tkIso()[h], weight);
		Fill(h1_hyp_lt_ee_relsusy_, hypType, iso_relsusy, weight);
	}

	if (abs(cms2.hyp_lt_p4()[h].eta()) < 1.5 && abs(cms2.hyp_lt_id()[h]) == 11) {
		Fill(h1_hyp_lt_eb_pt_, hypType, cms2.hyp_lt_p4()[h].Pt(), weight);
		Fill(h1_hyp_lt_eb_hoe_, hypType, cms2.els_hOverE()[cms2.hyp_lt_index()[h]], weight);
		Fill(h1_hyp_lt_eb_d0_, hypType, fabs(cms2.els_d0corr()[cms2.hyp_lt_index()[h]]), weight);
		Fill(h1_hyp_lt_eb_dPhiIn_, hypType, fabs(cms2.els_dPhiIn()[cms2.hyp_lt_index()[h]]), weight);
		Fill(h1_hyp_lt_eb_dEtaIn_, hypType, fabs(cms2.els_dEtaIn()[cms2.hyp_lt_index()[h]]), weight);
		Fill(h1_hyp_lt_eb_sigmaIEtaIEta_, hypType, cms2.els_sigmaIEtaIEta()[cms2.hyp_lt_index()[h]], weight);
		float E2x5MaxOver5x5 = cms2.els_e2x5Max()[cms2.hyp_lt_index()[h]] / cms2.els_e5x5()[cms2.hyp_lt_index()[h]];
		Fill(h1_hyp_lt_eb_E2x5MaxOver5x5_, hypType, E2x5MaxOver5x5, weight);
		Fill(h1_hyp_lt_eb_ecalIso_, hypType, cms2.els_ecalIso()[h], weight);
		Fill(h1_hyp_lt_eb_hcalIso_, hypType, cms2.els_hcalIso()[h], weight);
		Fill(h1_hyp_lt_eb_tkIso_, hypType, cms2.els_tkIso()[h], weight);
        Fill(h1_hyp_lt_eb_relsusy_, hypType, iso_relsusy, weight);
	}

	// find out what passed
	bool ltPassOld = false;
	bool ltPassNew = false;
	bool relSusyIso = false;
	bool relSusyIso_cand0 = false;
	bool relSusyIso_cand1 = false;
	bool isConv = false;
	if (abs(cms2.hyp_lt_id()[h]) == 11) {
		if (cms2.els_egamma_looseId().at(cms2.hyp_lt_index()[h])) ltPassOld = true;
		if (electronId_cand01(cms2.hyp_lt_index()[h])) ltPassNew = true;
		if (electronIsolation_relsusy(cms2.hyp_lt_index()[h], true) < 0.10) relSusyIso = true;
		if (electronIsolation_relsusy_cand0(cms2.hyp_lt_index()[h], true) < 0.10) relSusyIso_cand0 = true;
		if (electronIsolation_relsusy_cand1(cms2.hyp_lt_index()[h], true) < 0.10) relSusyIso_cand1 = true;
		if (isFromConversionPartnerTrack(cms2.hyp_lt_index()[h])) isConv = true;
	}

	//
	// fill numerator histograms
	//
	if (abs(cms2.hyp_lt_p4()[h].eta()) > 1.5) {
		if (ltPassOld) Fill(h1_hyp_lt_ee_pt_idold_, hypType, cms2.hyp_lt_p4()[h].Pt(), weight);
		if (ltPassNew) Fill(h1_hyp_lt_ee_pt_idnew_, hypType, cms2.hyp_lt_p4()[h].Pt(), weight);
		if (relSusyIso) Fill(h1_hyp_lt_ee_pt_isoold_, hypType, cms2.hyp_lt_p4()[h].Pt(), weight);
		if (relSusyIso_cand0) Fill(h1_hyp_lt_ee_pt_isonew_, hypType, cms2.hyp_lt_p4()[h].Pt(), weight);
		if (relSusyIso_cand1) Fill(h1_hyp_lt_ee_pt_isonew_cand1_, hypType, cms2.hyp_lt_p4()[h].Pt(), weight);
		if (relSusyIso_cand1 && ltPassNew && !isConv) Fill(h1_hyp_lt_ee_pt_id1_iso1_conv_, hypType, cms2.hyp_lt_p4()[h].Pt(), weight);
		if (!isConv) Fill(h1_hyp_lt_ee_pt_conv_, hypType, cms2.hyp_lt_p4()[h].Pt(), weight);
	}
	if (abs(cms2.hyp_lt_p4()[h].eta()) < 1.5)  {
		if (ltPassOld) Fill(h1_hyp_lt_eb_pt_idold_, hypType, cms2.hyp_lt_p4()[h].Pt(), weight);
		if (ltPassNew) Fill(h1_hyp_lt_eb_pt_idnew_, hypType, cms2.hyp_lt_p4()[h].Pt(), weight);
		if (relSusyIso) Fill(h1_hyp_lt_eb_pt_isoold_, hypType, cms2.hyp_lt_p4()[h].Pt(), weight);
		if (relSusyIso_cand0) Fill(h1_hyp_lt_eb_pt_isonew_, hypType, cms2.hyp_lt_p4()[h].Pt(), weight);
		if (relSusyIso_cand1) Fill(h1_hyp_lt_eb_pt_isonew_cand1_, hypType, cms2.hyp_lt_p4()[h].Pt(), weight);
		if (relSusyIso_cand1 && ltPassNew && !isConv) Fill(h1_hyp_lt_eb_pt_id1_iso1_conv_, hypType, cms2.hyp_lt_p4()[h].Pt(), weight);
		if (!isConv) Fill(h1_hyp_lt_eb_pt_conv_, hypType, cms2.hyp_lt_p4()[h].Pt(), weight);
	}

}

//
// Main function
//
int MyScanChain::ScanChain(bool isData, std::string sampleName, TChain *chain, int nEvents, std::string skimFilePrefix) {

	//
	// define counters
	//

	// count the (weighted and unweighted) number of candidates passing our cuts
	float             cands_passing[4];
	float             cands_passing_w2[4];
	unsigned int       cands_count[4];
	memset(cands_passing   , 0, sizeof(cands_passing       ));
	memset(cands_passing_w2        , 0, sizeof(cands_passing_w2    ));
	memset(cands_count             , 0, sizeof(cands_count         ));

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
	FormatAllAnaHistograms(sampleName);
	FormatAllDYEstHistograms(sampleName);

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
			// loop on hypothesis
			//
			std::vector<unsigned int> hyp_index_selected;
			hyp_index_selected.clear();
			for (size_t h = 0; h < cms2.hyp_type().size(); ++h) {

				//
				// fill basic electron ID histograms
				//
				FillAllEleIdHistograms(h, weight, sampleName);

				// apply lepton id 
				if (!looseLeptonSelectionNoIsoTTDil08(cms2.hyp_lt_id()[h], cms2.hyp_lt_index()[h])) continue;
				if (!looseLeptonSelectionNoIsoTTDil08(cms2.hyp_ll_id()[h], cms2.hyp_ll_index()[h])) continue;

				// apply isolation
				if (!passLeptonIsolationTTDil08(cms2.hyp_ll_id()[h], cms2.hyp_ll_index()[h])) continue;
				if (!passLeptonIsolationTTDil08(cms2.hyp_lt_id()[h], cms2.hyp_lt_index()[h])) continue;

				// opposite charge
				if (cms2.hyp_lt_id()[h] * cms2.hyp_ll_id()[h] > 0) continue;

				//
				// store indices of hypothesese that pass this preselection
				//	
				hyp_index_selected.push_back(h);

			} // end loop on hypothesis

			//
			// perform hypothesis disambiguation
			//

			int strasbourgDilType = -1;
			if (hyp_index_selected.size() == 0) continue;
			int hyp = 0;
			hyp = eventDilIndexByWeightTTDil08(hyp_index_selected, strasbourgDilType, false, false);

			//
			// make requirements of the selected hypothesis
			//

			// trigger
			if (!passTriggersMu9orLisoE15(cms2.hyp_type()[hyp])) continue;

			//
			// If we got to here then classify the event according to nJets
			//	

			std::vector<unsigned int> corCaloJets;
			corCaloJets.clear();
			for (size_t j = 0; j < cms2.jets_p4().size(); ++j) {
				if (!isGoodDilHypJet(j, hyp, 30.0, 2.4, 0.4, false)) continue;
				corCaloJets.push_back(j);
			}

            //
            // estimate the DY background before the z veto and MET cuts are applied
            //

            FillAllDYEstHistograms(hyp, weight, corCaloJets.size());

			//
			// apply remaining requirements
			//

            // met
            if (!passMet_OF20_SF30(hyp, false)) continue;

            // z mass window
            if (cms2.hyp_type()[hyp] == 0 || cms2.hyp_type()[hyp] == 3) {
                if (inZmassWindow(cms2.hyp_p4()[hyp].mass())) continue;
            }

			//
			// get hypothesis type and fill analysis results histograms
			//

			DileptonHypType hypType = hyp_typeToHypType(cms2.hyp_type()[hyp]);
			Fill(h1_hyp_njets_, hypType, corCaloJets.size(), weight);

			//
			// count...
			//

			if (corCaloJets.size() >= 2) {

				cands_passing[hypType] += weight;
				cands_passing_w2[hypType] += weight * weight;
				cands_count[hypType] ++;

				cands_passing[DILEPTON_ALL] += weight;
				cands_passing_w2[DILEPTON_ALL] += weight * weight;
				cands_count[DILEPTON_ALL] ++;

				//std::cout << "event number, hyp type: " << cms2.evt_event() << " \t" << cms2.hyp_type()[hyp] << std::endl;

			}


		} // end loop on files

	} // end loop on events

	if ( nEventsChain != nEventsTotal ) {
		std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
	}

	//
	// print table entry
	//

	std::cout.flush();
	std::cout << std::endl;
	for (unsigned int i = 0; i < 4; ++i) {
		std::string str = dilepton_hypo_names[i];
		std::cout << " & " << dilepton_hypo_names[i] << "\t";
	}
	std::cout << "\\\\ \\hline" << std::endl;
	std::cout << sampleName << "\t";
	for (unsigned int i = 0; i < 4; ++i) {
		std::cout << " & " << cands_passing[i] << " $\\pm$ " << sqrt(cands_passing_w2[i]) << "\t";
	}
	std::cout << "\\\\ \\hline" << std::endl;

	//
	// make sure we're back in the right root dir
	//

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

