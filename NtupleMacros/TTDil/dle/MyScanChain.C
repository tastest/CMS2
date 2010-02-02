
//
// ttbar -> ll
// Dave "the one but not the only" Evans 
//

// C++ includes
#include <iostream>
#include <vector>

// ROOT includes
#include "TChain.h"
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

void Fill(TH1F** hist, unsigned int hyp, float val, float weight)
{
	hist[hyp]->Fill(val, weight);
	hist[DILEPTON_ALL]->Fill(val, weight);
}

void FormatHist(TH1F** hist, std::string sampleName, std::string name, Int_t n, Float_t min, Float_t max)
{       
	// loop on EB, EE
	for (unsigned int i = 0; i < 4; ++i)
	{
		std::string str = dilepton_hypo_names[i];
		std::string title = name + "_" + str;
		hist[i] = new TH1F(Form("%s_%s_%s", sampleName.c_str(), name.c_str(), str.c_str()),
				title.c_str(), n, min, max);
		hist[i]->GetXaxis()->SetTitle(name.c_str());
	}
}    

//
// Main function
//
int ScanChain(bool isData, std::string sampleName, TChain *chain, int nEvents = -1, std::string skimFilePrefix="") {

	TObjArray *listOfFiles = chain->GetListOfFiles();

	unsigned int nEventsChain=0;
	if(nEvents==-1) 
		nEvents = chain->GetEntries();
	nEventsChain = nEvents;
	unsigned int nEventsTotal = 0;
	TDirectory *rootdir = gDirectory->GetDirectory("Rint:");

	//
	// define counters
	//

	// count the (weighted and unweighted) number of candidates passing our cuts
	float             cands_passing[4];
	float             cands_passing_w2[4];
	unsigned int       cands_count[4];
	unsigned int 	   hyp_count[4];
	memset(cands_passing   , 0, sizeof(cands_passing       ));
	memset(cands_passing_w2        , 0, sizeof(cands_passing_w2    ));
	memset(cands_count             , 0, sizeof(cands_count         ));
	memset(hyp_count             , 0, sizeof(hyp_count         ));

	//
	//
	//
	TH1F *h1_njets[4];
        FormatHist(h1_njets, sampleName, "hyp_njets", 10, -0.5, 9.5);

	// file loop
	//
	TIter fileIter(listOfFiles);
	TFile *currentFile = 0;
	while ( currentFile = (TFile*)fileIter.Next() ) {
		TFile f(currentFile->GetTitle());
		TTree *tree = (TTree*)f.Get("Events");
		cms2.Init(tree);

		//Event Loop
		unsigned int nEvents = tree->GetEntries();
		for( unsigned int event = 0; event < nEvents; ++event) {
			cms2.GetEntry(event);
			++nEventsTotal;

			// print out event being processed
			//if (nEventsTotal % 1000 == 0)
			//	std::cout << "Event: " << nEventsTotal << std::endl;

			// work out event weight
			float weight = cms2.evt_scale1fb()*0.01;

			// does this event contain two generated leptons with mother as W?

			//
			// loop on hypothesis
			//
			std::vector<unsigned int> hyp_index_selected;
			hyp_index_selected.clear();
			for (size_t h = 0; h < cms2.hyp_type().size(); ++h) {

				// apply lepton id 
				if (!looseLeptonSelectionNoIsoTTDil08(cms2.hyp_lt_id()[h], cms2.hyp_lt_index()[h])) continue;
				if (!looseLeptonSelectionNoIsoTTDil08(cms2.hyp_ll_id()[h], cms2.hyp_ll_index()[h])) continue;

				// apply isolation
				if (!passLeptonIsolationTTDil08(cms2.hyp_ll_id()[h], cms2.hyp_ll_index()[h])) continue;
				if (!passLeptonIsolationTTDil08(cms2.hyp_lt_id()[h], cms2.hyp_lt_index()[h])) continue;

				// opposite charge
				if (cms2.hyp_lt_id()[h] * cms2.hyp_ll_id()[h] > 0) continue;

				// z mass window
            			if (cms2.hyp_type()[h] == 0 || cms2.hyp_type()[h] == 3) {
					if (inZmassWindow(cms2.hyp_p4()[h].mass())) continue;
				}

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

			// met
			if (!passMet_OF20_SF30(hyp, false)) continue;

			//
			// If we got to here then classify the event according to nJets
			//	

			std::vector<unsigned int> corCaloJets;
			corCaloJets.clear();
			for (size_t j = 0; j < cms2.jets_p4().size(); ++j) {
				if (!isGoodDilHypJet(j, hyp, 30.0, 2.4, 0.4, false)) continue;
				corCaloJets.push_back(j);
			}

                        DileptonHypType hypType = hyp_typeToHypType(cms2.hyp_type()[hyp]);
			Fill(h1_njets, hypType, corCaloJets.size(), weight);

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

	std::cout << "sampleName \t" << dilepton_hypo_names[3] << "\t" << dilepton_hypo_names[1] << "\t" << dilepton_hypo_names[2] << "\t" << dilepton_hypo_names[0] << std::endl;
	//        for (unsigned int i = 0; i < 4; ++i) {
	//                std::string str = dilepton_hypo_names[i];
	//		std::cout << dilepton_hypo_names[i] << "\t";
	//	}
	//	std::cout << std::endl;
	//        std::cout << sampleName << "\t";


	std::cout << cands_passing[3] << " $\pm$ " << sqrt(cands_passing_w2[3]) << "\t & ";
	std::cout << cands_passing[1] << " $\pm$ " << sqrt(cands_passing_w2[1]) << "\t & ";
	std::cout << cands_passing[2] << " $\pm$ " << sqrt(cands_passing_w2[2]) << "\t & ";
	std::cout << cands_passing[0] << " $\pm$ " << sqrt(cands_passing_w2[0]) << "\t & ";
	std::cout << std::endl;

	//	for (unsigned int i = 0; i < 4; ++i) {

	//		std::cout << cands_passing[i] << " $\pm$ " << sqrt(cands_passing_w2[i]) << "\t &";

	//	}
	//	std::cout << "\\\\ \\hline" << std::endl;

	return 0;
}

