
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

// CMS2 includes
#include "CMS2.h"
#include "../CORE/electronSelections.h"
#include "../CORE/selections.h"

#include "../Tools/DileptonHypType.h"

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

//   jets_p4   
std::vector<LorentzVector> getCorCaloJets(int i_hyp) {
  std::vector<LorentzVector> calo_jets;
  calo_jets.clear();

  for (unsigned int jj=0; jj < cms2.jets_p4().size(); ++jj) {
    if ((dRbetweenVectors(cms2.hyp_lt_p4()[i_hyp],cms2.jets_p4()[jj]) < 0.4)||
        (dRbetweenVectors(cms2.hyp_ll_p4()[i_hyp],cms2.jets_p4()[jj]) < 0.4)
        ) continue;
    if (cms2.jets_p4()[jj].pt()*cms2.jets_cor()[jj] < 30) continue;
    if (fabs(cms2.jets_p4()[jj].Eta()) > 2.4) continue;
    //fkw July21 2009 if (cms2.jets_emFrac()[jj] < 0.1) continue;
    calo_jets.push_back(cms2.jets_p4()[jj]*cms2.jets_cor()[jj]);
  }

  if (calo_jets.size() > 1) {
       sort(calo_jets.begin(), calo_jets.end(),  comparePt);
  }
  return calo_jets;
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
	double             cands_passing[4];
	double             cands_passing_w2[4];
	unsigned int       cands_count[4];
	unsigned int 	   hyp_count[4];
        memset(cands_passing   , 0, sizeof(cands_passing       ));
        memset(cands_passing_w2        , 0, sizeof(cands_passing_w2    ));
        memset(cands_count             , 0, sizeof(cands_count         ));
        memset(hyp_count             , 0, sizeof(hyp_count         ));

	//
	// define histograms
	//
	TH1F *h1_nhyp[4];
	FormatHist(h1_nhyp, sampleName, "nhyp", 10, -0.5, 9.5);

	TH1F *h1_hyp_lt_pt[4];
	TH1F *h1_hyp_ll_pt[4];
        FormatHist(h1_hyp_lt_pt, sampleName, "hyp_lt_pt", 20, 0.0, 100.0);
        FormatHist(h1_hyp_ll_pt, sampleName, "hyp_ll_pt", 20, 0.0, 100.0);

        TH1F *h1_hyp_lt_pt_idnew[4];
        TH1F *h1_hyp_ll_pt_idnew[4];
        TH1F *h1_hyp_lt_pt_idold[4];
        TH1F *h1_hyp_ll_pt_idold[4];
        FormatHist(h1_hyp_lt_pt_idnew, sampleName, "hyp_lt_pt_idnew", 20, 0.0, 100.0);
        FormatHist(h1_hyp_ll_pt_idnew, sampleName, "hyp_ll_pt_idnew", 20, 0.0, 100.0);
        FormatHist(h1_hyp_lt_pt_idold, sampleName, "hyp_lt_pt_idold", 20, 0.0, 100.0);
        FormatHist(h1_hyp_ll_pt_idold, sampleName, "hyp_ll_pt_idold", 20, 0.0, 100.0);

	// file loop

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
			if (nEventsTotal % 1000 == 0)
				std::cout << "Event: " << nEventsTotal << std::endl;

			// work out event weight
		        float weight = cms2.evt_scale1fb();

			//
			// loop on hypothesis
			//
			std::vector<unsigned int> hyp_index_selected;
			for (size_t h = 0; h < cms2.hyp_type().size(); ++h) {

				// apply electron id
				if (!looseLeptonSelectionNoIsoTTDil08(cms2.hyp_lt_id()[h], cms2.hyp_lt_index()[h])) continue;
                                if (!looseLeptonSelectionNoIsoTTDil08(cms2.hyp_ll_id()[h], cms2.hyp_ll_index()[h])) continue;

				// apply isolation
				if (!passLeptonIsolationTTDil08(cms2.hyp_lt_id()[h], cms2.hyp_lt_index()[h])) continue;
                                if (!passLeptonIsolationTTDil08(cms2.hyp_ll_id()[h], cms2.hyp_ll_index()[h])) continue;

				// opposite charge
				if (cms2.hyp_lt_charge()[h] * cms2.hyp_ll_charge()[h] > 0) continue;

				// z mass window
				if (inZmassWindow(cms2.hyp_p4()[h].M())) continue;
			
				//
				// store indices of hypothesese that pass this preselection
				//	
				hyp_index_selected.push_back(h);

			} // end loop on hypothesis

			//
			// perform hypothesis disambiguation
			//

			int strasbourgDilType = 0;
			if (!hyp_index_selected.size()) continue;
			int hyp = eventDilIndexByWeightTTDil08(hyp_index_selected, strasbourgDilType, false, false);
		
			//
			// make requirements of the selected hypothesis
			//

			// trigger
			if (!passTriggersMu9orLisoE15(cms2.hyp_type()[hyp])) continue;

			// met
			if (!passMet_OF20_SF30(hyp, true)) continue;

			//
			// If we got to here then classify the event according to nJets
			//	

			std::vector<LorentzVector> corCaloJets = getCorCaloJets(hyp);
			std::cout << corCaloJets.size() << std::endl;


		} // end loop on files

	} // end loop on events

	if ( nEventsChain != nEventsTotal ) {
		std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
	}

	return 0;
}

