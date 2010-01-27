
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
#include "../Tools/DileptonHypType.h"

//
// Namespaces
//
using namespace tas;

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

			// apply event level cuts
			if (cms2.evt_tcmet() < 30.0) continue;

			// loop on hypothesis
			for (size_t h = 0; h < cms2.hyp_type().size(); ++h) {

				DileptonHypType hyp = hyp_typeToHypType(cms2.hyp_type()[h]);
				hyp_count[hyp] ++;
			
				Fill(h1_hyp_lt_pt, hyp, cms2.hyp_lt_p4()[h].Pt(), weight);
                               	Fill(h1_hyp_ll_pt, hyp, cms2.hyp_ll_p4()[h].Pt(), weight);

				if (hyp == DILEPTON_EE) {
					// new ID
					if (electronId_cand01(cms2.hyp_lt_index()[h])) Fill(h1_hyp_lt_pt_idnew, hyp, cms2.hyp_lt_p4()[h].Pt(), weight);
					if (electronId_cand01(cms2.hyp_lt_index()[h])) Fill(h1_hyp_ll_pt_idnew, hyp, cms2.hyp_ll_p4()[h].Pt(), weight);
					// old ID
                                        if (electronId_classBasedTight(cms2.hyp_lt_index()[h])) Fill(h1_hyp_lt_pt_idold, hyp, cms2.hyp_lt_p4()[h].Pt(), weight);
                                        if (electronId_classBasedTight(cms2.hyp_lt_index()[h])) Fill(h1_hyp_ll_pt_idold, hyp, cms2.hyp_ll_p4()[h].Pt(), weight);
				}


			} // end loop on hypothesis


			// fill histograms of number of hypothesis found
			for (unsigned int i = 0; i < 4; ++i) {
				Fill(h1_nhyp, i, hyp_count[i], 1);
				hyp_count[i] = 0;
			}
	

		} // end loop on files

	} // end loop on events

	if ( nEventsChain != nEventsTotal ) {
		std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
	}

	return 0;
}

