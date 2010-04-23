
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
#include "../../CORE/CMS2.h"

#include "../../CORE/eventSelections.h"
#include "../../CORE/electronSelections.h"
#include "../../CORE/metSelections.h"
#include "../../CORE/jetSelections.h"

#include "../../CORE/ttbarSelections.h"

#include "../../Tools/tools.cc"

//
// Namespaces
//
using namespace tas;

//
// definitions...
//

enum {

    PASS_TOP_PT,

};

//
// for jets
//

static const char jetbin_names[][128] = { "0j", "1j", "2j"};

//
// functions
//

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
    if (passMetAsIs_OF20_SF30(cms2.evt_pfmet(), hypType)) Fill(h1_dyest_mll_met_[jetbin], hypType, mass, weight);

    // fill the met histograms for "in" and "out" regions
    float mymet = cms2.evt_pfmet();
    if (inZmassWindow(mass)) {
        Fill(h1_dyest_met_in_[jetbin], hypType, mymet, weight);
    }
    else Fill(h1_dyest_met_out_[jetbin], hypType, mymet, weight);

}

//
// Main function
//

int MyScanChain::ScanChain(bool isData, std::string sampleName, TChain *chain, int nEvents, std::string skimFilePrefix) {

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

	FormatAllAnaHistograms(sampleName);
	FormatAllDYEstHistograms(sampleName);

    //
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

				// apply lepton id and iso
                if (!isGoodHypwIso(h)) continue;

				// opposite charge
				if (cms2.hyp_lt_id()[h] * cms2.hyp_ll_id()[h] > 0) continue;

				//
				// store indices of hypothesese that pass this preselection
				//	
				hyp_index_selected.push_back(h);

			} // end loop on hypothesis

            //
            // require a hypothesis was found
            //

            if (hyp_index_selected.size() == 0) continue;
            int hyp = 0;

            //
			// perform hypothesis disambiguation and get hyp type for selected hyp
			//

			int strasbourgDilType = -1;
			hyp = eventDilIndexByWeightTTDil08(hyp_index_selected, strasbourgDilType, false, false);
            DileptonHypType hypType = hyp_typeToHypType(cms2.hyp_type()[hyp]);

			//
			// make requirements of the selected hypothesis
			//

			// trigger
			//if (!passTriggersMu9orLisoE15(cms2.hyp_type()[hyp])) continue;

			//
			// If we got to here then classify the event according to nJets
			//	
std::cout << cms2.pfjets_p4().size() << std::endl;
std::cout << cms2.pfjets_cor().size() << std::endl;

unsigned int nJetsFound = 0;
            //unsigned int nJetsFound = nJets(hyp, JETS_TYPE_PF_CORR, JETS_CLEAN_HYP_E_MU, 0.4, 30.0, 2.4);

			//
			// estimate the DY background before the z veto and MET cuts are applied
			//

			FillAllDYEstHistograms(hyp, weight, nJetsFound);

			//
			// apply remaining requirements
			//

			// met
			if (!passMetAsIs_OF20_SF30(cms2.evt_pfmet(), hypType)) continue;

			// z mass window
			if (cms2.hyp_type()[hyp] == 0 || cms2.hyp_type()[hyp] == 3) {
				if (inZmassWindow(cms2.hyp_p4()[hyp].mass())) continue;
			}

			//
			// fill analysis results histograms
			//

			Fill(h1_hyp_njets_, hypType, nJetsFound, weight);

		} // end loop on files

	} // end loop on events

	if ( nEventsChain != nEventsTotal ) {
		std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
	}

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

