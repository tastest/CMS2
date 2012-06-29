//#include "SS2012AnalysisFilter.hpp"

#include <iostream>
#include <assert.h>
#include <string>
#include "TChain.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TTree.h"
#include "TTreeCache.h"
#include "Rtypes.h"
#ifdef __NON_ROOT_BUILD__
#include "CMS2.h"
#include "CORE/muonSelections.h"
#include "CORE/electronSelections.h"
#include "CORE/ssSelections.h"
#include "CORE/jetSelections.h"
#include "CORE/trackSelections.h"
#include "CORE/eventSelections.h"
#include "CORE/jetcorr/FactorizedJetCorrector.h"
#include "Tools/DileptonHypType.h"
#else
#include "CMS2.h"
#include "CORE/muonSelections.h"
#include "CORE/electronSelections.h"
#include "CORE/ssSelections.h"
#include "CORE/jetSelections.h"
#include "CORE/trackSelections.h"
#include "CORE/eventSelections.h"
#include "CORE/jetcorr/FactorizedJetCorrector.h"
#include "Tools/DileptonHypType.cc"
#endif

typedef ULong64_t uint64;

using namespace tas;

bool selectSameSign(bool do20_20, bool btag2, bool jets2, FactorizedJetCorrector* jet_pf_corrector=NULL)
{
    // loop over hypothesis
    for (unsigned int ihyp = 0; ihyp < cms2.hyp_p4().size(); ihyp++)
    {
		// eta < 2.5
        if (fabs(cms2.hyp_lt_p4().at(ihyp).eta()) > 2.5)
        {
            continue;
        }
        if (fabs(cms2.hyp_ll_p4().at(ihyp).eta()) > 2.5)
        {
            continue;
        }

        int lt_id   = cms2.hyp_lt_id().at(ihyp);
        int lt_idx  = cms2.hyp_lt_index().at(ihyp);
        float lt_pt = cms2.hyp_lt_p4().at(ihyp).pt();
        int ll_id   = cms2.hyp_ll_id().at(ihyp);
        int ll_idx  = cms2.hyp_ll_index().at(ihyp);
        float ll_pt = cms2.hyp_ll_p4().at(ihyp).pt();

        // require 20/10 GeV leptons
        if (min(lt_pt, ll_pt) < 10.0 || max(lt_pt, ll_pt) < 20.0)
        {
            continue;
        }
		if (do20_20)
		{
        	if (min(lt_pt, ll_pt) < 20.0)
        	{
        	    continue;
        	}
		}

        // jet type
        JetType jet_type = evt_isRealData() ? JETS_TYPE_PF_FAST_CORR_RESIDUAL : JETS_TYPE_PF_FAST_CORR;

        if (btag2)
        { 
            int nbtags = jet_pf_corrector ? samesign::nBtaggedJets(ihyp, jet_pf_corrector, jet_type, JETS_BTAG_CSVM, /*dR=*/0.4, /*jet_pt>*/40.0, /*|eta|<*/2.4, /*pt1>*/20.0, /*pt1>*/20.0) :
                                            samesign::nBtaggedJets(ihyp, jet_type, JETS_BTAG_CSVM,                   /*dR=*/0.4, /*jet_pt>*/40.0, /*|eta|<*/2.4, /*pt1>*/20.0, /*pt1>*/20.0);
            if (nbtags < 2)
            {
                continue;
            }
        }
        else if (jets2)
        {
            int njets = jet_pf_corrector  ? samesign::nJets(ihyp, jet_pf_corrector, jet_type, /*dR=*/0.4, /*jet_pt>*/40.0, /*|eta|<*/2.4, /*pt1>*/20.0, /*pt2>*/20.0) : 
                                            samesign::nJets(ihyp, jet_type,                   /*dR=*/0.4, /*jet_pt>*/40.0, /*|eta|<*/2.4, /*pt1>*/20.0, /*pt2>*/20.0);
            if (njets < 2)
            {
                continue;
            }
        }

		// need at least one good vertex
        int nGoodVertices = numberOfGoodVertices();
        if (nGoodVertices < 1)
        {
            continue;
        }

        // same sign
        if ((lt_id * ll_id) > 0)
        {
            // keep if FO/FO
            if (!samesign::isDenominatorLepton(lt_id, lt_idx) || !samesign::isDenominatorLepton(ll_id, ll_idx))
            {
                continue;
            }
            return true;
        }
        // oppposite sign and data
        else if ((lt_id * ll_id) < 0)
        {
            //if (jets2 || evt_isRealData())
            {
                // select specific events
                bool ee_or_em = ((hyp_typeToHypType(hyp_type().at(ihyp)) == DILEPTON_EE) || (hyp_typeToHypType(hyp_type().at(ihyp)) == DILEPTON_EMU));
                if (not (ee_or_em && samesign::isNumeratorHypothesis(ihyp)))
                {
                    continue;
                }
                return true;
            }        
        }
    }

    return false;
}

void ntupleFilterSameSign(const std::string &infile, const std::string &outfile, const std::string& jetcorr_path, bool do20_20, bool btag2, bool jets2, bool printPass)  
{
    // output file and tree
    TFile *output =TFile::Open(outfile.c_str(), "RECREATE");
    assert(output != 0);
    TTree *newtree = 0;

    const long long max_tree_size = 20000000000000000LL;
    TTree::SetMaxTreeSize(max_tree_size);

    FILE *log = 0; //for keeping any output desired on selection
    if( printPass ) {
        size_t pos = outfile.find(".root");
        assert( pos != string::npos );
        std::string outcpy = outfile;
        log = fopen( outcpy.replace(pos, 5, "_run_lumi_event").c_str(), "w" );
    }

    // set up on-the-fly residual JEC
    //std::vector<std::string> jetcorr_pf_filenames;
    //jetcorr_pf_filenames.push_back(Form("%s/jetcorr/data/GR_R_52_V7_L1FastJet_AK5PF.txt"   , jetcorr_path.empty() ? "." : jetcorr_path.c_str()));
    //jetcorr_pf_filenames.push_back(Form("%s/jetcorr/data/GR_R_52_V7_L2Relative_AK5PF.txt"  , jetcorr_path.empty() ? "." : jetcorr_path.c_str()));
    //jetcorr_pf_filenames.push_back(Form("%s/jetcorr/data/GR_R_52_V7_L3Absolute_AK5PF.txt"  , jetcorr_path.empty() ? "." : jetcorr_path.c_str()));
    //jetcorr_pf_filenames.push_back(Form("%s/jetcorr/data/GR_R_52_V7_L2L3Residual_AK5PF.txt", jetcorr_path.empty() ? "." : jetcorr_path.c_str()));
    //FactorizedJetCorrector* jet_pf_corrector = makeJetCorrector(jetcorr_pf_filenames); 

    TChain *chain = new TChain("Events");
    chain->Add(infile.c_str());
    if (chain->GetListOfBranches()->Contains("EventAuxiliary"))
    {
        chain->SetBranchStatus("EventAuxiliary", 0);
    }
    TObjArray *listOfFiles = chain->GetListOfFiles();
    const uint64 nEventsChain = chain->GetEntries();
    uint64 nEventsTotal = 0;
    uint64 nEventsSelected = 0;

    // file loop
    TIter fileIter(listOfFiles);
    TFile *currentFile = 0;
    bool first = true;
    int i_permille_old = 0;
    while ( (currentFile = (TFile*)fileIter.Next()) ) {
        TFile* f = TFile::Open(currentFile->GetTitle());
        TTree* tree = (TTree*)f->Get("Events");
    	TTreeCache::SetLearnEntries(10);
    	tree->SetCacheSize(128*1024*1024);
        tree->SetDirectory(f);

        // for the first file, clone the tree
        if ( first ) {
            output->cd();
            newtree = chain->CloneTree(0);
            newtree->SetDirectory(output);
            first = false;
        }

        // init
        cms2.Init(newtree);
        cms2.Init(tree);

        //f->cd();
        cout << "processing: " << f->GetName() << endl;

        // Event Loop
        const unsigned int nEvents = tree->GetEntries();
        for (unsigned int event = 0; event < nEvents; ++event, ++nEventsTotal)
        {
            int i_permille = (int)floor(10000 * nEventsTotal / float(nEventsChain));
            if (i_permille != i_permille_old) {
                // xterm magic from L. Vacavant and A. Cerri
                if (isatty(1))
                {
                    printf("\015\033[32m ---> \033[1m\033[31m%5.2f%%"
                            "\033[0m\033[32m <---\033[0m\015", i_permille/100.);
                    fflush(stdout);
                }
                i_permille_old = i_permille;
            }

            cms2.GetEntry(event);

            //cout << Form("run %d, ls %d, event %d", evt_run(), evt_lumiBlock(), evt_event()) << endl;

            //set condition to skip event
            if (not selectSameSign(do20_20, btag2, jets2)) 
            {
                continue;
            }

            ++nEventsSelected;
            if( printPass ) {
                fprintf(log, "%i %i %i\n", cms2.evt_run(), cms2.evt_lumiBlock(), cms2.evt_event());
            }

            // fill the new tree
            cms2.LoadAllBranches();
            output->cd();
            newtree->Fill();
        }

		delete f;
    }

    if (printPass)
    {
        fprintf(log, "\nTotal events run on: %llu\n", nEventsTotal);
        fprintf(log, "Num events selected: %llu\n", nEventsSelected ); //need two fprintf statements bc of some gcc bug
    }

    output->cd();
    newtree->Write();
    output->Close();
    delete output;
}

