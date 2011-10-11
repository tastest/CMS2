#include <assert.h>
#include <string>
#include "TChain.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TTree.h"

#include "CMS2.cc"

#include "CORE/utilities.cc"
#include "CORE/trackSelections.cc"
#include "CORE/eventSelections.cc"
#include "CORE/MITConversionUtilities.cc"
#include "CORE/electronSelectionsParameters.cc"
#include "CORE/electronSelections.cc"
#include "CORE/muonSelections.cc"
#include "CORE/mcSelections.cc"


#include "Rtypes.h"
typedef ULong64_t uint64;

using namespace tas;

// Used to filter an ntuple based on the 'select' function return
// cut is specified in the select() function

// bool electronpass_VBTF90_NOHOEEND(const unsigned int index, bool applyAlignementCorrection = false, bool removedEtaCutInEndcap = false)
// {
// 		// VBTF90 w/o no H/E in endcap
// 		double sigmaIEtaIEtaThresholds[2]       = {0.01 ,    0.03   };
// 		double dPhiInThresholds[2]              = {0.8  ,     0.7   };
// 		double dEtaInThresholds[2]              = {0.007,   0.009   };
// 		double hoeThresholds[2]                 = {0.12 ,    9999.  };

// 		//
// 	    // get corrected dEtaIn and dPhiIn
// 	    //
// 		float dEtaIn = cms2.els_dEtaIn()[index];
// 		float dPhiIn = cms2.els_dPhiIn()[index];
// 		if (applyAlignementCorrection) electronCorrection_pos(index, dEtaIn, dPhiIn);

// 		// barrel
// 		if (fabs(cms2.els_etaSC()[index]) < 1.479) 
// 		{
// 				if( fabs(dEtaIn) > dEtaInThresholds[0]) return false; 
// 				if(	fabs(dPhiIn) > dPhiInThresholds[0]) return false;
// 				if(	cms2.els_hOverE()[index] > hoeThresholds[0]) return false; 
// 			    if(	cms2.els_sigmaIEtaIEta()[index] > sigmaIEtaIEtaThresholds[0]) return false;

// 				return true;
// 		}

// 		// endcap 
// 		if (fabs(cms2.els_etaSC()[index]) > 1.479) 
// 		{
// 				bool passdEtaCut = fabs(dEtaIn) < dEtaInThresholds[1];
// 				if(removedEtaCutInEndcap) passdEtaCut = true;
		
// 				if( !passdEtaCut ) return false;	
// 				if( fabs(dEtaIn) > dEtaInThresholds[1]) return false; 
// 				if(	fabs(dPhiIn) > dPhiInThresholds[1]) return false;
// 				if(	cms2.els_hOverE()[index] > hoeThresholds[1]) return false; 
// 			    if(	cms2.els_sigmaIEtaIEta()[index] > sigmaIEtaIEtaThresholds[1]) return false;

// 				return true;
// 		}
// }

bool select (bool isData)
{
  int nel=0, nmu=0, ntau=0; //dummy variables since the function requires it
  if(!isData && leptonGenpCount(nel, nmu, ntau)>=2){return true;}
  if(cms2.hyp_p4().size()>=1){return true;}
  // //if contains 3 or more status==3 gen leptons, accept it
  // int nel=0, nmu=0, ntau=0; //dummy variables since the function requires it
  // if(!isData && leptonGenpCount(nel, nmu, ntau)>=3){return true;}


  // 		vector<int> v_idxgoodele;
  // 		vector<int> v_idxgoodmu;

  // 		// good electrons
  // 		for(unsigned int iels=0; iels < els_p4().size(); iels++)
  // 		{
  // 				if( abs(els_p4().at(iels).eta()) > 2.5 ) continue;
  // 				if( els_p4().at(iels).pt() < 10 ) continue;
  // 				if( electronIsolation_rel_v1(iels, true ) > 1.0 ) continue;	// NT RELISO < 1.0
  // 				if( !electronpass_VBTF90_NOHOEEND(iels) ) continue; // VBTF90 w/o H/E in endcap 

  // 				v_idxgoodele.push_back(iels);
  // 		}


  // 		// good muons 
  // 		for(unsigned int imus=0; imus < mus_p4().size(); imus++)
  // 		{
  // 				if( abs(mus_p4().at(imus).eta()) > 2.5 ) continue;
  // 				if( mus_p4().at(imus).pt() < 10 ) continue;
  // 				if( muonIsoValue(imus, false) > 1.0 ) continue;
  // 				if (((cms2.mus_type().at(imus)) & (1<<1)) == 0) continue; // global muon
  // 				if (((cms2.mus_type().at(imus)) & (1<<2)) == 0) continue; // tracker muon
				
  // 				v_idxgoodmu.push_back(imus);
  // 		}						
 

  // 		// more than 3 good leptons
  // 	    if( v_idxgoodele.size() + v_idxgoodmu.size() < 3) return false;


  // 		// max pt > 20 GeV ?
  // 		float maxpt = 0;
  // 		for(unsigned int ivels = 0; ivels<v_idxgoodele.size(); ivels++)
  // 				if( cms2.els_p4().at(v_idxgoodele[ivels]).pt() > maxpt ) maxpt	=	cms2.els_p4().at(v_idxgoodele[ivels]).pt();
  // 		for(unsigned int ivmus = 0; ivmus<v_idxgoodmu.size(); ivmus++)
  // 				if( cms2.mus_p4().at(v_idxgoodmu[ivmus]).pt() > maxpt ) maxpt	=	cms2.mus_p4().at(v_idxgoodmu[ivmus]).pt();
		
  // 		if( maxpt < 20 ) return false;

		return false;

}

void ntupleFilterDilepHyp(const std::string &infile, const std::string &outfile, bool printPass=false, bool isData = true, std::string runlist = "")  
{
    // set good run list
    //if (runlist != "") set_goodrun_file(runlist.c_str());


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

    TChain *chain = new TChain("Events");
    chain->Add(infile.c_str());
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
        TFile f(currentFile->GetTitle());
        //const char *name = f.GetName();
        TTree *tree = (TTree*)f.Get("Events");

        //chain->SetBranchStatus(chain->GetAlias("evt_scale1fb"     ), 0);
        //chain->SetBranchStatus(chain->GetAlias("evt_xsec_excl"    ), 0);
        //chain->SetBranchStatus(chain->GetAlias("evt_xsec_incl"    ), 0);
        //chain->SetBranchStatus(chain->GetAlias("evt_kfactor"      ), 0);
        //chain->SetBranchStatus(chain->GetAlias("evt_nEvts"        ), 0);
        //chain->SetBranchStatus(chain->GetAlias("evt_filt_eff"     ), 0);
        chain->SetBranchStatus("EventAuxiliary",0);
        // for the first file, clone the tree
        if ( first ) {
            newtree = chain->CloneTree(0);
            newtree->SetDirectory(output);
            first = false;

        }

        // init
        cms2.Init(newtree);
        cms2.Init(tree);


        // Event Loop
        const unsigned int nEvents = tree->GetEntries();
        for (unsigned int event = 0; event < nEvents; ++event, ++nEventsTotal) {
            int i_permille = (int)floor(10000 * nEventsTotal / float(nEventsChain));
            if (i_permille != i_permille_old) {
                // xterm magic from L. Vacavant and A. Cerri
                if (isatty(1)) {
                    printf("\015\033[32m ---> \033[1m\033[31m%5.2f%%"
                            "\033[0m\033[32m <---\033[0m\015", i_permille/100.);
                    fflush(stdout);
                }
                i_permille_old = i_permille;
            }

            cms2.GetEntry(event);
            //set condition to skip event
            if (not select(isData)) 
                continue;

            ++nEventsSelected;
            if( printPass ) {
                fprintf(log, "%i %i %i\n", cms2.evt_run(), cms2.evt_lumiBlock(), cms2.evt_event() );
            }

            cms2.LoadAllBranches();

            // fill the new tree
            newtree->Fill();
        }
    }

    if( printPass ) {
        fprintf(log, "\nTotal events run on: %llu\n", nEventsTotal);
        fprintf(log, "Num events selected: %llu\n", nEventsSelected ); //need two fprintf statements bc of some gcc bug
    }

    output->cd();
    newtree->Write();
    delete output;
}
