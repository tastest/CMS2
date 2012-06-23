#include <assert.h>
#include <string>
#include "TChain.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TTree.h"
#include "TDatabasePDG.h"

#include "CMS2.cc"
#include "CORE/utilities.cc"
#include "CORE/electronSelections.cc"
#include "CORE/electronSelectionsParameters.cc"
#include "CORE/MITConversionUtilities.cc"
#include "CORE/muonSelections.cc"
#include "CORE/eventSelections.cc"
#include "CORE/ttbarSelections.cc"
#include "CORE/trackSelections.cc"

#include "Rtypes.h"
typedef ULong64_t uint64;

using namespace tas;

// Used to filter an ntuple based on the 'select' function return
// infile is the path to the ntuple you want to filter (* allowed)
// outfile is the result
// cut is specified in the select() function

// eg:
// ntupleFilter("/data/tmp/cms2-V03-00-10/MinimumBias_BeamCommissioning09-rereco_FIRSTCOLL_v1/*.root","/data/tmp/wandrews/minbias/filtered_ntuple.root")
//infile and outfile must end in ".root"

// WARNING: all of the input files are put into 1 output file, so
// please be careful you don't create a file which is enormous (unless
// you can handle enormous files)


bool select (bool isData)
{

  //------------------------------------
  // check for electron T&P pair
  //------------------------------------

  // loop on tags
  for (unsigned int tag = 0; tag < cms2.els_p4().size(); ++tag) {

    if( !pass_electronSelection( tag , electronSelection_ssV5 , false , false ) ) continue; // SS ID/iso
    if( cms2.els_p4()[tag].Pt() < 20.)                                            continue; // pT > 20 GeV
    if( fabs(cms2.els_etaSC()[tag]) > 2.5)                                        continue; // |eta| < 2.5
    
    // loop on probes
    for (unsigned int probe = 0; probe < cms2.els_p4().size(); ++probe) {

      // check probe and tag do not overlap
      if (tag == probe) continue;

      // basic probe denominator
      if( cms2.els_p4()[probe].Pt() < 10.)      continue; // pT > 10 GeV
      if( fabs(cms2.els_etaSC()[probe]) > 2.5)  continue; // |eta| < 2.5

      float dilmass = (cms2.els_p4()[probe] + cms2.els_p4()[tag]).M();

      // return true if T&P mass is in Z window
      if( dilmass > 76 && dilmass < 106 ) return true;
    }
  }

  //------------------------------------
  // check for muon T&P pair
  //------------------------------------

  // loop on tags
  for (unsigned int tag = 0; tag < cms2.mus_p4().size(); ++tag) {

    if( !muonId( tag , OSGeneric_v3 )   )        continue; // OS ID/iso    
    if( cms2.mus_p4()[tag].Pt() < 20.0)          continue; // pT > 20 GeV
    if( fabs(cms2.mus_p4()[tag].Eta()) > 2.4)    continue; // |eta| < 2.4

    // loop on probes
    for (unsigned int probe = 0; probe < cms2.mus_p4().size(); ++probe) {

      // check probe and tag do not overlap
      if (tag == probe) continue;

      // basic probe denominator
      if (cms2.mus_p4()[probe].Pt() < 10.0)                continue; // pt cut
      if (fabs(cms2.mus_p4()[probe].Eta()) > 2.4)          continue; // eta cut

      float dilmass = (cms2.mus_p4()[probe] + cms2.mus_p4()[tag]).M();

      // return true if T&P mass is in Z window
      if( dilmass > 76 && dilmass < 106 ) return true;
    }
  }

  return false;
}

void ntupleFilter_StopTNPSkim(const std::string &infile, const std::string &outfile, bool printPass=false, bool isData = true, std::string runlist = "")  
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


