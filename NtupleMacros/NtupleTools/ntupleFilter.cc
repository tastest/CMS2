// $Id: ntupleFilter.cc,v 1.2 2009/12/06 04:10:56 jmuelmen Exp $

#include <assert.h>
#include <string>
#include "TChain.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TTree.h"

#include "../CORE/CMS2.h"
CMS2 cms2;

#include "Rtypes.h"
typedef ULong64_t uint64;

// Used to filter an ntuple based on a selection expression
// infile is the path to the ntuple you want to filter (* allowed)
// outfile is the result
// cut is specified in the select() function

// eg:
// ntupleFilter("/data/tmp/cms2-V03-00-10/MinimumBias_BeamCommissioning09-rereco_FIRSTCOLL_v1/*.root","/data/tmp/wandrews/minbias/filtered_ntuple.root")

// WARNING: all of the input files are put into 1 output file, so
// please be careful you don't create a file which is enormous (unless
// you can handle enormous files)

bool select ()
{
     // for example, require l1 tech bit 40 or 41
     return cms2.l1_techbits2() & (1 << 8) || cms2.l1_techbits2() & (1 << 9);
}

void ntupleFilter (const std::string &infile, const std::string &outfile)  
{
     // output file and tree
     TFile *output =TFile::Open(outfile.c_str(), "RECREATE");
     assert(output != 0);
     TTree *newtree = 0;

     const long long max_tree_size = 20000000000000000LL;
     TTree::SetMaxTreeSize(max_tree_size);

     TChain *chain = new TChain("Events");
     chain->Add(infile.c_str());
     TObjArray *listOfFiles = chain->GetListOfFiles();
     const uint64 nEventsChain = chain->GetEntries();
     uint64 nEventsTotal = 0;

     // file loop
     TIter fileIter(listOfFiles);
     TFile *currentFile = 0;
     bool first = true;
     int i_permille_old = 0;
     while ( currentFile = (TFile*)fileIter.Next() ) {
	  TFile f(currentFile->GetTitle());
	  const char *name = f.GetName();
	  TTree *tree = (TTree*)f.Get("Events");
	  
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
	       if (not select()) 
		    continue;

	       cms2.LoadAllBranches();
	       
	       // fill the new tree
	       newtree->Fill();
	  }
     }

     output->cd();
     newtree->Write();
     delete output;
}
