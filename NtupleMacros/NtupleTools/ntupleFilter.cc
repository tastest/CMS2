#include <string>
#include <assert.h>
#include "TFile.h"
#include "TChain.h"

//Used to filter an ntuple based on a branch
//infile is the path to the ntuple you want to filter (* allowed)
//outfile is the result
//cut is specified in the same way as you would in Events->Draw(...)

//eg:
//ntupleFilter("/data/tmp/cms2-V03-00-10/MinimumBias_BeamCommissioning09-rereco_FIRSTCOLL_v1/*.root","/data/tmp/wandrews/minbias/filtered_ntuple.root", "(l1_techbits2 & (1<<8)) || (l1_techbits2 & (1<<9))")

//WARNING: all of the input files are put into 1 output file, so please be careful you don't create a file which is enormous

void ntupleFilter( string infile, string outfile, string cut="" ) {
  
  long long max_tree_size = 20000000000000000LL;
  TTree::SetMaxTreeSize(max_tree_size);

  TChain *chain = new TChain("Events");
  chain->Add(infile.c_str());

  TFile *file = TFile::Open(outfile.c_str(), "RECREATE");

  assert( file != 0 );

  TTree* tree = chain->CopyTree( cut.c_str() );
  //TDirectory *old_gdir = gDirectory;

  //file->cd(); //go to the directory that actually contains the tree
  tree->Write();

  //now close whatever we can get our hands on
  // gDirectory = old_gdir;
  file->Close();
  //if (file != gDirectory)
  // 	gDirectory->Close();
}
