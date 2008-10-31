#include <assert.h>
#include <string>

#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TList.h"

int main (int argc, char **argv) 
{
//      TApplication app("crapo", argc, argv);
     long long max_tree_size = 200000000000LL;
     TTree::SetMaxTreeSize(max_tree_size);
     
     assert(argc >= 3);
     TList *list = new TList;
     TFile *file1 = TFile::Open(argv[2]);
     TTree *tree1 = dynamic_cast<TTree*>(file1->Get("Events"));
     for ( unsigned int counter = 3;
	   counter < argc;
	   ++counter ) {
       list->Add(dynamic_cast<TTree*>(TFile::Open(argv[counter])->Get("Events")));
     }
     TFile *output = new TFile(argv[1],"recreate");
     TTree *copy = tree1->CloneTree();
     if ( list->GetEntries() > 0 ) {
       copy->Merge(list);
     }
     copy->Write();
     output->Close();
     return 0;
}
