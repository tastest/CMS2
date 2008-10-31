#include <assert.h>
#include <string>

#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"

int main (int argc, char **argv) 
{
//      TApplication app("crapo", argc, argv);
     using std::string;
     assert(argc == 3);
     string infile = argv[1];
     string outfile = argv[2];
     TFile *f_in = TFile::Open(infile.c_str());
     TTree *events = dynamic_cast<TTree *>(f_in->Get("Events"));
     TFile *f_out = new TFile(outfile.c_str(), "recreate");
     TTree *t_copy = events->CloneTree();
     t_copy->Write();
     return 0;
}
