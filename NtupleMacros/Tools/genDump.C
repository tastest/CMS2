//------------------------------------------------------------------
// Dumps the GEN from the ntuple to the screen in interactive root
//
// Usage:
// root> .L genDump.C++
// root> TChain * ch = new TChain("Events")
// root> ch->Add("file.root")
// root> genDump(chain, first, nev)
//   will dump a number of events (nev) starting from event "first"
//
//-----------------------------------------------------------------
#include <iostream>
#include <vector>

#include "TChain.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TROOT.h"

#include "../CORE/CMS2.cc"
#include "../CORE/mcSelections.cc"
using namespace tas;

int genDump( TChain* chain, int first, int nev) {

  TObjArray *listOfFiles = chain->GetListOfFiles();

  //  unsigned int nEventsChain=0;
  int nEventsChain = chain->GetEntries();
  // nEventsChain = nEvents;
  unsigned int nEventsTotal = 0;
  int nDumped = 0;
  int nprocessed=0;

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
      // Progress feedback to the user
      if(nEventsTotal%1000 == 0) {
        // xterm magic from L. Vacavant and A. Cerri
        if (isatty(1)) {
          printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
          "\033[0m\033[32m <---\033[0m\015", (float)nEventsTotal/(nEventsChain*0.01));
          fflush(stdout);
        }
      }//if(nEventsTotal%20000 == 0) {
      nprocessed++;
      if (nprocessed < first) continue;
      cout << "  " <<endl;
      cout << " Dumping entry " << nprocessed-1 << endl;
      dumpDocLines();
      nDumped++;
      if (nDumped >= nev) return 0;

    }
  }

  return 0;
}
