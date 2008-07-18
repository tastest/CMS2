//now make the source file
#include <iostream>
#include <vector>

#include "TChain.h"
#include "TFile.h"

#include "TH1F.h"
#include "TH2F.h"

#include "CMS2.h"

int ScanChain( TChain* chain) {

  TObjArray *listOfFiles = chain->GetListOfFiles();

  unsigned int nEventsChain = chain->GetEntries();
  unsigned int nEventsTotal = 0;

  // file loop
  TIter fileIter(listOfFiles);
  TFile *currentFile = 0;
  while ( currentFile = (TFile*)fileIter.Next() ) {
    TFile f(currentFile->GetTitle());
    TTree *tree = (TTree*)f.Get("Events");
    Init(tree);
    
    //Event Loop
    unsigned int nEvents = tree->GetEntries();
    for( unsigned int event = 0; event < nEvents; ++event) {
      GetEntry(event);
      ++nEventsTotal;
      std::cout << "els size: " << els_p4.size() << " ";
      std::cout << "mus size: " << mus_p4.size() << std::endl;
      for (unsigned int hyp = 0;
           hyp < hyp_jets_p4.size();
           ++hyp) {
        std::cout << "hyp: " << hyp << "jet corrections:";
        for ( unsigned int jet = 0;
              jet < hyp_jets_p4[hyp].size();
              ++jet ) {
          std::cout << " " << hyp_jets_p4[hyp][jet].pt();
        }
        std::cout << endl;
      }
      if ( hyp_jets_p4.size() == 0 ) {
        std::cout << "no hypothesis!" << std::endl;
      }
    }
  }

  if ( nEventsChain != nEventsTotal ) {
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  }

  return 0;
}
