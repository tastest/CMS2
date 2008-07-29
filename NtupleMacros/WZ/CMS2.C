//now make the source file
#include <iostream>
#include <vector>

#include "TChain.h"
#include "TFile.h"

#include "TH1F.h"
#include "TH2F.h"

#include "CMS2.h"

#include "../Tools/selections.h"
#include "../Tools/selections.C"
#include "../Tools/utilities.C"

int ScanChain( TChain* chain) {

  TObjArray *listOfFiles = chain->GetListOfFiles();

  unsigned int nEventsChain = chain->GetEntries();
  unsigned int nEventsTotal = 0;

  const unsigned int allBuckets = 20;
  
  // trilepton candidate per bucket
  unsigned int trilepCounter[allBuckets];
  for ( unsigned int i = 0; i < allBuckets; ++i ) {
    trilepCounter[i] = 0;
  }

  // CUTS
  const int   goodLeptonsCut        = 3;   // good lepton cut
  const float triggerLeptonMinPtCut = 20.; // one of the leptons of the trilepton canidate has to have at least this pt
  const float leptonMinPtCut        = 10.; // all leptons of the trilepton canidate have to have at least this pt

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

      // Progress feedback to the user
      if ((nEventsTotal)%1000 == 0) std::cout << "Processing event: " << nEventsTotal << std::endl;

      // loop over trilepton candidates
      for ( unsigned int cand = 0; 
	    cand < hyp_trilep_bucket.size();
	    ++cand ) {
	
	// CUT: goodLeptonIsolated
	if ( !goodLeptonIsolated(hyp_trilep_bucket[cand],
				 hyp_trilep_first_index[cand],
				 hyp_trilep_second_index[cand],
				 hyp_trilep_third_index[cand]) ) continue;
	
	++trilepCounter[hyp_trilep_bucket[cand]];
      }
    }
  }

  for ( unsigned int i = 0; i < allBuckets; ++i ) {
    cout << "Bucket: " << i << " entries: " << trilepCounter[i] << endl << endl;
  }

  if ( nEventsChain != nEventsTotal ) {
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  }

  return 0;
}
