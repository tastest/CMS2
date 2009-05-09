#include <iostream>
#include <vector>

#include "TChain.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TEventList.h"
#include "TBranch.h"

#include "CMS2.h"
CMS2 cms2;
#include "../CORE/selections.cc"
#include "../CORE/utilities.cc"

using namespace tas;

int ScanChain( TChain* chain, const char* outputname) {

  // output file and tree
  TFile *output = new TFile(outputname,"RECREATE");
  TTree *newtree = 0;

  // additional branches: declaration
  TBranch *anotherBranch = 0;

  // additional variables to be filled into an additional branch into the new tree
  int anotherVariable;
  std::vector<int> *anotherVector = new std::vector<int>;

  // list of files
  TObjArray *listOfFiles = chain->GetListOfFiles();

  // events to process
  unsigned int nEventsChain = chain->GetEntries();
  unsigned int nEventsTotal = 0;

  // file loop
  TIter fileIter(listOfFiles);
  TFile *currentFile = 0;
  bool first = true;
  while ( currentFile = (TFile*)fileIter.Next() ) {
    TFile f(currentFile->GetTitle());
    TTree *tree = (TTree*)f.Get("Events");

    // for the first file, clone the tree
    if ( first ) {
      newtree = chain->CloneTree(0);
      newtree->SetDirectory(output);
      first = false;

      // book branches
      anotherBranch = newtree->Branch("anotherbranch",&anotherVariable, "another_branch/I");
      newtree->SetAlias("another_branch","anotherbranch");

      anotherBranch = newtree->Branch("anothervector",&anotherVector);
      newtree->SetAlias("another_vector","anothervector");

    }

    // init
    cms2.Init(newtree);
    cms2.Init(tree);

    // Event Loop
    unsigned int nEvents = tree->GetEntries();
    for( unsigned int event = 0; event < nEvents; ++event) {
      ++nEventsTotal;
      if ( nEventsTotal%10000 == 0 ) {
	cout << "Event: " << nEventsTotal << endl;
      }

      // reset additional variables
      anotherVariable = 0;
      anotherVector->clear();

      cms2.GetEntry(event);
      cms2.LoadAllBranches();

      // cut on njets >= 4
      if ( evt_njets() < 4 ) continue;

      // fill additional variables
      anotherVariable = evt_njets()-2;
      anotherVector->push_back(evt_njets()-4);
      anotherVector->push_back(evt_njets()-3);
      anotherVector->push_back(evt_njets()-2);
      anotherVector->push_back(evt_njets()-1);

      // fill the new tree
      newtree->Fill();

    }
  }

  output->cd();
  newtree->Write();
  delete output;

  if ( nEventsChain != nEventsTotal ) {
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  }
  
  return 0;
}

