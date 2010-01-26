
//
// ttbar -> ll
// Dave "the one but not the only" Evans 
//

// C++ includes
#include <iostream>
#include <vector>

// ROOT includes
#include "TChain.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TROOT.h"

// CMS2 includes
#include "CMS2.h"
#include "../CORE/electronSelections.h"

//
// CMS2 
//
CMS2 cms2;

//
// Namespaces
//
using namespace tas;

//
// Main function
//
int ScanChain(bool isData, TChain* chain, int nEvents = -1, std::string skimFilePrefix="") {

	TObjArray *listOfFiles = chain->GetListOfFiles();

	unsigned int nEventsChain=0;
	if(nEvents==-1) 
		nEvents = chain->GetEntries();
	nEventsChain = nEvents;
	unsigned int nEventsTotal = 0;
	TDirectory *rootdir = gDirectory->GetDirectory("Rint:");

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

			std::cout << "Event: " << event << std::endl;


		} // end loop on files
	} // end loop on events

	if ( nEventsChain != nEventsTotal ) {
		std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
	}

	return 0;
}

