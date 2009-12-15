/* Usage:
   root [0] .L ScanChain.C++
   root [1] TFile *_file0 = TFile::Open("merged_ntuple.root")
   root [2] TChain *chain = new TChain("Events")
   root [3] chain->Add("merged_ntuple.root")

   There are several places where one may create CMS2 cms2
   It can be done here (in a doAll.C script), i.e.:

   root [4] CMS2 cms2 

   It can be done in the source as is done below, or it can be
   ascertained by including CORE/CMS2.cc as is commented out
   below.  They are all the same, and everything will work so
   long as it is created somewhere globally.

   root [5] ScanChain(chain)
 */
#include <iostream>
#include <vector>

#include "TChain.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "../CORE/CMS2.h"
#include "myselections.cc"

CMS2 cms2;

/*
#include "CORE/CMS2.cc"
#include "CORE/selections.cc"
#include "CORE/utilities.cc"
 */

using namespace tas;

int ScanChain(bool isData, TChain* chain, int nEvents = -1, std::string skimFilePrefix="") {

	TObjArray *listOfFiles = chain->GetListOfFiles();

	unsigned int nEventsChain=0;
	if(nEvents==-1) 
		nEvents = chain->GetEntries();
	nEventsChain = nEvents;
	unsigned int nEventsTotal = 0;
	TDirectory *rootdir = gDirectory->GetDirectory("Rint:");

	TH1F *h1_l1_techbits2_pass = new TH1F("h1_l1_techbits2_pass", "l1_techbits2", 32, -0.5, 31.5);

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

			// how many good tracks
			//
			int nGoodTracksVz = 0;
			int nHighPurityTracks = 0;
			for (size_t t = 0; t < cms2.trks_ndof().size(); ++t) {

				// count high purity tracks
				if (isTrackQuality(t, (1<<highPurity))) nHighPurityTracks ++;

				// cut on track ndof
				if (cms2.trks_ndof()[t] < 8) continue;
				if (cms2.vtxs_position().size() > 0)
					if (fabs(cms2.vtxs_position()[0].z() - cms2.trks_vertex_p4()[t].z()) < 10.0) continue;
					else if (fabs(cms2.trks_vertex_p4()[t].z()) < 10.0) continue;
					nGoodTracksVz ++;
			}

			// count good vertexs
			//
			int nGoodVertex = 0;
			for (size_t v = 0; v < cms2.vtxs_position().size(); ++v) {

				if (cms2.vtxs_isFake()[v]) continue;
				if (fabs(cms2.vtxs_position()[v].z()) > 10.0
					|| fabs(cms2.vtxs_position()[v].x()) > 0.50
                                        || fabs(cms2.vtxs_position()[v].y()) > 0.50) continue;
				nGoodVertex ++;
			}


			// require bit 40 or 41 passed
			//
			if (!(cms2.l1_techbits2() & (1<<8) || cms2.l1_techbits2() & (1<<9))) continue;

			// require bits 36-39 DIDN't pass ???
			//
			if (cms2.l1_techbits2() & (1<<7) || cms2.l1_techbits2() & (1<<6) || 
				cms2.l1_techbits2() & (1<<5) || cms2.l1_techbits2() & (1<<4)) continue;

			// require bit zero for beams 
			//
                        if (isData && !(cms2.l1_techbits1() & (1<<0))) continue;

			// pixel digi requirement to get rid of monster events
			//
			if (nHighPurityTracks < 10 || cms2.trks_ndof().size() > 100) continue;
			if (cms2.trks_ndof().size()/float(nHighPurityTracks) < 0.20) continue;
			

			// fill histograms
			//



		} // end loop on files
	} // end loop on events

	if ( nEventsChain != nEventsTotal ) {
		std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
	}

	return 0;
}
