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
#include <map>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <set>

#include "TChain.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TROOT.h"

#include "CMS2.cc"

using namespace tas;

struct DorkyEventIdentifier {
// this is a workaround for not having unique event id's in MC
	unsigned long int run, event,lumi;
	bool operator < (const DorkyEventIdentifier &) const;
	bool operator == (const DorkyEventIdentifier &) const;
};

bool DorkyEventIdentifier::operator < (const DorkyEventIdentifier &other) const
{
	if (run != other.run)
		return run < other.run;
	if (event != other.event)
		return event < other.event;
	if(lumi != other.lumi)
		return lumi < other.lumi;
	return false;
}

bool DorkyEventIdentifier::operator == (const DorkyEventIdentifier &other) const
{
	if (run != other.run)
		return false;
	if (event != other.event)
		return false;
	return true;
}

std::set<DorkyEventIdentifier> already_seen;
bool is_duplicate (const DorkyEventIdentifier &id) {
	std::pair<std::set<DorkyEventIdentifier>::const_iterator, bool> ret =
		already_seen.insert(id);
	return !ret.second;
}

int ScanChain( TChain* chain, int nEvents = -1) {

	TObjArray *listOfFiles = chain->GetListOfFiles();

	// map to hold lumi per run
	std::map<int,std::vector<int> > runs;

	unsigned int nEventsChain=0;
	if(nEvents==-1) 
		nEvents = chain->GetEntries();
	nEventsChain = nEvents;
	unsigned int nEventsTotal = 0;
	// TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
// TH1F *samplehisto = new TH1F("samplehisto", "Example histogram", 200,0,200);
// samplehisto->SetDirectory(rootdir);
// file loop
	TIter fileIter(listOfFiles);
	TFile *currentFile = 0;
	while ( ( currentFile = (TFile*)fileIter.Next() ) ) {
		TFile f(currentFile->GetTitle());
		TTree *tree = (TTree*)f.Get("Events");
		cms2.Init(tree);

		//Event Loop
		unsigned int nEvents = tree->GetEntries();
		for( unsigned int event = 0; event < nEvents; ++event) {
			cms2.GetEntry(event);

			DorkyEventIdentifier id = { evt_run(),evt_event(), evt_lumiBlock() };
			if (is_duplicate(id) )
				continue;
			++nEventsTotal;
	// Progress feedback to the user
			if(nEventsTotal%1000 == 0) {
		// xterm magic from L. Vacavant and A. Cerri
				if (isatty(1)) {
					printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
						"\033[0m\033[32m <---\033[0m\015", (float)nEventsTotal/(nEventsChain*0.01));
					fflush(stdout);
				}
			}

			int run = evt_run();
			int lumi = evt_lumiBlock();

			if ( runs.find(run) == runs.end() ) {
				std::vector<int> tmp;
				tmp.push_back(lumi);
				runs[run] = tmp;
			} else {
				if ( std::find(runs.find(run)->second.begin(),runs.find(run)->second.end(),lumi) == runs.find(run)->second.end() ) {
					runs.find(run)->second.push_back(lumi);
				}
			}
		}
		delete tree;
		f.Close();
	}

	if ( nEventsChain != nEventsTotal ) {
		std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
	}

	// print runs/lumi
	std::ostringstream result("");
	result << "{";
	bool first = true;
	for ( std::map<int,std::vector<int> >::const_iterator run = runs.begin(); run != runs.end(); ++run ) {
		if ( first ) {
			result << "\"" << run->first << "\":[";
			first = false;
		} else {
			result << ",\"" << run->first << "\":[";			
		}
		std::vector<int> lumis = run->second;
		std::sort(lumis.begin(),lumis.end());
		result << "[" << lumis[0] << ",";
		int upperBound = lumis[0];
		for ( unsigned int counter = 1; counter < lumis.size(); ++counter ) {
			if ( upperBound + 1 == lumis[counter] ) {
				upperBound = lumis[counter];
			} else {
				result << upperBound << "],[" << lumis[counter] << ",";
				upperBound = lumis[counter];
			}
		}
		result << upperBound << "]]";
	}
	result << "}";
	// std::cout << result.str() << std::endl;
	std::ofstream output("runsNlumis.json");
	output << result.str() << std::endl;
	output.close();

// samplehisto->Draw();
	return 0;
}
