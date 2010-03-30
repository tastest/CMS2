
//
// ttbar -> ll
// Dave "the one but not the only" Evans 
//

#include "MyScanChain.h"

// ROOT includes
#include "TChain.h"
#include "TChainElement.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"

#include "Math/LorentzVector.h"

// CMS2 includes
#include "CMS2.h"
#include "../../CORE/electronSelections.h"
#include "../../CORE/eventSelections.h"
#include "../../CORE/jetSelections.h"
#include "../../Tools/DileptonHypType.h"

//
// Namespaces
//
using namespace tas;

//
//
//

double dRbetweenVectors(const LorentzVector &vec1, const LorentzVector &vec2 ){

	double dphi = std::min(::fabs(vec1.Phi() - vec2.Phi()), 2 * M_PI - fabs(vec1.Phi() - vec2.Phi()));
	double deta = vec1.Eta() - vec2.Eta();
	return sqrt(dphi*dphi + deta*deta);
} 

// to decdie if to fill the EB histo (zero in the array)
// or the EE histo (one in the array)
enum DetType MyScanChain::getSubdet(int eleIndex)
{
	// determine what detector the electron is in
	if (cms2.els_fiduciality()[eleIndex] & (1<<ISEB)) return DET_EB;
	else if (cms2.els_fiduciality()[eleIndex] & (1<<ISEE)) return DET_EE;
	std::cout << "ERROR! not in EB or EE" << std::endl;
	return DET_ALL;
}

void MyScanChain::Format2DHist(TH2F** hist, std::string name, Int_t nx, Float_t minx, Float_t maxx,
		Int_t ny, Float_t miny, Float_t maxy)
{
	// loop on EB, EE
	for (unsigned int i = 0; i < 3; ++i)
	{
		std::string det = det_names[i];
		hist[i] = new TH2F(Form("%s_%s_%s", sampleName_.c_str(), name.c_str(), det.c_str()),
				name.c_str(), nx, minx, maxx, ny, miny, maxy);
	}
}

void MyScanChain::FormatHist(TH1F** hist, std::string name, Int_t n, Float_t min, Float_t max)
{
	// loop on EB, EE
	for (unsigned int i = 0; i < 3; ++i)
	{
		std::string det = det_names[i];
		hist[i] = new TH1F(Form("%s_%s_%s", sampleName_.c_str(), name.c_str(), det.c_str()),
				name.c_str(), n, min, max);
		//hist[i]->SetFillColor(sampleName_.histo_color);
		hist[i]->GetXaxis()->SetTitle(name.c_str());
	}
}

bool MyScanChain::CheckCutsNM1(cuts_t apply, cuts_t remove, cuts_t passed)
{           
    if ((passed & (apply & (~remove))) == (apply & (~remove))) return true;
    return false;
}   
        
bool MyScanChain::CheckCuts(cuts_t apply, cuts_t passed)
{
    if ((apply & passed) == apply) return true;
    return false;
}


//
// Main function
//
int MyScanChain::ScanChain(bool isData, std::string sampleName, TChain *chain, int nEvents, std::string skimFilePrefix) {

	// set sampleName
	sampleName_ = sampleName;

	//
	//
	//
	TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
	if (rootdir == 0){
		std::cout<<"Head directory root: not found. Try Rint: ..."<<std::endl;
		rootdir = gROOT->GetDirectory("Rint:");
		if (rootdir){
			std::cout<<"OK: Got Rint:"<<std::endl;
		} else {
			std::cout<<"ERROR: no root: or Rint: found. Histograms will likely be lost"<<std::endl;
		}
	} 

	//
	// format histograms
	//
	// N-1
	FormatHist(h1_mass_, "mass", 200, 0, 200);
    FormatHist(h1_met_, "met", 400, 0, 400);
    FormatHist(h1_metin_, "metin", 400, 0, 400);

    FormatHist(h1_sumpt_, "sumpt", 1000, 0, 1000);
    FormatHist(h1_sumptin_, "sumptin", 1000, 0, 1000);


	// file loop
	//

	unsigned int nEventsChain=0;
	if(nEvents == -1) nEvents = chain->GetEntries();
	nEventsChain = nEvents;
	unsigned int nEventsTotal = 0;
	int i_permille_old = 0;

	TObjArray *listOfFiles = chain->GetListOfFiles();
	TIter fileIter(listOfFiles);
	TFile *currentFile = 0;
	while (TChainElement *currentFile = (TChainElement*)fileIter.Next()) {
		TFile *f = TFile::Open(currentFile->GetTitle());
		TTree *tree = (TTree*)f->Get("Events");
		cms2.Init(tree);

		//Event Loop
		ULong64_t nEvents = tree->GetEntries();
		for(ULong64_t event = 0; event < nEvents; ++event) {
			cms2.GetEntry(event);
			++nEventsTotal;

			// Progress feedback to the user
			int i_permille = (int)floor(1000 * nEventsTotal / float(nEventsChain));
			if (i_permille != i_permille_old) {
				// xterm magic from L. Vacavant and A. Cerri
				if (isatty(1)) {
					printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
							"\033[0m\033[32m <---\033[0m\015", i_permille/10.);
					fflush(stdout);
				}
				i_permille_old = i_permille;
			}

			// work out event weight
			//float weight = cms2.evt_scale1fb()*0.1;

            //
            // event level cuts
            //


            //
            // loop on hyps
            //


            for (size_t h = 0; h < cms2.hyp_p4().size(); ++h) {

                    //
                    // hyp level cuts
                    //

            }


		} // end loop on files

	} // end loop on events

	if ( nEventsChain != nEventsTotal ) {
		std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
	}

    rootdir = gROOT->GetDirectory("root:");
    if (rootdir) rootdir->cd();
    else{
        std::cout<<"Cant find root: . Current dir is "<<gDirectory->GetName()<<std::endl;
        rootdir = gROOT->GetDirectory("Rint:");
        if (rootdir){
            std::cout<<"OK, got Rint: "<<std::endl;
            rootdir->cd();
        } else {
            std::cout<<"Cant find Rint: either . Current dir is "<<gDirectory->GetName()<<std::endl;
        }
    }

	return 0;
}

