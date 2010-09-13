
//
// Dave "the one but not the only" Evans 
//

#include "MyScanChain.h"

// ROOT includes
#include "TChain.h"
#include "TChainElement.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"

// CMS2 includes
#include "../CORE/CMS2.h"
#include "../Tools/DileptonHypType.h"
#include "../Tools/tools.cc"
#include "../CORE/utilities.h"

// LHAPDF
#include "/tas07/disk00/jribnik/lhapdf/include/LHAPDF/LHAPDF.h"

//
// Namespaces
//
using namespace tas;

//
// names of analysis cuts
//

//
// for jets binning
//
enum JetBin {
    JET_0J,
    JET_1J,
    JET_GT2J,
    JET_ANY,
};

static const char jetbin_names[][128] = { "0j", "1j", "2j", "allj"};

//
// functions
//

void MyScanChain::Fill(TH1F** hist, const unsigned int hyp, const float &val, const float &weight)
{
    hist[hyp]->Fill(val, weight);
    hist[DILEPTON_ALL]->Fill(val, weight);
}

void MyScanChain::FormatHist(TH1F** hist, std::string sampleName, std::string name, int n, float min, float max)
{       
    // loop on EB, EE
    for (unsigned int i = 0; i < 4; ++i)
    {
        std::string str = dilepton_hypo_names[i];
        std::string title = name + "_" + str;
        hist[i] = new TH1F(Form("%s_%s_%s", sampleName.c_str(), name.c_str(), str.c_str()),
                title.c_str(), n, min, max);
        hist[i]->GetXaxis()->SetTitle(name.c_str());
        hist[i]->Sumw2();
    }
}    

//
// set the pdf to use
//

void MyScanChain::specifyPDF(std::string pdfName) 
{
    pdfName_ = pdfName;
}

//
// Main function
//

int MyScanChain::ScanChain(std::string sampleName, TChain *chain, float kFactor, int nEvents, std::string skimFilePrefix) {

    std::cout << "scanning " << sampleName << std::endl;
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
    rootdir->cd(); 

    //
    // set up LHAPDF
    //

    // initialise the pdf 
    const int pdfSubset = 0;
    LHAPDF::initPDFSet(pdfName_, LHAPDF::LHGRID, pdfSubset);

    // set up array to hold the weights for each event
    // size is set to maximum possible
    double pdfWeights[MAXWEIGHT];
    unsigned int nsets = 1;
    if (LHAPDF::numberPDF() > 1) nsets += LHAPDF::numberPDF();

    //
    // format histograms
    // and set up quantities to study
    //

    double acceptance[MAXWEIGHT];
    double nTotal[MAXWEIGHT];
    double nPass[MAXWEIGHT];
    for (unsigned int i = 0; i < MAXWEIGHT; ++i) {
        acceptance[i] = 0.0;
        nTotal[i] = 0.0;
        nPass[i] = 0.0;
    }

    //
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
            float weight = cms2.evt_scale1fb() * kFactor;

            //
            // does this event pass the analysis selection
            //

            int nGoodLep = 0;
            for (size_t i = 0; i < cms2.genps_id().size(); ++i) 
            {
                if (!(abs(cms2.genps_id()[i]) == 11 || abs(cms2.genps_id()[i]) == 13)) continue;
                if(cms2.genps_p4()[i].Pt() < 20.0 || abs(cms2.genps_p4()[i].Eta()) > 2.5) continue;
                //if (abs(cms2.genps_id_mother()[i]) != 24) continue;
                nGoodLep++;
            }

            bool passSelection = false;
            if (nGoodLep == 2) passSelection = true;

            //
            // do PDF analysis
            //

            pdfWeights[0] = 1.0;

            // calculate the central value
            LHAPDF::initPDF(0);
            double fx1Q0 = LHAPDF::xfx(cms2.pdfinfo_x1(), cms2.pdfinfo_scale(), cms2.pdfinfo_id1()) / cms2.pdfinfo_x1();
            double fx2Q0 = LHAPDF::xfx(cms2.pdfinfo_x2(), cms2.pdfinfo_scale(), cms2.pdfinfo_id2()) / cms2.pdfinfo_x2();

            // calculate the weight for the ith subset
            // for this event
            for (unsigned int subset = 1; subset < nsets; ++subset) 
            {
                LHAPDF::initPDF(subset);
                double fx1Qi = LHAPDF::xfx(cms2.pdfinfo_x1(), cms2.pdfinfo_scale(), cms2.pdfinfo_id1()) / cms2.pdfinfo_x1();
                double fx2Qi = LHAPDF::xfx(cms2.pdfinfo_x2(), cms2.pdfinfo_scale(), cms2.pdfinfo_id2()) / cms2.pdfinfo_x2();
                pdfWeights[subset] = (fx1Qi*fx2Qi)/(fx1Q0*fx2Q0);
            }

            //
            // histogram the observable quantities
            // or record in some other way
            //
            for (unsigned int subset = 0; subset < nsets; ++subset) 
            {
                nTotal[subset] += pdfWeights[subset];
                if (passSelection) nPass[subset] += pdfWeights[subset];
            }

        } // end loop on files

    } // end loop on events

    if ( nEventsChain != nEventsTotal ) {
        std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
    }

    //
    // Analyse the observables
    //

    for (unsigned int i = 0; i < nsets; ++i) {
        acceptance[i] = nPass[i] / nTotal[i];
    } 

    double plus_max = 0.0;
    double minus_max = 0.0;
    double X0 = acceptance[0];

    for (unsigned int subset = 0; subset < ((nsets - 1)/2); ++subset)
    {
        double Xi_up = acceptance[(subset*2) + 1];
        double Xi_down = acceptance[(subset*2) + 2];
        plus_max += pow(max(max(Xi_up - X0, Xi_down - X0), 0.0), 2);
        minus_max += pow(max(max(X0 - Xi_up, X0 - Xi_down), 0.0), 2);
    }

    plus_max = sqrt(plus_max);
    minus_max = sqrt(minus_max);

    std::cout << "[MyScanChain::ScanChain] " << sampleName << std::endl;
    std::cout << "[MyScanChain::ScanChain] Analysing PDF uncertainty on the acceptance" << std::endl;
    std::cout << "[MyScanChain::ScanChain] Central value is         : " << X0 << "$^{+" << plus_max << "}_{-" << minus_max << "}" << std::endl;
    std::cout << "[MyScanChain::ScanChain] Relative uncertainty is  : $^{+" << plus_max/X0 << "}_{-" << minus_max/X0 << "}" << std::endl;

    //
    // make sure we're back in the right root dir
    //

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

