
//
// Dave "the one but not the only" Evans 
//

#include "MyScanChain.h"
#include "../CORE/CMS2.h"

// ROOT includes
#include "TChain.h"
#include "TChainElement.h"
#include "TTreeCache.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"

#include <cmath>
#include <cassert>

// LHAPDF
#include "/tas/dlevans/HWW2012/CMSSW_5_2_3/src/LHAPDF-5.8.92b/include/LHAPDF/LHAPDF.h"

MyScanChain::MyScanChain(std::string genPdfName, unsigned int genPdfSubset)
{                   
    // gen set parameters
    genPdfName_     = genPdfName;
    genPdfSubset_   = genPdfSubset;

    // set up LHAPDF
    LHAPDF::setPDFPath("/tas/dlevans/HWW2012/CMSSW_5_2_3/src/LHAPDF-5.8.92b/PDFSets/");
    LHAPDF::initPDFSetM(genset_, genPdfName_);
    LHAPDF::initPDFM(genset_, genPdfSubset_);

}   

//
// Main function
//

int MyScanChain::ScanChain(std::string sampleName, TChain *chain, std::string pdfName)
{

    TObjArray *listOfFiles = chain->GetListOfFiles();
    if (listOfFiles->GetEntries() == 0) {
        std::cout << "[MyScanChain::ScanChain] " << sampleName << " is not defined" << std::endl;
        return 1;
    }
    else {
        std::cout << "[MyScanChain::ScanChain] " << sampleName << std::endl;
    }

    TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
    if (rootdir == 0){
        std::cout<<"Head directory root: not found. Try Rint: ..."<<std::endl;
        if (rootdir){
            std::cout<<"OK: Got Rint:"<<std::endl;
        } else {
            std::cout<<"ERROR: no root: or Rint: found. Histograms will likely be lost"<<std::endl;
        }
    } 
    rootdir->cd(); 

    //
    // setup pdf stuff
    //

    // cteq6ll is only available in LHpdf format
    if (pdfName != "cteq6ll")   LHAPDF::initPDFSetM(set_, pdfName + ".LHgrid");
    else                        LHAPDF::initPDFSetM(set_, pdfName + ".LHpdf");

    unsigned int nsets = 1 + LHAPDF::numberPDFM(set_);
    if (pdfName == "NNPDF20_as_0116_100" || pdfName == "NNPDF20_as_0122_100") nsets = 5;
    if (pdfName == "NNPDF20_as_0117_100" || pdfName == "NNPDF20_as_0121_100") nsets = 27;
    if (pdfName == "NNPDF20_as_0118_100" || pdfName == "NNPDF20_as_0120_100") nsets = 72;
    if (pdfName == "cteq6mE") nsets = 1;
    if (pdfName == "cteq6ll") nsets = 1;

    // are we calculating the central value of 
    // the observable? 
    // - e.g. no pdf re-weighting needed
    bool doingGenSet = false;
    if (pdfName == genPdfName_) doingGenSet = true;

    //
    // setup histograms
    //

    Int_t nbins = 40;
    Float_t min = 0.0;
    Float_t max = 200.0;
    std::vector <TH1F*> histArr;
    for (unsigned int i = 0; i < nsets; ++i) {
        histArr.push_back(new TH1F(Form("%s_%s_%i", sampleName.c_str(), pdfName.c_str(), i),
                Form("%s_%s_%i", sampleName.c_str(), pdfName.c_str(), i), nbins, min, max));
    }

    //
    // branch variables
    //

    float   x1          = 0.0;      // momentum fraction for parton1
    float   x2          = 0.0;
    int     id1         = 0;        // pdgid of parton1
    int     id2         = 0;
    float   Q           = 0.0;      // event momentum scale
    float   scale1fb    = 0.0;      // event weight
    float   var         = 0.0;      // variable to study pdf uncertainty wrt

    //
    // loop over pdf subsets
    //

    for (unsigned int subset = 0; subset < nsets; ++subset)
    {

        std::cout << "doing set, subset: " << set_ << ", " << subset << std::endl;
        LHAPDF::initPDFM(set_, subset);

        //
        // loop over content of sample
        //

        TIter fileIter(listOfFiles);
        while (TChainElement *currentFile = (TChainElement*)fileIter.Next()) {

            // get the tree
            TFile *f = TFile::Open(currentFile->GetTitle()); 
            assert(f);
            TTree *tree = (TTree*)f->Get("Events");
            assert(tree);
            TTreeCache::SetLearnEntries(10);
            tree->SetCacheSize(128*1024*1024);
            cms2.Init(tree);

            //
            // loop over events in file
            //

            ULong64_t nEvents = tree->GetEntries();
            for(ULong64_t event = 0; event < nEvents; ++event) {

                cms2.GetEntry(event);
                x1          = cms2.pdfinfo_x1();
                x2          = cms2.pdfinfo_x2();
                id1         = cms2.pdfinfo_id1();
                id2         = cms2.pdfinfo_id2();
                Q           = cms2.pdfinfo_scale();
                scale1fb    = cms2.evt_scale1fb();
                var         = cms2.evt_pfmet();

                // check if event passes kinematic selection
                if (!Cuts()) continue;

                // calculate the pdf weight
                double pdf_weight = 1.0;
                double experimental_weight = scale1fb;

                // if looper has been invoked to calculate the same pdf the sample was
                // generated with, then we know the pdf_weight will always be 1.0
                // so no need to actually do any work
                // e.g. this will represent the "central value" of the observable

                if (!doingGenSet) {
                    // generated pdf values
                    double fx1Q0gen = LHAPDF::xfxM(genset_, x1, Q, id1) / x1;
                    double fx2Q0gen = LHAPDF::xfxM(genset_, x2, Q, id2) / x2;
                    // subset pdf values
                    double fx1Qi = LHAPDF::xfxM(set_, x1, Q, id1) / x1;
                    double fx2Qi = LHAPDF::xfxM(set_, x2, Q, id2) / x2;
                    // calculate weight and fill histogram
                    pdf_weight = ((fx1Qi*fx2Qi)/(fx1Q0gen*fx2Q0gen));
                }

                // inclusive uncertainty, one fixed bin
                // but could equally easily fill with a physical variable...
                histArr[subset]->Fill(var, pdf_weight * experimental_weight);

            } // end loop on events

            delete tree;
            f->Close();
            delete f;

        } // end loop on files in chain

    } // end loop on subsets

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

// test event selection
bool MyScanChain::Cuts()
{
    return true;
}

