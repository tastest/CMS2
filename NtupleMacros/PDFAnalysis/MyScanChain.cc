
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
#include <cmath>

// SMURF
#include "SmurfTree.h"

// LHAPDF
#include "/tas03/home/dlevans/LHAPDF/include/LHAPDF/LHAPDF.h"
//#include "/afs/cern.ch/cms/slc5_amd64_gcc434/external/lhapdf/5.8.5/full/include/LHAPDF/LHAPDF.h"

MyScanChain::MyScanChain() {
}

//
// Main function
//

int MyScanChain::ScanChain(std::string sampleName, std::string pdfName, const char *file) {

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
    // file loop
    //

    TFile f(file, "READ");
    TTree *t = (TTree*)f.Get("tree");
    t->SetBranchAddress("Q", &Q_);   
    t->SetBranchAddress("x1", &x1_);                    
    t->SetBranchAddress("id1", &id1f_);
    t->SetBranchAddress("x2", &x2_); 
    t->SetBranchAddress("id2", &id2f_);
    t->SetBranchAddress("scale1fb", &scale1fb_);
    t->SetBranchAddress(sampleName.c_str(), &bdt_);    
    rootdir->cd();

    const unsigned int genset = 1;
    const unsigned int set = 2;
    std::string pdfDir = "/tas03/home/dlevans/lhapdf-5.8.6b2/PDFsets/";
    //std::string pdfDir = "/afs/cern.ch/cms/slc5_amd64_gcc434/external/lhapdf/5.8.5/share/lhapdf/PDFsets/";
    LHAPDF::initPDFSetM(genset, pdfDir +"cteq6mE.LHgrid");
    LHAPDF::initPDFSetM(set, pdfDir + pdfName + ".LHgrid");

    unsigned int nsets = 1 + LHAPDF::numberPDF(set);

    // number of replicas to sample related to alpha_S for NNPDF 
    if (pdfName == "NNPDF20_as_0116_100" || pdfName == "NNPDF20_as_0122_100") nsets = 5;
    if (pdfName == "NNPDF20_as_0117_100" || pdfName == "NNPDF20_as_0121_100") nsets = 27;
    if (pdfName == "NNPDF20_as_0118_100" || pdfName == "NNPDF20_as_0120_100") nsets = 72;

    // only do the central value of CTEQ6M
    if (pdfName == "cteq6mE") nsets = 1;

    Int_t nbins = 40;
    Float_t min = -1;
    Float_t max = 1.0;
    std::vector <TH1F*> histArr;
    for (unsigned int i = 0; i < nsets; ++i) {
        histArr.push_back(new TH1F(Form("%s_%s_%i", sampleName.c_str(), pdfName.c_str(), i),
                Form("%s_%s_%i", sampleName.c_str(), pdfName.c_str(), i), nbins, min, max));
    }

    ULong64_t nEventsTree = t->GetEntries();
    for(ULong64_t event = 0; event < nEventsTree; ++event) {

        t->GetEntry(event);
        id1_ = int(id1f_);
        id2_ = int(id2f_);

        // get the generated central value
        LHAPDF::initPDFM(genset, 0);
        double fx1Q0gen = LHAPDF::xfxM(genset, x1_, Q_, id1_) / x1_;
        double fx2Q0gen = LHAPDF::xfxM(genset, x2_, Q_, id2_) / x2_;

        for (unsigned int subset = 0; subset < nsets; ++subset)
        {

            // get the weight for this subset
            LHAPDF::initPDFM(set, subset);
            double fx1Qi = LHAPDF::xfxM(set, x1_, Q_, id1_) / x1_;
            double fx2Qi = LHAPDF::xfxM(set, x2_, Q_, id2_) / x2_;

            double weight = ((fx1Qi*fx2Qi)/(fx1Q0gen*fx2Q0gen));

            // full the histogram
            histArr[subset]->Fill(bdt_, weight * scale1fb_);
        }

    } // end loop on tree

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

    f.Close();
    return 0;

}

