
#ifndef PDFANALYSIS_H
#define PDFANALYSIS_H

// C++ includes
#include <iostream>
#include <vector>

#include "TROOT.h"

class TH1F;
class TChain;
class TDirectory;
class TTree;

#define MAXWEIGHT 110

// indices in LHAPDF for the genset and the 
// alternate set
static const unsigned int set_ = 2;
static const unsigned int genset_ = 1;

class MyScanChain {

    public:

        MyScanChain() {};
        MyScanChain(std::string genPdfName, unsigned int genPdfSubset);
        ~MyScanChain() {};

        int ScanChain(std::string sampleName, TChain *chain, std::string pdfName);

    private:

        // test event selection
        bool Cuts();

        // variables to read
        // from tree
        float           scale1fb_;
        float           Q_;
        float           x1_, x2_;
        int             id1_, id2_;

        // pdf parameters
        std::string     genPdfName_;
        unsigned int    genPdfSubset_;

};

#endif

