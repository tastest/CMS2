
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

class MyScanChain {

    public:

        MyScanChain();
        ~MyScanChain() {};

        int ScanChain(std::string sampleName, std::string pdfName, const char *file);

    private:

        // for the smurf tree
        float scale1fb_;
        float Q_;
        float x1_, x2_;
        float id1f_, id2f_;
        int id1_, id2_;
        float bdt_;

};

#endif

