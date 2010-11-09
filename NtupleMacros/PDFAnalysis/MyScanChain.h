
#ifndef PDFANALYSIS_H
#define PDFANALYSIS_H

// C++ includes
#include <iostream>
#include <vector>

#include "TROOT.h"

class TH1F;
class TChain;
class TDirectory;

#define MAXWEIGHT 101 

class MyScanChain {

    public:

        MyScanChain() {};
        ~MyScanChain() {};

        int ScanChain(std::string sampleName, TChain *chain, float kFactor = 1.0, int nEvents = -1, std::string skimFilePrefix="");
        void specifyPDF(std::string pdfName1, std::string pdfName2);

    private:

        float GetGenMeff();

        void Fill(TH1F** hist, const unsigned int hyp, const float &val, const float &weight);
        void FormatHist(TH1F** hist, std::string sampleName, std::string name, int n, float min, float max);
        std::string pdfName1_;
        std::string pdfName2_;

};

#endif

