#ifndef textFileMaker_h
#define textFileMaker_h

#include "TFile.h"
#include "TTree.h"

//class TChain;

class textFileMaker
{
    public:
        textFileMaker() {};
        ~textFileMaker() {};
        void ScanChain (TChain *);

	//private:

};

#endif
