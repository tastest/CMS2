#ifndef LEPTONTREELOOPER_H
#define LEPTONTREELOOPER_H

#include "TChain.h"
#include "TString.h"

#include <iostream>

class LeptonTreeLooper {

    public:
        LeptonTreeLooper();
        ~LeptonTreeLooper();

        void setGoodRunList(const char *runlist); 
        void loop(TChain *chain, TString name);

    private:

    bool runlistIsSet_;

};

#endif
