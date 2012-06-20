#ifndef LEPTONTREELOOPER_H
#define LEPTONTREELOOPER_H

#include "TChain.h"
#include "TString.h"

#include <iostream>

#include "Enums.h"
class LeptonTree;

class LeptonTreeLooper {

    public:
        LeptonTreeLooper();
        ~LeptonTreeLooper();

        void setGoodRunList(const char *runlist); 
        void unsetGoodRunList();
        void loop(TChain *chain, TString name, unsigned int plotBin);

    private:

        bool runlistIsSet_;
        bool testPlotBin(const unsigned int plotBin, const LeptonTree *tree);

bool passElectronFO2012(const LeptonTree *leptonTree);
bool passElectronIso2012(const LeptonTree *leptonTree);

};

#endif
