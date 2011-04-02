#include <iostream>

#include "TObjArray.h"
#include "TObjString.h"
#include "TPRegexp.h"
#include "TROOT.h"
#include "TString.h"

#include "CORE/CMS2.cc"
#include "CORE/electronSelections.cc"
#include "CORE/electronSelectionsParameters.cc"
#include "CORE/metSelections.cc"
#include "CORE/muonSelections.cc"
#include "CORE/trackSelections.cc"
#include "CORE/jetSelections.cc"
#include "CORE/mcSelections.cc"

#include "nlepskim.C"
#include "babymakercommon.C"
#include "emubabymaker.C"
#include "dilepbabymaker.C"
#include "trilepbabymaker.C"

void UpdateSkims(const char *sample, const char *inputFileName) 
{
    TPRegexp preg("\\S+/merged_ntuple_(\\d+).root");
    TString ident;
    if (preg.MatchB(TString(inputFileName)))
        ident = ((TObjString*)preg.MatchS(TString(inputFileName))->At(1))->GetString();
    else
        ident = "0";

    // emu

    TString babyFileName = Form("%s/emu_baby/emuskim_baby_%s.root", sample, ident.Data());

    std::cout << "Making a baby named " << babyFileName << std::endl;
    emubabymaker *emubaby = new emubabymaker();
    emubaby->ScanChain(inputFileName, babyFileName.Data());

    // dilep

    babyFileName = Form("%s/dilep_baby/dilepskim_baby_%s.root", sample, ident.Data());

    std::cout << "Making a baby named " << babyFileName << std::endl;
    dilepbabymaker *dilepbaby = new dilepbabymaker();
    dilepbaby->ScanChain(inputFileName, babyFileName.Data());

    // trilep

    babyFileName = Form("%s/trilep_baby/trilepskim_baby_%s.root", sample, ident.Data());

    std::cout << "Making a baby named " << babyFileName << std::endl;
    trilepbabymaker *trilepbaby = new trilepbabymaker();
    trilepbaby->ScanChain(inputFileName, babyFileName.Data());
}
