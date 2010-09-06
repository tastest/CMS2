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

    TString skimFileName = "";
    skimFileName.Append(TString(sample));
    skimFileName.Append("/emu_skim/emuskim_");
    skimFileName.Append(ident);
    skimFileName.Append(".root");

    nlepskim(sample, inputFileName, skimFileName.Data(), 1, true);

    TString babyFileName = skimFileName;
    babyFileName.ReplaceAll("emu_skim/emuskim_", "emu_baby/emuskim_baby_");

    std::cout << "Making a baby named " << babyFileName << std::endl;
    emubabymaker *emubaby = new emubabymaker();
    emubaby->ScanChain(skimFileName.Data(), babyFileName.Data());

    // dilep

    skimFileName = "";
    skimFileName.Append(TString(sample));
    skimFileName.Append("/dilep_skim/dilepskim_");
    skimFileName.Append(ident);
    skimFileName.Append(".root");

    nlepskim(sample, inputFileName, skimFileName.Data(), 2, false);

    babyFileName = skimFileName;
    babyFileName.ReplaceAll("dilep_skim/dilepskim_", "dilep_baby/dilepskim_baby_");

    std::cout << "Making a baby named " << babyFileName << std::endl;
    dilepbabymaker *dilepbaby = new dilepbabymaker();
    dilepbaby->ScanChain(skimFileName.Data(), babyFileName.Data());

    // trilep

    skimFileName = "";
    skimFileName.Append(TString(sample));
    skimFileName.Append("/trilep_skim/trilepskim_");
    skimFileName.Append(ident);
    skimFileName.Append(".root");

    nlepskim(sample, inputFileName, skimFileName.Data(), 3, false);

    babyFileName = skimFileName;
    babyFileName.ReplaceAll("trilep_skim/trilepskim_", "trilep_baby/trilepskim_baby_");

    std::cout << "Making a baby named " << babyFileName << std::endl;
    trilepbabymaker *trilepbaby = new trilepbabymaker();
    trilepbaby->ScanChain(skimFileName.Data(), babyFileName.Data());
}
