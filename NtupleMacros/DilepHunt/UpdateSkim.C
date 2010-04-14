#include <iostream>

#include "TObjArray.h"
#include "TObjString.h"
#include "TPRegexp.h"
#include "TROOT.h"
#include "TString.h"

#include "CORE/CMS2.cc"
#include "dilepskim.cc"
#include "twinmaker.C"

void UpdateSkim(const char *inputFileName) 
{
    TPRegexp preg("\\S+/merged_ntuple_(\\d+_\\d+).root");
    TString ident = ((TObjString*)preg.MatchS(TString(inputFileName))->At(1))->GetString();
    TString skimFileName = "";
    skimFileName.Append("skim/dilepskim_");
    skimFileName.Append(ident);
    skimFileName.Append(".root");

    dilepskim(inputFileName, skimFileName.Data());

    TString twinFileName = skimFileName;
    twinFileName.ReplaceAll("skim/dilepskim_", "twin/dilepskim_twin_");

    std::cout << "Making twins named " << twinFileName << std::endl;
    twinmaker *twin = new twinmaker();
    twin->ScanChain(skimFileName.Data(), twinFileName.Data());
}
