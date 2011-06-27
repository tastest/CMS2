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

    dilepskim(inputFileName, skimFileName.Data(), 1);

    TString babyFileName = skimFileName;
    babyFileName.ReplaceAll("skim/dilepskim_", "baby/dilepskim_baby_");

    std::cout << "Making a baby named " << babyFileName << std::endl;
    twinmaker *twin = new twinmaker();
    twin->ScanChain(skimFileName.Data(), babyFileName.Data());
}
