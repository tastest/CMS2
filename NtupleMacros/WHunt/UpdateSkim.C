#include <iostream>

#include "TObjArray.h"
#include "TObjString.h"
#include "TPRegexp.h"
#include "TROOT.h"
#include "TString.h"

#include "CORE/CMS2.cc"
#include "emuskim.cc"
#include "babymaker.C"

void UpdateSkim(const char *inputFileName) 
{
    TPRegexp preg("\\S+/merged_ntuple_(\\d+_\\d+).root");
    TString ident = ((TObjString*)preg.MatchS(TString(inputFileName))->At(1))->GetString();

    TString skimFileName = "";
    skimFileName.Append("skim/emuskim_");
    skimFileName.Append(ident);
    skimFileName.Append(".root");

    emuskim(inputFileName, skimFileName.Data());

    TString babyFileName = skimFileName;
    babyFileName.ReplaceAll("skim/emuskim_", "baby/emuskim_baby_");

    std::cout << "Making a baby named " << babyFileName << std::endl;
    babymaker *baby = new babymaker();
    baby->ScanChain(skimFileName.Data(), babyFileName.Data());
}
