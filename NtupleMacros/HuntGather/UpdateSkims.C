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

void UpdateSkims(const char *inputFileName) 
{
   // TPRegexp preg("\\S+/merged_ntuple_(\\d+_\\d+).root");
  TPRegexp preg("\\S+/skimmed_ntuple_(\\d+_\\d+).root");
  
  TString ident = ((TObjString*)preg.MatchS(TString(inputFileName))->At(1))->GetString();

    // emu

    TString skimFileName = "";
    skimFileName.Append("emu_skim/emuskim_");
    skimFileName.Append(ident);
    skimFileName.Append(".root");
    nlepskim(inputFileName, skimFileName.Data(), 1, true);
    TString babyFileName = skimFileName;
    babyFileName.ReplaceAll("emu_skim/emuskim_", "emu_baby/emuskim_baby_");
    std::cout << "Making a baby named " << babyFileName << std::endl;
    emubabymaker *emubaby = new emubabymaker();
    emubaby->ScanChain(skimFileName.Data(), babyFileName.Data());
    
    // dilep

    skimFileName = "";
    skimFileName.Append("dilep_skim/dilepskim_");
    skimFileName.Append(ident);
    skimFileName.Append(".root");

    nlepskim(inputFileName, skimFileName.Data(), 2, false);

    babyFileName = skimFileName;
    babyFileName.ReplaceAll("dilep_skim/dilepskim_", "dilep_baby/dilepskim_baby_");

    std::cout << "Making a baby named " << babyFileName << std::endl;
    dilepbabymaker *dilepbaby = new dilepbabymaker();
    dilepbaby->ScanChain(skimFileName.Data(), babyFileName.Data());

    // trilep

    skimFileName = "";
    skimFileName.Append("trilep_skim/trilepskim_");
    skimFileName.Append(ident);
    skimFileName.Append(".root");

    nlepskim(inputFileName, skimFileName.Data(), 3, false);

    babyFileName = skimFileName;
    babyFileName.ReplaceAll("trilep_skim/trilepskim_", "trilep_baby/trilepskim_baby_");

    std::cout << "Making a baby named " << babyFileName << std::endl;
    trilepbabymaker *trilepbaby = new trilepbabymaker();
    trilepbaby->ScanChain(skimFileName.Data(), babyFileName.Data());
}
