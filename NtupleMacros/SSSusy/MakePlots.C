#include <iostream>
void MakePlots(TString histofilename){
  gROOT->SetStyle("Plain");
  char * descr = getenv("CMS2_LOCATION");
  TString command=".x ";
  command.Append(descr);
  command.Append("NtupleMacros/Tools/setup.C");
  gROOT->ProcessLine(command);
  loadHist(histofilename);
  histofilename+="_out";
  browseStacks(true, false, histofilename);
}
