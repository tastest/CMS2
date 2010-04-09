// #include "TTimeStamp.h"
// #include "TString.h"
// #include <iostream>
//void UpdateSkim(void) 
{
  TTimeStamp * time = new TTimeStamp();
  std::cout<<"hello"<<std::endl;
  std::cout<<"And the date is: "<<time->GetDate()<<std::endl;
  std::cout<<"And the time is: "<<time->GetTime()<<std::endl;

  TString outputFileName = "";
  outputFileName.Append("skim/WHunt_update_");
  outputFileName+= time->GetDate();
  outputFileName.Append("_");
  outputFileName+= time->GetTime();
  outputFileName.Append(".root");

  gROOT->ProcessLine(".L emuskim.cc+");
  gROOT->ProcessLine("emuskim(\"RunsToProcess.txt\",outputFileName.Data())");

  TString babyoutputFileName = outputFileName;
  babyoutputFileName.ReplaceAll("skim/WHunt_update_", "baby/WHunt_update_baby_");

  TChain *c = new TChain("Events");
  c->Add(outputFileName);
  gROOT->ProcessLine(".L babymaker.C+");
  babymaker bah;
  bah.ScanChain(c, babyoutputFileName.Data());
}
