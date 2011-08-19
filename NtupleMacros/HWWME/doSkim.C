#include "TFile.h"
#include <iostream>
#include <fstream>


int doSkim() {

  gROOT->ProcessLine(".L SkimChain.C+");
  
  using namespace std;
  
  ProcessSample("data/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Spring11-PU_S1_START311_V1G1-v1_SingleLepton/V04-01-01/merged_ntuple_100.root","data/", "data/", "PartonSkim");

  return 1;
}
