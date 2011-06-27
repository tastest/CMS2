{


  using namespace std;

  // Load and compile something to allow proper treatment of vectors
  // Not clear that it is needed
  gROOT->LoadMacro("loader.C+");

  // Load various tools  
  gROOT->ProcessLine(".x setup.C(true)");
  gROOT->LoadMacro("ScanChain.C+");

  // variables
  TChain *chain = 0;

  // QCD pt 30
  chain = new TChain("Events");
  chain->Add("/data/tmp/cms2-V01-02-06/QCDpt30/merged_ntuple*.root");
  ScanChain(chain,"qcd_pt_30.root");
  delete chain;

}
