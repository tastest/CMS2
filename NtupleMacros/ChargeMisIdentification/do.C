{

  using namespace std;

  // Load and compile something to allow proper treatment of vectors
  // Not clear that it is needed
  gROOT->LoadMacro("loader.C+");

  // Load various tools  
  gROOT->ProcessLine(".x setup.C(false)");

  gROOT->LoadMacro("ScanChain.C+");

  TChain *chain = new TChain("Events");
  chain->Add("/uscms_data/d1/gutsche/ntuple/cms2-V01-02-06/SingleElectronFlatPt5To100/merged_ntuple*.root");
  ScanChain(chain);

  //save all the histograms
    
  const char* outFile = "myHist.root";
  hist::saveHist(outFile);
  hist::deleteHistos();
}
