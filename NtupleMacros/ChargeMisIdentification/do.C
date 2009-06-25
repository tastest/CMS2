{

  using namespace std;

  // Load various tools  
  gROOT->ProcessLine(".x setup.C(true)");
  gROOT->LoadMacro("ScanChain.C+");

  TChain *chain = new TChain("Events");
  chain->Add("/store/disk02/gutsche/cms2-V1-02-06/SingleElectronFlatPt5To100/merged_ntuple*.root");
  ScanChain(chain);

  //save all the histograms
    
  const char* outFile = "myHist.root";
  hist::saveHist(outFile);
  hist::deleteHistos();
}
