{

  using namespace std;

  // Load various tools  
  gROOT->ProcessLine(".x setup.C(true)");
  gROOT->LoadMacro("ScanChain.C+");

  TChain *chain = new TChain("Events");
  chain->Add("/nfs-3/userdata/cms2/SingleElectronFlatPt5To100_CMSSW_3_6_1_patch3/V03-04-25/merged_ntuple_1.root");
  //chain->Add("/nfs-3/userdata/cms2/SingleElectronFlatPt5To100_CMSSW_3_6_1_patch3/V03-04-25/*.root");


  ScanChain(chain);

  //save all the histograms
    
  const char* outFile = "myHist.root";
  hist::saveHist(outFile);
  hist::deleteHistos();
}
