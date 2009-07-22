{

  using namespace std;

  // Load various tools  
  gROOT->ProcessLine(".x setup.C(true)");
  gROOT->LoadMacro("ScanChain.C+");

  TChain *chain = new TChain("Events");
//  chain->Add("/store/disk02/gutsche/cms2-V1-02-06/SingleElectronFlatPt5To100/merged_ntuple*.root"); // cms-tas03
  chain->Add("/store/disk01/data/gutsche/cms2-V1-02-06/SingleElectronFlatPt5To100/merged_ntuple*.root"); // cms-tas01
//   chain->Add("/store/disk01/data/cms2-V01-03-01/WZ_incl_Summer08_IDEAL_V11_redigi_v1/*.root"); // cms-tas01
//   chain->Add("/store/disk01/data/cms2-V01-03-01/WW_Summer08_IDEAL_V11_redigi_v1/*.root"); // cms-tas01
//   chain->Add("/store/disk01/data/cms2-V01-03-01/ZZ_Summer08_IDEAL_V11_redigi_v1/*.root"); // cms-tas01
//   chain->Add("/store/disk01/data/cms2-V01-03-01/WJets-madgraph_Summer08_IDEAL_V11_redigi_v1/*.root"); // cms-tas01
//   chain->Add("/store/disk01/data/cms2-V01-03-01/ZJets-madgraph_Summer08_IDEAL_V11_redigi_v1/*.root"); // cms-tas01
//   chain->Add("/store/disk01/data/cms2-V01-03-01/TTJets-madgraph_Fall08_IDEAL_V11_redigi_v10/*.root"); // cms-tas01
//   chain->Add("/store/disk01/data/cms2-V01-03-01/SingleTop_tWChannel_Summer08_IDEAL_V11_redigi_v3/*.root"); // cms-tas01


  ScanChain(chain);

  //save all the histograms
    
  const char* outFile = "myHist.root";
  hist::saveHist(outFile);
  hist::deleteHistos();
}
