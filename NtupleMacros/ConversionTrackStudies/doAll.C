#include <vector>

void doAll() {
  
  //TString datapath = "/data/tmp/cms2/MinimumBias_BeamCommissioning09-Dec9thReReco-v1/V03-00-19/123734/*ntuple*.root";
  //TString datapath = "/data/tmp/cms2/ExpressPhysics_BeamCommissioning09-Express-v2/V03-00-20/124230/*ntuple*.root";
  TString datapath = "/data/tmp/cms2/MinimumBias_BeamCommissioning09-PromptReco-v2/V03-00-23/124230_124275/*ntuple*.root";
  //TString datapath = "/data/tmp/cms2/MinimumBias_BeamCommissioning09-Dec14thReReco-v1/V03-00-23/124120/*ntuple*.root";

  gROOT->ProcessLine(".L ScanChain.C++");
  gROOT->ProcessLine(".L histtools.C++");
  //  gROOT->ProcessLine(".L makePlots.C++");
  TChain *ch_data = new TChain("Events");
  ch_data->Add(datapath.Data());
  TChain *ch_mc = new TChain("Events");
  ch_mc->Add("/data/tmp/cms2-V03-00-19/MinBias_Summer09-STARTUP3X_V8I_900GeV-v2/merged_ntuple_1.root");
  //  ch_mc->Add("/data/tmp/cms2-V03-00-19/MinBias_Summer09-STARTUP3X_V8I_900GeV-v2/merged_ntuple_2.root");
  
  
  std::vector<unsigned int> goodRuns;
  //goodRuns.push_back(123734);
  //goodRuns.push_back(124120);
  goodRuns.push_back(124230);


  bool runningonGEN = false;
  bool requireTrackCuts  = true;
  TString description;
  description = ScanChain(ch_data, runningonGEN, requireTrackCuts, goodRuns);
  
  bool runningonGEN=true;
  //ScanChain(ch_mc, runningonGEN, requireTrackCuts);

  TString outfname = "Runs";
  for(unsigned int i = 0; i < goodRuns.size(); i++)
    outfname = outfname + TString(Form("_%d", goodRuns.at(i)));
  if(requireTrackCuts)
    outfname = outfname + "_withTrackingCuts_withTriggers";
  else
    outfname = outfname + "_withTriggers";


  outfname = outfname + ".root";
  //saveHist("Run124120_"+outfname);
  saveHist("Run124230_"+outfname);
  //  deleteHistos(); 
  //  makePlots("Run123734_"+outfname, description);
  
}  
  
