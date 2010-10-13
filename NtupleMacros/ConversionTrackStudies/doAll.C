#include <vector>

void doAll() {
  
  //TString datapath = "/data/tmp/cms2/MinimumBias_BeamCommissioning09-Dec9thReReco-v1/V03-00-19/123734/*ntuple*.root";
  //TString datapath = "/data/tmp/cms2/ExpressPhysics_BeamCommissioning09-Express-v2/V03-00-20/124230/*ntuple*.root";
  //TString datapath1 = "/data/tmp/cms2/MinimumBias_BeamCommissioning09-Dec14thReReco-v1/V03-00-23/all_900GeV/*ntuple*.root";
  //TString datapath2 = "/data/tmp/cms2/MinimumBias_BeamCommissioning09-PromptReco-v2/V03-00-23/124230_124275/*ntuple*.root";
  //TString datapath3 = "/data/tmp/cms2/MinimumBias_BeamCommissioning09-Dec14thReReco-v1/V03-00-23/124120/*ntuple*.root";
  //TString datapath1 = "/data/tmp/cms2/MinimumBias_BeamCommissioning09-Dec14thReReco-v1/V03-00-23/all_900GeV/skimGood/skim*.root";
  //TString datapath2 = "/data/tmp/cms2/MinimumBias_BeamCommissioning09-PromptReco-v2/V03-00-23/124230_124275/skimGood/skim*.root";
  //TString datapath3 = "/data/tmp/cms2/MinimumBias_BeamCommissioning09-Dec14thReReco-v1/V03-00-23/124120/skimGood/skim*.root";

  //fkw for the met pas

  //old rereco but new version of ntuple maker
  TString datapath1 = "/data/tmp/cms2/MinimumBias_BeamCommissioning09-Dec14thReReco-v1/V03-00-35/merged*.root"; 
  //new rereco new version of ntuple maker
  TString datapath2 = "/data/tmp/cms2/MinimumBias_BeamCommissioning09-SD_AllMinBias-Jan23Skim-v1/V03-00-35/merged*.root";
  //new mc new version of ntuplemaker
  TString datapath3 = "/data/tmp/cms2/MinBias_Summer09-STARTUP3X_V8P_900GeV-v1/merged_ntuple_1*.root";
  TString datapath4 = "/data/tmp/cms2/MinBias_Summer09-STARTUP3X_V8P_900GeV-v1/merged_ntuple_1*.root";

  gROOT->ProcessLine(".L ScanChain.C++");
  gROOT->ProcessLine(".L histtools.C++");
  //  gROOT->ProcessLine(".L makePlots.C++");
  TChain *ch_data = new TChain("Events");

  //run on Dec. 14th data
  ch_data->Add(datapath1.Data());

  //run on Jan 23rd data
  //fkw ch_data->Add(datapath2.Data());

  TChain *ch_mc = new TChain("Events");
  ch_mc->Add(datapath3.Data());
  ch_mc->Add(datapath4.Data());
  
  
  std::vector<unsigned int> goodRuns;
  //goodRuns.push_back(123734);
  //goodRuns.push_back(123591);
  //goodRuns.push_back(123592);
  goodRuns.push_back(123596);
  goodRuns.push_back(123615);
  goodRuns.push_back(123732);
  goodRuns.push_back(123815);
  goodRuns.push_back(123818);
  //goodRuns.push_back(123906);
  goodRuns.push_back(123908);
  //goodRuns.push_back(123909);
  //goodRuns.push_back(123970);
  //goodRuns.push_back(123977);
  //goodRuns.push_back(123978);
  //goodRuns.push_back(123985);
  //goodRuns.push_back(123987);
  goodRuns.push_back(124008);
  goodRuns.push_back(124009);
  goodRuns.push_back(124020);
  goodRuns.push_back(124022);
  goodRuns.push_back(124023);
  goodRuns.push_back(124024);
  goodRuns.push_back(124025);
  goodRuns.push_back(124027);
  goodRuns.push_back(124030);
  goodRuns.push_back(124120);
  //goodRuns.push_back(124230);
  //goodRuns.push_back(124275);

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
  //saveHist("Data_"+outfname);
  saveHist("Test900GeVGoodRun_"+outfname);
  //  deleteHistos(); 
  //  makePlots("Run123734_"+outfname, description);
  
}  
  
