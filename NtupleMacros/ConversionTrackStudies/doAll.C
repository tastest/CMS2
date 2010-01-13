#include <vector>

void doAll() {
  
  //TString datapath = "/data/tmp/cms2/MinimumBias_BeamCommissioning09-Dec9thReReco-v1/V03-00-19/123734/*ntuple*.root";
  //TString datapath = "/data/tmp/cms2/ExpressPhysics_BeamCommissioning09-Express-v2/V03-00-20/124230/*ntuple*.root";
  //TString datapath1 = "/data/tmp/cms2/MinimumBias_BeamCommissioning09-Dec14thReReco-v1/V03-00-23/all_900GeV/*ntuple*.root";
  //TString datapath2 = "/data/tmp/cms2/MinimumBias_BeamCommissioning09-PromptReco-v2/V03-00-23/124230_124275/*ntuple*.root";
  //TString datapath3 = "/data/tmp/cms2/MinimumBias_BeamCommissioning09-Dec14thReReco-v1/V03-00-23/124120/*ntuple*.root";
  TString datapath1 = "/data/tmp/cms2/MinimumBias_BeamCommissioning09-Dec14thReReco-v1/V03-00-23/all_900GeV/skimGood/skim*.root";
  TString datapath2 = "/data/tmp/cms2/MinimumBias_BeamCommissioning09-PromptReco-v2/V03-00-23/124230_124275/skimGood/skim*.root";
  TString datapath3 = "/data/tmp/cms2/MinimumBias_BeamCommissioning09-Dec14thReReco-v1/V03-00-23/124120/skimGood/skim*.root";

  gROOT->ProcessLine(".L ScanChain.C++");
  gROOT->ProcessLine(".L histtools.C++");
  //  gROOT->ProcessLine(".L makePlots.C++");
  TChain *ch_data = new TChain("Events");
  ch_data->Add(datapath1.Data());
  ch_data->Add(datapath2.Data());
  //fkw ch_data->Add(datapath3.Data());
  TChain *ch_mc = new TChain("Events");
  ch_mc->Add("/data/tmp/cms2-V03-00-19/MinBias_Summer09-STARTUP3X_V8I_900GeV-v2/merged_ntuple_1.root");
  ch_mc->Add("/data/tmp/cms2-V03-00-19/MinBias_Summer09-STARTUP3X_V8I_900GeV-v2/merged_ntuple_2.root");
  ch_mc->Add("/data/tmp/cms2-V03-00-19/MinBias_Summer09-STARTUP3X_V8I_900GeV-v2/merged_ntuple_3.root");
  ch_mc->Add("/data/tmp/cms2-V03-00-19/MinBias_Summer09-STARTUP3X_V8I_900GeV-v2/merged_ntuple_4.root");
  ch_mc->Add("/data/tmp/cms2-V03-00-19/MinBias_Summer09-STARTUP3X_V8I_900GeV-v2/merged_ntuple_5.root");
  ch_mc->Add("/data/tmp/cms2-V03-00-19/MinBias_Summer09-STARTUP3X_V8I_900GeV-v2/merged_ntuple_6.root");
  ch_mc->Add("/data/tmp/cms2-V03-00-19/MinBias_Summer09-STARTUP3X_V8I_900GeV-v2/merged_ntuple_7.root");
  ch_mc->Add("/data/tmp/cms2-V03-00-19/MinBias_Summer09-STARTUP3X_V8I_900GeV-v2/merged_ntuple_8.root");
  ch_mc->Add("/data/tmp/cms2-V03-00-19/MinBias_Summer09-STARTUP3X_V8I_900GeV-v2/merged_ntuple_9.root");
  ch_mc->Add("/data/tmp/cms2-V03-00-19/MinBias_Summer09-STARTUP3X_V8I_900GeV-v2/merged_ntuple_10.root");
  ch_mc->Add("/data/tmp/cms2-V03-00-19/MinBias_Summer09-STARTUP3X_V8I_900GeV-v2/merged_ntuple_11.root");
  ch_mc->Add("/data/tmp/cms2-V03-00-19/MinBias_Summer09-STARTUP3X_V8I_900GeV-v2/merged_ntuple_12.root");
  ch_mc->Add("/data/tmp/cms2-V03-00-19/MinBias_Summer09-STARTUP3X_V8I_900GeV-v2/merged_ntuple_13.root");
  ch_mc->Add("/data/tmp/cms2-V03-00-19/MinBias_Summer09-STARTUP3X_V8I_900GeV-v2/merged_ntuple_14.root");
  ch_mc->Add("/data/tmp/cms2-V03-00-19/MinBias_Summer09-STARTUP3X_V8I_900GeV-v2/merged_ntuple_15.root");
  ch_mc->Add("/data/tmp/cms2-V03-00-19/MinBias_Summer09-STARTUP3X_V8I_900GeV-v2/merged_ntuple_16.root");
  ch_mc->Add("/data/tmp/cms2-V03-00-19/MinBias_Summer09-STARTUP3X_V8I_900GeV-v2/merged_ntuple_17.root");
  ch_mc->Add("/data/tmp/cms2-V03-00-19/MinBias_Summer09-STARTUP3X_V8I_900GeV-v2/merged_ntuple_18.root");
  ch_mc->Add("/data/tmp/cms2-V03-00-19/MinBias_Summer09-STARTUP3X_V8I_900GeV-v2/merged_ntuple_19.root");
  ch_mc->Add("/data/tmp/cms2-V03-00-19/MinBias_Summer09-STARTUP3X_V8I_900GeV-v2/merged_ntuple_20.root");
  
  
  std::vector<unsigned int> goodRuns;
  //goodRuns.push_back(123734);
  goodRuns.push_back(123591);
  goodRuns.push_back(123592);
  goodRuns.push_back(123596);
  goodRuns.push_back(123615);
  goodRuns.push_back(123732);
  goodRuns.push_back(123815);
  goodRuns.push_back(123818);
  goodRuns.push_back(123906);
  goodRuns.push_back(123908);
  goodRuns.push_back(123909);
  goodRuns.push_back(123970);
  goodRuns.push_back(123977);
  goodRuns.push_back(123978);
  goodRuns.push_back(123985);
  goodRuns.push_back(123987);
  goodRuns.push_back(124009);
  goodRuns.push_back(124020);
  goodRuns.push_back(124022);
  goodRuns.push_back(124023);
  goodRuns.push_back(124024);
  goodRuns.push_back(124025);
  goodRuns.push_back(124027);
  goodRuns.push_back(124030);
  goodRuns.push_back(124120);
  goodRuns.push_back(124230);
  goodRuns.push_back(124275);

  bool runningonGEN = false;
  bool requireTrackCuts  = true;
  TString description;
  //description = ScanChain(ch_data, runningonGEN, requireTrackCuts, goodRuns);
  
  bool runningonGEN=true;
  ScanChain(ch_mc, runningonGEN, requireTrackCuts);

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
  saveHist("Test900GeVMC_"+outfname);
  //  deleteHistos(); 
  //  makePlots("Run123734_"+outfname, description);
  
}  
  
