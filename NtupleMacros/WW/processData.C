#if defined(__CINT__) && !defined(__MAKECINT__)
{
  gSystem->Load("libCMS2NtupleMacrosCORE");
  gSystem->Load("libCMS2NtupleMacrosLooper");
  gSystem->CompileMacro("processData.C","k");
  processData();
}
#endif 

#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TRegexp.h"
#include "TFile.h"

#ifndef __CINT__
 #include "wwtypes.h"
 #include "doAnalysis.h"
#endif

void processData()
{
  using namespace std;
  //
  // Output file
  //
  const char* outFile = "processed_data.root";

  //
  // Define how to handle complex Monte Carlo samples which may have 
  // more than one event type mixed in. Careful with this option. 
  // You don't need it for not mixed samples.
  //
  const bool identifyDYEvents = false;
  const bool zStudy = false;

  //
  // Flags for files to run over 
  // (0 and 1 are easier to modify)
  //
  bool runWW    = 1;
  bool runHWW   = 1;
  bool runGGWW  = 1;
  bool runWZ    = 1;
  bool runZZ    = 1;
  bool runWjets = 1;
  bool runDYee  = 1;
  bool runDYmm  = 1;
  bool runDYtt  = 1;
  bool runttbar = 1;
  bool runtW    = 1;
  bool runQCD   = 0; 
  bool runData  = 0;

  // 
  // Ntuple version
  //
  string version = "";
  if (gSystem->Getenv("VERSION")){
    version = gSystem->Getenv("VERSION");
    cout << "VERSION: " << version << endl;
  }

  //
  // ===================================================================================
  // 
  // Load various tools
  gROOT->ProcessLine(".x init.C");
  gSystem->Load("../Tools/MiniFWLite/libMiniFWLite");
  gROOT->LoadMacro("../Tools/getMyHistosNames.C");
  gROOT->LoadMacro("../Tools/histtools.C+");
  gROOT->LoadMacro("../Tools/browseStacks.C");
 
  // Define colors numbers:
  gStyle->SetPalette(1);
  enum EColor { kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan };

  // read dataset prefix
  string dataset = "data";
  if (gSystem->Getenv("DataDir")){
    dataset = gSystem->Getenv("DataDir");
    cout << "Dataset directory: " << dataset << endl;
  }
  
  const double integratedLumi = 1000.0; // pb^1  
  cout << "Integrated luminosity to scale to: " << integratedLumi << endl;

  if (runWW)
    ProcessSample(dataset+"/VVJetsTo4L_TuneD6T_7TeV-madgraph-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/"+version+"/*.root", WW, integratedLumi, 4.5*0.919, 963356*(682015./963356.), kRed, true);

  if (runGGWW)
    ProcessSample(dataset+"/GluGluToWWTo4L_TuneZ2_7TeV-gg2ww-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/*.root", ggWW, integratedLumi, 0.153, -1, kRed, true);

  if (runHWW){
    ProcessSample(dataset+"/GluGluToHToWWTo2L2Nu_M-120_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v2/V04-01-01/"+version+"/*.root", hWW120, integratedLumi, -1, -1, kBlue);
    ProcessSample(dataset+"/GluGluToHToWWTo2L2Nu_M-130_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/"+version+"/*.root", hWW130, integratedLumi, -1, -1, kBlue);
    ProcessSample(dataset+"/GluGluToHToWWTo2L2Nu_M-140_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/"+version+"/*.root", hWW140, integratedLumi, -1, -1, kBlue);
    ProcessSample(dataset+"/GluGluToHToWWTo2L2Nu_M-150_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/"+version+"/*.root", hWW150, integratedLumi, -1, -1, kBlue);
    ProcessSample(dataset+"/GluGluToHToWWTo2L2Nu_M-160_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/"+version+"/*.root", hWW160, integratedLumi, -1, -1, kBlue);
    ProcessSample(dataset+"/GluGluToHToWWTo2L2Nu_M-170_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/"+version+"/*.root", hWW170, integratedLumi, -1, -1, kBlue);
    ProcessSample(dataset+"/GluGluToHToWWTo2L2Nu_M-180_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/"+version+"/*.root", hWW180, integratedLumi, -1, -1, kBlue);
    ProcessSample(dataset+"/GluGluToHToWWTo2L2Nu_M-190_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/"+version+"/*.root", hWW190, integratedLumi, -1, -1, kBlue);
    ProcessSample(dataset+"/GluGluToHToWWTo2L2Nu_M-200_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/"+version+"/*.root", hWW200, integratedLumi, -1, -1, kBlue);
    ProcessSample(dataset+"/GluGluToHToWWTo2L2Nu_M-210_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/"+version+"/*.root", hWW210, integratedLumi, -1, -1, kBlue);
    ProcessSample(dataset+"/GluGluToHToWWTo2L2Nu_M-220_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/"+version+"/*.root", hWW220, integratedLumi, -1, -1, kBlue);
    ProcessSample(dataset+"/GluGluToHToWWTo2L2Nu_M-230_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/"+version+"/*.root", hWW230, integratedLumi, -1, -1, kBlue);
    ProcessSample(dataset+"/GluGluToHToWWTo2L2Nu_M-250_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/"+version+"/*.root", hWW250, integratedLumi, -1, -1, kBlue);
    ProcessSample(dataset+"/GluGluToHToWWTo2L2Nu_M-300_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/"+version+"/*.root", hWW300, integratedLumi, -1, -1, kBlue);
  }

  if (runWZ)
    ProcessSample(dataset+"/WZtoAnything_TuneZ2_7TeV-pythia6-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/"+version+"/*.root", WZ, integratedLumi, -1, -1, kBlue);
  
  if (runZZ)
    ProcessSample(dataset+"/ZZtoAnything_TuneZ2_7TeV-pythia6-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-01-01"+version+"/*.root", ZZ, integratedLumi, -1, -1, kGreen);
 
  if (runWjets){
    std::vector<string> wSamples;
    wSamples.push_back(dataset+"/WToENu_TuneZ2_7TeV-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/"+version+"/*.root");
    wSamples.push_back(dataset+"/WToMuNu_TuneZ2_7TeV-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/"+version+"/*.root");
    wSamples.push_back(dataset+"/WToTauNu_TuneZ2_7TeV-pythia6-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/"+version+"/*.root");
    ProcessSample(wSamples, Wjets,  integratedLumi , -1, -1, 40);
  }

  if (runDYee)
    ProcessSample(dataset+"/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Spring11-PU_S1_START311_V1G1-v1/V04-01-00/"+version+"/*root", DYee, integratedLumi, -1, -1, kMagenta, identifyDYEvents);
  
  if (runDYmm)
    ProcessSample(dataset+"/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Spring11-PU_S1_START311_V1G1-v1/V04-01-00/"+version+"/*root", DYmm, integratedLumi, -1, -1, kMagenta, identifyDYEvents);
  
  if (runDYtt)
    ProcessSample(dataset+"/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/"+version+"/*root", DYtt, integratedLumi, -1, -1, kMagenta, identifyDYEvents);
 
  if (runttbar)
    ProcessSample(dataset+"/TTJets_TuneZ2_7TeV-madgraph-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/"+version+"/*root", ttbar, integratedLumi, -1, -1, kYellow);

  if (runtW)
    ProcessSample(dataset+"/TToBLNu_TuneZ2_tW-channel_7TeV-madgraph_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/"+version+"/*.root", tW, integratedLumi, -1, -1, 63);
  
  std::vector<string> qcdSamples;
  qcdSamples.push_back(dataset+"/QCD_Pt30_Spring10-START3X_V26_S09-v1/V03-04-08/"+version+"/merged_ntuple*.root");
  qcdSamples.push_back(dataset+"/QCD_Pt80_Spring10-START3X_V26_S09-v1/V03-04-08/"+version+"/merged_ntuple*.root");
  if (runQCD)
    ProcessSample(qcdSamples, qcd, integratedLumi, -1, -1, 40, false, true);
  
  // RealData
  TString cms2_json_file = "files/merged_JsonReRecoSep17_JsonStreamExpressV2_35.49.txt";

  std::vector<string> dataSamples;
  dataSamples.push_back(dataset+"/Mu_Run2010A-Sep17ReReco_v2_RECO/V03-06-14/diLepPt1020Skim/"+version+"/*.root");
  dataSamples.push_back(dataset+"/Mu_Run2010B-PromptReco-v2_RECO/V03-06-14-00/diLepPt1020Skim/"+version+"/*.root");
  dataSamples.push_back(dataset+"/Mu_Run2010B-PromptReco-v2_RECO/V03-06-14/diLepPt1020Skim/"+version+"/*.root");
  dataSamples.push_back(dataset+"/EG_Run2010A-Sep17ReReco_v2_RECO/V03-06-14/diLepPt1020Skim/"+version+"/*.root");
  dataSamples.push_back(dataset+"/Electron_Run2010B-PromptReco-v2_RECO/V03-06-14-00/diLepPt1020Skim/"+version+"/*.root");
  dataSamples.push_back(dataset+"/Electron_Run2010B-PromptReco-v2_RECO/V03-06-14/diLepPt1020Skim/"+version+"/*.root");
  
  if (runData)
    ProcessSample(dataSamples, Data, 3.1, -1, -1, kBlack, false, false, zStudy, true, cms2_json_file);
  
  if (gSystem->Getenv("SkimSamples")) return;
  
  //
  // save all the histograms
  //
  //saveHist(outFile);

  TList* list = gDirectory->GetList() ;
  TIterator* iter = list->MakeIterator();
  
  TRegexp re("*",kTRUE) ;
  
  TFile outf(outFile,"RECREATE") ;
  while(TObject* obj = iter->Next()) {
    if (TString(obj->GetName()).Index(re)>=0) {
      obj->Write() ;
      cout << "." ;
      cout.flush() ;
    }
  }
  cout << endl ;
  outf.Close() ;
  
  delete iter ;
 
}
