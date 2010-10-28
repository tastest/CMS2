//
// DOCUMENT ME
//
#if defined(__CINT__) && !defined(__MAKECINT__)
{
  gSystem->Load("libCMS2NtupleMacrosCORE");
  gSystem->Load("libRooFit.so");
  gSystem->Load("libCMS2NtupleMacrosLooper");
  gSystem->CompileMacro("processData.C","k");
  processData();
}
#endif 

#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "RooDataSet.h"
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
  bool runWZ    = 0;
  bool runZZ    = 0;
  bool runWjets = 0;
  bool runDYee  = 0;
  bool runDYmm  = 0;
  bool runDYtt  = 0;
  bool runttbar = 0;
  bool runtW    = 0;
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

  RooDataSet *fullDataSet(0);
  
  // read dataset prefix
  string dataset = "data";
  if (gSystem->Getenv("DataDir")){
    dataset = gSystem->Getenv("DataDir");
    cout << "Dataset directory: " << dataset << endl;
  }
  
  const double integratedLumi = 100; // pb^1

  if (runWW)
    // ProcessSample(dataset+"/WW_Spring10-START3X_V26_S09-v1_DiLep/V03-04-08/"+version+"/*.root", WW, integratedLumi, 43, -1, fullDataSet, kRed);
    ProcessSample(dataset+"/WWTo2L2Nu_TuneZ2_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-14/"+version+"/*.root", WW, integratedLumi, 4.5, -1, fullDataSet, kRed);
  // ProcessSample("/tas/yygao/WWTo2L2Nu_TuneZ2_7TeV-pythia6_Fall10-START38_V12-v1/merged_ntuple*.root", WW, integratedLumi, 4.5, -1, fullDataSet, kRed);

  if (runWZ)
    ProcessSample(dataset+"/WZ_Spring10-START3X_V26_S09-v1/V03-04-08/"+version+"/*.root", WZ, integratedLumi, 18.2, -1, fullDataSet, kBlue);
  
  if (runZZ)
    ProcessSample(dataset+"/ZZ_Spring10-START3X_V26_S09-v1_DiLep/V03-04-08/"+version+"/*.root", ZZ, integratedLumi, 5.9, -1, fullDataSet, kGreen);
 
  if (runWjets){
    if ( 1 ) {
      // ProcessSample(dataset+"/WJets-madgraph_Spring10-START3X_V26_S09-v1_SingleLep/V03-04-08/"+version+"/dilep-skim.root", Wjets, 100.0, 28049, -1, fullDataSet, 40);
      // ProcessSample(dataset+"/WJets-madgraph_Spring10-START3X_V26_S09-v1_SingleLep/V03-04-13-07/"+version+"/merged_ntuple*root", Wjets, 100.0, 28049, -1, fullDataSet, 40);
      std::vector<string> wjetsSamples;
      wjetsSamples.push_back(dataset+"/WToENu_TuneZ2_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-14/"+version+"/*.root");
      wjetsSamples.push_back(dataset+"/WToMuNu_TuneZ2_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-14/"+version+"/*.root");
      wjetsSamples.push_back(dataset+"/WToTauNu_TuneZ2_7TeV-pythia6-tauola_Fall10-START38_V12-v1/V03-06-14/"+version+"/*.root");
      ProcessSample(wjetsSamples, Wjets, integratedLumi, 10438*0.742, -1, fullDataSet, 40);
    } else {
    // ProcessSample(dataset+"/W1Jets_Pt0to100-alpgenSpring10-START3X_V26_multilepton-v1/V03-04-08/"+version+"/*.root", Wjets, 100.0, 3790, 54751345, fullDataSet, 40);
    // ProcessSample(dataset+"/W0Jets_Pt0to100-alpgenSummer10-START36_V10_multilepton-v1/V03-04-08/"+version+"/*.root", Wjets, 100.0, 20727.9, 3586806, fullDataSet, 40);
    // ProcessSample(dataset+"/W2Jets_Pt0to100-alpgenSummer10-START36_V10_multilepton-v1/V03-04-08/"+version+"/*.root", Wjets, 100.0, 931, 659865, fullDataSet, 40);
    // ProcessSample(dataset+"/Wc1Jets-alpgenSummer10-START3X_V26_S09_multilepton-v1/V03-04-08/"+version+"/*.root", Wjets, 100.0, 420.18, 335056, fullDataSet, 40);
    // ProcessSample(dataset+"/WCJets_7TeV-madgraph_Spring10-START3X_V26-v1/V03-04-13-01/"+version+"/*.root", Wjets, 100.0, -1, -1, fullDataSet, 40);
    // ProcessSample(dataset+"/PhotonVJets-madgraph_Spring10-START3X_V26_S09-v1-bugfix/V03-04-08-01"+version+"/merged_ntuple*.root", Wjets, 100.0, -1, -1, fullDataSet, 40);

    // ProcessSample(dataset+"/W0Jets_Pt0to100-alpgenSummer10-START36_V10_multilepton-v1/V03-04-08/"+version+"/*.root", WZ, 100.0, 20727.9, 3586806, fullDataSet, 40);
    // ProcessSample(dataset+"/W1Jets_Pt0to100-alpgenSpring10-START3X_V26_multilepton-v1/V03-04-08/"+version+"/*.root", ZZ, 100.0, 3790, 54751345, fullDataSet, 40);
    // ProcessSample(dataset+"/W2Jets_Pt0to100-alpgenSummer10-START36_V10_multilepton-v1/V03-04-08/"+version+"/*.root", Wjets, 100.0, 931, 659865, fullDataSet, 40);
    // ProcessSample(dataset+"/Wc0Jets-alpgenSummer10-START3X_V26_S09_multilepton-v1/V03-04-08/"+version+"/*.root", DYee, 100.0, 420.18, 335056, fullDataSet, 40);
    // ProcessSample(dataset+"/Wc1Jets-alpgenSummer10-START3X_V26_S09_multilepton-v1/V03-04-08/"+version+"/*.root", DYmm, 100.0, 143.577, 103835, fullDataSet, 40);
    ProcessSample(dataset+"/PhotonVJets-madgraph_Spring10-START3X_V26_S09-v1-bugfix/V03-04-08-01"+version+"/merged_ntuple*.root", DYtt, 100.0, -1, -1, fullDataSet, 40);
    // ProcessSample(dataset+"/WCJets_7TeV-madgraph_Spring10-START3X_V26-v1/V03-04-13-01/"+version+"/*.root", ttbar, 100.0, -1, -1, fullDataSet, 40);
    // ProcessSample(dataset+"/WJets-madgraph_Spring10-START3X_V26_S09-v1_SingleLep/V03-04-08/"+version+"/dilep-skim.root", tW, 100.0, 28049, -1, fullDataSet, 40);
    }
  }

  if (runDYee)
    // ProcessSample(dataset+"/Zee_Spring10-START3X_V26_S09-v1/V03-04-08-01/"+version+"/*.root", DYee, integratedLumi, 1482.0, -1, fullDataSet, kMagenta, identifyDYEvents);
    ProcessSample(dataset+"/DYToEE_M-20_TuneZ2_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-14/"+version+"/*.root", DYee, integratedLumi, 1666, -1, fullDataSet, kMagenta, identifyDYEvents);
 
  if (runDYmm)
    // ProcessSample(dataset+"/Zmumu_Spring10-START3X_V26_S09-v1/V03-04-08-01/"+version+"/*.root", DYmm, integratedLumi, 1482.0, -1, fullDataSet, kMagenta, identifyDYEvents);
    ProcessSample(dataset+"/DYToMuMu_M-20_TuneZ2_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-14/"+version+"/*.root", DYmm, integratedLumi, 1666, -1, fullDataSet, kMagenta, identifyDYEvents);
 
  if (runDYtt)
    // ProcessSample(dataset+"/Ztautau_Spring10-START3X_V26_S09-v1/V03-04-08-01/"+version+"/*.root", DYtt, integratedLumi, 1482.0, -1, fullDataSet, kMagenta, identifyDYEvents);
    ProcessSample(dataset+"/DYToTauTau_M-20_TuneZ2_7TeV-pythia6-tauola_Fall10-START38_V12-v1/V03-06-14/"+version+"/*.root", DYtt, integratedLumi, 1666, -1, fullDataSet, kMagenta, identifyDYEvents);
 
  if (runttbar)
    //ProcessSample(dataset+"/TTbarJets-madgraph_Summer09-MC_31X_V3_7TeV-v2/"+version+"/merged_ntuple*.root", ttbar, integratedLumi, 165.0, -1, fullDataSet, kYellow);
    ProcessSample(dataset+"/TTbarJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-07/"+version+"/*.root", ttbar, integratedLumi, 165.0, -1, fullDataSet, kYellow);
    //ProcessSample(dataset+"/TTbar_Spring10-START3X_V26_S09-v1/V03-04-08"+version+"/merged_ntuple*.root", ttbar, integratedLumi, 165.0, -1, fullDataSet, kYellow);

  if (runtW)
    ProcessSample(dataset+"/SingleTop_tWChannel-madgraph_Spring10-START3X_V26_S09-v1/V03-04-07/"+version+"/*.root", tW, integratedLumi, 10.6, -1, fullDataSet, 63);
  
  std::vector<string> qcdSamples;
  //qcdSamples.push_back(dataset+"/QCD_Pt30_Summer09-MC_31X_V3_7TeV-v1/"+version+"/merged_ntuple*.root");
  //qcdSamples.push_back(dataset+"/QCD_Pt80_Summer09-MC_31X_V3_7TeV-v1/"+version+"/merged_ntuple*.root");
  qcdSamples.push_back(dataset+"/QCD_Pt30_Spring10-START3X_V26_S09-v1/V03-04-08/"+version+"/merged_ntuple*.root");
  qcdSamples.push_back(dataset+"/QCD_Pt80_Spring10-START3X_V26_S09-v1/V03-04-08/"+version+"/merged_ntuple*.root");
  if (runQCD)
    ProcessSample(qcdSamples, qcd, integratedLumi, -1, -1, fullDataSet, 40, false, true);
  
  /*
  //QCD file
  TChain *fqcd = new TChain("Events");
  if (runQCD) {
    fqcd->Add((dataset+"/cms2-V01-02-06/InclusiveMuPt15/merged_ntuple*.root").c_str());
    fqcd->Add((dataset+"/cms2-V01-02-06/InclusiveMu5Pt50/merged_ntuple*.root").c_str());
    fqcd->Add((dataset+"/cms2-V01-02-06/QCD_EMenriched_Pt20to30/merged_ntuple*.root").c_str());
    fqcd->Add((dataset+"/cms2-V01-02-06/QCD_EMenriched_Pt30to80/merged_ntuple*.root").c_str());
    fqcd->Add((dataset+"/cms2-V01-02-06/QCD_EMenriched_Pt80to170/merged_ntuple*.root").c_str());
    fqcd->Add((dataset+"/cms2-V01-02-06/QCD_BCtoE_Pt20to30/merged_ntuple*.root").c_str());
    fqcd->Add((dataset+"/cms2-V01-02-06/QCD_BCtoE_Pt30to80/merged_ntuple*.root").c_str());
    fqcd->Add((dataset+"/cms2-V01-02-06/QCD_BCtoE_Pt80to170/merged_ntuple*.root").c_str());
  }
  */
  /*
  if (runQCD) {
    cout << "Processing QCD.."<<endl;
    RooDataSet* data = ScanChain(fWjets, qcd, false);
    if ( data ){
      if ( fullDataSet )
	fullDataSet->append(*data);
      else
	fullDataSet=data;
    }
    hist::color("qcd", 40);
  }
  */
  
  // RealData
  //TString cms2_json_file = "Cert_TopSep11_Merged_135059-144114_allPVT.txt";
  TString cms2_json_file;// = "Cert_132440-147454_7TeV_StreamExpress_Collisions10_JSON.txt";
  // TString cms2_json_file = "Cert_132440-148058_7TeV_StreamExpress_Collisions10_JSON.txt";

  std::vector<string> dataSamples;
  //dataSamples.push_back("/tas/yygao/dataskims/diLepPt2020wwskim.root"); // this is the skim corresponds to the 3.1/pb with Json files at Sep11
  // dataSamples.push_back(dataset+"/Mu_Run2010A-Sep17ReReco_v2_RECO/V03-06-09/diLepPt1020Skim/*.root");
  // dataSamples.push_back(dataset+"/Mu_Run2010A-Sep17ReReco_v2_RECO/V03-06-09/singleLepPt10Skim/*.root");
  // dataSamples.push_back(dataset+"/Mu_Run2010A-Sep17ReReco_v2_RECO/V03-06-09/diLepPt1020Skim/skimmed_ntuple_142132*.root");
  // dataSamples.push_back(dataset+"/Mu_Run2010A-Sep17ReReco_v2_RECO/V03-06-09/diLepPt1020Skim/"+version+"/*.root");
  // dataSamples.push_back(dataset+"/Mu_Run2010B-PromptReco-v2_RECO/V03-06-09/diLepPt1020Skim/"+version+"/*.root");
  dataSamples.push_back(dataset+"/Mu_Run2010B-PromptReco-v2_RECO/V03-06-14/diLepPt1020Skim/"+version+"/*.root");
  // dataSamples.push_back(dataset+"/EG_Run2010A-Sep17ReReco_v2_RECO/V03-06-09/diLepPt1020Skim/"+version+"/*.root");
  // dataSamples.push_back(dataset+"/Electron_Run2010B-PromptReco-v2_RECO/V03-06-09/diLepPt1020Skim/"+version+"/*.root");
  dataSamples.push_back(dataset+"/Electron_Run2010B-PromptReco-v2_RECO/V03-06-14/diLepPt1020Skim/"+version+"/*.root");
  
  if (runData)
    ProcessSample(dataSamples, Data, 3.1, -1, -1, fullDataSet, kBlack, false, false, zStudy, true, cms2_json_file);
  
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
  if (fullDataSet) {
    std::string description("Full N-1 dataset:");
    if (runWW)     description+=" WW";
    if (runWZ)     description+=" WZ";
    if (runZZ)     description+=" ZZ";
    if (runWjets)  description+=" Wjets";
    if (runDYee)   description+=" DYee";
    if (runDYmm)   description+=" DYmm";
    if (runDYtt)   description+=" DYtt";
    if (runttbar)  description+=" ttbar";
    if (runtW)     description+=" tW";
    if (runQCD)    description+=" QCD";
    if (runData)   description+=" Data";
    fullDataSet->SetName("fulldataset");
    fullDataSet->SetTitle(description.c_str());
    fullDataSet->Write();
  }
  outf.Close() ;
  
  delete iter ;
 
}
