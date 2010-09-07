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
  
  //
  // Flags for files to run over 
  // (0 and 1 are easier to modify)
  //
  bool runWW    = 1;
  bool runWZ    = 1;
  bool runZZ    = 1;
  bool runWjets = 1;
  bool runDYee  = 1;
  bool runDYmm  = 1;
  bool runDYtt  = 1;
  bool runttbar = 1;
  bool runtW    = 1;
  bool runQCD   = 0; 

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

  if (runWW)
    ProcessSample(dataset+"/WW_Spring10-START3X_V26_S09-v1_DiLep/V03-04-08/"+version+"/merged_ntuple*.root", WW, 100.0, 43, fullDataSet, kRed);

  if (runWZ)
    ProcessSample(dataset+"/WZ_Spring10-START3X_V26_S09-v1/V03-04-08/"+version+"/merged_ntuple*.root", WZ, 100.0, 18.2, fullDataSet, kBlue);
  
  if (runZZ)
    ProcessSample(dataset+"/ZZ_Spring10-START3X_V26_S09-v1_DiLep/V03-04-08/"+version+"/merged_ntuple*.root", ZZ, 100.0, 5.9, fullDataSet, kGreen);
 
  if (runWjets)
    ProcessSample(dataset+"/WJets-madgraph_Spring10-START3X_V26_S09-v1_SingleLep/V03-04-08/"+version+"/dilep-skim.root", Wjets, 100.0, 28049, fullDataSet, 40);
  // ProcessSample(dataset+"/PhotonVJets-madgraph_Spring10-START3X_V26_S09-v1-bugfix/V03-04-08-01"+version+"/merged_ntuple*.root", Wjets, 100.0, -1, fullDataSet, 40);

  if (runDYee)
    ProcessSample(dataset+"/Zee_Spring10-START3X_V26_S09-v1/V03-04-08-01/"+version+"/merged_ntuple*.root", DYee, 100.0, 1482.0, fullDataSet, kMagenta, identifyDYEvents);
 
  if (runDYmm)
    ProcessSample(dataset+"/Zmumu_Spring10-START3X_V26_S09-v1/V03-04-08-01/"+version+"/merged_ntuple*.root", DYmm, 100.0, 1482.0, fullDataSet, kMagenta, identifyDYEvents);
 
  if (runDYtt)
    ProcessSample(dataset+"/Ztautau_Spring10-START3X_V26_S09-v1/V03-04-08-01/"+version+"/merged_ntuple*.root", DYtt, 100.0, 1482.0, fullDataSet, kMagenta, identifyDYEvents);

  if (runttbar)
    //ProcessSample(dataset+"/TTbarJets-madgraph_Summer09-MC_31X_V3_7TeV-v2/"+version+"/merged_ntuple*.root", ttbar, 100.0, 165.0, fullDataSet, kYellow);
    ProcessSample(dataset+"/TTbarJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-07/"+version+"/merged_ntuple*.root", ttbar, 100.0, 165.0, fullDataSet, kYellow);
    //ProcessSample(dataset+"/TTbar_Spring10-START3X_V26_S09-v1/V03-04-08"+version+"/merged_ntuple*.root", ttbar, 100.0, 165.0, fullDataSet, kYellow);

  if (runtW)
    ProcessSample(dataset+"/SingleTop_tWChannel-madgraph_Spring10-START3X_V26_S09-v1/V03-04-07/"+version+"/merged_ntuple*.root", tW, 100.0, 10.6, fullDataSet, 63);
  
  std::vector<string> qcdSamples;
  //qcdSamples.push_back(dataset+"/QCD_Pt30_Summer09-MC_31X_V3_7TeV-v1/"+version+"/merged_ntuple*.root");
  //qcdSamples.push_back(dataset+"/QCD_Pt80_Summer09-MC_31X_V3_7TeV-v1/"+version+"/merged_ntuple*.root");
  qcdSamples.push_back(dataset+"/QCD_Pt30_Spring10-START3X_V26_S09-v1/V03-04-08/"+version+"/merged_ntuple*.root");
  qcdSamples.push_back(dataset+"/QCD_Pt80_Spring10-START3X_V26_S09-v1/V03-04-08/"+version+"/merged_ntuple*.root");
  if (runQCD)
    ProcessSample(qcdSamples, qcd, 100.0, -1, fullDataSet, 40, false, true);
  
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
    fullDataSet->SetName("fulldataset");
    fullDataSet->SetTitle(description.c_str());
    fullDataSet->Write();
  }
  outf.Close() ;
  
  delete iter ;
 
}
