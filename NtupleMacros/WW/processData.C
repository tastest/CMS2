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
  const bool identifyDYEvents = true;
  const bool zStudy = false;

  //
  // Flags for files to run over 
  // (0 and 1 are easier to modify)
  //
  bool runTest   = 0;
  bool runWW     = 1;
  bool runHWW    = 0;
  bool runGGWW   = 1;
  bool runWZ     = 1;
  bool runZZ     = 1;
  bool runWjets  = 1;
  bool runWgamma = 0;
  bool runDYee   = 1;
  bool runDYmm   = 1;
  bool runDYtt   = 1;
  bool runttbar  = 1;
  bool runtW     = 1;
  bool runQCD    = 0; 
  bool runData   = 1;

  // 
  // Ntuple version
  //
  string version = "wwfilter";
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

  if (runTest){
    std::vector<string> samples;
    samples.push_back("/nfs-6/userdata/hww/GluGluToHToWWTo2L2Nu_M-160_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/merged_ntuple.root");
    // samples.push_back("/nfs-6/userdata/hww/GluGluToHToWWTo2L2Nu_M-160_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-26a/DiLeptonFilter/merged_ntuple.root");
    // samples.push_back("/nfs-4/userdata/cerati/GluGluToWWTo4L_TuneZ2_7TeV-gg2ww-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged_ntuple.root");
    ProcessSample(samples, SmurfTree::hww160, integratedLumi, -1, -1, kBlue);
  }

  if (runWW)
    ProcessSample(dataset+"/WWJetsTo2L2Nu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root", 
		  SmurfTree::qqww, integratedLumi, 4.5*0.919, 963356*(682015./963356.), kRed, true);

  if (runWgamma)
    ProcessSample(dataset+"/PhotonVJets_7TeV-madgraph_Spring11-PU_S1_START311_V1G1-v1/V04-01-12/"+version+"/*.root", 
		  SmurfTree::wgamma, integratedLumi, -1, -1, kRed, true);

  if (runGGWW)
    ProcessSample(dataset+"/GluGluToWWTo4L_TuneZ2_7TeV-gg2ww-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root", 
		  SmurfTree::ggww, integratedLumi, 0.153, -1, kRed, true);

  if (runHWW){
    std::vector<string> samples;

    samples.push_back(dataset+"/GluGluToHToWWTo2LAndTau2Nu_M-110_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    ProcessSample(samples, SmurfTree::hww110, integratedLumi, -1, -1, kBlue); samples.clear();

    samples.push_back(dataset+"/GluGluToHToWWTo2LAndTau2Nu_M-115_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    ProcessSample(samples, SmurfTree::hww115, integratedLumi, -1, -1, kBlue); samples.clear();

    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-120_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-120_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-120_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    ProcessSample(samples, SmurfTree::hww120, integratedLumi, -1, -1, kBlue); samples.clear();

    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-130_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-130_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-130_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    ProcessSample(samples, SmurfTree::hww130, integratedLumi, -1, -1, kBlue); samples.clear();

    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-140_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-140_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-140_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    ProcessSample(samples, SmurfTree::hww140, integratedLumi, -1, -1, kBlue); samples.clear();

    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-150_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-150_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-150_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    ProcessSample(samples, SmurfTree::hww150, integratedLumi, -1, -1, kBlue); samples.clear();

    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-160_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-160_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-160_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    ProcessSample(samples, SmurfTree::hww160, integratedLumi, -1, -1, kBlue); samples.clear();

    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-170_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-170_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-170_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    ProcessSample(samples, SmurfTree::hww170, integratedLumi, -1, -1, kBlue); samples.clear();

    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-180_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-180_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-180_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    ProcessSample(samples, SmurfTree::hww180, integratedLumi, -1, -1, kBlue); samples.clear();

    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-190_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-190_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-190_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    ProcessSample(samples, SmurfTree::hww190, integratedLumi, -1, -1, kBlue); samples.clear();

    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-200_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-200_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-200_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    ProcessSample(samples, SmurfTree::hww200, integratedLumi, -1, -1, kBlue); samples.clear();

    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-210_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-210_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-210_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    ProcessSample(samples, SmurfTree::hww210, integratedLumi, -1, -1, kBlue); samples.clear();

    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-220_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-220_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-220_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    ProcessSample(samples, SmurfTree::hww220, integratedLumi, -1, -1, kBlue); samples.clear();

    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-230_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-230_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-230_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    ProcessSample(samples, SmurfTree::hww230, integratedLumi, -1, -1, kBlue); samples.clear();

    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-250_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-250_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-250_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    ProcessSample(samples, SmurfTree::hww250, integratedLumi, -1, -1, kBlue); samples.clear();

    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-300_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-300_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-300_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    ProcessSample(samples, SmurfTree::hww300, integratedLumi, -1, -1, kBlue); samples.clear();

    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-350_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-350_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-350_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    ProcessSample(samples, SmurfTree::hww350, integratedLumi, -1, -1, kBlue); samples.clear();

    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-400_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-400_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-400_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    ProcessSample(samples, SmurfTree::hww400, integratedLumi, -1, -1, kBlue); samples.clear();

    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-450_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-450_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-450_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    ProcessSample(samples, SmurfTree::hww450, integratedLumi, -1, -1, kBlue); samples.clear();

    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-500_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-500_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-500_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    ProcessSample(samples, SmurfTree::hww500, integratedLumi, -1, -1, kBlue); samples.clear();

    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-550_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-550_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-550_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    ProcessSample(samples, SmurfTree::hww550, integratedLumi, -1, -1, kBlue); samples.clear();

    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-600_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-600_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-600_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    ProcessSample(samples, SmurfTree::hww600, integratedLumi, -1, -1, kBlue); samples.clear();

  }

  if (runWZ)
    ProcessSample(dataset+"/WZJetsTo3LNu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root", SmurfTree::wz, integratedLumi, -1, -1, kBlue);
  
  if (runZZ)
    ProcessSample(dataset+"/ZZJetsTo2L2Nu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root", SmurfTree::zz, integratedLumi, 0.283, -1, kGreen);
  
  if (runWjets){
    std::vector<string> wSamples;
    wSamples.push_back(dataset+"/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(wSamples, SmurfTree::wjets,  integratedLumi , -1, -1, 40);
  }

  if (runttbar)
    ProcessSample(dataset+"/TTJets_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v2/V04-02-36/vvfilter/*root", 
		  SmurfTree::ttbar, integratedLumi, -1, -1, kYellow);

  if (runtW) {
    std::vector<string> samples;
    samples.push_back(dataset+"/T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::tw, integratedLumi, -1, -1, 63);
  }

  if (runDYee){
    std::vector<string> samples;
    samples.push_back(dataset+"/DYToEE_M-10To20_CT10_TuneZ2_7TeV-powheg-pythia_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::dyee, integratedLumi, -1, -1, kMagenta, identifyDYEvents);
  }

  if (runDYmm){
    std::vector<string> samples;
    samples.push_back(dataset+"/DYToMuMu_M-10To20_CT10_TuneZ2_7TeV-powheg-pythia_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::dymm, integratedLumi, -1, -1, kMagenta, identifyDYEvents);
  }
  if (runDYtt){
    std::vector<string> samples;
    samples.push_back(dataset+"/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::dytt, integratedLumi, -1, -1, kMagenta, identifyDYEvents);
  }
 
  std::vector<string> qcdSamples;
  qcdSamples.push_back(dataset+"/QCD_Pt30_Spring10-START3X_V26_S09-v1/V03-04-08/"+version+"/merged_ntuple*.root");
  qcdSamples.push_back(dataset+"/QCD_Pt80_Spring10-START3X_V26_S09-v1/V03-04-08/"+version+"/merged_ntuple*.root");
  if (runQCD)
    ProcessSample(qcdSamples, SmurfTree::qcd, integratedLumi, -1, -1, 40, false, true);
  
  // RealData
  // TString cms2_json_file = "files/Cert_160404-166861_7TeV_PromptReco_Collisions11_JSON_715ipb.txt";
  // TString cms2_json_file = "files/Cert_EPSFINAL_May10ReReco_v2_PromptReco_160404_167913_JSON.txt";
  // TString cms2_json_file = "files/HWW.conservativeCertificationLP11.json";
  TString cms2_json_file = "files/Cert_160404-178677_7TeV_PromptReco_Collisions11_JSON.txt";

  std::vector<string> dataSamples;
  dataSamples.push_back(dataset+"/DoubleElectron_Run2011A-May10ReReco-v1/V04-02-36/vvfilter/*.root");
  dataSamples.push_back(dataset+"/DoubleMu_Run2011A-May10ReReco-v1/V04-02-36/vvfilter/*.root");
  dataSamples.push_back(dataset+"/MuEG_Run2011A-May10ReReco-v1/V04-02-36/vvfilter/*.root");
  dataSamples.push_back(dataset+"/SingleElectron_Run2011A-May10ReReco-v1/V04-02-36/vvfilter/*.root");
  dataSamples.push_back(dataset+"/SingleMu_Run2011A-May10ReReco-v1/V04-02-36/vvfilter/*.root");

  dataSamples.push_back(dataset+"/DoubleElectron_Run2011A-PromptReco-v4/V04-02-36/vvfilter/*.root");
  dataSamples.push_back(dataset+"/MuEG_Run2011A-PromptReco-v4/V04-02-36/vvfilter/*.root");
  dataSamples.push_back(dataset+"/SingleElectron_Run2011A-PromptReco-v4/V04-02-36/vvfilter/*.root");
  dataSamples.push_back(dataset+"/SingleMu_Run2011A-PromptReco-v4/V04-02-36/vvfilter/*.root");
  dataSamples.push_back(dataset+"/DoubleMu_Run2011A-PromptReco-v4/V04-02-36/vvfilter/*.root");

  dataSamples.push_back(dataset+"/DoubleElectron_Run2011A-05Aug2011-v1/V04-02-36/vvfilter/*.root");
  dataSamples.push_back(dataset+"/MuEG_Run2011A-05Aug2011-v1/V04-02-36/vvfilter/*.root");
  dataSamples.push_back(dataset+"/SingleElectron_Run2011A-05Aug2011-v1/V04-02-36/vvfilter/*.root");
  dataSamples.push_back(dataset+"/SingleMu_Run2011A-05Aug2011-v1/V04-02-36/vvfilter/*.root");
  dataSamples.push_back(dataset+"/DoubleMu_Run2011A-05Aug2011-v1/V04-02-36/vvfilter/*.root");

  dataSamples.push_back(dataset+"/DoubleElectron_Run2011A-PromptReco-v6/V04-02-36/vvfilter/*.root");
  dataSamples.push_back(dataset+"/MuEG_Run2011A-PromptReco-v6/V04-02-36/vvfilter/*.root");
  // dataSamples.push_back(dataset+"/SingleElectron_Run2011A-PromptReco-v6/V04-02-36/vvfilter/*.root");
  // dataSamples.push_back(dataset+"/SingleMu_Run2011A-PromptReco-v6/V04-02-36/vvfilter/*.root");
  dataSamples.push_back(dataset+"/DoubleMu_Run2011A-PromptReco-v6/V04-02-36/vvfilter/*.root");

  dataSamples.push_back(dataset+"/DoubleElectron_Run2011B-PromptReco-v1/V04-02-36/vvfilter/*.root");
  dataSamples.push_back(dataset+"/MuEG_Run2011B-PromptReco-v1/V04-02-36/vvfilter/*.root");
  // dataSamples.push_back(dataset+"/SingleElectron_Run2011B-PromptReco-v1/V04-02-36/vvfilter/*.root");
  // dataSamples.push_back(dataset+"/SingleMu_Run2011B-PromptReco-v1/V04-02-36/vvfilter/*.root");
  dataSamples.push_back(dataset+"/DoubleMu_Run2011B-PromptReco-v1/V04-02-36/vvfilter/*.root");


  if (runData)
    ProcessSample(dataSamples, SmurfTree::data, 3.1, -1, -1, kBlack, false, false, zStudy, true, cms2_json_file);
  
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
