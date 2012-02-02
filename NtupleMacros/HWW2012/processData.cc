#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TRegexp.h"
#include "TFile.h" 
#include "wwtypes.h"
#include "doAnalysis.h"

int main()
{
  using namespace std;

  //
  // Define how to handle complex Monte Carlo samples which may have 
  // more than one event type mixed in. Careful with this option. 
  // You don't need it for not mixed samples.
  //
  const bool identifyDYEvents = true;

  //
  // Flags for files to run over 
  // (0 and 1 are easier to modify)
  //
  bool runTest   = 1;
  bool runWW     = 0;
  bool runHWW    = 0;
  bool runGGWW   = 0;
  bool runWZ     = 0;
  bool runZZ     = 0;
  bool runWjets  = 0;
  bool runWgamma = 0;
  bool runDYmgee = 0;
  bool runDYmgmm = 0;
  bool runDYee   = 0;
  bool runDYmm   = 0;
  bool runDYtt   = 0;
  bool runttbar  = 0;
  bool runtW     = 0;
  bool runData   = 0;

  // 
  // Ntuple version
  //
  string version = "wwfilter";
  if (gSystem->Getenv("VERSION")){
    version = gSystem->Getenv("VERSION");
    cout << "VERSION: " << version << endl;
  }
 
  // read dataset prefix
  string dataset = "data";
  if (gSystem->Getenv("DataDir")){
    dataset = gSystem->Getenv("DataDir");
    cout << "Dataset directory: " << dataset << endl;
  }
  
  const double integratedLumi = 1000.0; // pb^1  
  cout << "Integrated luminosity to scale to: " << integratedLumi << endl;

  TString cms2_json_file = "";
  if (runTest){
    std::vector<string> samples;
    //samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-160_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/merged_ntuple.root");
    //ProcessSample(samples, SmurfTree::hww160, integratedLumi, -1, -1);
    samples.push_back("/home/users/cerati/UpdateCMS2ForLeptonMVA/CMSSW_4_2_7_patch1_V04-02-36/src/CMS2/NtupleMaker/test/ntuple.root");
    ProcessSample(samples, SmurfTree::data, 3.1, -1, -1, false, true, cms2_json_file);
  }

  if (runWW)
    ProcessSample(dataset+"/WWJetsTo2L2Nu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root", 
		  SmurfTree::qqww, integratedLumi, -1, -1, true);

  if (runWgamma) {
    std::vector<string> samples;
    samples.push_back(dataset+"/WGToMuNuG_TuneZ2_7TeV-madgraph_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/WGToENuG_TuneZ2_7TeV-madgraph_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/WGToTauNuG_TuneZ2_7TeV-madgraph_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples,SmurfTree::wgamma, integratedLumi, -1, -1, true);
  }

  if (runGGWW)
    ProcessSample(dataset+"/GluGluToWWTo4L_TuneZ2_7TeV-gg2ww-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root", 
		  SmurfTree::ggww, integratedLumi, -1, -1, true);

  if (runHWW){
    std::vector<string> samples;

    samples.push_back(dataset+"/GluGluToHToWWTo2LAndTau2Nu_M-110_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    ProcessSample(samples, SmurfTree::hww110, integratedLumi, -1, -1); samples.clear();

    samples.push_back(dataset+"/GluGluToHToWWTo2LAndTau2Nu_M-115_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    ProcessSample(samples, SmurfTree::hww115, integratedLumi, -1, -1); samples.clear();

    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-120_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-120_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-120_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::hww120, integratedLumi, -1, -1); samples.clear();

    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-130_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-130_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-130_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::hww130, integratedLumi, -1, -1); samples.clear();

    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-140_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-140_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-140_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::hww140, integratedLumi, -1, -1); samples.clear();

    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-150_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-150_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-150_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::hww150, integratedLumi, -1, -1); samples.clear();

    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-160_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-160_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-160_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::hww160, integratedLumi, -1, -1); samples.clear();

    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-170_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-170_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-170_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::hww170, integratedLumi, -1, -1); samples.clear();

    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-180_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-180_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-180_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::hww180, integratedLumi, -1, -1); samples.clear();

    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-190_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-190_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-190_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::hww190, integratedLumi, -1, -1); samples.clear();

    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-200_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-200_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-200_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::hww200, integratedLumi, -1, -1); samples.clear();

    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-210_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-210_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-210_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::hww210, integratedLumi, -1, -1); samples.clear();

    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-220_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-220_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-220_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::hww220, integratedLumi, -1, -1); samples.clear();

    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-230_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-230_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-230_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::hww230, integratedLumi, -1, -1); samples.clear();

    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-250_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-250_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-250_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::hww250, integratedLumi, -1, -1); samples.clear();

    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-300_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-300_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-300_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::hww300, integratedLumi, -1, -1); samples.clear();

    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-350_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-350_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-350_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::hww350, integratedLumi, -1, -1); samples.clear();

    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-400_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-400_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-400_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::hww400, integratedLumi, -1, -1); samples.clear();

    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-450_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-450_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-450_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::hww450, integratedLumi, -1, -1); samples.clear();

    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-500_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-500_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-500_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::hww500, integratedLumi, -1, -1); samples.clear();

    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-550_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-550_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-550_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::hww550, integratedLumi, -1, -1); samples.clear();

    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-600_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-600_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-600_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::hww600, integratedLumi, -1, -1); samples.clear();

  }

  if (runWZ)
    ProcessSample(dataset+"/WZJetsTo3LNu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root", SmurfTree::wz, integratedLumi, -1, -1);
  
  if (runZZ)
    ProcessSample(dataset+"/ZZJetsTo2L2Nu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root", SmurfTree::zz, integratedLumi, -1, -1);
  
  if (runWjets){
    std::vector<string> wSamples;
    wSamples.push_back(dataset+"/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(wSamples, SmurfTree::wjets,  integratedLumi , -1, -1);
  }

  if (runttbar)
    ProcessSample(dataset+"/TTTo2L2Nu2B_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*root", 
		  SmurfTree::ttbar, integratedLumi, -1, -1);

  if (runtW) {
    std::vector<string> samples;
    samples.push_back(dataset+"/T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::tw, integratedLumi, -1, -1);
  }

  if (runDYmgee){
    std::vector<string> samples;
    samples.push_back("/nfs-7/userdata/cms2/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged_ntuple*.root");
    ProcessSample(samples, SmurfTree::dyee, integratedLumi, -1, -1, identifyDYEvents);
  }

  if (runDYmgmm){
    std::vector<string> samples;
    samples.push_back("/nfs-7/userdata/cms2/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged_ntuple*.root");
    ProcessSample(samples, SmurfTree::dymm, integratedLumi, -1, -1, identifyDYEvents);
  }

  if (runDYee){
    std::vector<string> samples;
    samples.push_back(dataset+"/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::dyee, integratedLumi, -1, -1, identifyDYEvents);
  }

  if (runDYmm){
    std::vector<string> samples;
    samples.push_back(dataset+"/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/merged_ntuple*.root");
    ProcessSample(samples, SmurfTree::dymm, integratedLumi, -1, -1, identifyDYEvents);
  }

  if (runDYtt){
    std::vector<string> samples;
    samples.push_back(dataset+"/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::dytt, integratedLumi, -1, -1, identifyDYEvents);
  }
 
  // RealData
  //TString cms2_json_file = "files/Cert_160404-178677_7TeV_PromptReco_Collisions11_JSON.txt";

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
  dataSamples.push_back(dataset+"/SingleMu_Run2011B-PromptReco-v1/V04-02-36/vvfilter/*.root");
  dataSamples.push_back(dataset+"/DoubleMu_Run2011B-PromptReco-v1/V04-02-36/vvfilter/*.root");

  if (runData)
    ProcessSample(dataSamples, SmurfTree::data, 3.1, -1, -1, false, true, cms2_json_file);
  
  return 0; 
}
