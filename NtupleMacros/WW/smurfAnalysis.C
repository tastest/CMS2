// 
// It is assumed that fake inputs are calculated as number of events
// when one of the lepton fails nominal selection, but satisfy the fakable
// and the othe lepton satisfies the nominal selection
//
// The fake rate ratio is assumed to be a straight ratio of number of
// leptons that passed the nominal selection devided by number of those
// that passed the fakable object requirements.
// 

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TSystem.h"
#include <iostream>
#include <iomanip>
#include <string>
#include "../../../Smurf/Core/SmurfTree.h"
#include "../Tools/goodrun.cc"

struct Result{
  Result(){ 
    yield[0]=0; yield[1]=0; yield[2]=0; yield[3]=0;
     err2[0]=0;  err2[1]=0;  err2[2]=0;  err2[3]=0;
  }
  void add(unsigned int i, double value){
    yield[i] += value;
    err2[i] += value*value;
  }
  double yield[4];
  double err2[4];
  double total_yield(){ return yield[0]+yield[1]+yield[2]+yield[3]; }
  double total_error(){ return sqrt(err2[0]+err2[1]+err2[2]+err2[3]); }
};
enum SelectionType{
  WW0Jet,
  WW1Jet,
  WW2Jets,
  HWWCutBased0Jet,
  HWWCutBased1Jet,
  HWWCutBased2Jets,
  HWWMVABased0Jet,
  HWWMVABased1Jet,
  HWWMVABased2Jets
};

struct FakeBkg{
  // data estimate not corrected for WW contamination
  Result data_mu;
  Result data_el;
  // MC based expectation for WW contamination in the prediction
  Result ww_mu; 
  Result ww_el;
  // MC prediction
  Result mc;
};

struct Measurement{
  FakeBkg fakeBkgData;
  Result fakeBkgMC;
  SelectionType type;
  SmurfTree::DataType sig_type;
};

struct SmurfAnalysis{
  struct Entry {
    SmurfTree::DataType sample;
    // electron fake histograms by types
    TH2F* hist_el[4]; 
    TH2F* hist_el_ss[4]; 
    // muon fake histograms by types
    TH2F* hist_mu[4]; 
    TH2F* hist_mu_ss[4]; 
    // final selection yield
    TH2F* hist_final[4]; 
    TH2F* hist_final_ss[4]; 
    Result nTopTaggedEvents;
    Result nTop1BJetEvents;
    Result nTop1BJetTopTaggedEvents;
    Result finalSelection;
    Result nZpeak;
  };
  SmurfAnalysis(double lumi,
		const char* dir = "/smurf/dmytro/samples/smurf",
		const char* el_fakerate_file = "files/ww_el_fr.root",
		const char* el_fakerate_name = "el_fr_v4",
		const char* mu_fakerate_file = "files/ww_mu_fr.root",
		const char* mu_fakerate_name = "mu_fr_m2");
  bool passedExtraCuts(SmurfTree& tree);
  void showYields(double lumi = -1, bool tex_format=false); // -1 for lumi means internal lumi
  void estimateFakeBackground(bool verbose = true);
  void estimateTopBackground();
  void estimateDYBackground();
  void setGoodRunList(const char* file){
    json_ = file;
    set_goodrun_file_json(json_);
  }
  void processSamples();
  bool isStandardModel(SmurfTree::DataType sample){
    switch (sample){
    case SmurfTree::qqww: case SmurfTree::ggww: case SmurfTree::ttbar: 
    case SmurfTree::tw: case SmurfTree::dyee: case SmurfTree::dymm: case SmurfTree::dytt:
    case SmurfTree::wjets: case SmurfTree::wz: case SmurfTree::zz: case SmurfTree::wgamma:
    case SmurfTree::qcd:
      return true;
    default:
      return false;
    }
  }
  bool isGGHiggs(SmurfTree::DataType sample){
    switch (sample){
    case SmurfTree::hww115: case SmurfTree::hww120: case SmurfTree::hww130: 
    case SmurfTree::hww140: case SmurfTree::hww150: case SmurfTree::hww160: case SmurfTree::hww170:
    case SmurfTree::hww180: case SmurfTree::hww190: case SmurfTree::hww200: case SmurfTree::hww210:
    case SmurfTree::hww220: case SmurfTree::hww230: case SmurfTree::hww250: case SmurfTree::hww300:
    case SmurfTree::hww350: case SmurfTree::hww400: case SmurfTree::hww450: case SmurfTree::hww500:
    case SmurfTree::hww550: case SmurfTree::hww600:
      return true;
    default:
      return false;
    }
  }
  bool isRelevantToFakeBkg(SmurfTree::DataType sample){
    return sample==SmurfTree::wjets || sample==SmurfTree::data || sample==SmurfTree::qqww ||
      sample==SmurfTree::ttbar;
  }
  bool isRelevantToTopBkg(SmurfTree::DataType sample){
    return sample==SmurfTree::data || sample==SmurfTree::ttbar || sample==SmurfTree::tw || sample==SmurfTree::qqww;
  }
  bool isRelevantToBeProcessed(SmurfTree::DataType sample){
    switch (sample){
    case SmurfTree::qqww: case SmurfTree::ggww: case SmurfTree::ttbar: 
    case SmurfTree::tw: // case SmurfTree::dyee: case SmurfTree::dymm: case SmurfTree::dytt:
    case SmurfTree::wjets: //case SmurfTree::wz: case SmurfTree::zz: case SmurfTree::wgamma: case SmurfTree::qcd:
    case SmurfTree::hww130: case SmurfTree::hww160: case SmurfTree::data:
      return true;
    default:
      return false;
    }
    // return sample==SmurfTree::data;
    // return sample==SmurfTree::data || isStandardModel(sample) || isGGHiggs(sample);
  }

  double lumi_;
  bool processOnlyImportantSamples;
  TH2F *elFakeRate0_;
  TH2F *elFakeRate1_;
  TH2F *elFakeRate2_;
  TH2F *muFakeRate0_;
  TH2F *muFakeRate1_;
  TH2F *muFakeRate2_;
  std::vector<Entry> entries_;
  const char* dir_;
  Measurement measurement_;
  const char* json_;
  double rInOutEl_; // R in/out for Zee
  double rInOutMu_; // R in/out for Zmm
  double kElMu_; // yield difference between els and mus k=sqrt(Nee/Nmm) in Z peak
      
  static const unsigned int cut_base  = SmurfTree::BaseLine|SmurfTree::FullMET|SmurfTree::ZVeto|SmurfTree::TopVeto|SmurfTree::ExtraLeptonVeto;
  static const unsigned int cut_el_1  = SmurfTree::Lep2FullSelection|SmurfTree::Lep1LooseEleV4|cut_base;
  static const unsigned int cut_el_2  = SmurfTree::Lep1FullSelection|SmurfTree::Lep2LooseEleV4|cut_base;
  static const unsigned int cut_mu_1  = SmurfTree::Lep2FullSelection|SmurfTree::Lep1LooseMuV2|cut_base;
  static const unsigned int cut_mu_2  = SmurfTree::Lep1FullSelection|SmurfTree::Lep2LooseMuV2|cut_base;
  static const unsigned int cut_final = SmurfTree::Lep1FullSelection|SmurfTree::Lep2FullSelection|cut_base;

  static const unsigned int cut_final_complete = SmurfTree::ChargeMatch|cut_base;

  static const unsigned int cut_top1B = SmurfTree::Lep1FullSelection|SmurfTree::Lep2FullSelection|SmurfTree::BaseLine|SmurfTree::ChargeMatch|SmurfTree::FullMET|SmurfTree::ZVeto|SmurfTree::OneBJet;
  static const unsigned int cut_topTagged = SmurfTree::Lep1FullSelection|SmurfTree::Lep2FullSelection|SmurfTree::BaseLine|SmurfTree::ChargeMatch|SmurfTree::FullMET|SmurfTree::ZVeto|SmurfTree::TopTag|SmurfTree::ExtraLeptonVeto;
  static const unsigned int cut_final_nomass    = SmurfTree::Lep1FullSelection|SmurfTree::Lep2FullSelection|SmurfTree::BaseLine|SmurfTree::ChargeMatch|SmurfTree::FullMET|SmurfTree::TopVeto|SmurfTree::ExtraLeptonVeto;
  
  // same binning for electrons and muons
  static const unsigned int nPtPoints = 6;
  static const unsigned int nEtaPoints = 5;
  static const Double_t ptBins[nPtPoints];
  static const Double_t etaBins[nEtaPoints];
  static const char*    types[4];
  static const char*  pm;

private:
  bool MCheckBinLimits(const TArrayD* h1Array, const TArrayD* h2Array);
  bool CheckConsistency(const TH1* h1, const TH1* h2);
  void getIntegral(const TH2F* h, double& integral, double& error, const TH2F* fakeRate, 
		   int etaBin=-1, int ptBin=-1); // last two options for slicing in eta and pt bins
  void addSample(SmurfTree::DataType sample);

};

bool HWWCuts_SmurfV5(SmurfTree& tree, SmurfTree::DataType sig_type){
  if (tree.njets_!=0) return false;
  switch (sig_type){
  case SmurfTree::hww120: 
  case SmurfTree::hww115: 
    return tree.lep1_.pt()>20 && tree.lep2_.pt()>10 &&
      tree.dilep_.mass()<40 && fabs(tree.dPhi_)<2.0 && 
      tree.mt_>70 && tree.mt_<120;
  case SmurfTree::hww130: 
    return tree.lep1_.pt()>25 && tree.lep2_.pt()>10 &&
      tree.dilep_.mass()<45 && fabs(tree.dPhi_)<1.5 && 
      tree.mt_>75 && tree.mt_<125;
  case SmurfTree::hww140:
    return tree.lep1_.pt()>25 && tree.lep2_.pt()>15 &&
      tree.dilep_.mass()<45 && fabs(tree.dPhi_)<1.5 &&
      tree.mt_>80 && tree.mt_<130;
  case SmurfTree::hww150:
    return tree.lep1_.pt()>27 && tree.lep2_.pt()>25 &&
      tree.dilep_.mass()<50 && fabs(tree.dPhi_)<1.5 &&
      tree.mt_>80 && tree.mt_<150;
  case SmurfTree::hww160:
    return tree.lep1_.pt()>30 && tree.lep2_.pt()>25 &&
      tree.dilep_.mass()<50 && fabs(tree.dPhi_)<1.0 &&
      tree.mt_>90 && tree.mt_<160;;
  case SmurfTree::hww170:
    return tree.lep1_.pt()>34 && tree.lep2_.pt()>25 &&
      tree.dilep_.mass()<50 && fabs(tree.dPhi_)<M_PI/180*60;
  case SmurfTree::hww180:
    return tree.lep1_.pt()>36 && tree.lep2_.pt()>25 &&
      tree.dilep_.mass()<60 && fabs(tree.dPhi_)<M_PI/180*70;
  case SmurfTree::hww190:
    return tree.lep1_.pt()>38 && tree.lep2_.pt()>25 &&
      tree.dilep_.mass()<80 && fabs(tree.dPhi_)<M_PI/180*90;
  case SmurfTree::hww200:
    return tree.lep1_.pt()>40 && tree.lep2_.pt()>25 &&
      tree.dilep_.mass()<90 && fabs(tree.dPhi_)<M_PI/180*100;
  case SmurfTree::hww210:
    return tree.lep1_.pt()>44 && tree.lep2_.pt()>25 &&
      tree.dilep_.mass()<110 && fabs(tree.dPhi_)<M_PI/180*110;
  case SmurfTree::hww220:
    return tree.lep1_.pt()>48 && tree.lep2_.pt()>25 &&
      tree.dilep_.mass()<120 && fabs(tree.dPhi_)<M_PI/180*120;
  case SmurfTree::hww230:
    return tree.lep1_.pt()>52 && tree.lep2_.pt()>25 &&
      tree.dilep_.mass()<130 && fabs(tree.dPhi_)<M_PI/180*130;
  case SmurfTree::hww250:
    return tree.lep1_.pt()>55 && tree.lep2_.pt()>25 &&
      tree.dilep_.mass()<150 && fabs(tree.dPhi_)<M_PI/180*140;
  default: return false;
  }
}

bool SmurfAnalysis::passedExtraCuts(SmurfTree& tree){
  switch (measurement_.type){
  case WW0Jet:
    return tree.njets_==0;
  case WW1Jet:
    return tree.njets_==1 && 
      (tree.type_==1 || tree.type_==2 || fabs(tree.dPhiDiLepJet1_)<M_PI*165/180);
  case WW2Jets:
    return tree.njets_==2;
  case HWWCutBased0Jet:
    return tree.njets_==0 && HWWCuts_SmurfV5(tree,measurement_.sig_type);
  default:
    return false;
  }
  return false;
}
const Double_t SmurfAnalysis::ptBins[SmurfAnalysis::nPtPoints] = {10., 15., 20., 25., 30., 35.};
const Double_t SmurfAnalysis::etaBins[SmurfAnalysis::nEtaPoints] = {0.0, 1.0, 1.479, 2.0, 2.5};
const char*    SmurfAnalysis::types[4] = {"mm","me","em","ee"};
const char*    SmurfAnalysis::pm = "+/-";


//___________________________________________________________________________

SmurfAnalysis::SmurfAnalysis(double lumi, const char* dir, 
			     const char* el_fakerate_file,
			     const char* el_fakerate_name,
			     const char* mu_fakerate_file,
			     const char* mu_fakerate_name):
  lumi_(lumi),processOnlyImportantSamples(true),dir_(dir),json_(0),
  rInOutEl_(0.16), rInOutMu_(0.20), kElMu_(0.823352)
{
  measurement_.type = WW0Jet;
  measurement_.sig_type = SmurfTree::qqww;
  TFile *el_fakeRateFile = TFile::Open(el_fakerate_file);
  assert(el_fakeRateFile);
  elFakeRate0_ = dynamic_cast<TH2F*>( el_fakeRateFile->Get(Form("%s_15",el_fakerate_name)) );
  assert(elFakeRate0_);
  elFakeRate1_ = dynamic_cast<TH2F*>( el_fakeRateFile->Get(Form("%s_35",el_fakerate_name)) );
  assert(elFakeRate1_);
  elFakeRate2_ = dynamic_cast<TH2F*>( el_fakeRateFile->Get(Form("%s_50",el_fakerate_name)) );
  assert(elFakeRate2_);
  TFile *mu_fakeRateFile = TFile::Open(mu_fakerate_file);
  assert(mu_fakeRateFile);
  muFakeRate0_ = dynamic_cast<TH2F*>( mu_fakeRateFile->Get(Form("%s_5",mu_fakerate_name)) );
  assert(muFakeRate0_);
  muFakeRate1_ = dynamic_cast<TH2F*>( mu_fakeRateFile->Get(Form("%s_15",mu_fakerate_name)) );
  assert(muFakeRate1_);
  muFakeRate2_ = dynamic_cast<TH2F*>( mu_fakeRateFile->Get(Form("%s_30",mu_fakerate_name)) );
  assert(muFakeRate2_);
}

void SmurfAnalysis::processSamples()
{
  using namespace std;
  entries_.clear();
  for (unsigned int i=0; i<58; ++i) //58 is the current number of samples defined
    addSample(SmurfTree::DataType(i));
}

void SmurfAnalysis::estimateDYBackground()
{
  const bool binomialError = false;
  if (entries_.empty()){
    printf("No samples are processed.\n");
    return;
  }
  printf("\nDrell-Yan background estimation inputs :\n\b");;    
  printf(" \t     %s      \t      %s      \t      %s      \t      %s      \t    Total     \n", types[0], types[1], types[2], types[3]);
  for (unsigned int i=0; i<entries_.size(); ++i){
    if ( !isRelevantToTopBkg(entries_.at(i).sample) ) continue;
    printf("%s", SmurfTree::name(entries_.at(i).sample).c_str());
    for (unsigned int j=0; j<4; ++j){
      printf(" \t%4.1f+/-%3.1f", entries_.at(i).nZpeak.yield[j], sqrt(entries_.at(i).nZpeak.err2[j]));
    }
    printf(" \t%4.1f+/-%3.1f\n", entries_.at(i).nZpeak.total_yield(),
	   entries_.at(i).nZpeak.total_error());
  }
  /*
    double A    = entries_.at(i).nTop1BJetTopTaggedEvents.total_yield();
    double errA = entries_.at(i).nTop1BJetTopTaggedEvents.total_error();
    double B    = entries_.at(i).nTop1BJetEvents.total_yield();
    double errB = entries_.at(i).nTop1BJetEvents.total_error();
    double C    = entries_.at(i).nTopTaggedEvents.total_yield();
    double errC = entries_.at(i).nTopTaggedEvents.total_error();
    
    // compute efficiencies
    double eff_1b = B>0 ? A/B : 0;
    double eff_1b_err(0);
    if ( binomialError ){
      eff_1b_err = sqrt(eff_1b*(1-eff_1b)/B);
    } else {
      eff_1b_err = eff_1b>0 ? eff_1b * sqrt(pow(errA/A,2)+pow(errB/B,2)): 0;
    }
    double eff = 1-pow(1-eff_1b,2);
    double eff_err = 2*(1-eff_1b)*eff_1b_err;

    double scale = eff > 0 ? (1-eff)/eff : 0;
    printf("\ttag eff: %0.3f+/-%0.3f \ttop bkg: %0.2f+/-%0.2f\n", 
	   eff, eff_err, C*scale, sqrt(C)*scale);
  }  

  //
  // MC prediction
  //
  printf("\nTop yield MC estimation:\n\b");;    
  printf(" \t     %s      \t      %s      \t      %s      \t      %s      \t    Total     \n", types[0], types[1], types[2], types[3]);

  for (unsigned int i=0; i<entries_.size(); ++i){
    if (entries_.at(i).sample!=SmurfTree::ttbar && entries_.at(i).sample!=SmurfTree::tw) continue;
    printf("%s", SmurfTree::name(entries_.at(i).sample).c_str());
    for (unsigned int j=0; j<4; ++j)
      printf(" \t%5.2f%s%4.2f", entries_.at(i).finalSelection.yield[j], pm, sqrt(entries_.at(i).finalSelection.err2[j]));
    printf(" \t%5.2f%s%4.2f\n", entries_.at(i).finalSelection.total_yield(), pm, entries_.at(i).finalSelection.total_error());
  }
  */
}

void SmurfAnalysis::estimateTopBackground()
{
  const bool binomialError = false;
  printf("\nTop background estimation inputs (nTopTaggedEvents/nTop1BJetEvents/nTop1BJetTopTaggedEvents):\n\b");;    
  printf(" \t     %s      \t      %s      \t      %s      \t      %s      \t    Total     \n", types[0], types[1], types[2], types[3]);
  if (entries_.empty()){
    printf("No samples are processed.\n");
    return;
  }
  for (unsigned int i=0; i<entries_.size(); ++i){
    if ( !isRelevantToTopBkg(entries_.at(i).sample) ) continue;
    printf("%s", SmurfTree::name(entries_.at(i).sample).c_str());
    for (unsigned int j=0; j<4; ++j){
      printf(" \t%4.1f/%4.1f/%4.1f", 
	     entries_.at(i).nTopTaggedEvents.yield[j], 
	     entries_.at(i).nTop1BJetEvents.yield[j],
	     entries_.at(i).nTop1BJetTopTaggedEvents.yield[j]);
    }
    printf(" \t%4.1f/%4.1f/%4.1f", 
	   entries_.at(i).nTopTaggedEvents.total_yield(), 
	   entries_.at(i).nTop1BJetEvents.total_yield(), 
	   entries_.at(i).nTop1BJetTopTaggedEvents.total_yield());
    double A    = entries_.at(i).nTop1BJetTopTaggedEvents.total_yield();
    double errA = entries_.at(i).nTop1BJetTopTaggedEvents.total_error();
    double B    = entries_.at(i).nTop1BJetEvents.total_yield();
    double errB = entries_.at(i).nTop1BJetEvents.total_error();
    double C    = entries_.at(i).nTopTaggedEvents.total_yield();
    double errC = entries_.at(i).nTopTaggedEvents.total_error();
    
    // compute efficiencies
    double eff_1b = B>0 ? A/B : 0;
    double eff_1b_err(0);
    if ( binomialError ){
      eff_1b_err = sqrt(eff_1b*(1-eff_1b)/B);
    } else {
      eff_1b_err = eff_1b>0 ? eff_1b * sqrt(pow(errA/A,2)+pow(errB/B,2)): 0;
    }
    double eff = 1-pow(1-eff_1b,2);
    double eff_err = 2*(1-eff_1b)*eff_1b_err;

    double scale = eff > 0 ? (1-eff)/eff : 0;
    printf("\ttag eff: %0.3f+/-%0.3f \ttop bkg: %0.2f+/-%0.2f\n", 
	   eff, eff_err, C*scale, sqrt(C)*scale);
  }  

  //
  // MC prediction
  //
  printf("\nTop yield MC estimation:\n\b");;    
  printf(" \t     %s      \t      %s      \t      %s      \t      %s      \t    Total     \n", types[0], types[1], types[2], types[3]);

  for (unsigned int i=0; i<entries_.size(); ++i){
    if (entries_.at(i).sample!=SmurfTree::ttbar && entries_.at(i).sample!=SmurfTree::tw && entries_.at(i).sample!=SmurfTree::qqww) continue;
    printf("%s", SmurfTree::name(entries_.at(i).sample).c_str());
    for (unsigned int j=0; j<4; ++j)
      printf(" \t%5.2f%s%4.2f", entries_.at(i).finalSelection.yield[j], pm, sqrt(entries_.at(i).finalSelection.err2[j]));
    printf(" \t%5.2f%s%4.2f\n", entries_.at(i).finalSelection.total_yield(), pm, entries_.at(i).finalSelection.total_error());
  }
}

void SmurfAnalysis::showYields(double lumi, bool tex_format)
{
  printf("\nEvent yields:\n\b");;    
  double weight = 1;
  if (lumi>0) weight = lumi/lumi_;
  for (unsigned int i=0; i<entries_.size(); ++i){
    printf("%s", SmurfTree::name(entries_.at(i).sample).c_str());
    for (unsigned int j=0; j<4; ++j)
      printf(" \t%5.2f%s%4.2f", entries_.at(i).finalSelection.yield[j]*weight, pm, sqrt(entries_.at(i).finalSelection.err2[j])*weight);
    printf(" \t%5.2f%s%4.2f\n", entries_.at(i).finalSelection.total_yield()*weight, pm, entries_.at(i).finalSelection.total_error()*weight);
  }
}

void SmurfAnalysis::estimateFakeBackground(bool verbose)
{
  if (entries_.empty()){
    printf("No samples are processed.\n");
    return;
  }

  //
  // Count Pass-Failed pairs
  //
  TFile* debug = TFile::Open("debug.root","RECREATE");
  assert(debug);

  if (verbose){
    printf("\nPass-Fail selectionyields :\n\b");;    
    printf(" \t  EleFake-%s  \t  EleFake-%s  \t  EleFake-%s ", types[1], types[2], types[3]);
    printf(" \t  MuFake-%s  \t  MuFake-%s  \t  MuFake-%s  \n", types[0], types[1], types[2]);
    for (unsigned int i=0; i<entries_.size(); ++i){
      if ( !isRelevantToFakeBkg(entries_.at(i).sample) ) continue;
      printf("%s", SmurfTree::name(entries_.at(i).sample).c_str());
      for (unsigned int j=0; j<4; ++j){
	double err_el(0);
	double yield_el = entries_.at(i).hist_el[j]->IntegralAndError(0,nEtaPoints,0,nPtPoints,err_el);
	if (j!=0) {
	  printf(" \t%5.2f%s%4.2f", yield_el, pm, err_el);
	  entries_.at(i).hist_el[j]->Write();
	}
      }
      for (unsigned int j=0; j<4; ++j){
	double err_mu(0);
	double yield_mu = entries_.at(i).hist_mu[j]->IntegralAndError(0,nEtaPoints,0,nPtPoints,err_mu);
	if (j!=3) {
	  printf(" \t%5.2f%s%4.2f", yield_mu, pm, err_mu);
	  entries_.at(i).hist_mu[j]->Write();
	}
      }
      printf("\n");
    }
  }
  debug->Close();

  //
  // Apply fake rates 
  //

  printf("\nJet induced fake background estimation :\n\b");;    
  printf("      \t  EleFake-%s       \t  EleFake-%s       \t  EleFake-%s ", types[1], types[2], types[3]);
  printf("      \t  MuFake-%s       \t  MuFake-%s       \t  MuFake-%s       \t  Total-El       \t  Total-Mu       \t  Total\n", types[0], types[1], types[2]);
	 
  for (unsigned int i=0; i<entries_.size(); ++i){
    if ( !isRelevantToFakeBkg(entries_.at(i).sample) ) continue;
    printf("%s", SmurfTree::name(entries_.at(i).sample).c_str());
    TH2F* hist_el(0);
    TH2F* hist_mu(0);
    double err_el[3];
    double yield_el[3];
    double err_mu[3];
    double yield_mu[3];
    for (unsigned int j=0; j<4; ++j){
      getIntegral(entries_.at(i).hist_el[j], yield_el[0], err_el[0], elFakeRate0_);
      getIntegral(entries_.at(i).hist_el[j], yield_el[1], err_el[1], elFakeRate1_);
      getIntegral(entries_.at(i).hist_el[j], yield_el[2], err_el[2], elFakeRate2_);
      if (j!=0) printf(" \t%4.1f%s%3.1f+%3.1f-%3.1f", yield_el[1], pm, err_el[1],
		       std::max(yield_el[2]-yield_el[1],yield_el[0]-yield_el[1]),
		       fabs(std::min(yield_el[2]-yield_el[1],yield_el[0]-yield_el[1])) );
      if (!hist_el){
	hist_el = (TH2F*)entries_.at(i).hist_el[j]->Clone();
	hist_el->SetDirectory(0);
      } else {
	hist_el->Add(entries_.at(i).hist_el[j]);
      }
    }
    for (unsigned int j=0; j<4; ++j){
      getIntegral(entries_.at(i).hist_mu[j], yield_mu[0], err_mu[0], muFakeRate0_);
      getIntegral(entries_.at(i).hist_mu[j], yield_mu[1], err_mu[1], muFakeRate1_);
      getIntegral(entries_.at(i).hist_mu[j], yield_mu[2], err_mu[2], muFakeRate2_);
      if (j!=3) printf(" \t%4.1f%s%3.1f+%3.1f-%3.1f", yield_mu[1], pm, err_mu[1],
		       std::max(yield_mu[2]-yield_mu[1],yield_mu[0]-yield_mu[1]),
		       fabs(std::min(yield_mu[2]-yield_mu[1],yield_mu[0]-yield_mu[1])) );
      if (!hist_mu){
	hist_mu = (TH2F*)entries_.at(i).hist_mu[j]->Clone();
	hist_mu->SetDirectory(0);
      } else {
	hist_mu->Add(entries_.at(i).hist_mu[j]);
      }
    }
    getIntegral(hist_el, yield_el[0], err_el[0], elFakeRate0_);
    getIntegral(hist_el, yield_el[1], err_el[1], elFakeRate1_);
    getIntegral(hist_el, yield_el[2], err_el[2], elFakeRate2_);
    getIntegral(hist_mu, yield_mu[0], err_mu[0], muFakeRate0_);
    getIntegral(hist_mu, yield_mu[1], err_mu[1], muFakeRate1_);
    getIntegral(hist_mu, yield_mu[2], err_mu[2], muFakeRate2_);
    double err_el_p = std::max(yield_el[2]-yield_el[1],yield_el[0]-yield_el[1]);
    double err_el_m = fabs(std::min(yield_el[2]-yield_el[1],yield_el[0]-yield_el[1]));
    double err_mu_p = std::max(yield_mu[2]-yield_mu[1],yield_mu[0]-yield_mu[1]);
    double err_mu_m = fabs(std::min(yield_mu[2]-yield_mu[1],yield_mu[0]-yield_mu[1]));
    printf(" \t%4.1f%s%3.1f+%3.1f-%3.1f", yield_el[1], pm, err_el[1], err_el_p, err_el_m);
    printf(" \t%4.1f%s%3.1f+%3.1f-%3.1f", yield_mu[1], pm, err_mu[1], err_mu_p, err_mu_m);
    printf(" \t%4.1f%s%3.1f+%3.1f-%3.1f", yield_el[1]+yield_mu[1], pm, sqrt(err_el[1]*err_el[1]+err_mu[1]*err_mu[1]),
	   sqrt(err_el_p*err_el_p+err_mu_p*err_mu_p),sqrt(err_el_m*err_el_m+err_mu_m*err_mu_m));
    printf("\n");
  }

  //
  // MC prediction
  //
  printf("\nWjets yield MC estimation:\n\b");;    
  printf(" \t     %s      \t      %s      \t      %s      \t      %s      \t    Total     \n", types[0], types[1], types[2], types[3]);
  for (unsigned int i=0; i<entries_.size(); ++i){
    if (entries_.at(i).sample!=SmurfTree::wjets) continue;
    printf("%s", SmurfTree::name(entries_.at(i).sample).c_str());
    for (unsigned int j=0; j<4; ++j)
      printf(" \t%5.2f%s%4.2f", entries_.at(i).finalSelection.yield[j], pm, sqrt(entries_.at(i).finalSelection.err2[j]));
    printf(" \t%5.2f%s%4.2f\n", entries_.at(i).finalSelection.total_yield(), pm, entries_.at(i).finalSelection.total_error());
  }

  //
  // SS closure test
  //

  printf("\nJet induced fake background estimation (same-sign):\n\b");;    
  printf("      \t  EleFake-%s       \t  EleFake-%s       \t  EleFake-%s ", types[1], types[2], types[3]);
  printf("      \t  MuFake-%s       \t  MuFake-%s       \t  MuFake-%s       \t  Total-El       \t  Total-Mu       \t  Total\n", types[0], types[1], types[2]);
	 
  for (unsigned int i=0; i<entries_.size(); ++i){
    if ( !isRelevantToFakeBkg(entries_.at(i).sample) && 
	 !isRelevantToBeProcessed(entries_.at(i).sample) ) continue;
    printf("%s", SmurfTree::name(entries_.at(i).sample).c_str());
    TH2F* hist_el(0);
    TH2F* hist_mu(0);
    double err_el[3];
    double yield_el[3];
    double err_mu[3];
    double yield_mu[3];
    for (unsigned int j=0; j<4; ++j){
      getIntegral(entries_.at(i).hist_el_ss[j], yield_el[0], err_el[0], elFakeRate0_);
      getIntegral(entries_.at(i).hist_el_ss[j], yield_el[1], err_el[1], elFakeRate1_);
      getIntegral(entries_.at(i).hist_el_ss[j], yield_el[2], err_el[2], elFakeRate2_);
      if (j!=0) printf(" \t%4.1f%s%3.1f+%3.1f-%3.1f", yield_el[1], pm, err_el[1],
		       std::max(yield_el[2]-yield_el[1],yield_el[0]-yield_el[1]),
		       fabs(std::min(yield_el[2]-yield_el[1],yield_el[0]-yield_el[1])) );
      if (!hist_el){
	hist_el = (TH2F*)entries_.at(i).hist_el_ss[j]->Clone();
	hist_el->SetDirectory(0);
      } else {
	hist_el->Add(entries_.at(i).hist_el_ss[j]);
      }
    }
    for (unsigned int j=0; j<4; ++j){
      getIntegral(entries_.at(i).hist_mu_ss[j], yield_mu[0], err_mu[0], muFakeRate0_);
      getIntegral(entries_.at(i).hist_mu_ss[j], yield_mu[1], err_mu[1], muFakeRate1_);
      getIntegral(entries_.at(i).hist_mu_ss[j], yield_mu[2], err_mu[2], muFakeRate2_);
      if (j!=3) printf(" \t%4.1f%s%3.1f+%3.1f-%3.1f", yield_mu[1], pm, err_mu[1],
		       std::max(yield_mu[2]-yield_mu[1],yield_mu[0]-yield_mu[1]),
		       fabs(std::min(yield_mu[2]-yield_mu[1],yield_mu[0]-yield_mu[1])) );
      if (!hist_mu){
	hist_mu = (TH2F*)entries_.at(i).hist_mu_ss[j]->Clone();
	hist_mu->SetDirectory(0);
      } else {
	hist_mu->Add(entries_.at(i).hist_mu_ss[j]);
      }
    }
    getIntegral(hist_el, yield_el[0], err_el[0], elFakeRate0_);
    getIntegral(hist_el, yield_el[1], err_el[1], elFakeRate1_);
    getIntegral(hist_el, yield_el[2], err_el[2], elFakeRate2_);
    getIntegral(hist_mu, yield_mu[0], err_mu[0], muFakeRate0_);
    getIntegral(hist_mu, yield_mu[1], err_mu[1], muFakeRate1_);
    getIntegral(hist_mu, yield_mu[2], err_mu[2], muFakeRate2_);
    double err_el_p = std::max(yield_el[2]-yield_el[1],yield_el[0]-yield_el[1]);
    double err_el_m = fabs(std::min(yield_el[2]-yield_el[1],yield_el[0]-yield_el[1]));
    double err_mu_p = std::max(yield_mu[2]-yield_mu[1],yield_mu[0]-yield_mu[1]);
    double err_mu_m = fabs(std::min(yield_mu[2]-yield_mu[1],yield_mu[0]-yield_mu[1]));
    printf(" \t%4.1f%s%3.1f+%3.1f-%3.1f", yield_el[1], pm, err_el[1], err_el_p, err_el_m);
    printf(" \t%4.1f%s%3.1f+%3.1f-%3.1f", yield_mu[1], pm, err_mu[1], err_mu_p, err_mu_m);
    printf(" \t%4.1f%s%3.1f+%3.1f-%3.1f", yield_el[1]+yield_mu[1], pm, sqrt(err_el[1]*err_el[1]+err_mu[1]*err_mu[1]),
	   sqrt(err_el_p*err_el_p+err_mu_p*err_mu_p),sqrt(err_el_m*err_el_m+err_mu_m*err_mu_m));
    printf("\n");
  }

  //
  // SS yeilds
  //
  printf("\nSame-sign yeilds:\n\b");;    
  printf(" \t     %s      \t      %s      \t      %s      \t      %s      \t    Total     \n", types[0], types[1], types[2], types[3]);

  for (unsigned int i=0; i<entries_.size(); ++i){
    if (!isRelevantToBeProcessed(entries_.at(i).sample) && 
	entries_.at(i).sample!=SmurfTree::data) continue;
    printf("%s", SmurfTree::name(entries_.at(i).sample).c_str());
    double t(0);
    double t2(0);
    for (unsigned int j=0; j<4; ++j){
      double err(0);
      double yield = entries_.at(i).hist_final_ss[j]->IntegralAndError(0,nEtaPoints,0,nPtPoints,err);
      printf(" \t%5.2f%s%4.2f", yield, pm, err);
      t+=yield;
      t2+=err*err;
    }
    printf(" \t%5.2f%s%4.2f", t, pm, sqrt(t2));
    printf("\n");
  }
//   for (unsigned int i=0; i<entries_.size(); ++i){
//     if (entries_.at(i).sample!=SmurfTree::wjets) continue;
//     printf("%s", SmurfTree::name(entries_.at(i).sample).c_str());
//     for (unsigned int j=0; j<4; ++j)
//       printf(" \t%5.2f%s%4.2f", entries_.at(i).finalSelection.yield[j], pm, sqrt(entries_.at(i).finalSelection.err2[j]));
//     printf(" \t%5.2f%s%4.2f\n", entries_.at(i).finalSelection.total_yield(), pm, entries_.at(i).finalSelection.total_error());
//   }

}

void SmurfAnalysis::addSample(SmurfTree::DataType sample)
{
  if ( !isRelevantToBeProcessed(sample) ) return;
  FileStat_t buf;
  std::string file_name = std::string(dir_)+"/"+SmurfTree::name(sample)+".root";
  // check if file exists
  if (gSystem->GetPathInfo(file_name.c_str(), buf)) return;
  // std::cout << "processing file " << file_name << std::endl;
  Entry entry;
  entry.sample = sample;
  for (unsigned int i=0; i<4; ++i){
    entry.hist_el[i] = new TH2F(Form("%s_el_%s", SmurfTree::name(sample).c_str(), types[i]), "",
				nEtaPoints-1, etaBins, nPtPoints-1, ptBins);
    entry.hist_el[i]->Sumw2();
    entry.hist_el[i]->SetDirectory(0);

    entry.hist_el_ss[i] = new TH2F(Form("%s_el_ss_%s", SmurfTree::name(sample).c_str(), types[i]), "",
				   nEtaPoints-1, etaBins, nPtPoints-1, ptBins);
    entry.hist_el_ss[i]->Sumw2();
    entry.hist_el_ss[i]->SetDirectory(0);

    entry.hist_mu[i] = new TH2F(Form("%s_mu_%s", SmurfTree::name(sample).c_str(), types[i]), "",
				nEtaPoints-1, etaBins, nPtPoints-1, ptBins);
    entry.hist_mu[i]->Sumw2();
    entry.hist_mu[i]->SetDirectory(0);

    entry.hist_mu_ss[i] = new TH2F(Form("%s_mu_ss_%s", SmurfTree::name(sample).c_str(), types[i]), "",
				   nEtaPoints-1, etaBins, nPtPoints-1, ptBins);
    entry.hist_mu_ss[i]->Sumw2();
    entry.hist_mu_ss[i]->SetDirectory(0);
    
    entry.hist_final[i] = new TH2F(Form("%s_final_%s", SmurfTree::name(sample).c_str(), types[i]), "",
				   nEtaPoints-1, etaBins, nPtPoints-1, ptBins);
    entry.hist_final[i]->Sumw2();
    entry.hist_final[i]->SetDirectory(0);

    entry.hist_final_ss[i] = new TH2F(Form("%s_final_ss_%s", SmurfTree::name(sample).c_str(), types[i]), "",
				      nEtaPoints-1, etaBins, nPtPoints-1, ptBins);
    entry.hist_final_ss[i]->Sumw2();
    entry.hist_final_ss[i]->SetDirectory(0);
  }
  
  const double max_eta = etaBins[nEtaPoints-1];
  const double max_pt  = ptBins[nPtPoints-1];

  SmurfTree tree;
  tree.LoadTree(file_name.c_str());
  tree.InitTree();
  Long64_t nentries = tree.tree_->GetEntries();
  unsigned int nEvents(0);
  for (Long64_t i = 0; i < nentries; i++){
    tree.tree_->GetEntry(i);
    if ( sample == SmurfTree::data && json_ )
      if( !goodrun(tree.run_, tree.lumi_) ) continue;
    double weight = 1.0;
    if ( sample != SmurfTree::data ) {
      weight = tree.scale1fb_*lumi_;
      if ( isGGHiggs(sample) ){
	switch (measurement_.type){
	case WW0Jet:
	case HWWCutBased0Jet:
	case HWWMVABased0Jet:
	  weight = weight*1.13;
	  break;
	case WW1Jet:
	case HWWCutBased1Jet:
	case HWWMVABased1Jet:
	  weight = weight*0.86;
	  break;
	case WW2Jets:
	case HWWCutBased2Jets:
	case HWWMVABased2Jets:
	  weight = weight*0.73;
	  break;
	}
      }
    }

    // Control samples
    if ( (tree.cuts_ & cut_top1B) == cut_top1B ) {
      entry.nTop1BJetEvents.add(tree.type_,weight);
      if ( (tree.cuts_ & SmurfTree::TopTagNotInJets) > 0 ) entry.nTop1BJetTopTaggedEvents.add(tree.type_,weight);
    }

    // Proper event selection
    if ( !passedExtraCuts(tree) ) continue;

    // overflow goes into last bin
    double eta1 = fabs(tree.lep1_.eta());
    if (eta1 > max_eta) eta1 = max_eta-.1;
    double pt1 = tree.lep1_.pt();
    if (pt1 > max_pt) pt1 = max_pt-.1;
    double eta2 = fabs(tree.lep2_.eta());
    if (eta2 > max_eta) eta2 = max_eta-.1;
    double pt2 = tree.lep2_.pt();
    if (pt2 > max_pt) pt2 = max_pt-.1;
    // first electron is fake
    if ( (tree.cuts_ & cut_el_1)==cut_el_1 && (tree.cuts_ & cut_final)!=cut_final ){
      if ( tree.lq1_*tree.lq2_<0 )
	entry.hist_el[tree.type_]->Fill(eta1,pt1,weight);
      else
	entry.hist_el_ss[tree.type_]->Fill(eta1,pt1,weight);
    }	
    // second electron is fake
    if ( (tree.cuts_ & cut_el_2)==cut_el_2 && (tree.cuts_ & cut_final)!=cut_final ){
      if ( tree.lq1_*tree.lq2_<0 )
	entry.hist_el[tree.type_]->Fill(eta2,pt2,weight);
      else
	entry.hist_el_ss[tree.type_]->Fill(eta2,pt2,weight);
      // cout << "Second electron is fake. Type: " << tree.type_ << "\t integral: " << entry.hist_el[tree.type_]->Integral() << endl;
    }
    // first muon is fake
    if ( (tree.cuts_ & cut_mu_1)==cut_mu_1 && (tree.cuts_ & cut_final)!=cut_final ){
      if ( tree.lq1_*tree.lq2_<0 )
	entry.hist_mu[tree.type_]->Fill(eta1,pt1,weight);
      else
	entry.hist_mu_ss[tree.type_]->Fill(eta2,pt2,weight);
    }	
    // second muon is fake
    if ( (tree.cuts_ & cut_mu_2)==cut_mu_2 && (tree.cuts_ & cut_final)!=cut_final ){
      if ( tree.lq1_*tree.lq2_<0 )
	entry.hist_mu[tree.type_]->Fill(eta2,pt2,weight);
      else
	entry.hist_mu_ss[tree.type_]->Fill(eta2,pt2,weight);
      // cout << "Second muon is fake. Type: " << tree.type_ << "\t integral: " << entry.hist_mu[tree.type_]->Integral() << endl;
    }	
    if ( (tree.cuts_ & cut_final)==cut_final ) {
      if ( tree.lq1_*tree.lq2_<0 ){
	entry.hist_final[tree.type_]->Fill(eta1,pt1,weight);
	entry.finalSelection.add(tree.type_,weight);
	++nEvents;
      } else {
	entry.hist_final_ss[tree.type_]->Fill(eta1,pt1,weight);
      }	
    }
    if ( (tree.cuts_ & cut_topTagged) == cut_topTagged ) entry.nTopTaggedEvents.add(tree.type_,weight);
    if ( (tree.cuts_ & cut_final_nomass) == cut_final_nomass ){
      if (std::min(tree.pmet_,tree.pTrackMet_)>35){ // use the same MET cut in all final states
	if ( fabs(tree.dilep_.mass() - 91.1876) < 15 ){
	  entry.nZpeak.add(tree.type_,weight);
	  // std::cout << "Found event in Z window (run/lumi/event): " 
	  // << tree.run_ << " / " << tree.lumi_ << " / " << tree.event_ << std::endl; 
	}
      }
    }
  }
  // std::cout << "Events passed full selection: " << nEvents << std::endl;
  // fill histograms
  entries_.push_back(entry);
}

// ==================================================================================== //

bool SmurfAnalysis::MCheckBinLimits(const TArrayD* h1Array, const TArrayD* h2Array)
{
  Int_t fN = h1Array->fN;
  if ( fN != 0 ) {
    if ( h2Array->fN != fN ) {
      printf("DifferentBinLimits\n");
      return false;
    }
    else {
      for ( int i = 0; i < fN; ++i ) {
	if ( fabs(h1Array->GetAt(i)-h2Array->GetAt(i)) > 1E-3 ) {
	  printf("DifferentBinLimits\n");
	  return false;
	}
      }
    }
  }

  return true;
}

bool SmurfAnalysis::CheckConsistency(const TH1* h1, const TH1* h2)
{
  // Check histogram compatibility   Int_t nbinsx = h1->GetNbinsX();
  Int_t nbinsx = h1->GetNbinsX();
  Int_t nbinsy = h1->GetNbinsY();
  Int_t nbinsz = h1->GetNbinsZ();

  // Check whether the histograms have the same number of bins.
  if (nbinsx != h2->GetNbinsX() || nbinsy != h2->GetNbinsY() || nbinsz != h2->GetNbinsZ()) {
    printf("DifferentNumberOfBins\n");
    return false;
  }
  // Check that the axis limits of the histograms are the same
  if (h1->GetXaxis()->GetXmin() != h2->GetXaxis()->GetXmin() ||
      h1->GetXaxis()->GetXmax() != h2->GetXaxis()->GetXmax() ||
      h1->GetYaxis()->GetXmin() != h2->GetYaxis()->GetXmin() ||
      h1->GetYaxis()->GetXmax() != h2->GetYaxis()->GetXmax() ||
      h1->GetZaxis()->GetXmin() != h2->GetZaxis()->GetXmin() ||
      h1->GetZaxis()->GetXmax() != h2->GetZaxis()->GetXmax()) {
    printf("DifferentAxisLimits\n");
    return false;
  }

  bool ret = true;

  ret &= MCheckBinLimits(h1->GetXaxis()->GetXbins(), h2->GetXaxis()->GetXbins());
  ret &= MCheckBinLimits(h1->GetYaxis()->GetXbins(), h2->GetYaxis()->GetXbins());
  ret &= MCheckBinLimits(h1->GetZaxis()->GetXbins(), h2->GetZaxis()->GetXbins());

  return ret;
}

void SmurfAnalysis::getIntegral(const TH2F* h, double& integral, double& error, const TH2F* fakeRate, int etaBin, int ptBin)
{
  integral = 0;
  error = 0;

  // first check if histogram have the same binning
  if ( fakeRate && ! CheckConsistency(h,fakeRate) ){
    printf("Error: histograms are not consistent. Abort");
    return;
  }
  if ( fakeRate == 0 ){
    for ( Int_t x=0; x<h->GetNbinsX(); ++x )
      for ( Int_t y=0; y<h->GetNbinsY(); ++y ){
	integral += h->GetBinContent(x+1,y+1);
	error += h->GetBinError(x+1,y+1)*h->GetBinError(x+1,y+1);
      }
  } else {
    for ( Int_t x=0; x<h->GetNbinsX(); ++x ){
      if (etaBin>=0 && x+1!= etaBin ) continue;
      for ( Int_t y=0; y<h->GetNbinsY(); ++y ){
	if (ptBin>=0 && y+1!= ptBin ) continue;
	double bin    = h->GetBinContent(x+1,y+1);
	double binerr = h->GetBinError(x+1,y+1);
	double fr     = fakeRate->GetBinContent(x+1,y+1);
	double frerr  = fakeRate->GetBinError(x+1,y+1);
	integral += bin*fr/(1-fr);
	error += pow(binerr*fr/(1-fr),2) + pow(bin*frerr/(1-fr)/(1-fr), 2);
      }
    }
  }
  error = sqrt(error);
}

