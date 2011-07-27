#ifndef DOANALYSIS_H
#define DOANALYSIS_H

#include <vector>
#include <set>
#include "TFile.h"
#include "TTree.h"
#include "TBitSet.hh"
#include "TH1F.h"
#include "TH2F.h"
#include "TROOT.h"
#include "TRint.h"
#include "TChain.h"
#include "Math/LorentzVector.h"
#include "TLorentzVector.h"
#include "Math/VectorUtil.h"


using namespace std;

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

inline double fast_sign(double f) {
    if (f > 0) return 1;
    return (f == 0) ? 0 : -1;
}


enum hyp_selection {
  kcut_Trigger,
  kcut_OS,
  kcut_GoodVertex,
  kcut_LT,
  kcut_LL,
  kcut_Mll,
  kcut_zsel,
  kcut_met,
  kcut_jetveto,
  kcut_softmuonveto,
  kcut_extraleptonveto,
  kcut_toptag,
  kNCuts
};



enum dilepton_flav {
  mumu,
  emu,
  ee,
  all,
  kNdilepFlav
};

enum HypTypeInNtuples {
  MuMu, 
  MuEl , 
  ElMu, 
  ElEl
};

enum jet_types {
  jptJet,
  pfJet,
  GenJet,
  CaloJet,
  TrkJet,
  kNJetTypes
};


struct EventIdentifier {
  unsigned long int run, event, lumi;
  float trks_d0;
  float hyp_lt_pt, hyp_lt_eta, hyp_lt_phi;
  bool operator < (const EventIdentifier &) const;
  bool operator == (const EventIdentifier &) const;
};


// main functions
void ScanChain(const char* process, 
	       TChain *chain, 
	       TFile *utilFile_,  
	       int nEvents = -1, 
	       double IntLumi=100, 
	       double Xsect=1.0, 
	       int nProcessedEvents=-1, 
	       std::string skimFilePrefix="", 
	       bool realData=false, 
	       bool identifyEvents=false);

void ProcessSample(const char* process, 
		   std::vector<std::string> file_patterns,
		   TFile *utilFile_,  
		   int nEvents = -1, 
		   double IntLumi=100, 
		   double Xsect=1.0, 
		   int nProcessedEvents=-1, 
		   std::string skimFilePrefix="",
		   bool realData=false, 
		   bool identifyEvents=false);

void ProcessSample(const char* process, 
		   std::string file_pattern, 
		   TFile *utilFile_, 
		   int nEvents = -1, 
		   double IntLumi=100, 
		   double Xsect=1.0, 
		   int nProcessedEvents=-1, 
		   std::string skimFilePrefix="", 
		   bool realData=false, 
		   bool identifyEvents=false);

//Utils

bool is_duplicate (const EventIdentifier &id);
int ApplyEventSelection(unsigned int, bool realData);
void CalculateFakeRateProb();
void progress( int nEventsTotal, int nEventsChain );
int getHypothesisType( int NTtype );

double nearestDeltaPhi(double Phi, int i_hyp);
bool comparePt(LorentzVector lv1, LorentzVector lv2);
std::vector<LorentzVector> getJets(int type, int i_hyp, double etThreshold, double etaMax, bool sortJets, bool btag);
bool isGoodVertex(size_t ivtx);
unsigned int nGoodVertex();
double dzPV(const LorentzVector& vtx, const LorentzVector& p4, const LorentzVector& pv);
double mt(double pt1, double pt2, double dphi);

// met
double metValue();
double metPhiValue();
double projectedMet(unsigned int, double, double);
bool passedMetRequirements(unsigned int i_hyp);

//trigger
bool  passedTrigger(TString trigName);
bool passedTriggerRequirements();
bool defaultBTag(int type, unsigned int iJet); 
double BTag(int type, unsigned int iJet); 

//muon ID
bool goodMuonIsolated(unsigned int i);
bool fakableMuon(unsigned int i);
bool ww_muBase(unsigned int index);
bool ww_mud0PV(unsigned int index);
double ww_mud0ValuePV(unsigned int index);
bool ww_mudZPV(unsigned int index);
bool ww_muId(unsigned int index);
double ww_muIsoVal(unsigned int index);
bool ww_muIso(unsigned int index);

//electron ID
bool goodElectronIsolated(unsigned int i);
bool fakableElectron(unsigned int i);
bool ww_elBase(unsigned int index);
bool ww_eld0PV(unsigned int index);
bool ww_eldZPV(unsigned int index);
bool ww_elId(unsigned int index);
double ww_elIsoVal(unsigned int index);
bool ww_elIso(unsigned int index);


unsigned int numberOfSoftMuons(int i_hyp, bool nonisolated, const std::vector<LorentzVector>& = std::vector<LorentzVector>());
unsigned int numberOfExtraLeptons(int i_hyp, double minPt);
bool toptag(int type, int i_hyp, double minPt);






// Functions to fill the effiiency / fakerate related histograms
void fillEffHist(const char* process, double weight);
void fillKtHist(const char* process, double weight);
void fillFOHist();
void findClosestEleFO(LorentzVector v_parton, double& minDR, int& idx_minDR);
void findClosestMuFO(LorentzVector v_parton, double& minDR, int& idx_minDR);
void findClosestGenPs(LorentzVector v_parton, double& minDR, int& idx_minDR);

// Utility Functions
bool isIdentified(const char* process);
void getEff(double & numer, double & denom, double & eff, double & efferr );

// 
// This function assumes the two 1D histograms are filled with event counts rather 
// than weighted event counts
// 
void InitMCUtilHist(const char* process, TFile *utilFile_);
void fill1DEffHist(TH1F* hist_numer, TH1F* hist_denom, TH1F* hist_eff);
void fill2DEffHist(TH2F* hist_numer, TH2F* hist_denom, TH2F* hist_eff);
void saveMCUtilOutput(const char* process, TFile *utilFile_);

#endif
