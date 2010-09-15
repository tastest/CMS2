#ifndef WW_doAnalysis_h
#define WW_doAnalysis_h
#include "Math/LorentzVector.h"
#include "Rtypes.h"
#include <vector>
#include <set>
#include "wwtypes.h"
#include "TChain.h"
#include <fstream>

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

//
// Electron Id
//

bool   ww_elBase(unsigned int i);
bool   ww_elId(unsigned int i);
bool   ww_eld0(unsigned int i);
bool   ww_eld0PV(unsigned int i);
bool   ww_elIso(unsigned int i);
double ww_elIsoVal(unsigned int i);

// combined analysis selectors
bool goodElectronWithoutIsolation(unsigned int i);
bool goodElectronIsolated(unsigned int i);

//
// Muon Id
//

bool   ww_muBase(unsigned int i);
bool   ww_muId(unsigned int i);
bool   ww_muIso(unsigned int i);
double ww_muIsoVal(unsigned int i);
bool   ww_mud0(unsigned int i);
bool   ww_mud0PV(unsigned int i);

unsigned int numberOfSoftMuons(int i_hyp, bool nonisolated);

// combined analysis selectors
bool goodMuonWithoutIsolation(unsigned int i);
bool goodMuonIsolated(unsigned int i);

//
// trigger
//

bool passedTriggerRequirements(HypTypeInNtuples type);

//
// Met
//

// the following two functions set met values used in 
// all other selectors and functions
double metValue();
double metPhiValue();

bool ww2009_met(unsigned int i_hyp);

// analysis MET requirements
bool passedMetRequirements(unsigned int i_hyp);

double projectedMet(unsigned int i_hyp);

bool metBalance (unsigned int i_hyp);

double nearestDeltaPhi(double Phi, int i_hyp);

//
// Jets
//
enum JetType { CaloJet, jptJet, pfJet, TrkJet, GenJet };
std::vector<LorentzVector> getJets(JetType type, 
				   int i_hyp, 
				   double etThreshold,
				   double maxEta,
				   bool sorted = false);

// analysis jet type is set here.
JetType jetType();

unsigned int numberOfJets(unsigned int i_hyp);

//
// Various Cuts
//

//
// Datasets for fakes
// 


//
// Tools
//

HypTypeInNtuples hypType(unsigned int i_hyp);

struct EventIdentifier {
  unsigned long int run, event, lumi;
  float trks_d0;
  float hyp_lt_pt, hyp_lt_eta, hyp_lt_phi;
  bool operator < (const EventIdentifier &) const;
  bool operator == (const EventIdentifier &) const;
};

unsigned int getDrellYanType();
unsigned int getVVType();
bool isDYee();
bool isDYmm();
bool isDYtt();
bool isWW();
bool isWZ();
bool isZZ();

//
// Not Classified
// 

 
// filter events by process
bool filterByProcess( enum Sample sample );
bool isIdentified( enum Sample sample );

void checkIsolation(int i_hyp, double weight);

class RooDataSet;
void getIsolationSidebandsAfterSelections(int i_hyp, 
					  double weight, 
					  RooDataSet* dataset, 
					  bool passedAllLeptonRequirements);

void find_most_energetic_jets(int i_hyp, double weight);
void hypo (int i_hyp, double kFactor, RooDataSet* dataset = 0); 

RooDataSet* MakeNewDataset(const char* name);

void AddIsoSignalControlSample( int i_hyp, double kFactor, RooDataSet* dataset = 0 );
class TChain;
RooDataSet* ScanChain( TChain* chain, 
		       Sample sample, 
		       double integratedLumi,
		       double xsec,
		       int nProcessedEvents,
		       bool identifyEvents,
		       bool qcdBackground = false);
void SkimChain(TChain* chain);
bool passedSkimSelection();

void ProcessSample( std::string file_pattern, 
		    Sample sample, 
		    double integratedLumi,
		    double xsec,
		    int nProcessedEvents,
		    RooDataSet* output_dataset, 
		    Color_t color, 
		    bool identifyEvents = false,
		    bool qcdBackground = false);
void ProcessSample( std::vector<std::string> file_patterns, 
		    Sample sample,
		    double integratedLumi,
		    double xsec,
		    int nProcessedEvents,
		    RooDataSet* output_dataset, 
		    Color_t color, 
		    bool identifyEvents = false,
		    bool qcdBackground = false);

#endif
