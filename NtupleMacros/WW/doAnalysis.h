#ifndef WW_doAnalysis_h
#define WW_doAnalysis_h
#include "Math/LorentzVector.h"
#include "Rtypes.h"
#include <vector>
#include <set>
#include "wwtypes.h"
#include "CORE/CMS2.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

//
// Electron Id
//

bool   ww_elId(unsigned int i);
bool   ww_eld0(unsigned int i);
double ww_elIso(unsigned int i);

bool   ww2009_elId(unsigned int i);
bool   ww2009_eld0(unsigned int i);
double ww2009_elIso(unsigned int i);

// combined analysis selectors
bool goodElectronWithoutIsolation(unsigned int i);
bool goodElectronIsolated(unsigned int i);

//
// Muon Id
//

bool   ww_muId(unsigned int i);
bool   ww_mud0(unsigned int i);
double ww_muIso(unsigned int i);

bool   ww2009_muId(unsigned int i);
bool   ww2009_mud0(unsigned int i);
double ww2009_muIso(unsigned int i);

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

struct EventIdentifier {
  unsigned long int run, event, lumi;
  float trks_d0;
  float hyp_lt_pt, hyp_lt_eta, hyp_lt_phi;
  bool operator < (const EventIdentifier &) const;
  bool operator == (const EventIdentifier &) const;
};

struct hypo_monitor{
  std::vector<std::pair<std::string,unsigned int> > counters;
  void count(unsigned int index, const char* name){
    unsigned int current_size = counters.size();
    for ( unsigned int i=current_size; i<=index; ++i ) 
      counters.push_back( std::pair<std::string,unsigned int>("",0) );
    counters[index].first = name;
    counters[index].second++;
  }
  void print(){
    for ( unsigned int i=0; i<counters.size(); ++i ) 
      std::cout << counters[i].first << "\t" << counters[i].second << std::endl;
  }
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
RooDataSet* ScanChain( TChain* chain, Sample sample, bool identifyEvents );

#endif
