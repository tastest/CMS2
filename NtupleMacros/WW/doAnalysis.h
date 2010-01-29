#ifndef WW_doAnalysis_h
#define WW_doAnalysis_h
#include "Math/LorentzVector.h"
#include "Rtypes.h"
#include <vector>
#include <set>
#include "CORE/selections.h"
#include "wwtypes.h"

struct DorkyEventIdentifier {
     // this is a workaround for not having unique event id's in MC
     unsigned long int run, event, lumi;
     float trks_d0;
     float hyp_lt_pt, hyp_lt_eta, hyp_lt_phi;
     bool operator < (const DorkyEventIdentifier &) const;
     bool operator == (const DorkyEventIdentifier &) const;
};

// filter events by process
bool filterByProcess( enum Sample sample );
bool isIdentified( enum Sample sample );

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
