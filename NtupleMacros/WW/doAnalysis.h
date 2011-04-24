#ifndef WW_doAnalysis_h
#define WW_doAnalysis_h
#include "Math/LorentzVector.h"
#include "Rtypes.h"
#include <vector>
#include <set>
#include "wwtypes.h"
#include "TChain.h"
#include <fstream>
#include <vector>
#include "../../../Smurf/Core/SmurfTree.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;
typedef UInt_t wwcuts_t; // 32 bits only!

//
// Electron Id
//

bool   ww_elBase(unsigned int i);
bool   ww_elId(unsigned int i);
bool   ww_eld0(unsigned int i);
bool   ww_eld0PV(unsigned int i);
bool   ww_eldZPV(unsigned int i);
bool   ww_elIso(unsigned int i);
double ww_elIsoVal(unsigned int i);

// combined analysis selectors
bool goodElectronWithoutIsolation(unsigned int i);
bool goodElectronIsolated(unsigned int i);
enum EleFOTypes { EleFOV1, EleFOV2, EleFOV3, EleFOV4};
bool fakableElectron(unsigned int i,EleFOTypes);

//
// Muon Id
//

bool   ww_muBase(unsigned int i);
bool   ww_muId(unsigned int i);
bool   ww_muIso(unsigned int i);
double ww_muIsoVal(unsigned int i);
bool   ww_mud0(unsigned int i);
bool   ww_mud0PV(unsigned int i);
bool   ww_mudZPV(unsigned int i);

unsigned int numberOfSoftMuons(int i_hyp, bool nonisolated,
			       const std::vector<LorentzVector>& = std::vector<LorentzVector>());

// combined analysis selectors
bool goodMuonWithoutIsolation(unsigned int i);
bool goodMuonIsolated(unsigned int i);
enum MuFOTypes { MuFOV1, MuFOV2 };
bool fakableMuon(unsigned int i, MuFOTypes);

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

double projectedMet(unsigned int i_hyp, double met, double phi);

bool metBalance (unsigned int i_hyp);

double nearestDeltaPhi(double Phi, int i_hyp);

//
// Jets
//
enum WWJetType { CaloJet, jptJet, pfJet, TrkJet, GenJet };
std::vector<LorentzVector> getJets(WWJetType type, 
				   int i_hyp, 
				   double etThreshold,
				   double maxEta,
				   bool sorted = false,
				   bool btag = false);

// analysis jet type is set here.
WWJetType jetType();

unsigned int numberOfJets(unsigned int i_hyp);

//
// Various Cuts
//

bool isGoodVertex(size_t ivtx);

//
// Datasets for fakes
// 


//
// Tools
//
bool defaultBTag(WWJetType type, unsigned int iJet);

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
// N-1
// return true if the cuts to apply - the cuts to remove were passed in the cuts that "passed"
bool CheckCutsNM1(wwcuts_t apply, wwcuts_t remove, wwcuts_t passed);
// Simple check if the desired cuts to apply are set in the cuts that "passed"
bool CheckCuts(wwcuts_t apply, wwcuts_t passed);

//
// Not Classified
// 

 
// filter events by process
bool filterByProcess( enum SmurfTree::DataType sample );
bool isIdentified( enum SmurfTree::DataType sample );

// void checkIsolation(int i_hyp, double weight);

void find_leading_genjet(int i_hyp, double etaMin, double etaMax, double vetoCone, double & genJetMax);
void find_leading_jptjet(int i_hyp, double etaMin, double etaMax, double vetoCone, double & jptMax, int &jptMaxIndex);
void find_leading_calojet(int i_hyp, double etaMin, double etaMax, double vetoCone, double & caloJetMax);
void find_leading_trkjet(int i_hyp, double etaMin, double etaMax, double vetoCone, double & trkJetMax);
void find_leading_pfjet(int i_hyp, double etaMin, double etaMax, double vetoCone, double & pfJetMax);
void find_most_energetic_jets(int i_hyp, double weight, bool realData, double etaMin, double etaMax);
void getJetResponseFromZBalance(int i_hyp, double weight, bool realData, double etaMin, double etaMax);
void fill_val_plots(int i_hyp, wwcuts_t cuts_passed, double weight);
void fill_dyest_histograms(int i_hyp, float weight);
				
unsigned int bestZHyp();
bool hypo (int i_hyp, double weight, bool zStudy = false, bool realData = false ); 

class TChain;
void ScanChain( TChain* chain, 
	   SmurfTree::DataType sample, 
	   double integratedLumi,
	   double xsec,
	   int nProcessedEvents,
	   bool identifyEvents,
	   bool qcdBackground = false,
	   bool zStudy = false,
	   bool realData = false,
	   TString cms2_json_file = "");
void SkimChain(TChain* chain, bool mergeFiles=false);
bool passedSkimSelection();

void ProcessSample( std::string file_pattern, 
		    SmurfTree::DataType sample, 
		    double integratedLumi,
		    double xsec,
		    int nProcessedEvents,
		    Color_t color, 
		    bool identifyEvents = false,
		    bool qcdBackground = false,
		    bool zStudy = false,
		    bool realData = false,
		    TString cms2_json_file = "");
void ProcessSample( std::vector<std::string> file_patterns, 
		    SmurfTree::DataType sample,
		    double integratedLumi,
		    double xsec,
		    int nProcessedEvents,
		    Color_t color, 
		    bool identifyEvents = false,
		    bool qcdBackground = false,
		    bool zStudy = false,
		    bool realData = false,
		    TString cms2_json_file = "");

#endif
