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
typedef std::pair<LorentzVector,unsigned int> JetPair;

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
bool goodElectronTMVA(unsigned int i);
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
bool   ww_mudZPV(unsigned int i, float cut=0.1);

unsigned int numberOfSoftMuons(int i_hyp, bool nonisolated,
			       const std::vector<JetPair>& = std::vector<JetPair>());

// combined analysis selectors
bool goodMuonWithoutIsolation(unsigned int i);
bool goodMuonIsolated(unsigned int i);
enum MuFOTypes { MuFOV1, MuFOV2 };
bool fakableMuon(unsigned int i, MuFOTypes);

//
// trigger
//

bool passedTriggerRequirements();

//
// Met
//

// the following two functions set met values used in 
// all other selectors and functions
double metValue();
double metPhiValue();

// analysis MET requirements
bool passedMetRequirements(unsigned int i_hyp);

double projectedMet(unsigned int i_hyp, double met, double phi);

double nearestDeltaPhi(double Phi, int i_hyp);

//
// Jets
//
enum WWJetType { CaloJet, jptJet, pfJet, TrkJet, GenJet };
std::vector<JetPair> getJets(WWJetType type, 
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
// Tools
//
bool defaultBTag(WWJetType type, unsigned int iJet);

HypTypeInNtuples hypType(unsigned int i_hyp);
class CMS2;
class EventIdentifier {
  unsigned long int run, event, lumi;
  float trks_d0;
  bool data;
 public:
  EventIdentifier(CMS2& cms2, bool isData);
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

float getHiggsPt();
float getHiggsPtWeight(float pt, int mH);
float getZPt();
float getZRapidity();
float getZMass();
float getDYNNLOWeight(float pt, float rap, float mass);
				
bool hypo (int i_hyp, double weight, bool realData = false ); 

double mt(double pt1, double pt2, double dphi);

class TChain;
void ScanChain( TChain* chain, 
	   SmurfTree::DataType sample, 
	   double integratedLumi,
	   double xsec,
	   int nProcessedEvents,
	   bool identifyEvents,
	   bool realData = false,
	   TString cms2_json_file = "");

void ProcessSample( std::string file_pattern, 
		    SmurfTree::DataType sample, 
		    double integratedLumi,
		    double xsec,
		    int nProcessedEvents,
		    bool identifyEvents = false,
		    bool realData = false,
		    TString cms2_json_file = "");
void ProcessSample( std::vector<std::string> file_patterns, 
		    SmurfTree::DataType sample,
		    double integratedLumi,
		    double xsec,
		    int nProcessedEvents,
		    bool identifyEvents = false,
		    bool realData = false,
		    TString cms2_json_file = "");

#endif
