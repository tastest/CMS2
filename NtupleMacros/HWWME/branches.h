#ifndef BRANCHES_H
#define BRANCHES_H
#include <vector>
#include <set>
#include "TFile.h"
#include "TTree.h"
#include "TBitSet.hh"

#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"

#include "../../../Smurf/Core/SmurfTree.h"


inline double fast_sign(double f) {
    if (f > 0) return 1;
    return (f == 0) ? 0 : -1;
}


typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

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

TBitSet* _cutWord;
TBitSet* _cutMask;

int eventCount[kNdilepFlav];
double eventYield[kNdilepFlav];


const TString evInfo_format =
  "ientry/I:"
  "runNumber/I:"
  "lumiblockNumber/I:"
  "eventNumber/I:"
  "weight/D:"
  "dilflavor/I:"
  "Met/D:"
  "MetSpec/D:"
  "MetX/D:"
  "MetY/D:"
  "Njets/I:"
  "Qprod/I:"
  "dimass/D:"
  "dPhiLeptons/D:"
  "dEtaLeptons/D:"
  "dRLeptons/D:"

  "lep1_Type/I:"
  "lep1_FakeType/I:"
  "lep1_Px/D:"
  "lep1_Py/D:"
  "lep1_Pz/D:"
  "lep1_E/D:"
  "lep1_Charge/I:"

  "lep2_Type/I:"
  "lep2_FakeType/I:"
  "lep2_Px/D:"
  "lep2_Py/D:"
  "lep2_Pz/D:"
  "lep2_E/D:"
  "lep2_Charge/I";

struct evInfo_t {
    int    ientry;
    int    runNumber;
    int    lumiblockNumber;
    int    eventNumber;
    double weight;
    int    dilflavor;
    double tcMet;
    double tcMetSpec;
    double tcMetX;
    double tcMetY;
    int    Njets;
    int    Qprod;
    double dimass;
    double dPhiLeptons;
    double dEtaLeptons;
    double dRLeptons;

    int    lep1_Type;
    int    lep1_FakeType;
    double lep1_Px;
    double lep1_Py;
    double lep1_Pz;
    double lep1_E;
    int    lep1_Charge;

    int    lep2_Type;
    int    lep2_FakeType;
    double lep2_Px;
    double lep2_Py;
    double lep2_Pz;
    double lep2_E;
    int    lep2_Charge;

};

evInfo_t _evInfo;

struct EventIdentifier {
  unsigned long int run, event, lumi;
  float trks_d0;
  float hyp_lt_pt, hyp_lt_eta, hyp_lt_phi;
  bool operator < (const EventIdentifier &) const;
  bool operator == (const EventIdentifier &) const;
};


static std::set<EventIdentifier> already_seen;

bool is_duplicate (const EventIdentifier &id);


bool ApplyPreSelection(unsigned int);
int ApplyEventSelection(unsigned int, bool realData);
void FillSmurfNtuple(SmurfTree& tree, unsigned int i_hyp, double weight, const char *process);

void CalculateFakeRateProb();


//Utils
bool inZmassWindow(float);
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
double projectedMet(int);

//trigger
bool  passedTrigger(TString trigName);
bool passedTriggerRequirements();
bool defaultBTag(int type, unsigned int iJet); 
double BTag(int type, unsigned int iJet); 

//muon ID
bool goodMuonIsolated(unsigned int i);
bool fakableMuon(unsigned int i);
bool ww_mud0PV(unsigned int index);
bool ww_muId(unsigned int index);
double ww_muIsoVal(unsigned int index);
bool ww_muIso(unsigned int index);

//electron ID
bool goodElectronIsolated(unsigned int i);
bool fakableElectron(unsigned int i);
bool ww_eld0PV(unsigned int index);
bool ww_elId(unsigned int index);
double ww_elIsoVal(unsigned int index);
bool ww_elIso(unsigned int index);


unsigned int numberOfSoftMuons(int i_hyp, bool nonisolated, const std::vector<LorentzVector>& = std::vector<LorentzVector>());
unsigned int numberOfExtraLeptons(int i_hyp, double minPt);
bool toptag(int type, int i_hyp, double minPt);

double weight;
double Integral;


unsigned int evtevent_;
unsigned int evtrun_;
unsigned int evtlumiBlock_;
float evtscale1fb_;
float evtmet_;
float evtmetPhi_;
float evtmetMuonCorr_;
float evtmetMuonCorrPhi_;
float evttcmet_;
float evttcmetPhi_;
float evttcsumet_;
float genmet_;
float genmetPhi_;
float evtmetMuonJESCorr_;
float evtmetMuonJESCorrPhi_;
float evtsumet_;
float evtsumetMuonCorr_;
float evtpfmet_;
float evtpfmetPhi_;
std::vector<int> *genpsid_;
std::vector<int> *genpsidmother_;
std::vector<float> *jetscor_;
std::vector<float> *jetsemFrac_;
int hyptype_;
int hypltid_;
int hypllid_;
ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *hypp4_;
ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *hypltp4_;
ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *hypllp4_;
std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > *genpsp4_;
std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > *hypjetsp4_;
std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > *genjetsp4_;
std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > *jetsp4_;
std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > *jptsp4_;
int ngenps_;
int nhypjets_;
int ngenjets_;
int njets_;
int njpts_;

#endif
