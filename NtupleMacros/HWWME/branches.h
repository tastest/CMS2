#ifndef BRANCHES_H
#define BRANCHES_H
#include <vector>
#include <set>
#include "TFile.h"
#include "TTree.h"
#include "TBitSet.hh"

#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"

inline double fast_sign(double f) {
    if (f > 0) return 1;
    return (f == 0) ? 0 : -1;
}


typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

enum hyp_selection {
  kcut_Trigger,
  kcut_OS,
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
int ApplyEventSelection(unsigned int);
int FillNeededVariables( unsigned int);
void CalculateFakeRateProb();


//Utils
bool inZmassWindow(float);
double projectedMet(int);
double nearestDeltaPhi(double Phi, int i_hyp);
bool comparePt(LorentzVector lv1, LorentzVector lv2);
std::vector<LorentzVector> getJets(int type, int i_hyp, double etThreshold, double etaMax, bool sortJets, bool btag);
bool isGoodVertex(size_t ivtx);
double dzPV(const LorentzVector& vtx, const LorentzVector& p4, const LorentzVector& pv);

//trigger
bool passedTriggerRequirements();
bool defaultBTag(int type, unsigned int iJet); 
double BTag(int type, unsigned int iJet); 

//muon ID
bool goodMuonIsolated(unsigned int i);
bool ww_mud0PV(unsigned int index);
bool ww_muId(unsigned int index);
double ww_muIsoVal(unsigned int index);
bool ww_muIso(unsigned int index);

//electron ID
bool goodElectronIsolated(unsigned int i);
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
TFile *outFile_;
TTree *outTree_;
TTree *outGlobalTree_;

void InitSkimmedTree(std::string skimFilePrefix="") {

   if(skimFilePrefix != "")
      outFile_ = TFile::Open(string(skimFilePrefix + "_skimmednTuple.root").c_str(),"RECREATE");
   else 
   outFile_->cd();

   outTree_ = new TTree("Events", "");
   outGlobalTree_ = new TTree("Global", "");
   
   outTree_->Branch("evInfo", &_evInfo, evInfo_format);
   outGlobalTree_->Branch("Integral", &Integral,"Integral/D");


   //book the branches
   outTree_->Branch("evtevent_",   &evtevent_,   "evtevent_/i");
   outTree_->SetAlias("evt_event",   "evtevent_");
   outTree_->Branch("evtrun_",   &evtrun_,   "evtrun_/i");
   outTree_->SetAlias("evt_run",   "evtrun_");
   outTree_->Branch("evtlumiBlock_",   &evtlumiBlock_,   "evtlumiBlock_/i");
   outTree_->SetAlias("evt_lumiBlock",   "evtlumiBlock_");

   outTree_->Branch("evtscale1fb_",   &evtscale1fb_,   "evtscale1fb_/F");
   outTree_->SetAlias("evt_scale1fb",   "evtscale1fb_");
   outTree_->Branch("evtmet_",   &evtmet_,   "evtmet_/F");
   outTree_->SetAlias("evt_met",   "evtmet_");
   outTree_->Branch("evtmetPhi_",   &evtmetPhi_,   "evtmetPhi_/F");
   outTree_->SetAlias("evt_metPhi",   "evtmetPhi_");
   outTree_->Branch("evtmetMuonCorr_",   &evtmetMuonCorr_,   "evtmetMuonCorr_/F");
   outTree_->SetAlias("evt_metMuonCorr",   "evtmetMuonCorr_");
   outTree_->Branch("evtmetMuonCorrPhi_",   &evtmetMuonCorrPhi_,   "evtmetMuonCorrPhi_/F");
   outTree_->SetAlias("evt_metMuonCorrPhi",   "evtmetMuonCorrPhi_");
   outTree_->Branch("evttcmet_",   &evttcmet_,   "evttcmet_/F");
   outTree_->SetAlias("evt_tcmet",   "evttcmet_");
   outTree_->Branch("evttcmetPhi_",   &evttcmetPhi_,   "evttcmetPhi_/F");
   outTree_->SetAlias("evt_tcmetPhi",   "evttcmetPhi_");
   outTree_->Branch("evttcsumet_",   &evttcsumet_,   "evttcsumet_/F");
   outTree_->SetAlias("evt_tcsumet",   "evttcsumet_");
   outTree_->Branch("genmet_",   &genmet_,   "genmet_/F");
   outTree_->SetAlias("gen_met",   "genmet_");
   outTree_->Branch("genmetPhi_",   &genmetPhi_,   "genmetPhi_/F");
   outTree_->SetAlias("gen_metPhi",   "genmetPhi_");
   outTree_->Branch("evtmetMuonJESCorr_",   &evtmetMuonJESCorr_,   "evtmetMuonJESCorr_/F");
   outTree_->SetAlias("evt_metMuonJESCorr",   "evtmetMuonJESCorr_");
   outTree_->Branch("evtmetMuonJESCorrPhi_",   &evtmetMuonJESCorrPhi_,   "evtmetMuonJESCorrPhi_/F");
   outTree_->SetAlias("evt_metMuonJESCorrPhi",   "evtmetMuonJESCorrPhi_");
   outTree_->Branch("evtsumet_",   &evtsumet_,   "evtsumet_/F");
   outTree_->SetAlias("evt_sumet",   "evtsumet_");
   outTree_->Branch("evtsumetMuonCorr_",   &evtsumetMuonCorr_,   "evtsumetMuonCorr_/F");
   outTree_->SetAlias("evt_sumetMuonCorr",   "evtsumetMuonCorr_");
   outTree_->Branch("evtpfmet_",   &evtpfmet_,   "evtpfmet_/F");
   outTree_->SetAlias("evt_pfmet",   "evtpfmet_");
   outTree_->Branch("evtpfmetPhi_",   &evtpfmetPhi_,   "evtpfmetPhi_/F");
   outTree_->SetAlias("evt_pfmetPhi",   "evtpfmetPhi_");
   outTree_->Branch("ints_genpsid_",   "std::vector<int>",   &genpsid_);
   outTree_->SetAlias("genps_id",   "ints_genpsid_");
   outTree_->Branch("ints_genpsidmother_",   "std::vector<int>",   &genpsidmother_);
   outTree_->SetAlias("genps_id_mother",   "ints_genpsidmother_");
   outTree_->Branch("floats_jetscor_",   "std::vector<float>",   &jetscor_);
   outTree_->SetAlias("jets_cor",   "floats_jetscor_");
   outTree_->Branch("floats_jetsemFrac_",   "std::vector<float>",   &jetsemFrac_);
   outTree_->SetAlias("jets_emFrac",   "floats_jetsemFrac_");
   outTree_->Branch("hyptype_",   &hyptype_,   "hyptype_/I");
   outTree_->SetAlias("hyp_type",   "hyptype_");
   outTree_->Branch("hypltid_",   &hypltid_,   "hypltid_/I");
   outTree_->SetAlias("hyp_lt_id",   "hypltid_");
   outTree_->Branch("hypllid_",   &hypllid_,   "hypllid_/I");
   outTree_->SetAlias("hyp_ll_id",   "hypllid_");
   outTree_->Branch("ngenps_",   &ngenps_,   "ngenps_/I");
   outTree_->SetAlias("n_genps",   "ngenps_");
   outTree_->Branch("nhypjets_",   &nhypjets_,   "nhypjets_/I");
   outTree_->SetAlias("n_hyp_jets",   "nhypjets_");
   outTree_->Branch("ngenjets_",   &ngenjets_,   "ngenjets_/I");
   outTree_->SetAlias("n_genjets",   "ngenjets_");
   outTree_->Branch("njets_",   &njets_,   "njets_/I");
   outTree_->SetAlias("n_jets",   "njets_");
   outTree_->Branch("njpts_",   &njpts_,   "njpts_/I");
   outTree_->SetAlias("n_jpts",   "njpts_");
} 
#endif
