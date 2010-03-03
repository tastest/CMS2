#ifndef SELECTION_H
#define SELECTION_H
#include "TLorentzVector.h"
#include "Math/LorentzVector.h"
#include <vector>
#include <iostream>
#include <utility>

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

enum SubDetector
{
  EBp         =     0,
  EBm         =     1,
  EEp         =     2,
  EEm         =     3,
  HBHEa       =     5,
  HBHEb       =     6,
  HBHEc       =     7,
  HF          =     8,
  HO          =     9,
  RPC         =    12,
  DT0         =    13,
  DTp         =    14,
  DTm         =    15,
  CSCp        =    16,
  CSCm        =    17,
  CASTOR      =    20,
  TIBTID      =    24,
  TOB         =    25,
  TECp        =    26,
  TECm        =    27,
  BPIX        =    28,
  FPIX        =    29,
  ESp         =    30,
  ESm         =    31,
  nPartitions =    24
};

bool isSubDetectorGood( unsigned int cuts );

bool isTrackQuality( int index, int cuts);

bool inZmassWindow (float mass);
bool supertightElectron (int index);
bool deltaPhiInElectron (int index);
bool goodElectronWithoutIsolation(int index);
bool goodElectronWithoutIsolationWithoutd0(int index);
bool goodMuonWithoutIsolation(int index);
bool passElectronIsolation(int index, bool use_calo_iso);
bool passElectronIsolationLoose(int index, bool use_calo_iso);
bool passElectronIsolationLoose2(int index, bool use_calo_iso);
bool passElectronIsolation_1_6(int index, bool use_calo_iso);
bool passElectronIsolationLoose_1_6(int index, bool use_calo_iso);
bool passElectronIsolationLoose2_1_6(int index, bool use_calo_iso);
bool passMuonIsolation(int index);
bool passMuonIsolationLoose(int index);
bool goodMuonIsolated(int index);
bool goodElectronIsolated(int index, bool use_calo_iso = false);
bool goodLooseElectronWithoutIsolation(int index);
double nearestDeltaPhi(double Phi, int index);

bool isChargeFlip(int elIndex);
//old cuts on conversions
bool conversionElectron(int electron);
//new conversion stuff...cut on shared fraction of hits implemented
//when looking for electron's ctf track
bool isconversionElectron09(int elIdx);
std::pair<float, float> getConversionInfo(LorentzVector trk1_p4, 
					  int trk1_q, float trk1_d0, 
					  LorentzVector trk2_p4,
					  int trk2_q, float trk2_d0,
					  float bField);


int numberOfExtraMuons(int i_hyp, bool nonisolated);
bool passMuonBVeto_1_6 (int i_dilep, bool soft_nonisolated);

bool passTriLepVeto (int i_dilep);
int tagMuonIdx (int i_dilep);
double tagMuonPt (int i_dilep);
double tagMuonRelIso (int i_dilep);
int additionalZcounter();
bool additionalZveto();
int nTrkJets(int i_hyp);
bool passTrkJetVeto(int i_hyp);
std::vector<LorentzVector> JPTs(int i_hyp, double etThreshold);
unsigned int nJPTs(int i_hyp, double etThreshod);
unsigned int nJPTs(int i_hyp); // with default threshold = 20
// bool passCaloTrkjetCombo ();

double mu_rel_iso (int index);
double el_rel_iso (int index, bool use_calo_iso);
double reliso_lt (int i_hyp, bool use_calo_iso = false);
double reliso_ll (int i_hyp, bool use_calo_iso = false);

int conversionPartner (int i_el);
double conversionDeltaPhi (int i_conv, int i_el);

bool passTrackIsolation(int index);
int passTrackZVeto(int hyp_index);

int NjetVeto(std::vector<TLorentzVector>& Jet, double min_et);
bool trueElectron(int index);
bool trueMuon(int index);
bool trueMuonFromW(int index);
bool trueElectronFromW(int index);

// muon cleaning by ratio of track and global Pt
bool muonReconstructionCleaning(int i_hyp, float threshold);

bool passLeptonIsolation(int id, int index, bool use_ele_calo_iso);
double muonTrkIsolationPAT(int index); 
double muonCalIsolationPAT(int index); 
bool electron20Eta2p4(int index);
bool looseElectronSelectionNoIsoTTDil08(int index) ;
double electronTrkIsolationPAT(int index); 
double electronCalIsolationPAT(int index); 
float electronTrkIsolationTTDil08(int index);
float electronCalIsolationTTDil08(int index);
bool looseElectronSelectionTTDil08(int index) ;
bool passElectronIsolationTTDil08(int index);
bool muon20Eta2p4(int index);
bool looseMuonSelectionNoIsoTTDil08(int index) ;
bool lepton20Eta2p4(int id, int index);
float muonTrkIsolationTTDil08(int index);
float muonCalIsolationTTDil08(int index);
float leptonTrkIsolationTTDil08(int id, int index);
float leptonCalIsolationTTDil08(int id, int index);
bool looseMuonSelectionTTDil08(int index) ;
bool passMuonIsolationTTDil08(int index) ;
bool passLeptonIsolationTTDil08(int id, int index);
bool looseLeptonSelectionNoIsoTTDil08(int id, int index);
bool looseLeptonSelectionTTDil08(int id, int index);
bool additionalZvetoTTDil08() ;
bool haveExtraMuon(int hypIdx);
bool haveExtraMuon5(int hypIdx);
bool passTriggersMu9orLisoE15(int dilType) ;
bool passTriggersTTDil08JanTrial(int dilType) ;
int genpCountPDGId(int id0, int id1=-1, int id2=-1); 
int genpCountPDGId_Pt20h24(int id0, int id1=-1, int id2=-1);  
int genpDileptonType();
int eventDilIndexByWeightTTDil08(const std::vector<unsigned int>& goodHyps, int& strasbourgDilType, bool printDebug = false,
				 bool usePtOnlyForWeighting = false);

// SUSY Dilepton group cuts (VJets09)

double inv_mu_rel_iso(int index);
double inv_el_rel_iso(int index, bool use_calo_iso);
bool passMuonIsolationVJets09(int index);
bool passElectronIsolationVJets09(int index, bool use_calo_iso);
bool passLeptonIsolationVJets09(int id, int index);
bool looseElectronSelectionVJets09(int index);
bool looseMuonSelectionVJets09(int index);
bool passLeptonIDVJets09(int id, int index);
int numberOfExtraMuonsVJets09(int i_hyp);
int numberOfExtraElectronsVJets09(int i_hyp);

//SUSY dilepton selection TAS group

bool comparePt(const LorentzVector &lv1,  
	       const LorentzVector &lv2);
bool GoodSusyElectronWithoutIsolation(int index);
bool GoodSusyElectronWithoutIsolationNoD0(int index);
bool PassSusyElectronIsolationLoose(int index, bool use_calo_iso);
bool GoodSusyMuonWithIsolation(int index); 
bool GoodSusyMuonWithoutIsolation(int index);
double inv_mu_relsusy_iso(int index);
double inv_el_relsusy_iso(int index, bool use_calo_iso);
bool GoodSusyElectronWithIsolation(int index, bool use_calo_iso); 
bool GoodSusyLeptonWithIsolation(int id, int index);
bool GoodSusyLeptonID(int id, int index);
bool PassSusyMuonIsolation(int index);
bool PassSusyElectronIsolation(int index, bool use_calo_iso);
bool PassSusyLeptonIsolation(int id, int index);
bool GoodSusyLeptonID(int id, int index);
bool GoodSusyTrigger(int dilType);
int numberOfExtraElectronsSUSY(int i_hyp);
int numberOfExtraMuonsSUSY(int i_hyp);
std::vector<LorentzVector> getCaloJets(int i_hyp); 
std::vector<LorentzVector> getJPTJets(int i_hyp); 
bool additionalZvetoSUSY09(int i_hyp);
bool isFakeableElSUSY09(int iEl);
bool isFakeableMuSUSY09(int iMu);
bool isNumElSUSY09(int iEl);
bool isNumMuSUSY09(int iMu);

//TTDil08 fake rate functions
bool isNumElTTDil08(int iEl);
bool isFakeableElTTDil08(int iEl);
bool isNumMuTTDil08(int iMu);
bool isFakeableMuTTDil08(int iMu);

int findPrimTrilepZ(int i_hyp, double &mass);
bool vetoAddZ(int i_hyp, int unusedLepton, double &mass);
std::vector<LorentzVector> JPTsTrilep(int i_hyp, double etThreshold);
unsigned int nJPTsTrilep(int i_hyp, double etThreshold);
double nearestDeltaPhiTrilep(double Phi, int i_hyp);
double nearestDeltaPhiJet(double Phi, int i_hyp);
#endif
