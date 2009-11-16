#ifndef SELECTION_H
#define SELECTION_H
#include "TLorentzVector.h"
#include "Math/LorentzVector.h"
#include <vector>
#include <iostream>
#include <utility>

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

enum TrackAlgorithm { undefAlgorithm=0, ctf=1, rs=2, cosmics=3, iter0=4, 
                      iter1=5, iter2=6, iter3=7, iter4=8, iter5=9, iter6=10, iter7=11, iter8=12, iter9=13,iter10=14,
                      outInEcalSeededConv=15, inOutEcalSeededConv=16, 
                      nuclInter=17,
                      standAloneMuon=18,globalMuon=19,cosmicStandAloneMuon=20,cosmicGlobalMuon=21,
                      iter1LargeD0=22,iter2LargeD0=23,iter3LargeD0=24,iter4LargeD0=25,iter5LargeD0=26,
                      bTagGhostTracks=27,
                      beamhalo=28, 
                      gsf=29,
                      algoSize=30 };

enum TrackQuality { undefQuality=-1, loose=0, tight=1, highPurity=2, confirmed=3, goodIterative=4, qualitySize=5};

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
bool pass2Met(int index, const class TVector3& corr);
double nearestDeltaPhi(double Phi, int index);
double MetSpecial(double MET, double MetPhi, int index);
bool pass4Met(int index, const class TVector3& corr);
bool met10(int index, const class TVector3& corr);
bool met1(int index, const class TVector3& corr);
bool sumEt10(double sumEt);
bool sumEt1(double sumEt);

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
bool isDYee();
bool isDYmm();
bool isDYtt();
bool isWe();
bool isWm();
bool isWt();
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

int getDrellYanType();
int getZZType ();
void dumpDocLines();
int NjetVeto(std::vector<TLorentzVector>& Jet, double min_et);
bool trueElectron(int index);
bool trueMuon(int index);
bool trueMuonFromW(int index);
bool trueMuonFromW_WJets(int index);
bool trueElectronFromW(int index);
bool trueElectronFromW_WJets(int index);
//bool isFakeDenominatorElectron(int index);
//bool isFakeNumeratorElectron(int index, int type); 

// muon cleaning by ratio of track and global Pt
bool muonReconstructionCleaning(int i_hyp, float threshold);

// new MET factorization
bool metSimple (float threshold, const TVector3& corr);
bool metBalance (int i_hyp, const TVector3& corr);
bool metProjected (int i_hyp, const TVector3& corr);
bool pass5Met (int i_hyp, const TVector3& corr);
//


bool passLeptonIsolation(int id, int index, bool use_ele_calo_iso);
void leptonGenpCount(int& nele, int& nmuon, int& ntau) ; 
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
bool passMetVJets09(float value, bool useTcMet);
int numberOfExtraMuonsVJets09(int i_hyp);
int numberOfExtraElectronsVJets09(int i_hyp);

//SUSY dilepton selection TAS group

bool comparePt(const LorentzVector &lv1,  
	       const LorentzVector &lv2);
bool GoodSusyElectronWithoutIsolation(int index);
bool GoodSusyElectronWithoutIsolationNoD0(int index);
bool GoodSusyElectronWithIsolationLoose(int index, bool use_calo_iso); 
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
int ttbarconstituents(int i_hyp);
bool additionalZvetoSUSY09(int i_hyp);
bool idIsBeauty(int id);
bool idIsCharm(int id);
int leptonIsFromW(int idx, int id, LorentzVector v);
bool isFakeableElSUSY09(int iEl);
bool isFakeableMuSUSY09(int iMu);
bool isNumElSUSY09(int iEl);
bool isNumMuSUSY09(int iMu);

//TTDil08 fake rate functions
bool isNumElTTDil08(int iEl);
bool isFakeableElTTDil08(int iEl);
bool isNumMuTTDil08(int iMu);
bool isFakeableMuTTDil08(int iMu);

bool trueGammaFromMuon(int electron);

int findPrimTrilepZ(int i_hyp, double &mass);
bool vetoAddZ(int i_hyp, int unusedLepton, double &mass);
std::vector<LorentzVector> JPTsTrilep(int i_hyp, double etThreshold);
unsigned int nJPTsTrilep(int i_hyp, double etThreshold);
double nearestDeltaPhiTrilep(double Phi, int i_hyp);
double nearestDeltaPhiJet(double Phi, int i_hyp);
double MetSpecialTrilep(double MET, double MetPhi, int index);
#endif
