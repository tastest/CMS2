#ifndef WW_analysisSelections_h
#define WW_analysisSelections_h

#include "Math/LorentzVector.h"
#include "Rtypes.h"
#include "TPRegexp.h"

#include <set>
#include <vector>

#include "../Tools/ElectronIDMVA.h"
#include "../Tools/MuonIDMVA.h"
#include "../Tools/EGammaMvaEleEstimator.h"
#include "../Tools/MuonMVAEstimator.h"

#include "../../../Smurf/Core/SmurfTree.h"
#include "wwtypes.h"
#include "analysisEnums.h"
#include "analysisObjects.h"

//
// Electron Id
//

bool   ww_elBase(unsigned int i);
bool   ww_elId(unsigned int i, bool useLHeleId, int useMVAeleId, EGammaMvaEleEstimator* egammaMvaEleEstimator);
bool   ww_eld0(unsigned int i);
bool   ww_eld0PV(unsigned int i);
bool   ww_eldZPV(unsigned int i);
bool   ww_elIso(unsigned int i);
double ww_elIsoVal(unsigned int i);

// combined analysis selectors
bool goodElectronTMVA(EGammaMvaEleEstimator* egammaMvaEleEstimator, int useMVAeleId, unsigned int i); 
bool goodElectronWithoutIsolation(unsigned int i, bool useLHeleId, int useMVAeleId, EGammaMvaEleEstimator* egammaMvaEleEstimator);
bool goodElectronIsolated(unsigned int i, bool useLHeleId, int useMVAeleId, EGammaMvaEleEstimator* egammaMvaEleEstimator, bool lockToCoreSelectors);
bool fakableElectron(unsigned int i,EleFOTypes);
bool ElectronFOV4(unsigned int i);

//
// Muon Id
//

bool   ww_muBase(unsigned int i);
bool   ww_muId(unsigned int i, bool useMVAmuId, MuonIDMVA *mva);
bool   ww_muIso(unsigned int i);
bool   ww_muIso(unsigned int i, MuonMVAEstimator* muonMVAEstimator, std::vector<Int_t> IdentifiedMu, std::vector<Int_t> IdentifiedEle);
double ww_muIsoVal(unsigned int i);
bool   ww_mud0(unsigned int i);
bool   ww_mud0PV(unsigned int i);
bool   ww_mudZPV(unsigned int i, float cut=0.1);

// combined analysis selectors
bool goodMuonTMVA(MuonIDMVA *mva, unsigned int i);
bool goodMuonWithoutIsolation(unsigned int i, bool useMVAmuId, MuonIDMVA *mva);
bool goodMuonIsolated( unsigned int i, bool lockToCoreSelectors, bool useMVAmuId, MuonIDMVA *mva, 
					   MuonMVAEstimator* muonMVAEstimator, std::vector<Int_t> IdentifiedMu, std::vector<Int_t> IdentifiedEle );
bool fakableMuon(unsigned int i, MuFOTypes, MuonMVAEstimator* muonMVAEstimator,
                std::vector<Int_t> IdentifiedMu, std::vector<Int_t> IdentifiedEle);
bool passMuonRingsMVA(unsigned int mu, MuonMVAEstimator* muonMVAEstimator, std::vector<Int_t> IdentifiedMu, std::vector<Int_t> IdentifiedEle);

//
// extra leptons
//

//
// leptons
//

std::vector<LeptonPair> getExtraLeptons(int i_hyp, double minPt,  bool useLHeleId, int useMVAeleId, 
										EGammaMvaEleEstimator* egammaMvaEleEstimator, bool useMVAmuId, MuonIDMVA *mumva);
unsigned int numberOfExtraLeptons(int i_hyp, double minPt, bool useLHeleId, int useMVAeleId, 
								  EGammaMvaEleEstimator* egammaMvaEleEstimator, bool useMVAmuId, MuonIDMVA *mumva);
unsigned int numberOfSoftMuons(int i_hyp, bool nonisolated, const std::vector<JetPair>& = std::vector<JetPair>());

//
// trigger
//

bool passedTriggerRegExp(TPMERegexp trigName);
bool passedTrigger(TString trigName, unsigned int minRun = 0, unsigned int maxRun = 999999999);
bool passedTriggerRequirements();
bool passedTriggerRequirementsWithRuns();

//
// Met
//

// analysis MET requirements
double  minmet(unsigned int i_hyp);
bool 	passedMetRequirements(unsigned int i_hyp, FactorizedJetCorrector *jet_corrector_pfL1FastJetL2L3);
bool 	passedMetRequirements(unsigned int i_hyp, unsigned int njets, std::vector<JetPair> jets );
double 	projectedMet(unsigned int i_hyp, double met, double phi);
double 	nearestDeltaPhi(double Phi, int i_hyp);

//
// Vertex selections
//

bool isGoodVertex(size_t ivtx);
unsigned int nGoodVertex();
int primaryVertex();

//
// Jets 
//
bool passMVAJetId(double corjetpt, double jeteta, double mvavalue, unsigned int tightness);        

//
// Other
//

bool inZmassWindow(float mass, double delta=15.0);

//
// event top tagging
//

bool toptag(WWJetType type, int i_hyp, double minPt, FactorizedJetCorrector *jet_corrector_pfL1FastJetL2L3,
        	std::vector<JetPair> ignoreJets=std::vector<JetPair>());

#endif
