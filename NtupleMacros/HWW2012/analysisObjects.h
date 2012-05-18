#ifndef WW_analysisObjects_h
#define WW_analysisObjects_h

#include "analysisEnums.h"
#include "wwtypes.h"

#ifndef __CINT__
#include "CORE/CMS2.h"
#include "jetcorr/FactorizedJetCorrector.h"
#include "CORE/jetSelections.h"
#endif

//
// hyps
//

HypTypeInNtuples hypType(unsigned int i_hyp);

//
// met
//

double metValue();
double metPhiValue();
double sumetValue();

//
// jets
//

WWJetType jetType();
Bool_t comparePt(JetPair lv1, JetPair lv2);

std::vector<JetPair> getJets(WWJetType type,
                 LorentzVector &lt,
                 LorentzVector &ll,
                 double etThreshold,
                 double maxEta,
                 bool applyJEC,
                 FactorizedJetCorrector *jet_corrector_pfL1FastJetL2L3,
                 bool sorted = false,
                 bool btag = false);

std::vector<JetPair> getJets(WWJetType type,
                 int i_hyp,
                 double etThreshold,
                 double maxEta,
                 bool applyJEC,
                 FactorizedJetCorrector *jet_corrector_pfL1FastJetL2L3,
                 bool sorted = false,
                 bool btag = false);

unsigned int numberOfJets(unsigned int i_hyp, bool applyJEC, FactorizedJetCorrector *jet_corrector_pfL1FastJetL2L3);
std::vector<JetPair> getDefaultJets(unsigned int i_hyp, bool applyJEC, FactorizedJetCorrector *jet_corrector_pfL1FastJetL2L3, bool btagged=false);

//
// Btagging
//

double BTag(LorentzVector jetP4);
double BTag(WWJetType type, unsigned int iJet);
double BTag(WWJetType type, unsigned int iJet, float corjetpt);
bool   defaultBTag(WWJetType type, unsigned int iJet, float jec=1.);

#endif

