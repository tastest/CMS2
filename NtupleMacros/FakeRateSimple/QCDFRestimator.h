// -*- C++ -*-
#ifndef QCDFRestimator_H
#define QCDFRestimator_H
#include "CMS2.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"

class QCDFRestimator {

public: 
  
  int ScanChainQCD ( TChain* chain, TString prefix="", float kFactor=1.0, 
		     int prescale=1, float pthatmin = 0, float pthatmax = 999999999.9);
  int ScanChainWJets( TChain* chain, TString prefix="", float kFactor=1.0, 
		      int prescale=1);
  
  Float_t GetValueTH2F(Float_t x, Float_t y, TH2F* hist);
  void bookHistos(const char *prefix);
    
  bool isTrueMuFromW(int iMu);
  bool isTrueElFromW(int iEl);
  bool isTrueLeptonfromW(int pid);
  bool isElFromMu(int iEl);
  bool passFakeJetTrigger(float unCorrJetPtCut);
  //bool isFakeableEl(int iEl);
  //bool isNumEl(int iEl);
  //bool isFakeableMu(int iMu);
  //bool isNumMu(int iMu);
  bool testJetsForElectrons(vector<LorentzVector>&, const LorentzVector&);
  //FO
  TH2F *h_FOptvseta[2];
  TH1F *h_FOpt[2];
  TH1F *h_FOeta[2];
  TH1F *h_FOhfpt[2];
  TH1F *h_FOhfeta[2];
  TH1F *h_FOlqpt[2];
  TH1F *h_FOlqeta[2];
  TH1F *h_FOgpt[2];
  TH1F *h_FOgeta[2];

  //FO split PDG id:
  TH1F *h_FOmc3Id[2];
  TH1F *h_FOmc3dR[2];
  
  //tight selection info
  TH2F *h_numptvseta[2];
  TH1F *h_numpt[2];
  TH1F *h_numeta[2];
  TH1F *h_numhfpt[2];
  TH1F *h_numhfeta[2];
  TH1F *h_numlqpt[2];
  TH1F *h_numlqeta[2];
  TH1F *h_numgpt[2];
  TH1F *h_numgeta[2];
  

  TH1F *h_nummc3Id[2];
  TH1F *h_nummc3dR[2];
    
  //FR
  TH2F *h_FRptvseta[2];
  TH2F *h_FRErrptvseta[2];
  TH1F *h_FRpt[2];
  TH1F *h_FReta[2];
  TH1F *h_FRmc3Id[2];
  TH1F *h_FRmc3Id_largedR[2];
  TH1F *h_FRhfpt[2];
  TH1F *h_FRhfeta[2];
  TH1F *h_FRlqpt[2];
  TH1F *h_FRlqeta[2];
  TH1F *h_FRgpt[2];
  TH1F *h_FRgeta[2];

  //nJets
  TH1F *h_predictednJets[2];
  TH1F *h_actualnJets[2];  
  TH3F *h_nJets3D[2];

  //nJets for FO object
  TH1F *h_FOnJets[2];
};
#endif
