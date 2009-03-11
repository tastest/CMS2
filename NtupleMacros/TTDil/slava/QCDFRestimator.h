// -*- C++ -*-
#ifndef QCDFRestimator_H
#define QCDFRestimator_H
#include "CORE/CMS2.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"

class QCDFRestimator {

public: 
  
  int ScanChainQCD ( TChain* chain, char * prefix="", float kFactor=1.0, 
		     int prescale=1, float pthatmin = 0, float pthatmax = 999999999.9);
  int ScanChainWJets( TChain* chain, char * prefix="", float kFactor=1.0, 
		      int prescale=1);
  
  Float_t GetValueTH2F(Float_t x, Float_t y, TH2F* hist);
  void bookHistos(char *prefix);
    
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
  
  //tight selection info
  TH2F *h_numptvseta[2];
  TH1F *h_numpt[2];
  TH1F *h_numeta[2];
  
  //FR
  TH2F *h_FRptvseta[2];
  TH2F *h_FRErrptvseta[2];
  TH1F *h_FRpt[2];
  TH1F *h_FReta[2];
  
  //nJets
  TH1F *h_predictednJets[2];
  TH1F *h_actualnJets[2];  
  TH3F *h_nJets3D[2];
};
#endif
