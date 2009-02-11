// -*- C++ -*-
#ifndef ttDilCounts_looper_H
#define ttDilCounts_looper_H
// switch to V6 later
#include "CMS2V010201.h"
#include "TH1F.h"

class ttDilCounts_looper {

	 public: 

//   static const unsigned int LOOSEDIL = 1<<10;
//   static const unsigned int LOOSEDIL_OS = (1<<10) + (1<<13);
//   static const unsigned int LOOSEDIL_OS_ONEWEIGHT = (1<<10) + (1<<13) + (1<<14);

  

  int ScanChain ( TChain* chain, char * prefix="", float kFactor=1.0, int prescale=1, bool oldjets=true, unsigned int cutsMask=31);

  TH1F* hnJet[4];       // Njet distributions
  TH1F* helePt[4][3];      // electron Pt
  TH1F* hmuPt[4][3];       // muon Pt
  TH1F* hmuPtFromSilicon[4][3];    // muon Pt (from tracker)
  TH1F* hminLepPt[4][3];   // minimum lepton Pt
  TH1F* hmaxLepPt[4][3];   // maximum lepton Pt
  TH1F* helePhi[4][3];     // electron phi
  TH1F* hmuPhi[4][3];      // muon phi
  TH1F* hdphiLep[4][3];    // delta phi between leptons
  TH1F* heleEta[4][3];     // electron eta
  TH1F* hmuEta[4][3];      // muon eta
  TH1F* hdilMass[4][3];    // dilepton mass
  TH1F* hdilMassTightWindow[4][3]; // dilepton mass, but zooming around Z
  TH1F* hdilPt[4][3];       // dilepton Pt
  TH1F* hmet[4][3];       // MET
  TH1F* hmetPhi[4][3];       // MET phi
  TH2F* hmetVsDilepPt[4][3];  // MET vs dilepton Pt

  TH2F* hmetOverPtVsDphi[4][3]; // MET/Lepton Pt vs DeltaPhi between MET and Lepton Pt
  TH2F* hdphillvsmll[4][3]; // delta phi between leptons vs dilepton mass
  TH1F* hptJet1[4][3];   // Pt of 1st jet
  TH1F* hptJet2[4][3];   // Pt of 2nd jet
  TH1F* hptJet3[4][3];   // Pt of 3rd jet
  TH1F* hptJet4[4][3];   // Pt of 4th jet
  TH1F* hetaJet1[4][3];   // eta of 1st jet
  TH1F* hetaJet2[4][3];   // eta of 2nd jet
  TH1F* hetaJet3[4][3];   // eta of 3rd jet
  TH1F* hetaJet4[4][3];   // eta of 4th jet
  TH1F* numTightLep[4][3]; // number of tight leptons per event.
  TH1F* heleSumPt[4][3];   // sumPt for electron isolation
  TH1F* hmuSumPt[4][3];   // sumPt for muon isolation
  TH1F* hmuSumIso[4][3];  // sum of trk pt, em et, had et in cone of 0.3
  TH1F* heleRelIso[4][3]; //  Iso variable defined as pt/(pt+sum) for electron
  TH1F* hmuRelIso[4][3]; //  Iso variable defined as pt/(pt+sum) for muons

  // Unfortunately, our ntuple has no info for electron isolation other than candidate electrons.
  // When counting electrons, we thus can not apply an isolation criteria at this point !!!
  // For muons we only count good isolated muons.
  TH1F* hnJetLepVeto[4]; //njet distribution after requiring numTightLep < 3.

  // The statement below should work but does not work due to bug in root when TH2 are also used
  // Rene Brun promised a fix.
  //TH1::SetDefaultSumw2(kTRUE); // do errors properly based on weights

  void bookHistos(char *prefix);
};
#endif
