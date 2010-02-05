// -*- C++ -*-
#ifndef ttDilCounts_looper_H
#define ttDilCounts_looper_H
//#include "CORE/CMS2.h"
#include "TH1F.h"
#include "TH2F.h"
#include <vector>
#include "TChain.h"

#ifndef ProcDSChain_H
#define ProcDSChain_H
struct ProcDSChain {
  ProcDSChain(TChain* ch, float sc = 1, bool doW = true, bool chDup = false): events(ch), scale1fb(sc),  useWeigtFromBranch(doW), 
									      checkDuplicates(chDup) {}
  TChain* events;
  float scale1fb;
  bool useWeigtFromBranch;
  bool checkDuplicates;
};
#endif


class ttDilCounts_looper {

	 public: 

//   static const unsigned int LOOSEDIL = 1<<10;
//   static const unsigned int LOOSEDIL_OS = (1<<10) + (1<<13);
//   static const unsigned int LOOSEDIL_OS_ONEWEIGHT = (1<<10) + (1<<13) + (1<<14);

  
  struct EIDiif {
    EIDiif():i0(0),i1(0),f0(0) {}
    EIDiif(int ai0, int ai1, float af0):i0(ai0),i1(ai1),f0(af0) {}
    bool operator==(const EIDiif& rhs) const {return (i0==rhs.i0 && i1==rhs.i1 && f0==rhs.f0);}
    bool operator<(const EIDiif& rhs) const {
      return (i0 != rhs.i0 ? i0 < rhs.i0 : (i1 != rhs.i1 ? i1 < rhs.i1 : (f0 != rhs.f0 ? f0 < rhs.f0 : false)));
    }
    int i0;
    int i1;
    float f0;
  };

  int ScanChain ( std::string fName, std::string prefix, float kFactor=1.0, int prescale=1, unsigned int cutsMask=31);
  int ScanChain ( TChain* chain, std::string prefix, float kFactor=1.0, int prescale=1, unsigned int cutsMask=31);
  int ScanChain ( std::vector<ProcDSChain>& pds, std::string prefix, float kFactor=1.0, int prescale=1, unsigned int cutsMask=31);
  void fill1D(TH1F* h, double val, double weight);

  TH1F* hnJet[4];       // Njet distributions
  TH1F* hnJetinZwindow[4];  //usefull for DY estimate
  TH1F* hnJetoutZwindow[4]; //usefull for DY estimate
  TH1F* helePt[4][6];      // electron Pt
  TH1F* hmuPt[4][6];       // muon Pt
  TH1F* hmuPtFromSilicon[4][6];    // muon Pt (from tracker)
  TH1F* hminLepPt[4][6];   // minimum lepton Pt
  TH1F* hmaxLepPt[4][6];   // maximum lepton Pt
  TH1F* helePhi[4][6];     // electron phi
  TH1F* hmuPhi[4][6];      // muon phi
  TH1F* hdphiLep[4][6];    // delta phi between leptons
  TH1F* heleEta[4][6];     // electron eta
  TH1F* hmuEta[4][6];      // muon eta
  TH1F* hdilMass[4][6];    // dilepton mass
  TH1F* hdilMassTightWindow[4][6]; // dilepton mass, but zooming around Z
  TH1F* hdilPt[4][6];       // dilepton Pt
  TH1F* hmet[4][6];       // MET
  TH1F* hmetPhi[4][6];       // MET phi
  TH1F* hpatmet[4][6];       // pat MET
  TH1F* hpatmetPhi[4][6];       // pat MET phi
  TH1F* htcmet[4][6];       // tc MET
  TH1F* htcmetPhi[4][6];       // tc MET phi
  TH1F* hpfmet[4][6];       // pf MET
  TH1F* hpfmetPhi[4][6];       // pf MET phi
  TH1F* hpfmetSpec[4][6];       // pf MET special (projected)

  TH2F* hmetVsDilepPt[4][6];  // MET vs dilepton Pt
  TH2F* hmetOverPtVsDphi[4][6]; // MET/Lepton Pt vs DeltaPhi between MET and Lepton Pt

  TH2F* hpatmetVsDilepPt[4][6];  // PAT MET vs dilepton Pt
  TH2F* hpatmetOverPtVsDphi[4][6]; // PAT MET/Lepton Pt vs DeltaPhi between MET and Lepton Pt

  TH2F* htcmetVsDilepPt[4][6];  // tc MET vs dilepton Pt
  TH2F* htcmetOverPtVsDphi[4][6]; // tc MET/Lepton Pt vs DeltaPhi between MET and Lepton Pt
  TH2F* hpfmetVsDilepPt[4][6];  // pf MET vs dilepton Pt
  TH2F* hpfmetOverPtVsDphi[4][6]; // pf MET/Lepton Pt vs DeltaPhi between MET and Lepton Pt

  TH2F* hdphillvsmll[4][6]; // delta phi between leptons vs dilepton mass
  TH1F* hptJet1[4][6];   // Pt of 1st jet
  TH1F* hptJet2[4][6];   // Pt of 2nd jet
  TH1F* hptJet3[4][6];   // Pt of 3rd jet
  TH1F* hptJet4[4][6];   // Pt of 4th jet
  TH1F* hetaJet1[4][6];   // eta of 1st jet
  TH1F* hetaJet2[4][6];   // eta of 2nd jet
  TH1F* hetaJet3[4][6];   // eta of 3rd jet
  TH1F* hetaJet4[4][6];   // eta of 4th jet
  
  TH1F* hSumJSpt[4][6]; 
  TH1F* hSumJSMTpt[4][6]; 
  TH1F* hSumJStcMTpt[4][6]; 
  TH1F* hvecSumJSpt[4][6]; 
  TH1F* hvecSumJSmLLpt[4][6]; 
  TH2F* hvecSumJSmLLptVspatmet[4][6]; 
  TH2F* hvecSumJSmLLptVstcmet[4][6]; 
  TH2F* hvecSumJSmLLptVspfmet[4][6]; 
 

  TH1F* numTightLep[4][6]; // number of tight leptons per event.
  TH1F* hmuSumIso[4][6];  // sum of trk pt, em et, had et in cone of 0.3
  TH1F* helSumIso[4][6];  // sum of trk pt, em et, had et in cone of 0.3
  TH1F* helRelIso[4][6]; //  Iso variable defined as pt/(pt+sum) for electron
  TH1F* hmuRelIso[4][6]; //  Iso variable defined as pt/(pt+sum) for muons
  TH1F* helRelIsoTrack[4][6]; //  Iso variable defined as pt/(pt+sum) for electron
  TH1F* hmuRelIsoTrack[4][6]; //  Iso variable defined as pt/(pt+sum) for muons
  TH1F* helRelIsoCalo[4][6]; //  Iso variable defined as pt/(pt+sum) for electron
  TH1F* hmuRelIsoCalo[4][6]; //  Iso variable defined as pt/(pt+sum) for muons

  // Unfortunately, our ntuple has no info for electron isolation other than candidate electrons.
  // When counting electrons, we thus can not apply an isolation criteria at this point !!!
  // For muons we only count good isolated muons.
  TH1F* hnJetLepVeto[4]; //njet distribution after requiring numTightLep < 3.

  int bigBlob[8000000];
  // The statement below should work but does not work due to bug in root when TH2 are also used
  // Rene Brun promised a fix.
  //TH1::SetDefaultSumw2(kTRUE); // do errors properly based on weights

  std::string compactConfig;
  void bookHistos(std::string& prefix);
};
#endif
