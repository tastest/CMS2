// -*- C++ -*-
#ifndef Ana_looper_H
#define Ana_looper_H
#include "TH1F.h"
#define NCHANNELS 1
#define NHISTS 1

using namespace std;

class Ana_looper{
 public:
 
  Ana_looper();
  ~Ana_looper();

  int ScanChain( TChain* chain, int nEvents ,char* sample, float kFactor, int prescale,std::string  skimFilePrefix ); 

  void bookHistos(char* sample, int nchannels, int nhistsets);
 
  //  TH1F* els_pt[NCHANNELS][NHISTS];
  //TH1F* njets[NCHANNELS][NHISTS];
  TH1F* trkIso03[NCHANNELS][NHISTS];

  std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > *ran_trksp4_;
  std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > *trks_trk_p4_;
  std::vector<float> *ran_isoTrk03_mu_;
  TFile *outFile_;
  TTree *outTree_;

};
#endif
