// -*- C++ -*-
#ifndef Ana_looper_H
#define Ana_looper_H
#include "TH1F.h"
#define NCHANNELS 2
#define NHISTS 2

using namespace std;

class Ana_looper{
 public:
 
 
  int ScanChain( TChain* chain, int nEvents ,char* sample, float kFactor, int prescale); 

  void bookHistos(char* sample, int nchannels, int nhistsets);
 
  TH1F* els_pt[NCHANNELS][NHISTS];
  TH1F* njets[NCHANNELS][NHISTS];

};
#endif
