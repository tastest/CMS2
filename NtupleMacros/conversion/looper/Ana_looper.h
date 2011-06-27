// -*- C++ -*-
#ifndef Ana_looper_H
#define Ana_looper_H
#include "TH1F.h"
#define NCHANNELS 1
#define NHISTS 7

using namespace std;

class Ana_looper{
 public:
 
 
  int ScanChain( TChain* chain, int nEvents ,char* sample, float kFactor, int prescale); 

  void bookHistos(char* sample, int nchannels, int nhistsets);
  
  bool isconversionElectron09(int elIdx);
  
  bool isconversionElectron_PIXHIT(int ele_index); 
  
  double dRBetweenVectors(LorentzVector v1, LorentzVector v2);

  std::pair<float, float> getConversionInfo(LorentzVector trk1_p4, 
					  int trk1_q, float trk1_d0, 
					  LorentzVector trk2_p4,
					  int trk2_q, float trk2_d0,
							float bField);

  TH1F* els_pt[NCHANNELS][NHISTS];
  TH1F* els_eta[NCHANNELS][NHISTS];
  TH1F* njets[NCHANNELS][NHISTS];

};
#endif
