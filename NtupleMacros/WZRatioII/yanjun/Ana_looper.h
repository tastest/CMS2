// -*- C++ -*-
#ifndef Ana_looper_H
#define Ana_looper_H
#include "TH1F.h"
#define NCHANNELS 4
#define NHISTS 2

using namespace std;

enum {
	kOneEleBit              = 0x1 <<  0, // 0x00000001                                                                                                          
	kOneMuonBit             = 0x1 <<  1, // 0x00000002                                                                                                          
	kTwoEleBit              = 0x1 <<  2, // 0x00000004                                                                                                          
	kTwoMuonBit             = 0x1 <<  3, // 0x00000008                                                                                                          
	
	kPassWMETBit            = 0x1 <<  4, // 0x00000010                                                                                                          
	kPassWMassBit           = 0x1 <<  5, // 0x00000020                                                                                                          
	kPassZVetoBit           = 0x1 <<  6,   // 0x00000040                                                                                                          
	
	kPassZMETBit            = 0x1 <<  7,   // 0x00000080                                                                                                          
	kPassZMassBit           = 0x1 <<  8,   // 0x00000100  
	kOppChargeBit           = 0x1 <<  9,   // 0x00000200 
	kZeeGenBit              = 0x1 <<  10,  // 0x00000400 
	kZmmGenBit              = 0x1 <<  11,  // 0x00000800 
	kWenuGenBit             = 0x1 <<  12,  // 0x00001000 
	kWmnuGenBit             = 0x1 <<  13,  // 0x00002000 
	kZttGenBit              = 0x1 <<  14  // 0x00004000 
  };   
class Ana_looper{
 public:
 
 
  int ScanChain( TChain* chain, int nEvents ,char* sample, float kFactor, int prescale); 

  void bookHistos(char* sample, int nchannels, int nhistsets);
  
  double Trans_W_Mass(TVector3& tcMET, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > & p4);
 
  TH1F* els_pt[NCHANNELS][NHISTS];
  TH1F* njets[NCHANNELS][NHISTS];

};
#endif
