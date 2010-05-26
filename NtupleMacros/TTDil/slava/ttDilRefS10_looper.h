// -*- C++ -*-
#ifndef ttDilRefS10_looper_H
#define ttDilRefS10_looper_H
//#include "CORE/CMS2.h"
#include "TH1F.h"
#include "TH2F.h"
#include <vector>
#include "TChain.h"

#include "ProcDSS.h"

class ttDilRefS10_looper {
  
public: 
  
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
  
  int ScanChain ( std::string fName, std::string prefix, float kFactor=1.0, int prescale=1, 
		  unsigned long long int cutsMask=0);
  int ScanChain ( TChain* chain, std::string prefix, float kFactor=1.0, int prescale=1, 
		  unsigned long long int cutsMask=0);
  int ScanChain ( ProcDSS& pds, unsigned long long int cutsMask);
  void fill1D(TH1F* h, double val, double weight);

  TH1F* hnJet[4];       // Njet distributions

  std::string compactConfig;
  void bookHistos(std::string& prefix);
};
#endif
