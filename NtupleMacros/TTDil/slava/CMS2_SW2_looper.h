// -*- C++ -*-
#ifndef CMS2_SW2_looper_H
#define CMS2_SW2_looper_H
#include "CMS2_SW2.h"
#include "TH1F.h"

class CMS2_SW2_looper {

	 public: 
int ScanChain( TChain* chain, int nEvents=-1);
	 TH1F *samplehisto;
};
#endif
