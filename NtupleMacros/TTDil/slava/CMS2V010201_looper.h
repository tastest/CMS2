// -*- C++ -*-
#ifndef CMS2V010201_looper_H
#define CMS2V010201_looper_H
#include "CMS2V010201.h"
#include "TH1F.h"

class CMS2V010201_looper {

	 public: 
int ScanChain( TChain* chain, int nEvents=-1);
	 TH1F *samplehisto;
};
#endif
