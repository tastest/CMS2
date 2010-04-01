
#ifndef MYSCANCHAIN_H
#define MYSCANCHAIN_H

// C++ includes
#include <iostream>
#include <vector>

// looper includes
#include "EffMulti.h"
#include "Math/VectorUtil.h"


// typedefs
typedef UInt_t      uint32;
typedef ULong64_t   uint64;
typedef uint64  	cuts_t;

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

//
// for det types
//
enum DetType {
    DET_EE = 0,
    DET_EB = 1,
    DET_ALL = 2,
};
static const char det_names[][128] = { "ee", "eb", "all" };

// forward definitions
class TH1F;
class TH2F;
class TChain;
class TDirectory;

class MyScanChain {

	public:

		MyScanChain() {};
		~MyScanChain() {};

		int ScanChain(bool isData, std::string sampleName, TChain *chain, int nEvents = -1, std::string skimFilePrefix="");

	private:

		// get subdetector for histogram filling
		enum DetType getSubdet(int eleIndex);

		// N-1
		// return true if the cuts to apply - the cuts to remove were passed
		// in the cuts that "passed"
		bool CheckCutsNM1(cuts_t apply, cuts_t remove, cuts_t passed);

		// Simple check if the desired cuts to apply are set in the 
		// cuts that "passed"
		bool CheckCuts(cuts_t apply, cuts_t passed);

		// do stuff with histogram
		void Format2DHist(TH2F** hist, std::string name, Int_t nx, Float_t minx, Float_t maxx, Int_t ny, Float_t miny, Float_t maxy);
		void FormatHist(TH1F** hist, std::string name, Int_t n, Float_t min, Float_t max);
		void FormatEffHist(EffMulti** hist, bool lessThan, float thresholdEB, float ThresholdEE, std::string name);

		// sample name
		std::string sampleName_;

		// general  
		//
		TH1F    *h1_pt_[3];
		TH1F    *h1_eta_[3];
		TH1F    *h1_phi_[3];

		TH1F 	*h1_nm1_tcmet_[3];
        TH1F    *h1_nm1_pfmet_[3];
        TH1F    *h1_nm1_jetveto_[3];
        TH1F    *h1_nm1_iso_[3];
        TH1F    *h1_nm1_secondpt_[3];

};

#endif

