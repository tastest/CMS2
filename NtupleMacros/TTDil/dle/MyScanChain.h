
#ifndef MYSCANCHAIN_H
#define MYSCANCHAIN_H

// C++ includes
#include <iostream>
#include <vector>



class TH1F;
class TChain;
class TDirectory;

class MyScanChain {

	public:

		MyScanChain() {};
		~MyScanChain() {};

		int ScanChain(bool isData, std::string sampleName, TChain *chain, int nEvents = -1, std::string skimFilePrefix="");

	private:

    	TH1F *h1_njets_[4];

		void Fill(TH1F** hist, unsigned int hyp, float val, float weight);
		void FormatHist(TDirectory *rootdir, TH1F** hist, std::string sampleName, std::string name, int n, float min, float max);

};

#endif

