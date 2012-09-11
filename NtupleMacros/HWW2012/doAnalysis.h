#ifndef WW_doAnalysis_h
#define WW_doAnalysis_h
#include "Math/LorentzVector.h"
#include "Rtypes.h"
#include <vector>
#include <set>
#include "../HWW2012CORE/wwtypes.h"
#include "TChain.h"
#include <fstream>
#include <vector>
#include "../../../Smurf/Core/SmurfTree.h"

#include "../HWW2012CORE/analysisEnums.h"

bool hypo (int i_hyp, double weight, bool realData = false ); 

class TChain;
void ScanChain( TChain* chain, 
	   SmurfTree::DataType sample, 
	   double integratedLumi,
	   double xsec,
	   int nProcessedEvents,
	   bool identifyEvents,
	   bool realData = false,
	   TString cms2_json_file = "",
	   int beginrun=0, 
	   int endrun=999999);

void ProcessSample( std::string file_pattern, 
		    SmurfTree::DataType sample, 
		    double integratedLumi,
		    double xsec,
		    int nProcessedEvents,
		    bool identifyEvents = false,
		    bool realData = false,
		    TString cms2_json_file = "",
			int beginrun=0, 
			int endrun=999999);

void ProcessSample( std::vector<std::string> file_patterns, 
		    SmurfTree::DataType sample,
		    double integratedLumi,
		    double xsec,
		    int nProcessedEvents,
		    bool identifyEvents = false,
		    bool realData = false,
		    TString cms2_json_file = "",
			int beginrun=0, 
			int endrun=999999);

#endif
