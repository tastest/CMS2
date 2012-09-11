#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TRegexp.h"
#include "TFile.h" 
#include "../HWW2012CORE/wwtypes.h"
#include "LeptonTreeMaker.h"
#include "SmurfDataTypes.h"
#include "processLeptonTree.h"
#include <stdlib.h>

int main(int argc, char *argv[])
{
	// ------------------------------------------------------------------------
	// 		For BatchSubmission
	// ------------------------------------------------------------------------

	if (argc == 3) {
		SmurfTree::DataType dataType = (SmurfTree::DataType)atoi(argv[1]);
		std::string dataFile = argv[2];
		bool realData = false;
		if (dataType == 0) realData = true;
		const std::string cms2_json_file = "./files/Cert_190456-195396_8TeV_PromptReco_Collisions12_JSON_v2.jmu";
		processLeptonTree("test", dataType, dataFile, realData, "");
		return 0;
	}
	
	// ------------------------------------------------------------------------
	//  	For local submission
	// ------------------------------------------------------------------------
	bool dodata=1;
	bool domc=1;

	// 
	// Data
	//
	if(dodata)
	{
		TString goodrunlist = "Cert_190456-194479_8TeV_PromptReco_Collisions12_JSON_cms2.txt";
		//processLeptonTree("test", SmurfTree::data, "/hadoop/cms/store/user/yanjuntu/CMSSW_5_2_3_patch4_V05-02-27/SingleElectron_Run2012B-PromptReco-v1_AOD/merged/merged_ntuple_195552_0.root", true, "");
		processLeptonTree("test", SmurfTree::data, "/hadoop/cms/store/user/yanjuntu/CMSSW_5_2_3_patch4_V05-02-27/SingleElectron_Run2012B-PromptReco-v1_AOD/merged/merged_ntuple_195774_0.root", true, "");
	}

	// 
	// MC
	// 
	if(domc)
	{
		processLeptonTree("test1", SmurfTree::dyee, "/hadoop/cms/store/group/snt/papers2012/Summer12MC/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12-PU_S7_START52_V9-v2/V05-02-27/merged_ntuple_300.root", false, "");
		//processLeptonTree("test", SmurfTree::dyee, "/home/users/jaehyeok/HWW/MakeCMS2ntuples/ICHEP2012/CMSSW_5_2_3_patch4_V05-02-28/src/post_ntuple.root", false, "");
	}
	
	return 0; 
}
