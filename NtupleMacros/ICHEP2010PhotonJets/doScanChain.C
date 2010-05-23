#include "TChain.h"
#include "ScanChain.C"
#include "histtools.C"

//void doAll() {
void doScanChain(string name, bool runonData=true, bool runonGEN=false) {

  //set_goodrun_file("goodruns_134987.txt");
  set_goodrun_file("goodruns.txt");

  //////DATA
  string dataEGMay6 = "/nfs-3/userdata/cms2/MinimumBias_Commissioning10-May6thPDSkim2_SD_EG-v1/V03-04-09/merged_ntuple*.root";
  string dataJetMETMay6 = "/nfs-3/userdata/cms2/MinimumBias_Commissioning10-May6thPDSkim2_SD_JetMETTau-v1/V03-04-09/merged_ntuple*.root";
  string dataJetMETMonMay6 = "/nfs-3/userdata/cms2/MinimumBias_Commissioning10-May6thPDSkim2_SD_JetMETTauMonitor-v1/V03-04-09/merged_ntuple*.root";
  string dataJetMETMonSomeMay6 = "/nfs-3/userdata/cms2/MinimumBias_Commissioning10-May6thPDSkim2_SD_JetMETTauMonitor-v1/V03-04-09/merged_ntuple_10*.root";

  string datastring = dataJetMETMonSomeMay6;

  vector<string> v_data;
  //v_data.push_back(datapunapr1);
  //v_data.push_back(dataEGMay6);

  TChain *ch_data = new TChain("Events");
  if( v_data.size() == 0 )
	ch_data->Add(datastring.c_str());
  else {
	datastring = "";
	for( unsigned int i=0; i<v_data.size(); i++ ) {
	  ch_data->Add(v_data[i].c_str());
	  datastring += v_data[i]+",";
	}
  }

  //////MC
  string v25b_all    = "/nfs-3/userdata/cms2/MinBias_Spring10-START3X_V25B_356ReReco-v1/V03-03-07/merged_ntuple*.root";
  string v25b_11     = "/nfs-3/userdata/cms2/MinBias_Spring10-START3X_V25B_356ReReco-v1/V03-03-07/merged_ntuple_2*.root";
  string v25b_single = "/nfs-3/userdata/cms2/MinBias_Spring10-START3X_V25B_356ReReco-v1/V03-03-07/merged_ntuple_22.root";

  string mcstring = v25b_all;

  TChain *ch_mc  = new TChain("Events");
  ch_mc->Add(mcstring.c_str());

  bool requireTrackCuts = true;

  //////RUN data
  if( runonData ) {
	cout << "\nScanChaining data:\n" << datastring << endl << endl;
	ScanChain(ch_data, false, requireTrackCuts); //2nd arg is run on mc
  }

  /////RUN (D)MC
  if( runonGEN ) {
	cout << "\nScanChaining MC:\n" << mcstring << endl << endl;
	ScanChain(ch_mc, true, requireTrackCuts);
  }

  saveHist((name+".root").c_str());

}  
  
