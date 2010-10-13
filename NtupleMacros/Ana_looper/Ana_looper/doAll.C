//void doAll(unsigned int bitmask, bool skipFWLite = false){
doAll(){
 float kWW       = 1.;
  float kWZ       = 1.;
  float kZZ       = 1.;
  float kWjets    = 1.; 
  float kWcharm   = 1;
  float kDYee     = 1.;  
  float kDYmm     = 1.;  
  float kDYtautau = 1.; 
  float kppMuX    = 1.; 
  float kEM       = 1.;
  float ktW       = 1.; 
  float kVQQ      = 1;

  // Prescales
  int prettdil    = 1;
  int prettotr    = 1;
  int preWW       = 1;
  int preWZ       = 1;
  int preZZ       = 1;
  int preWjets    = 1;
  int preWcharm   = 1;
  int preDYee     = 1;
  int preDYmm     = 1;
  int preDYtautau = 1;
  int preppMuX    = 1;
  int preEM       = 1;
  int pretW       = 1;
  int preVQQ      = 1;

  bool runttdil    = true;
  bool runttotr    = true;
  bool runWW       = true;
  bool runWZ       = true;
  bool runZZ       = true;
  bool runWjets    = true;
  bool runWcharm   = true;
  bool runDYee     = true;
  bool runDYmm     = true;
  bool runDYtautau = true;
  bool runppMuX    = true;
  bool runEM       = true;
  bool runtW       = true;
  bool runVQQ      = true;

  
using namespace std;
 gROOT->ProcessLine(".x setup.C(true)");
 gSystem->CompileMacro("Ana_looper.C", "++k", "libAna_looper");
 
 string dataset;
 if ( ! gSystem->Getenv("CMS2_NTUPLE_LOCATION") ){
   cout << "ERROR: Dataset location is not set. Please set CMS2_NTUPLE_LOCATION." <<endl;
   return;
 }
 dataset = gSystem->Getenv("CMS2_NTUPLE_LOCATION");

 TChain *MinBias_skimdata = new TChain("Events");
 TChain *MinBias_skimdata_2360 = new TChain("Events");
 TChain *MinBias_mc_2360 = new TChain("Events");
 TChain *MinBias_mc = new TChain("Events");
 
 MinBias_skimdata->Add("/data/tmp/yanjuntu/cms2/cms2-V03-00-23/MinimumBias_BeamCommissioning09-BSCNOBEAMHALO-Dec14thSkim_v1/filtered*900.root"); 
 MinBias_mc->Add("/data/tmp/yanjuntu/cms2/cms2-V03-00-20/MinBias_Summer09-STARTUP3X_V8I_900GeV-v2/filtered*.root"); 
 MinBias_mc_2360->Add("/data/tmp/yanjuntu/cms2/cms2-V03-00-23/MinBias_Summer09-STARTUP3X_V8D_2360GeV-v2/filtered*.root"); 
 MinBias_skimdata_2360->Add("/data/tmp/yanjuntu/cms2/cms2-V03-00-23/MinimumBias_BeamCommissioning09-BSCNOBEAMHALO-Dec14thSkim_v1/filtered*2360.root"); 

 Ana_looper* looper = new Ana_looper();

 cout << "Processing MinBias_skimdata" << endl;
 // looper->ScanChain(MinBias_skimdata, -1, "Data",1,  1, "DataCorrZ");
 hist::color("Data", kBlack);
 
 cout << "Processing MinBias_skimdata at 2360 GeV" << endl;
 looper->ScanChain(MinBias_skimdata_2360, -1, "Data",1,  1, "DataCorrZ2360");
 hist::color("Data_2360", kBlack);
 
 cout << "Processing MinBias_mc" << endl;
 //looper->ScanChain(MinBias_mc, -1, "MC",  1,  1, "MCCorrZ");
 hist::color("MinBias_mc", kMagenta);

 cout << "Processing MinBias_mc at 2360 GeV" << endl;
 //looper->ScanChain(MinBias_mc_2360, -1, "MC",  1,  1, "MCCorrZ2360");
 hist::color("MinBias_mc_2360", kMagenta);
 
 const char* outFile = "myHist.root";
 hist::saveHist(outFile);
 hist::deleteHistos();
}
