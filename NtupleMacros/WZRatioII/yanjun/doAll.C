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

 
 TChain *chDYtautau = new TChain("Events");
 chDYtautau->Add((dataset+"/cms2-V01-03-01/Ztautau_M20_Summer08_IDEAL_V11_redigi_v1/merged_ntuple*.root").c_str());
 TChain *chDYee = new TChain("Events");
 chDYee->Add((dataset+"/cms2-V01-03-01/Zee_M20_Summer08_IDEAL_V11_redigi_v1-SingleLepton/merged_ntuple*.root").c_str());

 Ana_looper* looper = new Ana_looper();
 if (runDYtautau) {
    cout << "Processing DY->tautau" << endl;
    looper->ScanChain(chDYtautau, -1, "DYtautau",kDYtautau,  preDYtautau);
    hist::color("DYtautau", kBlack);
  }
  if (runDYee) {
    cout << "Processing DY->ee" << endl;
    looper->ScanChain(chDYee, -1, "DYee",  kDYee,  preDYee);
    hist::color("DYee", kMagenta);
  }
 const char* outFile = "myHist.root";
 hist::saveHist(outFile);
 hist::deleteHistos();
}
