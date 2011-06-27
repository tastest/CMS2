//void doAll(unsigned int bitmask, bool skipFWLite = false){
doAll(){

  float kElectron       = 1.; 
  float kGamma          = 1;

  // Prescales
 
  int preElectron       = 1;
  int preGamma          = 1;

 
  bool runElectron       = true;
  bool runGamma          = true;

  
using namespace std;
 gROOT->ProcessLine(".x setup.C(true)");
 gSystem->CompileMacro("Ana_looper.C", "++k", "libAna_looper");
 
 string dataset;
 if ( ! gSystem->Getenv("CMS2_NTUPLE_LOCATION") ){
   cout << "ERROR: Dataset location is not set. Please set CMS2_NTUPLE_LOCATION." <<endl;
   return;
 }
 dataset = gSystem->Getenv("CMS2_NTUPLE_LOCATION");

 
 TChain *chElectron = new TChain("Events");
 chElectron->Add("/hadoop/cms/store/user/yanjuntu/conversion_sample_cms2-V01-03-01/SingleElectron83055667ceda31da68d5aa532f72919d/preprocessing/ntuple*.root");
 TChain *chGamma = new TChain("Events");
 chGamma->Add("/hadoop/cms/store/user/yanjuntu/conversion_sample_cms2-V01-03-01/SingleGamma4016bf049eded18e8e9dac4c09773936/preprocessing/ntuple*.root");

 Ana_looper* looper = new Ana_looper();
 if (runElectron) {
    cout << "Processing single electron gun sample" << endl;
    looper->ScanChain(chElectron, -1, "Electron",kElectron,  preElectron);
    hist::color("Electron", kBlack);
  }
  if (runGamma) {
    cout << "Processing single gamma gun sample" << endl;
    looper->ScanChain(chGamma, -1, "Gamma",  kGamma,  preGamma);
    hist::color("Gamma", kMagenta);
  }
 const char* outFile = "myHist.root";
 hist::saveHist(outFile);
 hist::deleteHistos();
}
